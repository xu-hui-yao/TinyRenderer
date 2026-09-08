// Real-time path tracing kernel: the front end of the interactive renderer.
//
// Same integrator as the offline megakernel (path_trace.slang) - both call
// common/scene_common.slang's trace_path() - but structured for a persistent
// session that renders ONE frame at a time:
//
//   - RNG seeded by a frame counter instead of an spp index, so successive
//     frames get decorrelated subpixel jitter (free TAA-quality AA once the
//     temporal accumulation pass blends them).
//   - `spp_per_frame` samples accumulated in-kernel, so the host issues
//     exactly one dispatch per frame regardless of samples/frame.
//   - Radiance clamping, the single most effective firefly suppression for
//     1spp rendering (a single sample is a wildly high-variance estimator,
//     so one unlucky light path would otherwise detonate a whole neighbourhood
//     once the spatial filter starts gathering across it).
//
// ON TOP of the offline kernel it also writes what the denoiser needs (P2):
//
//   - a G-Buffer (albedo / normal+depth / world position+material id) from the
//     PRIMARY hit only - the surface the pixel actually sees
//   - DEMODULATED illumination, ill um = radiance / albedo, plus the first two
//     luminance moments
//
// Demodulation is what lets the filter smooth illumination without smearing
// texture: the high-frequency albedo detail is divided out before filtering and
// multiplied back at the end. Albedo is deliberately left at (1,1,1) for
// specular/delta materials and emitters, whose radiance carries no surface
// texture (dividing by a near-black albedo there would blow up instead).
//
// The (sum, count) accumulation buffer below is kept in parallel so that
// demodulation can be disabled for validation: with demodulate == 0 the kernel
// writes radiance straight through and must stay bit-identical to the offline
// megakernel, which is what src/gpu/interactive_verify_main.cpp checks.

#include "common/scene_common.slang"

// ---- Storage images: the denoiser's screen-sized buffers. See
// include/gpu/gpu_image.h for why these are images and not buffers. ----
[[vk::binding(16, 0)]] RWTexture2D<float4> img_albedo;       // RGBA16F: diffuse albedo (1,1,1) if not demodulable
[[vk::binding(17, 0)]] RWTexture2D<float4> img_normal_depth; // RGBA32F: (oct_normal.xy, linear_z, depth_gradient)
[[vk::binding(18, 0)]] RWTexture2D<float4> img_pos_id;       // RGBA32F: (world_pos.xyz, asfloat(material_id))
[[vk::binding(19, 0)]] RWTexture2D<float4> img_illum;        // RGBA16F: demodulated illumination for this frame
[[vk::binding(20, 0)]] RWTexture2D<float2> img_moments;      // RG32F:   (luminance, luminance^2)

[[vk::binding(0, 0)]] RWStructuredBuffer<float4> accum_buffer; // running sum of radiance (.xyz) + sample count (.w)

struct PushConstants {
    uint frame_index;   // decorrelates RNG / jitter between frames
    uint max_depth;
    uint rr_depth;
    uint num_lights;
    uint environment_light_id; // INVALID_INDEX if none
    uint accumulate;    // 0 => start a fresh average (camera moved), 1 => add to it
    float radiance_clamp;
    uint spp_per_frame;
    uint demodulate;    // 0 => write raw radiance, albedo = 1 (validation / debug)
};
[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

[shader("compute")]
[numthreads(8, 8, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= camera.width || tid.y >= camera.height) return;

    uint pixel_index = tid.y * camera.width + tid.x;

    PathConfig cfg;
    cfg.max_depth = pc.max_depth;
    cfg.rr_depth = pc.rr_depth;
    cfg.num_lights = pc.num_lights;
    cfg.environment_light_id = pc.environment_light_id;

    float3 sum = float3(0, 0, 0);
    // G-Buffer features come from the FIRST sample's primary hit. With
    // spp_per_frame > 1 the samples jitter within the pixel and can land on
    // different surfaces; any of them is a valid representative, and jittering
    // the choice per-sample would make the G-Buffer itself noisy - which is
    // worse, since the filters' edge-stopping functions depend on it.
    PrimaryHit primary;
    bool primary_valid = false;

    for (uint s = 0u; s < pc.spp_per_frame; ++s) {
        // `s` is folded into the RNG salt so the samples within one frame are
        // decorrelated from each other as well as from other frames.
        Rng rng = rng_init(pixel_index, pc.frame_index, s);
        float2 jitter = rng_next2(rng);

        float3 ray_o, ray_d;
        float ray_tmin, ray_tmax;
        camera_ray(float(tid.x) + jitter.x, float(tid.y) + jitter.y, ray_o, ray_d, ray_tmin, ray_tmax);

        PrimaryHit ph;
        float3 result = trace_path(ray_o, ray_d, ray_tmin, ray_tmax, rng, cfg, ph);

        // Clamp per-sample radiance. Applied before accumulation so a single
        // firefly cannot poison the running sum that history reuse will keep
        // reading for many subsequent frames.
        if (pc.radiance_clamp > 0.0) {
            float lum = max(max(result.x, result.y), result.z);
            if (lum > pc.radiance_clamp) result *= pc.radiance_clamp / lum;
        }
        sum += result;

        if (s == 0u) { primary = ph; primary_valid = ph.valid; }
    }

    float4 prev = (pc.accumulate != 0u) ? accum_buffer[pixel_index] : float4(0, 0, 0, 0);
    accum_buffer[pixel_index] = prev + float4(sum, float(pc.spp_per_frame));

    // ---- G-Buffer + demodulated illumination ----
    float3 albedo = float3(1, 1, 1);
    if (primary_valid && pc.demodulate != 0u)
        albedo = max(primary.albedo, float3(0.01, 0.01, 0.01)); // floor avoids 1/albedo blowups

    float3 illum = sum / float(pc.spp_per_frame);
    illum /= albedo;

    float lum = dot(illum, float3(0.212671, 0.715160, 0.072169));

    // NOTE: the depth gradient the temporal pass needs is NOT computed here.
    // ddx()/ddy() are fragment-shader derivatives and require quad
    // neighbourhoods; Metal (via MoltenVK) rejects them outright in a compute
    // kernel. It is derived by finite differences in svgf_temporal.slang
    // instead, which is where SVGF computes it anyway.
    img_albedo[tid.xy]       = float4(albedo, 1.0);
    img_normal_depth[tid.xy] = float4(oct_encode(primary.normal), primary.linear_depth, 1.0);
    img_pos_id[tid.xy]       = float4(primary.position, asfloat(primary_valid ? primary.material_id : INVALID_INDEX));
    img_illum[tid.xy]        = float4(illum, 1.0);
    img_moments[tid.xy]      = float2(lum, lum * lum);
}
