// Stage 1, checkpoint 3: full megakernel path tracer.
//
// This is a line-by-line port of src/integrators/path.cpp's li() loop. All of
// the machinery it needs - GPU-side re-implementations of every BSDF in
// src/bsdfs/*.cpp, every texture in src/textures/*.cpp, area-light sampling
// matching src/components/mesh.cpp's Mesh::sample_direction()/pdf_direction(),
// camera ray generation and BVH traversal - now lives in
// common/scene_common.slang, which is shared verbatim with the wavefront
// backend (common/wavefront_common.slang) and the real-time backend
// (rt_path_trace.slang). See that file's header for the CPU-parity notes
// (RNG choice, envmap sampling simplification, the deliberately un-"fixed"
// square_to_cosine_hemisphere_pdf, etc.) - they all still apply here, because
// the math was moved, not changed.
//
// What remains in THIS file is only what is specific to the megakernel
// backend: its output buffer (binding 0), its push-constant layout, and the
// full bounce loop in main().

#include "common/scene_common.slang"

// ============================================================================
// Megakernel-specific bindings + push constants
// ============================================================================
// Bindings 1..15 (the shared flattened scene) come from scene_common.slang.
// Binding 0 is this kernel's own output.

[[vk::binding(0, 0)]] RWStructuredBuffer<float4> accum_buffer; // running sum of radiance (.xyz) + sample count (.w), divided on readback

struct PushConstants {
    uint spp_index;    // which sample-per-pixel pass this dispatch computes (used to seed the RNG)
    uint max_depth;
    uint rr_depth;
    uint num_lights;
    uint environment_light_id; // INVALID_INDEX if none
    uint accumulate;   // 0 on the very first dispatch (clear), 1 otherwise (add)
};
[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

// ============================================================================
// One dispatch = one sample per pixel. The bounce loop itself lives in
// common/scene_common.slang's trace_path() (shared with rt_path_trace.slang);
// all this kernel does is generate a jittered primary ray, seed the RNG from
// its spp index, and accumulate the result.
// ============================================================================

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= camera.width || tid.y >= camera.height) return;

    uint pixel_index = tid.y * camera.width + tid.x;
    Rng rng = rng_init(pixel_index, pc.spp_index, 0u);

    // Random subpixel jitter for antialiasing across spp - the CPU integrator
    // gets this from its reconstruction filter's sampler-driven pixel offset.
    float2 jitter = rng_next2(rng);

    float3 ray_o, ray_d;
    float ray_tmin, ray_tmax;
    camera_ray(float(tid.x) + jitter.x, float(tid.y) + jitter.y, ray_o, ray_d, ray_tmin, ray_tmax);

    PathConfig cfg;
    cfg.max_depth = pc.max_depth;
    cfg.rr_depth = pc.rr_depth;
    cfg.num_lights = pc.num_lights;
    cfg.environment_light_id = pc.environment_light_id;

    float3 result = trace_path(ray_o, ray_d, ray_tmin, ray_tmax, rng, cfg);

    // ---- Accumulate (mirrors ImageBlock's running-average accumulation) ----
    float4 prev = (pc.accumulate != 0u) ? accum_buffer[pixel_index] : float4(0, 0, 0, 0);
    accum_buffer[pixel_index] = prev + float4(result, 1.0);
}
