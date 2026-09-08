// Display pass: remodulate the filtered illumination, tonemap, and pack to
// RGBA8 for the host to blit.
//
// Two jobs:
//   1. Undo the demodulation rt_path_trace.slang applied (multiply the filtered
//      illumination by the albedo that was divided out). This is what keeps
//      texture detail crisp: the filters only ever smoothed illumination, so
//      nothing here reintroduces blur.
//   2. Tonemap + encode to sRGB, on the GPU. The host then uploads the result
//      as-is with zero per-pixel work - the offline preview path's CPU
//      convert_to_ldr() is what would otherwise blow the frame budget.
//
// Changing anything in the Display group (exposure / tonemap / gamma) only
// affects THIS pass, so the accumulation is never invalidated: dragging the
// exposure slider at 1000 accumulated spp costs one cheap dispatch.
//
// The debug views are not a nicety - SVGF has half a dozen coupled parameters
// and the only practical way to tune them (or to tell "wrong filter settings"
// apart from "wrong G-Buffer") is to look at each intermediate buffer.

// math_common has no bindings; scene_common must NOT be included here (it binds
// the scene on 1..7, colliding with the textures below).
#include "common/math_common.slang"

[[vk::binding(0, 0)]] RWStructuredBuffer<uint> display_buffer;

[[vk::binding(1, 0)]] Texture2D<float4> img_illum;        // filtered, still demodulated
[[vk::binding(2, 0)]] Texture2D<float4> img_albedo;
[[vk::binding(3, 0)]] Texture2D<float4> img_normal_depth; // (oct normal.xy, linear_z, depth_grad)
[[vk::binding(4, 0)]] Texture2D<float2> img_moments;
[[vk::binding(5, 0)]] Texture2D<float> img_hist_len;
[[vk::binding(6, 0)]] Texture2D<float4> img_pos_id;
[[vk::binding(7, 0)]] Texture2D<float4> img_illum_raw;    // unfiltered (this frame's new sample)

// The output is a flat RWStructuredBuffer, so the kernel CANNOT recover the
// image dimensions itself: StructuredBuffer::GetDimensions() returns
// (element count, stride) = (width*height, 4), not (width, height). Using it
// as if it were a texture silently computes pixel_index = y*921600 + x, which
// writes far past the end of the buffer and corrupts whatever allocation
// follows it - which is exactly what happened to the accumulation buffer here.
// Dimensions must therefore come in as push constants.
struct PushConstants {
    float exposure;
    uint tonemap;        // 0 = clamp, 1 = Reinhard, 2 = ACES
    uint use_srgb;       // 1 = sRGB transfer function (matches save_png), 0 = plain pow(1/gamma)
    float gamma;         // only used when use_srgb == 0
    uint modulate;       // 1 = multiply by albedo (normal path); 0 = show raw illumination
    uint debug_view;     // see DEBUG_* below; 0 is the final image
    uint width;
    uint height;
    float variance_scale; // display gain for the variance debug view
    // 1 = read this frame's unfiltered sample instead of the filtered history.
    // This is the "denoiser off" path - note it reads the RAW buffer rather
    // than skipping the temporal dispatch, because when the temporal pass is
    // disabled nothing is ever written to the history pair the display pass
    // would otherwise read.
    uint use_raw_illum;
};

static const uint DEBUG_FINAL          = 0;
static const uint DEBUG_ILLUM_RAW      = 1;
static const uint DEBUG_ILLUM_FILTERED = 2;
static const uint DEBUG_ALBEDO         = 3;
static const uint DEBUG_NORMAL         = 4;
static const uint DEBUG_DEPTH          = 5;
static const uint DEBUG_VARIANCE       = 6;
static const uint DEBUG_HIST_LEN       = 7;
static const uint DEBUG_POSITION       = 8;

[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

// ACES filmic approximation (Narkowicz 2015), the same curve the CPU-side
// Bitmap::save_png() uses for its ACES tonemap mode.
float3 tonemap_aces(float3 x)
{
    const float a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
    return clamp((x * (a * x + b)) / (x * (c * x + d) + e), 0.0, 1.0);
}

float linear_to_srgb(float c)
{
    // Matches Bitmap::save_png()'s transfer function.
    return c <= 0.0031308 ? c * 12.92 : 1.055 * pow(c, 1.0 / 2.4) - 0.055;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= pc.width || tid.y >= pc.height) return;

    uint pixel_index = tid.y * pc.width + tid.x;
    const int2 ip = int2(tid.xy);

    float3 color;

    if (pc.debug_view == DEBUG_ILLUM_RAW) {
        color = img_illum_raw[ip].rgb;
    } else if (pc.debug_view == DEBUG_ILLUM_FILTERED) {
        color = img_illum[ip].rgb;
    } else if (pc.debug_view == DEBUG_ALBEDO) {
        color = img_albedo[ip].rgb;
    } else if (pc.debug_view == DEBUG_NORMAL) {
        color = oct_decode(img_normal_depth[ip].xy) * 0.5 + 0.5;
    } else if (pc.debug_view == DEBUG_DEPTH) {
        // Normalised to a fixed range; the absolute scale is scene-dependent,
        // so this view is for judging SHAPES and discontinuities, not distances.
        color = float3(img_normal_depth[ip].z * 0.05);
    } else if (pc.debug_view == DEBUG_VARIANCE) {
        float2 m = img_moments[ip];
        float var = max(m.y - m.x * m.x, 0.0);
        color = float3(sqrt(var) * pc.variance_scale);
    } else if (pc.debug_view == DEBUG_HIST_LEN) {
        // Ramp: dark = freshly reset, bright = long converged history.
        float h = img_hist_len[ip] / 64.0;
        color = float3(h, h, h);
    } else if (pc.debug_view == DEBUG_POSITION) {
        color = frac(img_pos_id[ip].xyz * 0.25);
    } else {
        // ---- Final image: remodulate, then tonemap ----
        color = (pc.use_raw_illum != 0u ? img_illum_raw : img_illum)[ip].rgb;
        if (pc.modulate != 0u)
            color *= img_albedo[ip].rgb;
    }

    // Debug views bypass tone mapping: they are meant to be inspected raw, and
    // running them through the filmic curve would hide exactly the
    // out-of-range values one is usually looking for.
    if (pc.debug_view == DEBUG_FINAL) {
        color *= exp2(pc.exposure);

        if (pc.tonemap == 1) {
            color = color / (1.0 + color); // Reinhard
        } else if (pc.tonemap == 2) {
            color = tonemap_aces(color);
        } else {
            color = clamp(color, 0.0, 1.0);
        }

        if (pc.use_srgb != 0u) {
            color = float3(linear_to_srgb(color.r), linear_to_srgb(color.g), linear_to_srgb(color.b));
        } else {
            float inv_gamma = 1.0 / max(pc.gamma, 1e-3);
            color = pow(max(color, float3(0, 0, 0)), float3(inv_gamma, inv_gamma, inv_gamma));
        }
    }

    color = clamp(color, 0.0, 1.0);

    uint r = uint(color.r * 255.0 + 0.5);
    uint g = uint(color.g * 255.0 + 0.5);
    uint b = uint(color.b * 255.0 + 0.5);
    uint a = 255u;

    display_buffer[pixel_index] = (a << 24) | (b << 16) | (g << 8) | r;
}
