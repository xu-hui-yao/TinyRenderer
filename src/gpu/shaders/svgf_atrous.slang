// Spatiotemporal Variance-Guided Filtering: SPATIAL pass (a-trous edge-avoiding
// filter). See Schied et al. 2017, "Spatiotemporal Variance-Guided Filtering".
//
// The temporal pass (svgf_temporal.slang) already cut noise by accumulating
// history, but at interactive sample counts the result is still far too noisy to
// look at, and every pixel only has its own history to average with. This pass
// gathers across NEIGHBOURING pixels to trade the residual noise for a little
// blur - which is only acceptable if the blur does not cross surface
// boundaries, hence the three edge-stopping functions below.
//
// The a-trous wavelet structure is what makes it affordable: instead of a
// growing kernel, each iteration keeps a fixed 5x5 stencil and doubles the
// STRIDE. Iteration i therefore gathers over a (4*2^i + 1)^2 neighbourhood for
// O(iterations * 25) taps total, rather than the O(N^2) of one big kernel.
//
// Every weight is guided by the per-pixel variance the temporal pass estimated:
// where the signal is already converged, luminance differences are trusted to
// mean "real detail" and the filter keeps its hands off; where variance is high,
// differences are assumed to be noise and smoothing proceeds. That is the whole
// point of computing the second moment back in rt_path_trace.slang.
//
// Deliberately includes math_common (no bindings) rather than scene_common
// (binds the scene on 1..15), which would collide with the textures below.

#include "common/math_common.slang"

[[vk::binding(0, 0)]] RWTexture2D<float4> img_illum_out; // filtered illumination
[[vk::binding(1, 0)]] Texture2D<float4> img_illum_in;   // previous iteration (or the temporal output)
[[vk::binding(2, 0)]] Texture2D<float2> img_moments;    // (luminance, luminance^2) for variance guidance
[[vk::binding(3, 0)]] Texture2D<float4> img_normal_depth;
[[vk::binding(4, 0)]] Texture2D<float4> img_pos_id;

struct PushConstants {
    uint width;
    uint height;
    uint stride;       // 2^iteration_index; the host passes it explicitly
    float phi_color;   // luminance edge-stopping scale (sigma_l)
    float phi_normal;  // normal edge-stopping EXPONENT (sigma_n)
    float phi_depth;   // depth edge-stopping scale (sigma_z)
    uint iteration;    // 0-based; used only for diagnostics
    uint debug_mode;   // reserved
};
[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

// Separable 5-tap binomial kernel: [1/16, 1/4, 3/8, 1/4, 1/16].
// Indexed by |offset| (0, 1, 2).
float kernel1d(int d)
{
    if (d == 0) return 3.0 / 8.0;
    if (d == 1) return 1.0 / 4.0;
    return 1.0 / 16.0;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= pc.width || tid.y >= pc.height) return;

    const int2 ip = int2(tid.xy);
    const int2 dim = int2(pc.width, pc.height);
    const int stride = int(pc.stride);

    const float4 center_illum = img_illum_in[ip];
    const float4 center_nd    = img_normal_depth[ip];
    const float4 center_pos   = img_pos_id[ip];

    const float center_z     = center_nd.z;
    const float3 center_n    = oct_decode(center_nd.xy);
    const float center_lum   = luminance(center_illum.rgb);
    const float center_matid = center_pos.w;

    // ---- Variance from the moments, smoothed with a 3x3 Gaussian ----
    // Filtering the moments BEFORE taking E[x^2] - E[x]^2 matters: the raw
    // per-pixel variance is itself extremely noisy at 1spp, and using it
    // unsmoothed makes the filter's behaviour flicker pixel to pixel.
    float2 mom_sum = float2(0, 0);
    float mom_wsum = 0.0;
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            int2 q = ip + int2(dx, dy);
            if (any(q < int2(0, 0)) || any(q >= dim)) continue;
            float w = (dx == 0 ? 2.0 : 1.0) * (dy == 0 ? 2.0 : 1.0); // separable [1 2 1]
            mom_sum += img_moments[q] * w;
            mom_wsum += w;
        }
    }
    mom_sum /= max(mom_wsum, 1e-6);
    float variance = max(mom_sum.y - mom_sum.x * mom_sum.x, 0.0);
    float sigma_l = sqrt(variance);

    // ---- A-trous gather ----
    float3 color_sum = float3(0, 0, 0);
    float wsum = 0.0;

    for (int dy = -2; dy <= 2; ++dy) {
        for (int dx = -2; dx <= 2; ++dx) {
            int2 q = ip + int2(dx, dy) * stride;

            // Clamp to the border rather than skipping: dropping off-screen taps
            // would darken edge pixels, since the weights would no longer sum to 1.
            q = clamp(q, int2(0, 0), dim - int2(1, 1));

            const float4 q_nd  = img_normal_depth[q];
            const float4 q_pos = img_pos_id[q];
            const float3 q_col = img_illum_in[q].rgb;

            float w = kernel1d(abs(dx)) * kernel1d(abs(dy));
            if (w <= 0.0) continue;

            // -- Normal stopping --
            // A large exponent makes the falloff very sharp, so only genuinely
            // coplanar neighbours contribute. This is the strongest of the three
            // cues and what keeps geometric edges crisp.
            float3 q_n = oct_decode(q_nd.xy);
            float ndot = max(dot(center_n, q_n), 0.0);
            w *= pow(ndot, pc.phi_normal);

            // -- Depth stopping --
            // Tolerance scales with the surface's depth gradient, so a sloped
            // surface (where depth changes a lot across a pixel) is not rejected
            // just because the neighbour is genuinely further away.
            float grad_z = length(float2(dx, dy)) * stride;
            float z_tol = pc.phi_depth * grad_z + 1e-3;
            w *= exp(-abs(center_z - q_nd.z) / z_tol);

            // -- Luminance stopping (variance guided) --
            // Divide by the local standard deviation: in a high-variance region
            // a given luminance step is probably noise and gets smoothed, while
            // in a low-variance region the same step is probably real detail and
            // is preserved.
            float l_diff = abs(center_lum - luminance(q_col));
            w *= exp(-l_diff / max(pc.phi_color * sigma_l, 1e-4));

            // -- Material id --
            // Cheap extra guard: never blend across different materials, which
            // would otherwise show up as colour bleeding at object silhouettes.
            if (abs(q_pos.w - center_matid) > 0.5) w *= 0.0;

            color_sum += q_col * w;
            wsum += w;
        }
    }

    float3 result = wsum > 1e-6 ? color_sum / wsum : center_illum.rgb;

    img_illum_out[ip] = float4(result, 1.0);
}
