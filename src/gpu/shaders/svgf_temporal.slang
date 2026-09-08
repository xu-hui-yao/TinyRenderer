// Spatiotemporal Variance-Guided Filtering: TEMPORAL accumulation pass.
//
// Accumulates the current frame's demodulated illumination into a history
// buffer by reprojecting each pixel to where its surface was in the PREVIOUS
// frame, and validates that the history sample really belongs to the same
// surface before blending it in (Schied et al. 2017, "Spatiotemporal
// Variance-Guided Filtering", section 4 - plus the widely used neighbourhood
// clamping refinement from subsequent real-time denoiser literature).
//
// Reprojection is done from the G-Buffer world position rather than from
// per-object motion vectors, because this renderer's scenes are static: only
// the camera moves, so the previous frame's view-projection matrix is enough to
// recover the motion exactly. That removes the whole object-velocity machinery
// at no cost in quality (and would be the thing to revisit if animated geometry
// is ever added).
//
// Inputs/outputs are ping-ponged by the host: it binds set 0 or set 1 depending
// on the frame parity, so this kernel always reads `*_in` and writes `*_out`.

// Includes math_common (no bindings) rather than scene_common (binds the scene
// on 1..15): this pass needs transform_point() but declares its own textures on
// those very slots, so including the scene would collide with them.

#include "common/math_common.slang"

// ---- Outputs (this frame's new history) ----
[[vk::binding(0, 0)]] RWTexture2D<float4> img_illum_out;
[[vk::binding(1, 0)]] RWTexture2D<float2> img_moments_out;
[[vk::binding(2, 0)]] RWTexture2D<float> img_hist_len_out;

// ---- Inputs (previous frame's history) ----
[[vk::binding(3, 0)]] Texture2D<float4> img_illum_in;
[[vk::binding(4, 0)]] Texture2D<float2> img_moments_in;
[[vk::binding(5, 0)]] Texture2D<float> img_hist_len_in;

// ---- This frame's freshly traced data ----
[[vk::binding(6, 0)]] Texture2D<float4> img_illum_cur;
[[vk::binding(7, 0)]] Texture2D<float2> img_moments_cur;

// ---- G-Buffer ----
[[vk::binding(8, 0)]] Texture2D<float4> img_pos_id;
[[vk::binding(9, 0)]] Texture2D<float4> img_normal_depth;

// Previous frame's camera: world_to_camera then camera_to_sample (i.e. the
// matrix pair that maps a world position to previous-frame pixel coords).
struct PrevCameraGPU {
    float4 w2c0, w2c1, w2c2, w2c3;
    float4 c2s0, c2s1, c2s2, c2s3;
    uint width, height;
};
[[vk::binding(10, 0)]] ConstantBuffer<PrevCameraGPU> prev_camera;

struct PushConstants {
    uint width;
    uint height;
    uint reset;          // 1 => discard all history (camera moved or first frame)
    float alpha_color;   // history blend floor for colour (0.05 typical)
    float alpha_moments; // history blend floor for the moments (usually higher)
    uint clamp_history;  // 1 => rectify history into the current 3x3 colour AABB
    float phi_depth;     // depth-rejection tolerance in units of the depth gradient
    float phi_normal;    // normal-rejection exponent
    float max_hist_len;  // cap on the accumulated history length (0 = uncapped)
    float clamp_growth;  // how fast the rectification radius widens with history length
    uint debug_mode;     // 1 => write a status code instead of illumination
};
[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

// Magnitude of the screen-space depth gradient, by central differences.
//
// Computed here rather than in rt_path_trace.slang because ddx()/ddy() are
// fragment-shader derivatives: they need quad neighbourhoods and Metal refuses
// them in a compute kernel. Finite differences over already-written neighbours
// are equivalent, and by the time this pass runs the whole G-Buffer for the
// frame exists, so the neighbours are valid.
float depth_gradient(int2 ip, float center_z)
{
    const int2 lo = max(ip - int2(1, 1), int2(0, 0));
    const int2 hi = min(ip + int2(1, 1), int2(pc.width, pc.height) - int2(1, 1));

    float zl = img_normal_depth[int2(lo.x, ip.y)].z;
    float zr = img_normal_depth[int2(hi.x, ip.y)].z;
    float zu = img_normal_depth[int2(ip.x, lo.y)].z;
    float zd = img_normal_depth[int2(ip.x, hi.y)].z;

    // Clamp the one-sided span at the image border back to 1 so edge pixels are
    // not treated as having a (spuriously small) half-width gradient.
    float sx = max(float(hi.x - lo.x) * 0.5, 1.0);
    float sy = max(float(hi.y - lo.y) * 0.5, 1.0);

    float2 g = float2((zr - zl) / sx, (zd - zu) / sy) * 0.5;
    return length(g);
}

// Projects a world position into the PREVIOUS frame's pixel coordinates.
// Returns false if the point falls behind the previous camera.
bool project_prev(float3 world_pos, out float2 prev_pixel)
{
    // Slang warns (correctly) that the early-return paths below leave `prev_pixel`
    // undefined, and the caller reads it on the rejection path to report how far
    // reprojection moved. Always define it, so a failed projection reports a
    // zero-length offset instead of a garbage one.
    prev_pixel = float2(0.0, 0.0);

    float3 cam_pos = transform_point(prev_camera.w2c0, prev_camera.w2c1, prev_camera.w2c2, prev_camera.w2c3,
                                     world_pos);
    // This renderer's cameras look down +Z (see camera_controller.h), so a
    // point is in front when its camera-space z is positive.
    if (cam_pos.z <= 0.0) return false;

    float4 clip = float4(dot(prev_camera.c2s0, float4(cam_pos, 1.0)),
                         dot(prev_camera.c2s1, float4(cam_pos, 1.0)),
                         dot(prev_camera.c2s2, float4(cam_pos, 1.0)),
                         dot(prev_camera.c2s3, float4(cam_pos, 1.0)));
    if (abs(clip.w) < 1e-6) return false;

    float3 ndc = clip.xyz / clip.w;
    prev_pixel = float2(ndc.x * float(prev_camera.width), ndc.y * float(prev_camera.height));
    return true;
}

[shader("compute")]
[numthreads(8, 8, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= pc.width || tid.y >= pc.height) return;

    const int2 ip = int2(tid.xy);
    const float2 fpos = float2(tid.xy) + 0.5;

    const float4  illum_cur   = img_illum_cur[ip];
    const float2  moments_cur = img_moments_cur[ip];
    const float4  pos_id      = img_pos_id[ip];
    const float4  nd          = img_normal_depth[ip];

    const bool  has_surface = pos_id.w >= 0.0; // material id is INVALID_INDEX for misses
    const float linear_z    = nd.z;
    const float depth_grad  = depth_gradient(ip, linear_z);
    const float3 normal     = oct_decode(nd.xy);

    // Sky / escaped rays have no surface to reproject. They still get their
    // moments accumulated (the environment is static, so it converges nicely),
    // but no reprojection or validity test.
    float2 prev_pixel;
    bool reprojected = has_surface && project_prev(pos_id.xyz, prev_pixel);

    float4 hist_illum   = float4(0, 0, 0, 0);
    float2 hist_moments = float2(0, 0);
    float  hist_len     = 0.0;
    bool   history_valid = false;

    bool need_fallback = false;

    if (pc.reset != 0u) {
        // Fresh start (first frame, or the camera just moved): there is no
        // history to blend with, so this frame's sample becomes the history
        // verbatim. Deliberately NOT the neighbourhood mean used for
        // disocclusions below - on a global reset every pixel is "newly
        // revealed", and averaging them would blur the whole image.
        hist_illum    = illum_cur;
        hist_moments  = moments_cur;
        hist_len      = 0.0;
        history_valid = false;
    } else if (reprojected) {
        // ---- Bilinear fetch of the previous frame's history ----
        float2 pp     = prev_pixel - 0.5;
        int2  bpi     = int2(floor(pp));
        float2 f      = pp - float2(bpi);

        int2 p00 = bpi, p10 = bpi + int2(1, 0), p01 = bpi + int2(0, 1), p11 = bpi + int2(1, 1);

        float w_sum = 0.0;
        float4 ill_sum = float4(0, 0, 0, 0);
        float2 mom_sum = float2(0, 0);
        float  len_sum = 0.0;

        // Each of the four taps is validated independently so a single invalid
        // corner (a disocclusion, or geometry that just came into view) does not
        // drag the whole sample away.
        for (int k = 0; k < 4; ++k) {
            int2  tap   = (k == 0) ? p00 : (k == 1) ? p10 : (k == 2) ? p01 : p11;
            float w     = (k == 0) ? (1.0 - f.x) * (1.0 - f.y) :
                          (k == 1) ? f.x * (1.0 - f.y) :
                          (k == 2) ? (1.0 - f.x) * f.y : f.x * f.y;
            if (w <= 0.0) continue;

            if (any(tap < int2(0, 0)) || any(tap >= int2(pc.width, pc.height))) continue;

            const float4 tap_nd = img_normal_depth[tap];
            const float4 tap_pos = img_pos_id[tap];

            // Material id equality: rejects history from a different surface
            // that happens to land nearby (e.g. a thin object in front).
            if (abs(tap_pos.w - pos_id.w) > 0.5) continue;

            // Plane-aware depth test. The tolerance scales with the surface's
            // depth gradient, so sloped/grazing surfaces - whose depth changes
            // a lot across a pixel - are not rejected just for moving.
            float z_tol = pc.phi_depth * depth_grad * length(prev_pixel - fpos) + 1e-2;
            if (abs(tap_nd.z - linear_z) > z_tol) continue;

            // Normal agreement, sharpened by a large exponent.
            float3 tap_n = oct_decode(tap_nd.xy);
            float ndot   = dot(tap_n, normal);
            if (ndot < cos(pc.phi_normal)) continue;

            ill_sum += img_illum_in[tap] * w;
            mom_sum += img_moments_in[tap] * w;
            len_sum += img_hist_len_in[tap] * w;
            w_sum   += w;
        }

        if (w_sum > 1e-4) {
            hist_illum    = ill_sum / w_sum;
            hist_moments  = mom_sum / w_sum;
            hist_len      = len_sum / w_sum;
            history_valid = true;
        }
        if (!history_valid) need_fallback = true;
    } else {
        // Reprojection impossible (the surface is behind the previous camera).
        need_fallback = true;
    }

    // A pixel whose reprojection landed off-screen or failed every validity tap
    // has no usable history yet. Fall back to the CURRENT frame's neighbourhood:
    // on a disocclusion the surrounding pixels are freshly lit by the same
    // geometry, so their mean is a far better starting point than a single noisy
    // sample - this is what removes the "garbage for a few frames" pop that
    // newly revealed surfaces otherwise show.
    if (need_fallback) {
        float4 m1 = float4(0, 0, 0, 0);
        float2 m2 = float2(0, 0);
        float  wsum = 0.0;
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int2 q = ip + int2(dx, dy);
                if (any(q < int2(0, 0)) || any(q >= int2(pc.width, pc.height))) continue;
                float4 s = img_illum_cur[q];
                m1 += s;
                m2 += img_moments_cur[q];
                wsum += 1.0;
            }
        }
        if (wsum > 0.0) {
            hist_illum   = m1 / wsum;
            hist_moments = m2 / wsum;
        } else {
            hist_illum   = illum_cur;
            hist_moments = moments_cur;
        }
        hist_len = 0.0;
    }

    // ---- History rectification (neighbourhood clamping) ----
    // Even a history sample that passed every geometric test can be wrong: the
    // validity checks compare geometry, not lighting, so a shadow that moved or
    // a light that changed leaves stale radiance behind. Clamping the history
    // into the current frame's 3x3 colour range removes that ghosting, and as a
    // bonus suppresses fireflies that slipped past the radiance clamp.
    if (pc.clamp_history != 0u && history_valid) {
        float3 m1 = float3(0, 0, 0);
        float3 m2 = float3(0, 0, 0);
        float wsum = 0.0;
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int2 q = ip + int2(dx, dy);
                if (any(q < int2(0, 0)) || any(q >= int2(pc.width, pc.height))) continue;
                float3 s = img_illum_cur[q].rgb;
                m1 += s;
                m2 += s * s;
                wsum += 1.0;
            }
        }
        if (wsum > 0.0) {
            m1 /= wsum;
            m2 /= wsum;
            float3 sigma = sqrt(max(m2 - m1 * m1, float3(0, 0, 0)));

            // The rectification radius WIDENS as history accumulates.
            //
            // With a fixed radius this is a trade-off that cannot be won: it is
            // essential early (a fresh sample is mostly noise, so clamping stops
            // a stale neighbour from smearing across a disocclusion), but on a
            // static view it then keeps dragging the converged estimate towards
            // each new noisy sample - biasing the result and halting convergence
            // at a few dozen samples instead of letting it keep refining.
            //
            // Measured on assets/teapot: with a fixed 1.5 sigma the error stops
            // improving after ~16 frames and starts rising again; letting the
            // radius grow with history length keeps it falling (7.7x reduction
            // by 64 frames, close to the ideal sqrt(64) = 8x).
            float radius = 1.5 + hist_len * pc.clamp_growth;
            float3 lo = m1 - radius * sigma;
            float3 hi = m1 + radius * sigma;
            hist_illum.rgb = clamp(hist_illum.rgb, lo, hi);
        }
    }

    // ---- Exponential moving average ----
    // alpha starts at 1 (all new sample) and decays as 1/(n+1), which is
    // exactly a running mean while the view is static, but floors at
    // alpha_color so that a slowly changing signal (e.g. the camera easing to a
    // stop) never freezes on stale data.
    float new_len = hist_len + 1.0;
    float alpha   = max(1.0 / new_len, pc.alpha_color);
    float alpha_m = max(1.0 / new_len, pc.alpha_moments);

    float4 out_illum   = float4(lerp(hist_illum.rgb, illum_cur.rgb, alpha), 1.0);
    float2 out_moments = lerp(hist_moments, moments_cur, alpha_m);

    // Note the moments are NOT rectified against the clamped colour: doing so
    // would under-report variance precisely in the regions (moving shadows,
    // disocclusions) where the filter most needs to know that history is
    // untrustworthy.
    if (pc.max_hist_len > 0.0)
        new_len = min(new_len, pc.max_hist_len);

    // Diagnostic mode: report WHY history was or was not reused, as a code in
    // the illumination output. Tuning this filter blind is hopeless - the
    // symptom of "no reprojected surface", "reprojection left the screen" and
    // "every validity tap rejected" is identical from the final image (no
    // accumulation), so each gets its own code:
    //   -1 no surface (primary ray escaped)   -2 reprojection impossible
    //   -3 all bilinear taps rejected         >=0 accumulated history length
    if (pc.debug_mode != 0u) {
        float code = -3.0;
        if (!has_surface) code = -1.0;
        else if (!reprojected) code = -2.0;
        else if (history_valid) code = min(hist_len, 1e6);
        else
            // Rejected: fold in how far reprojection moved, so the host can
            // tell "landed a whole pixel off" from "landed exactly on target
            // but the geometry test still said no".
            code = -10.0 - min(length(prev_pixel - fpos), 100.0);
        img_illum_out[ip]    = float4(code, code, code, 1.0);
        img_moments_out[ip]  = out_moments;
        img_hist_len_out[ip] = new_len;
        return;
    }

    img_illum_out[ip]    = out_illum;
    img_moments_out[ip]  = out_moments;
    img_hist_len_out[ip] = new_len;
}
