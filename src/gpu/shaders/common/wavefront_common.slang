// Stage 2: shared code for the wavefront path tracer.
//
// This is #include'd (not compiled standalone - it has no [shader] entry
// point, so it must live outside shaders/*.slang's non-recursive glob; see
// src/gpu/CMakeLists.txt) by each of the five wavefront_*.slang stage
// kernels. It contains everything that was previously inlined into
// path_trace.slang's single main() - GPU-side struct layouts, math/BSDF/
// texture/light helper functions, and BVH traversal - UNCHANGED, plus the
// additional bindings/push-constants needed to carry per-path state and
// work queues BETWEEN separate kernel dispatches instead of within one
// big loop. See path_trace.slang's own header comment for the CPU-parity
// notes (RNG choice, envmap sampling simplification, etc.) - they apply
// here identically, since the actual shading math is copied verbatim.
//
// ============================================================================
// Wavefront architecture overview
// ============================================================================
//
// Per-pixel path state (ray, throughput, rng, etc.) is persistent across
// bounces in device buffers indexed by PIXEL index (not by "queue slot"),
// so each kernel can look up/update a path's state without needing it
// passed through a queue payload. A pixel's path is represented in the
// CURRENT bounce's *work queue* only while it is still alive; queues store
// just the pixel index (a uint), making them cheap to build/compact.
//
// Two `active_indices_*` buffers are used as ping-pong queues: on bounce D,
// one holds the list of pixels still alive entering this bounce (the
// "read" queue) and the other is atomically appended to by wavefront_shade
// to build the list of pixels alive AFTER this bounce (the "write" queue).
// `pc.read_is_a` (flipped by the host every bounce) selects which physical
// buffer plays which role, avoiding a copy between bounces.
//
// Per-bounce sequence (see src/gpu/wavefront_main.cpp for the exact
// dispatch/barrier orchestration):
//   1. wavefront_extend   - BVH-traverse each active path's current ray,
//                           write the hit result (or miss) to a slot-indexed
//                           hit buffer.
//   2. wavefront_shade    - resolve material, accumulate emission (with
//                           MIS), evaluate+sample the BSDF, and for paths
//                           that survive: atomically append to the NEXT
//                           active queue AND (if applicable) atomically
//                           append a shadow ray to the NEE queue.
//   3. wavefront_prepare_shadow_indirect - tiny (1-thread) kernel: snapshot
//                           the shadow queue's size into an indirect
//                           dispatch-args buffer for step 4. Must run
//                           strictly between shade and shadow (see its own
//                           file for why the ordering here is load-bearing).
//   4. wavefront_shadow   - trace each queued shadow ray; if unoccluded,
//                           add its pre-computed NEE contribution into
//                           accum_buffer.
//   5. wavefront_advance  - tiny (1-thread) kernel: snapshot the next-active
//                           queue's size into the indirect dispatch-args
//                           buffer for bounce D+1's extend/shade, then reset
//                           the atomic counters for D+1's shade to build into.
//
// accum_buffer[pixel] is updated via plain (non-atomic) `+=` in both
// wavefront_shade (emission) and wavefront_shadow (NEE), which is safe
// WITHOUT atomics because within a single bounce, each pixel appears in
// at most one of {a terminating emission-hit in shade, a queued shadow ray}
// - never both - and every queue (by construction) contains each currently-
// active pixel at most once. See the conversation history / PR description
// for the full data-race argument.

// ============================================================================
// Data structures (byte layouts matching the host-side Gpu*Gpu packed structs
// in src/gpu/wavefront_main.cpp)
// ============================================================================

struct GPUCameraGPU {
    float4 c2w0, c2w1, c2w2, c2w3;
    float4 s2c0, s2c1, s2c2, s2c3;
    uint width, height;
    float near_clip, far_clip;
};

struct GPUBVHNodeGPU {
    float3 bmin;
    uint offset;
    float3 bmax;
    uint meta;
};

struct GPUTextureGPU {
    uint type; // 0=Constant 1=Checkerboard 2=Bitmap
    uint channels;
    uint width;
    uint height;
    uint pixel_offset;
    float color0_r, color0_g, color0_b;
    float color1_r, color1_g, color1_b;
    float scale_u, scale_v;
};

struct GPUMaterialGPU {
    uint type; // see MAT_* constants below
    uint flags;
    uint tex_reflectance;
    uint tex_specular_reflectance;
    uint tex_specular_transmittance;
    uint tex_eta;
    uint tex_k;
    uint tex_alpha;
    uint tex_opacity;
    uint tex_bump;
    uint nested_bsdf;
    uint front_bsdf;
    uint back_bsdf;
    float eta;
    float inv_eta;
    float alpha;
    float scale;
    uint has_reflection;
    uint has_transmission;
    uint nonlinear;
    float fdr_int;
    float fdr_ext;
    float inv_eta_2;
    float specular_sampling_weight;
    uint transmittance_lut_offset;
    uint transmittance_lut_count;
    float internal_reflectance;
};

struct GPULightGPU {
    uint type; // 0=Area 1=Envmap
    uint mesh_id;
    uint radiance_tex;
    uint cdf_offset;
    uint cdf_count;
    float inv_total_area;
    uint base_triangle; // global triangle id of this light's owning mesh's local triangle 0
    uint spatial_varying;
    float4 to_world0, to_world1, to_world2, to_world3;
    float4 to_world_inv0, to_world_inv1, to_world_inv2, to_world_inv3;
    float bsc_x, bsc_y, bsc_z;
    float bounding_sphere_radius;
};

static const uint MAT_DIFFUSE = 0;
static const uint MAT_CONDUCTOR = 1;
static const uint MAT_DIELECTRIC = 2;
static const uint MAT_THINDIELECTRIC = 3;
static const uint MAT_PLASTIC = 4;
static const uint MAT_ROUGHCONDUCTOR = 5;
static const uint MAT_ROUGHDIELECTRIC = 6;
static const uint MAT_ROUGHPLASTIC = 7;
static const uint MAT_MASK = 8;
static const uint MAT_BUMPMAP = 9;
static const uint MAT_TWOSIDED = 10;

static const uint LIGHT_AREA = 0;
static const uint LIGHT_ENVMAP = 1;

static const uint INVALID_INDEX = 0xFFFFFFFFu;
static const float M_PI_F = 3.14159265358979323846;
static const float M_INV_PI_F = 0.31830988618379067154;
static const float M_INV_TWOPI_F = 0.15915494309189533577;
static const float M_INV_FOUR_PI_F = 0.07957747154594766788;
static const float EPS = 1e-4;

// Threads per 1D workgroup for every wavefront kernel except raygen (which
// dispatches over pixels directly) and the two 1-thread housekeeping
// kernels. Must match the group size the host uses to compute indirect
// dispatch arguments (see wavefront_main.cpp's GROUP_SIZE).
static const uint GROUP_SIZE = 256u;

// ============================================================================
// Bindings
// ============================================================================
// Bindings 0-15 mirror path_trace.slang's scene data exactly (same order),
// so GPUScene upload code is shared between the megakernel and wavefront
// executables (see wavefront_main.cpp).

[[vk::binding(0, 0)]] RWStructuredBuffer<float4> accum_buffer; // running sum of radiance (.xyz) + sample count (.w)
[[vk::binding(1, 0)]] ConstantBuffer<GPUCameraGPU> camera;
[[vk::binding(2, 0)]] StructuredBuffer<float4> positions;
[[vk::binding(3, 0)]] StructuredBuffer<float4> normals;
[[vk::binding(4, 0)]] StructuredBuffer<float2> uvs;
[[vk::binding(5, 0)]] StructuredBuffer<uint> indices;
[[vk::binding(6, 0)]] StructuredBuffer<GPUBVHNodeGPU> bvh_nodes;
[[vk::binding(7, 0)]] StructuredBuffer<uint> bvh_primitives;
[[vk::binding(8, 0)]] StructuredBuffer<uint> triangle_material_id;
[[vk::binding(9, 0)]] StructuredBuffer<uint> triangle_light_id;
[[vk::binding(10, 0)]] StructuredBuffer<GPUMaterialGPU> materials;
[[vk::binding(11, 0)]] StructuredBuffer<GPUTextureGPU> textures;
[[vk::binding(12, 0)]] StructuredBuffer<float> texture_pixels;
[[vk::binding(13, 0)]] StructuredBuffer<float> material_lut;
[[vk::binding(14, 0)]] StructuredBuffer<GPULightGPU> lights;
[[vk::binding(15, 0)]] StructuredBuffer<float> light_triangle_cdf;

// ---- Per-pixel persistent path state (indexed by PIXEL index) ----
[[vk::binding(16, 0)]] RWStructuredBuffer<float4> path_ray_o;       // .xyz = origin,    .w = t_min
[[vk::binding(17, 0)]] RWStructuredBuffer<float4> path_ray_d;       // .xyz = direction, .w = t_max
[[vk::binding(18, 0)]] RWStructuredBuffer<float4> path_throughput;  // .xyz = throughput,.w = eta
[[vk::binding(19, 0)]] RWStructuredBuffer<uint> path_rng;           // Rng.state
[[vk::binding(20, 0)]] RWStructuredBuffer<float2> path_prev_bsdf;   // .x = pdf, .y = delta (0/1)
[[vk::binding(21, 0)]] RWStructuredBuffer<float4> path_prev_p;      // .xyz = previous vertex position

// ---- Hit results (indexed by QUEUE SLOT, aligned with the active-index
// list read by that bounce's extend/shade dispatch) ----
[[vk::binding(22, 0)]] RWStructuredBuffer<uint> hit_tri;    // INVALID_INDEX => miss
[[vk::binding(23, 0)]] RWStructuredBuffer<float4> hit_tuv;  // .x=t .y=u .z=v

// ---- Active-path queues (ping-pong; see file header) ----
[[vk::binding(24, 0)]] RWStructuredBuffer<uint> active_indices_a;
[[vk::binding(25, 0)]] RWStructuredBuffer<uint> active_indices_b;

// counters[0] = size of the CURRENT bounce's active queue (read by extend/shade)
// counters[1] = size of the NEXT bounce's active queue (atomically built by shade)
// counters[2] = size of this bounce's shadow-ray queue (atomically built by shade)
[[vk::binding(26, 0)]] RWStructuredBuffer<uint> counters;

// ---- Shadow-ray (NEE) queue, indexed by QUEUE SLOT ----
[[vk::binding(27, 0)]] RWStructuredBuffer<float4> shadow_o;      // .xyz = origin,    .w = t_max
[[vk::binding(28, 0)]] RWStructuredBuffer<float4> shadow_d;      // .xyz = direction
[[vk::binding(29, 0)]] RWStructuredBuffer<uint> shadow_pixel;    // owning pixel index
[[vk::binding(30, 0)]] RWStructuredBuffer<float4> shadow_contrib; // .xyz = radiance to add if unoccluded

// Indirect dispatch args: [0..2] = VkDispatchIndirectCommand for the
// active-queue-sized dispatches (extend/shade), [3..5] = for the
// shadow-queue-sized dispatch. Written by wavefront_advance /
// wavefront_prepare_shadow_indirect respectively.
[[vk::binding(31, 0)]] RWStructuredBuffer<uint> indirect_args;

struct PushConstants {
    uint spp_index;             // seeds the RNG in raygen
    uint max_depth;
    uint rr_depth;
    uint num_lights;
    uint environment_light_id;  // INVALID_INDEX if none
    uint depth;                 // current bounce index (0-based) for this extend/shade dispatch
    uint read_is_a;             // 1 => read active_indices_a/write _b this bounce, 0 => vice versa
    uint total_pixels;          // width*height, used by raygen's bounds check
};
[[vk::push_constant]] ConstantBuffer<PushConstants> pc;

// Ping-pong queue accessors - centralizing the `read_is_a` branch here means
// every stage kernel just calls these two functions instead of repeating the
// branch (and risking the two ever falling out of sync with each other).
uint read_active(uint slot) { return pc.read_is_a != 0u ? active_indices_a[slot] : active_indices_b[slot]; }
void write_next_active(uint slot, uint pixel)
{
    if (pc.read_is_a != 0u) active_indices_b[slot] = pixel; else active_indices_a[slot] = pixel;
}

// ============================================================================
// RNG - see path_trace.slang's file header note (1): a GPU-friendly 32-bit
// PCG hash, NOT the CPU's 64-bit pcg32. Statistically valid, different noise
// pattern; only the converged image is expected to match the CPU renderer.
// ============================================================================

struct Rng {
    uint state;
};

uint pcg_hash(uint v)
{
    uint state = v * 747796405u + 2891336453u;
    uint word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

Rng rng_init(uint pixel_index, uint spp_index, uint stream_salt)
{
    Rng r;
    r.state = pcg_hash(pixel_index * 9781u + spp_index * 6271u + stream_salt * 26699u + 1u);
    return r;
}

float rng_next(inout Rng r)
{
    r.state = r.state * 747796405u + 2891336453u;
    uint word = ((r.state >> ((r.state >> 28u) + 4u)) ^ r.state) * 277803737u;
    uint result = (word >> 22u) ^ word;
    return float(result >> 8) / float(1u << 24);
}

float2 rng_next2(inout Rng r)
{
    return float2(rng_next(r), rng_next(r));
}

// ============================================================================
// Small math helpers (mirrors core/common.h, core/warp.h/.cpp, core/fresnel.h)
// ============================================================================

float safe_sqrt_(float x) { return sqrt(max(x, 0.0)); }
float mulsign_(float value, float sign) { return sign >= 0.0 ? value : -value; }

float2 square_to_uniform_disk_concentric(float2 s)
{
    float x = 2.0 * s.x - 1.0;
    float y = 2.0 * s.y - 1.0;

    bool is_zero = (x == 0.0 && y == 0.0);
    bool quad13 = abs(x) < abs(y);

    float r = quad13 ? y : x;
    float rp = quad13 ? x : y;

    float phi = 0.25 * M_PI_F * rp / r;
    if (quad13) phi = 0.5 * M_PI_F - phi;
    if (is_zero) phi = 0.0;

    return float2(r * cos(phi), r * sin(phi));
}

float3 square_to_cosine_hemisphere(float2 s)
{
    float2 p = square_to_uniform_disk_concentric(s);
    float z = sqrt(max(1.0 - dot(p, p), 0.0));
    return float3(p.x, p.y, z);
}

float square_to_cosine_hemisphere_pdf(float3 v) { return M_PI_F * v.z; }

float3 square_to_uniform_sphere(float2 s)
{
    float z = 1.0 - 2.0 * s.y;
    float r = safe_sqrt_(1.0 - z * z);
    // Both sin and cos must be taken of the SAME angle (2*pi*s.x) - see the
    // matching comment in src/core/warp.cpp.
    float phi = 2.0 * M_PI_F * s.x;
    float sn = sin(phi);
    float cn = cos(phi);
    return float3(r * cn, r * sn, z);
}
float square_to_uniform_sphere_pdf() { return M_INV_FOUR_PI_F; }

float2 square_to_uniform_triangle(float2 s)
{
    float t = safe_sqrt_(1.0 - s.x);
    return float2(1.0 - t, t * s.y);
}

void fresnel_dielectric(float cos_theta_i, float eta, out float r, out float cos_theta_t, out float eta_it, out float eta_ti)
{
    bool outside = cos_theta_i >= 0.0;
    float rcp_eta = 1.0 / eta;
    eta_it = outside ? eta : rcp_eta;
    eta_ti = outside ? rcp_eta : eta;

    float cos_theta_t_sqr = 1.0 - eta_ti * eta_ti * (1.0 - cos_theta_i * cos_theta_i);
    float cos_theta_i_abs = abs(cos_theta_i);
    float cos_theta_t_abs = safe_sqrt_(cos_theta_t_sqr);

    bool index_matched = (eta == 1.0);
    bool special_case = index_matched || cos_theta_i_abs == 0.0;
    float r_sc = index_matched ? 0.0 : 1.0;

    float a_s = (cos_theta_i_abs - eta_it * cos_theta_t_abs) / (eta_it * cos_theta_t_abs + cos_theta_i_abs);
    float a_p = (cos_theta_t_abs - eta_it * cos_theta_i_abs) / (eta_it * cos_theta_i_abs + cos_theta_t_abs);

    float rr = 0.5 * (a_s * a_s + a_p * a_p);
    r = special_case ? r_sc : rr;
    cos_theta_t = -mulsign_(cos_theta_t_abs, cos_theta_i);
}

float fresnel_conductor(float cos_theta_i, float eta, float k)
{
    float cos_theta_i_2 = cos_theta_i * cos_theta_i;
    float sin_theta_i_2 = 1.0 - cos_theta_i_2;
    float sin_theta_i_4 = sin_theta_i_2 * sin_theta_i_2;

    float eta_r = eta, eta_i = k;
    float temp1 = eta_r * eta_r - eta_i * eta_i - sin_theta_i_2;
    float a_2_pb_2 = safe_sqrt_(temp1 * temp1 + 4.0 * eta_i * eta_i * eta_r * eta_r);
    float a = safe_sqrt_(0.5 * (a_2_pb_2 + temp1));

    float term1 = a_2_pb_2 + cos_theta_i_2;
    float term2 = 2.0 * cos_theta_i * a;
    float r_s = (term1 - term2) / (term1 + term2);

    float term3 = a_2_pb_2 * cos_theta_i_2 + sin_theta_i_4;
    float term4 = term2 * sin_theta_i_2;
    float r_p = r_s * (term3 - term4) / (term3 + term4);

    return 0.5 * (r_s + r_p);
}

float3 reflect_local(float3 wi) { return float3(-wi.x, -wi.y, wi.z); }
float3 reflect_about(float3 wi, float3 m) { return m * 2.0 * dot(wi, m) - wi; }
float3 refract_local(float3 wi, float cos_theta_t, float eta_ti) { return float3(-eta_ti * wi.x, -eta_ti * wi.y, cos_theta_t); }
float3 refract_about(float3 wi, float3 m, float cos_theta_t, float eta_ti) { return m * (dot(wi, m) * eta_ti + cos_theta_t) - wi * eta_ti; }

float mis_weight(float pdf_a, float pdf_b)
{
    pdf_a *= pdf_a;
    pdf_b *= pdf_b;
    float denom = pdf_a + pdf_b;
    return denom == 0.0 ? 0.0 : pdf_a / denom;
}

float lerp_gather(uint offset, uint size, float x)
{
    x *= float(size - 1);
    int index = min(int(x), int(size) - 2);
    float v0 = material_lut[offset + uint(index)];
    float v1 = material_lut[offset + uint(index) + 1u];
    return lerp(v0, v1, x - float(index));
}

// ============================================================================
// Microfacet distribution (mirrors core/microfacet.h's TMicrofacetDistribution)
// ============================================================================

float microfacet_eval(float3 m, float alpha)
{
    float alpha_uv = alpha * alpha;
    float cos_theta = m.z;
    float t1 = m.x / alpha * m.x / alpha;
    float t2 = m.y / alpha * m.y / alpha;
    float t3 = t1 + t2 + m.z * m.z;
    float result = 1.0 / (M_PI_F * alpha_uv * t3 * t3);
    return (result * cos_theta > 0.0) ? result : 0.0;
}

float microfacet_smith_g1(float3 v, float3 m, float alpha)
{
    float xy_alpha_2 = alpha * v.x * alpha * v.x + alpha * v.y * alpha * v.y;
    float tan_theta_alpha_2 = xy_alpha_2 / (v.z * v.z);
    float result = 2.0 / (1.0 + sqrt(1.0 + tan_theta_alpha_2));
    if (xy_alpha_2 == 0.0) result = 1.0;
    if (dot(v, m) * v.z <= 0.0) result = 0.0;
    return result;
}

float microfacet_g(float3 wi, float3 wo, float3 m, float alpha) { return microfacet_smith_g1(wi, m, alpha) * microfacet_smith_g1(wo, m, alpha); }

float microfacet_pdf(float3 wi, float3 m, float alpha)
{
    float result = microfacet_eval(m, alpha);
    result *= microfacet_smith_g1(wi, m, alpha) * abs(dot(wi, m)) / wi.z;
    return result;
}

float3 microfacet_sample(float3 wi, float2 s, float alpha, out float pdf_out)
{
    float3 wi_p = normalize(float3(alpha * wi.x, alpha * wi.y, wi.z));

    float sin_theta_2 = wi_p.x * wi_p.x + wi_p.y * wi_p.y;
    float inv_sin_theta = 1.0 / sqrt(max(sin_theta_2, 1e-12));
    float2 sc = float2(wi_p.x, wi_p.y) * inv_sin_theta;
    if (abs(sin_theta_2) <= 4.0 * 1.1920929e-7) sc = float2(1.0, 0.0);
    float cos_phi = sc.x;
    float sin_phi = sc.y;

    float cos_theta = wi_p.z;

    float2 p = square_to_uniform_disk_concentric(s);
    float sh = 0.5 * (1.0 + cos_theta);
    p.y = lerp(safe_sqrt_(1.0 - p.x * p.x), p.y, sh);
    float xx = p.x, yy = p.y;
    float zz = safe_sqrt_(1.0 - dot(p, p));
    float sin_theta_i = safe_sqrt_(1.0 - cos_theta * cos_theta);
    float norm = 1.0 / (sin_theta_i * yy + cos_theta * zz);
    float2 slope = float2(cos_theta * yy - sin_theta_i * zz, xx) * norm;

    float2 slope1 = float2((cos_phi * slope.x - sin_phi * slope.y) * alpha,
                          (sin_phi * slope.x + cos_phi * slope.y) * alpha);

    float3 m = normalize(float3(-slope1.x, -slope1.y, 1.0));
    pdf_out = microfacet_eval(m, alpha) * microfacet_smith_g1(wi, m, alpha) * abs(dot(wi, m)) / wi.z;
    return m;
}

// ============================================================================
// Texture evaluation (mirrors src/textures/{constant,checkerboard,bitmap}.cpp)
// ============================================================================

float3 tex_color0(GPUTextureGPU t) { return float3(t.color0_r, t.color0_g, t.color0_b); }
float3 tex_color1(GPUTextureGPU t) { return float3(t.color1_r, t.color1_g, t.color1_b); }

float3 bitmap_sample(GPUTextureGPU t, float2 uv)
{
    uv.x = clamp(uv.x, 0.0, 1.0);
    uv.y = clamp(uv.y, 0.0, 1.0);
    float x = uv.x * float(t.width - 1);
    float y = uv.y * float(t.height - 1);
    int x0 = int(floor(x));
    int y0 = int(floor(y));
    int x1 = clamp(x0 + 1, 0, int(t.width) - 1);
    int y1 = clamp(y0 + 1, 0, int(t.height) - 1);
    float fx = x - float(x0);
    float fy = y - float(y0);

    uint ch = t.channels;
    uint base = t.pixel_offset;
    uint w = t.width;

    float3 c00, c10, c01, c11;
    if (ch == 3u) {
        uint i00 = base + (uint(y0) * w + uint(x0)) * 3u;
        uint i10 = base + (uint(y0) * w + uint(x1)) * 3u;
        uint i01 = base + (uint(y1) * w + uint(x0)) * 3u;
        uint i11 = base + (uint(y1) * w + uint(x1)) * 3u;
        c00 = float3(texture_pixels[i00], texture_pixels[i00 + 1], texture_pixels[i00 + 2]);
        c10 = float3(texture_pixels[i10], texture_pixels[i10 + 1], texture_pixels[i10 + 2]);
        c01 = float3(texture_pixels[i01], texture_pixels[i01 + 1], texture_pixels[i01 + 2]);
        c11 = float3(texture_pixels[i11], texture_pixels[i11 + 1], texture_pixels[i11 + 2]);
    } else {
        uint i00 = base + (uint(y0) * w + uint(x0));
        uint i10 = base + (uint(y0) * w + uint(x1));
        uint i01 = base + (uint(y1) * w + uint(x0));
        uint i11 = base + (uint(y1) * w + uint(x1));
        c00 = float3(texture_pixels[i00], texture_pixels[i00], texture_pixels[i00]);
        c10 = float3(texture_pixels[i10], texture_pixels[i10], texture_pixels[i10]);
        c01 = float3(texture_pixels[i01], texture_pixels[i01], texture_pixels[i01]);
        c11 = float3(texture_pixels[i11], texture_pixels[i11], texture_pixels[i11]);
    }

    float3 c0 = c00 * (1.0 - fx) + c10 * fx;
    float3 c1 = c01 * (1.0 - fx) + c11 * fx;
    return c0 * (1.0 - fy) + c1 * fy;
}

float3 tex_eval(uint tex_id, float2 uv)
{
    if (tex_id == INVALID_INDEX) return float3(0.0, 0.0, 0.0);
    GPUTextureGPU t = textures[tex_id];
    if (t.type == 0u) {
        return tex_color0(t);
    } else if (t.type == 1u) {
        float2 uvs2 = float2(uv.x * t.scale_u, uv.y * t.scale_v);
        float2 fr = uvs2 - floor(uvs2);
        bool m0 = (fr.x > 0.5) == (fr.y > 0.5);
        return m0 ? tex_color0(t) : tex_color1(t);
    } else {
        return bitmap_sample(t, uv);
    }
}

float tex_eval_1(uint tex_id, float2 uv)
{
    float3 c = tex_eval(tex_id, uv);
    return c.x;
}

// ============================================================================
// Surface hit (mirrors SurfaceIntersection3f fields used by BSDFs)
// ============================================================================

struct SurfaceHit {
    float2 uv;
    float3 shading_s, shading_t, shading_n;
    float3 wi;
    float3 dp_du, dp_dv;
};

float3 frame_to_world(SurfaceHit h, float3 v) { return h.shading_s * v.x + h.shading_t * v.y + h.shading_n * v.z; }
float3 frame_to_local(SurfaceHit h, float3 v) { return float3(dot(v, h.shading_s), dot(v, h.shading_t), dot(v, h.shading_n)); }

SurfaceHit flip_hit(SurfaceHit h)
{
    SurfaceHit r = h;
    r.wi.z *= -1.0;
    return r;
}
float3 flip_vec(float3 v) { float3 r = v; r.z *= -1.0; return r; }

SurfaceHit bump_perturbed_hit(SurfaceHit h, uint bump_tex, float scale)
{
    float2 grad_uv = float2(0.0, 0.0);
    grad_uv *= scale;

    float3 dp_du = h.dp_du + h.shading_n * (grad_uv.x - dot(h.shading_n, h.dp_du));
    float3 dp_dv = h.dp_dv + h.shading_n * (grad_uv.y - dot(h.shading_n, h.dp_dv));

    SurfaceHit r = h;
    float3 n = normalize(cross(dp_du, dp_dv));
    float3 s = normalize(dp_du - n * dot(dp_du, n));
    float3 t = cross(n, s);
    r.shading_n = n;
    r.shading_s = s;
    r.shading_t = t;
    r.wi = frame_to_local(r, frame_to_world(h, h.wi));
    return r;
}

// ============================================================================
// BSDF sample result
// ============================================================================

struct BSDFSample {
    float3 wo;
    float pdf;
    float eta;
    bool delta;
};

BSDFSample make_bsdf_sample(float3 wo) { BSDFSample b; b.wo = wo; b.pdf = 0.0; b.eta = 1.0; b.delta = false; return b; }

// ============================================================================
// Per-material eval/sample/pdf - mirrors src/bsdfs/*.cpp exactly.
// ============================================================================

float3 eval_diffuse(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (!(cti > 0.0 && cto > 0.0)) return float3(0, 0, 0);
    return tex_eval(m.tex_reflectance, h.uv) * M_INV_PI_F * cto;
}
float pdf_diffuse(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    float p = square_to_cosine_hemisphere_pdf(wo);
    return (cti > 0.0 && cto > 0.0) ? p : 0.0;
}
BSDFSample sample_diffuse(GPUMaterialGPU m, SurfaceHit h, float2 s2, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    if (h.wi.z <= 0.0) return bs;
    bs.wo = square_to_cosine_hemisphere(s2);
    bs.pdf = square_to_cosine_hemisphere_pdf(bs.wo);
    bs.eta = 1.0;
    bs.delta = false;
    if (bs.pdf > 0.0) value = tex_eval(m.tex_reflectance, h.uv);
    return bs;
}

BSDFSample sample_conductor(GPUMaterialGPU m, SurfaceHit h, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    float cti = h.wi.z;
    if (cti <= 0.0) return bs;
    bs.wo = reflect_local(h.wi);
    bs.eta = 1.0;
    bs.pdf = 1.0;
    bs.delta = true;
    float3 eta = tex_eval(m.tex_eta, h.uv);
    float3 k = tex_eval(m.tex_k, h.uv);
    float3 f = float3(fresnel_conductor(cti, eta.x, k.x), fresnel_conductor(cti, eta.y, k.y), fresnel_conductor(cti, eta.z, k.z));
    if (m.tex_specular_reflectance != INVALID_INDEX) f *= tex_eval(m.tex_specular_reflectance, h.uv);
    value = f;
    return bs;
}

BSDFSample sample_dielectric(GPUMaterialGPU m, SurfaceHit h, float sample1, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    float cti = h.wi.z;

    float r_i, cos_theta_t, eta_it, eta_ti;
    fresnel_dielectric(cti, m.eta, r_i, cos_theta_t, eta_it, eta_ti);
    float t_i = 1.0 - r_i;

    bool selected_r; float weight;
    bool has_r = m.has_reflection != 0u, has_t = m.has_transmission != 0u;
    if (has_r && has_t) {
        selected_r = sample1 <= r_i;
        bs.pdf = selected_r ? r_i : t_i;
        weight = 1.0;
    } else if (has_r || has_t) {
        selected_r = has_r;
        bs.pdf = 1.0;
        weight = has_r ? r_i : t_i;
    } else {
        return bs;
    }
    bool selected_t = !selected_r;

    bs.wo = selected_r ? reflect_local(h.wi) : refract_local(h.wi, cos_theta_t, eta_ti);
    bs.eta = selected_r ? 1.0 : eta_it;
    bs.delta = true;

    float3 result = float3(weight, weight, weight);
    if (selected_r) {
        if (m.tex_specular_reflectance != INVALID_INDEX) result *= tex_eval(m.tex_specular_reflectance, h.uv);
    }
    if (selected_t) {
        float factor = eta_ti;
        result *= factor * factor;
        if (m.tex_specular_transmittance != INVALID_INDEX) result *= tex_eval(m.tex_specular_transmittance, h.uv);
    }
    value = result;
    return bs;
}

BSDFSample sample_thindielectric(GPUMaterialGPU m, SurfaceHit h, float sample1, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);

    float rr, cos_theta_t, eta_it, eta_ti;
    fresnel_dielectric(abs(h.wi.z), m.eta, rr, cos_theta_t, eta_it, eta_ti);
    rr *= 2.0 / (1.0 + rr);
    float t = 1.0 - rr;

    bool selected_r; float weight;
    bool has_r = m.has_reflection != 0u, has_t = m.has_transmission != 0u;
    if (has_r && has_t) {
        selected_r = sample1 <= rr;
        bs.pdf = selected_r ? rr : t;
        weight = 1.0;
    } else if (has_r || has_t) {
        selected_r = has_r;
        bs.pdf = 1.0;
        weight = has_r ? rr : t;
    } else {
        return bs;
    }

    float3 result = float3(weight, weight, weight);
    if (selected_r) {
        if (m.tex_specular_reflectance != INVALID_INDEX) result *= tex_eval(m.tex_specular_reflectance, h.uv);
    } else {
        if (m.tex_specular_transmittance != INVALID_INDEX) result *= tex_eval(m.tex_specular_transmittance, h.uv);
    }

    bs.wo = selected_r ? reflect_local(h.wi) : -h.wi;
    bs.eta = 1.0;
    bs.delta = true;
    value = result;
    return bs;
}

float3 eval_plastic(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (!(cti > 0.0 && cto > 0.0)) return float3(0, 0, 0);
    float f_i, ct, ei, et; fresnel_dielectric(cti, m.eta, f_i, ct, ei, et);
    float f_o, ct2, ei2, et2; fresnel_dielectric(cto, m.eta, f_o, ct2, ei2, et2);

    float3 diffuse = tex_eval(m.tex_reflectance, h.uv);
    float3 denom = m.nonlinear != 0u ? (float3(1, 1, 1) - diffuse * m.fdr_int) : float3(1.0 - m.fdr_int, 1.0 - m.fdr_int, 1.0 - m.fdr_int);
    diffuse /= denom;
    diffuse *= square_to_cosine_hemisphere_pdf(wo) * m.inv_eta_2 * (1.0 - f_i) * (1.0 - f_o);
    return diffuse;
}
float pdf_plastic(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (!(cti > 0.0 && cto > 0.0)) return 0.0;
    float f_i, ct, ei, et; fresnel_dielectric(cti, m.eta, f_i, ct, ei, et);
    float prob_spec = f_i * m.specular_sampling_weight;
    float prob_diff = (1.0 - f_i) * (1.0 - m.specular_sampling_weight);
    prob_diff = prob_diff / (prob_diff + prob_spec);
    return square_to_cosine_hemisphere_pdf(wo) * prob_diff;
}
BSDFSample sample_plastic(GPUMaterialGPU m, SurfaceHit h, float sample1, float2 sample2, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    float cti = h.wi.z;
    if (cti <= 0.0) return bs;

    float f_i, ct, ei, et; fresnel_dielectric(cti, m.eta, f_i, ct, ei, et);
    float prob_spec = f_i * m.specular_sampling_weight;
    float prob_diff = (1.0 - f_i) * (1.0 - m.specular_sampling_weight);
    prob_spec = prob_spec / (prob_spec + prob_diff);
    prob_diff = 1.0 - prob_spec;

    bool sample_spec = sample1 < prob_spec;
    bs.eta = 1.0; bs.pdf = 0.0; bs.delta = true;

    if (sample_spec) {
        bs.wo = reflect_local(h.wi);
        bs.pdf = prob_spec;
        float3 v = float3(f_i / bs.pdf, f_i / bs.pdf, f_i / bs.pdf);
        if (m.tex_specular_reflectance != INVALID_INDEX) v *= tex_eval(m.tex_specular_reflectance, h.uv);
        value = v;
    } else {
        bs.wo = square_to_cosine_hemisphere(sample2);
        bs.pdf = prob_diff * square_to_cosine_hemisphere_pdf(bs.wo);
        // Smooth lobe - must clear the `delta` flag set above, otherwise the
        // next bounce skips the emitter pdf and double-counts this vertex's
        // diffuse contribution (mirrors src/bsdfs/plastic.cpp).
        bs.delta = false;
        float f_o, ct2, ei2, et2; fresnel_dielectric(bs.wo.z, m.eta, f_o, ct2, ei2, et2);
        float3 v = tex_eval(m.tex_reflectance, h.uv);
        float3 denom = m.nonlinear != 0u ? (float3(1, 1, 1) - v * m.fdr_int) : float3(1.0 - m.fdr_int, 1.0 - m.fdr_int, 1.0 - m.fdr_int);
        v /= denom;
        v *= m.inv_eta_2 * (1.0 - f_i) * (1.0 - f_o) / prob_diff;
        value = v;
    }
    return bs;
}

float3 eval_roughconductor(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (!(cti > 0.0 && cto > 0.0)) return float3(0, 0, 0);
    float3 hv = normalize(wo + h.wi);
    float alpha = tex_eval(m.tex_alpha, h.uv).x;
    float d = microfacet_eval(hv, alpha);
    if (d == 0.0) return float3(0, 0, 0);
    float g = microfacet_g(h.wi, wo, hv, alpha);
    float temp = d * g / (4.0 * cti);
    float3 eta = tex_eval(m.tex_eta, h.uv);
    float3 k = tex_eval(m.tex_k, h.uv);
    float cosih = dot(h.wi, hv);
    float3 f = float3(fresnel_conductor(cosih, eta.x, k.x), fresnel_conductor(cosih, eta.y, k.y), fresnel_conductor(cosih, eta.z, k.z));
    float3 result = f * temp;
    if (m.tex_specular_reflectance != INVALID_INDEX) result *= tex_eval(m.tex_specular_reflectance, h.uv);
    return result;
}
float pdf_roughconductor(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    float3 mv = normalize(wo + h.wi);
    if (!(cti > 0.0 && cto > 0.0 && dot(h.wi, mv) > 0.0 && dot(wo, mv) > 0.0)) return 0.0;
    float alpha = tex_eval(m.tex_alpha, h.uv).x;
    return microfacet_pdf(h.wi, mv, alpha) / (4.0 * dot(wo, mv));
}
BSDFSample sample_roughconductor(GPUMaterialGPU m, SurfaceHit h, float2 sample2, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    if (h.wi.z <= 0.0) return bs;

    float alpha = tex_eval(m.tex_alpha, h.uv).x;
    float pdf_m;
    float3 mv = microfacet_sample(h.wi, sample2, alpha, pdf_m);
    bs.pdf = pdf_m;
    bs.wo = reflect_about(h.wi, mv);
    bs.eta = 1.0;
    bs.delta = false;

    if (!(bs.pdf != 0.0 && bs.wo.z > 0.0)) return bs;

    float weight = microfacet_smith_g1(bs.wo, mv, alpha);
    bs.pdf /= 4.0 * dot(bs.wo, mv);

    float3 eta = tex_eval(m.tex_eta, h.uv);
    float3 k = tex_eval(m.tex_k, h.uv);
    float cosih = dot(h.wi, mv);
    float3 f = float3(fresnel_conductor(cosih, eta.x, k.x), fresnel_conductor(cosih, eta.y, k.y), fresnel_conductor(cosih, eta.z, k.z));
    float3 result = f * weight;
    if (m.tex_specular_reflectance != INVALID_INDEX) result *= tex_eval(m.tex_specular_reflectance, h.uv);
    value = result;
    return bs;
}

float3 eval_roughdielectric(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (cti == 0.0) return float3(0, 0, 0);
    bool refl = cti * cto > 0.0;
    float eta = cti > 0.0 ? m.eta : m.inv_eta;
    float inv_eta = cti > 0.0 ? m.inv_eta : m.eta;

    float3 mv = normalize(h.wi + wo * (refl ? 1.0 : eta));
    mv = mulsign_(1.0, mv.z) * mv;

    float alpha = tex_eval(m.tex_alpha, h.uv).x;
    float d = microfacet_eval(mv, alpha);
    float fr, ct, ei, et; fresnel_dielectric(dot(h.wi, mv), m.eta, fr, ct, ei, et);
    float g = microfacet_g(h.wi, wo, mv, alpha);

    float3 result = float3(0, 0, 0);
    bool eval_r = (m.has_reflection != 0u) && refl;
    bool eval_t = (m.has_transmission != 0u) && !refl;

    if (eval_r) {
        float3 v = float3(fr * d * g / (4.0 * abs(cti)), fr * d * g / (4.0 * abs(cti)), fr * d * g / (4.0 * abs(cti)));
        if (m.tex_specular_reflectance != INVALID_INDEX) v *= tex_eval(m.tex_specular_reflectance, h.uv);
        result = v;
    }
    if (eval_t) {
        float scale = inv_eta * inv_eta;
        float temp = dot(h.wi, mv) + eta * dot(wo, mv);
        float val = abs(scale * (1.0 - fr) * d * g * eta * eta * dot(h.wi, mv) * dot(wo, mv) / (cti * temp * temp));
        float3 v = float3(val, val, val);
        if (m.tex_specular_transmittance != INVALID_INDEX) v *= tex_eval(m.tex_specular_transmittance, h.uv);
        result = v;
    }
    return result;
}
float pdf_roughdielectric(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (cti == 0.0) return 0.0;
    bool refl = cti * cto > 0.0;
    if (!((m.has_reflection != 0u && refl) || (m.has_transmission != 0u && !refl))) return 0.0;

    float eta = cti > 0.0 ? m.eta : m.inv_eta;
    float3 mv = normalize(h.wi + wo * (refl ? 1.0 : eta));
    mv = mulsign_(1.0, mv.z) * mv;

    if (!(dot(h.wi, mv) * cti > 0.0 && dot(wo, mv) * cto > 0.0)) return 0.0;

    float temp = dot(h.wi, mv) + eta * dot(wo, mv);
    float dwh_dwo = refl ? 1.0 / (4.0 * dot(wo, mv)) : eta * eta * dot(wo, mv) / (temp * temp);

    float alpha = tex_eval(m.tex_alpha, h.uv).x;
    float alpha_scaled = alpha * (1.2 - 0.2 * sqrt(abs(cti)));
    float prob = microfacet_pdf(mulsign_(1.0, cti) * h.wi, mv, alpha_scaled);

    if (m.has_transmission != 0u && m.has_reflection != 0u) {
        float fr, ct, ei, et; fresnel_dielectric(dot(h.wi, mv), m.eta, fr, ct, ei, et);
        prob *= refl ? fr : (1.0 - fr);
    }
    return prob * abs(dwh_dwo);
}
BSDFSample sample_roughdielectric(GPUMaterialGPU m, SurfaceHit h, float sample1, float2 sample2, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    float cti = h.wi.z;
    if (cti == 0.0) return bs;

    float alpha = tex_eval(m.tex_alpha, h.uv).x;
    float alpha_scaled = alpha * (1.2 - 0.2 * sqrt(abs(cti)));

    float pdf_m;
    float3 mv = microfacet_sample(mulsign_(1.0, cti) * h.wi, sample2, alpha_scaled, pdf_m);
    bs.pdf = pdf_m;
    if (bs.pdf == 0.0) return bs;

    float fr, cos_theta_t, eta_it, eta_ti; fresnel_dielectric(dot(h.wi, mv), m.eta, fr, cos_theta_t, eta_it, eta_ti);

    bool selected_r; float weight;
    bool has_r = m.has_reflection != 0u, has_t = m.has_transmission != 0u;
    if (has_r && has_t) {
        selected_r = sample1 <= fr;
        weight = 1.0;
        bs.pdf *= selected_r ? fr : (1.0 - fr);
    } else if (has_r || has_t) {
        selected_r = has_r;
        weight = has_r ? fr : (1.0 - fr);
    } else {
        return bs;
    }
    bool selected_t = !selected_r;

    bs.eta = selected_r ? 1.0 : eta_it;
    float dwh_dwo = 0.0;
    float3 result = float3(weight, weight, weight);

    if (selected_r) {
        bs.wo = reflect_about(h.wi, mv);
        dwh_dwo = 1.0 / (4.0 * dot(bs.wo, mv));
        if (m.tex_specular_reflectance != INVALID_INDEX) result *= tex_eval(m.tex_specular_reflectance, h.uv);
    }
    if (selected_t) {
        bs.wo = refract_about(h.wi, mv, cos_theta_t, eta_ti);
        float temp = dot(h.wi, mv) + bs.eta * dot(bs.wo, mv);
        dwh_dwo = bs.eta * bs.eta * dot(bs.wo, mv) / (temp * temp);
        float factor = eta_ti;
        result *= factor * factor;
        if (m.tex_specular_transmittance != INVALID_INDEX) result *= tex_eval(m.tex_specular_transmittance, h.uv);
    }

    result *= microfacet_smith_g1(bs.wo, mv, alpha);
    bs.pdf *= abs(dwh_dwo);
    bs.delta = false;
    value = result;
    return bs;
}

float3 eval_roughplastic(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (!(cti > 0.0 && cto > 0.0)) return float3(0, 0, 0);

    float3 hv = normalize(wo + h.wi);
    float d = microfacet_eval(hv, m.alpha);
    float fr, ct, ei, et; fresnel_dielectric(dot(h.wi, hv), m.eta, fr, ct, ei, et);
    float g = microfacet_g(h.wi, wo, hv, m.alpha);
    float3 value = float3(fr * d * g / (4.0 * cti), fr * d * g / (4.0 * cti), fr * d * g / (4.0 * cti));
    if (m.tex_specular_reflectance != INVALID_INDEX) value *= tex_eval(m.tex_specular_reflectance, h.uv);

    uint lut = m.transmittance_lut_offset, lutn = m.transmittance_lut_count;
    float t_i = lerp_gather(lut, lutn, cti);
    float t_o = lerp_gather(lut, lutn, cto);
    float3 diffuse = tex_eval(m.tex_reflectance, h.uv);
    float3 denom = m.nonlinear != 0u ? (float3(1, 1, 1) - diffuse * m.internal_reflectance) : float3(1.0 - m.internal_reflectance, 1.0 - m.internal_reflectance, 1.0 - m.internal_reflectance);
    diffuse /= denom;
    value += diffuse * (M_INV_PI_F * m.inv_eta_2 * cto * t_i * t_o);
    return value;
}
float pdf_roughplastic(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (!(cti > 0.0 && cto > 0.0)) return 0.0;

    uint lut = m.transmittance_lut_offset, lutn = m.transmittance_lut_count;
    float t_i = lerp_gather(lut, lutn, cti);
    float prob_spec = (1.0 - t_i) * m.specular_sampling_weight;
    float prob_diff = t_i * (1.0 - m.specular_sampling_weight);
    prob_spec = prob_spec / (prob_spec + prob_diff);
    prob_diff = 1.0 - prob_spec;

    float3 hv = normalize(wo + h.wi);
    float result = microfacet_eval(hv, m.alpha) * microfacet_smith_g1(h.wi, hv, m.alpha) / (4.0 * cti);
    result *= prob_spec;
    result += prob_diff * square_to_cosine_hemisphere_pdf(wo);
    return result;
}
BSDFSample sample_roughplastic(GPUMaterialGPU m, SurfaceHit h, float sample1, float2 sample2, out float3 value)
{
    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    value = float3(0, 0, 0);
    float cti = h.wi.z;
    if (cti <= 0.0) return bs;

    uint lut = m.transmittance_lut_offset, lutn = m.transmittance_lut_count;
    float t_i = lerp_gather(lut, lutn, cti);
    float prob_spec = (1.0 - t_i) * m.specular_sampling_weight;
    float prob_diff = t_i * (1.0 - m.specular_sampling_weight);
    prob_spec /= (prob_spec + prob_diff);

    bool sample_spec = sample1 < prob_spec;
    bs.eta = 1.0;

    if (sample_spec) {
        float pdf_m;
        float3 mv = microfacet_sample(h.wi, sample2, m.alpha, pdf_m);
        bs.wo = reflect_about(h.wi, mv);
    } else {
        bs.wo = square_to_cosine_hemisphere(sample2);
    }

    bs.pdf = pdf_roughplastic(m, h, bs.wo);
    bs.delta = false;
    if (bs.pdf <= 0.0) return bs;
    float3 result = eval_roughplastic(m, h, bs.wo);
    value = result / bs.pdf;
    return bs;
}

// ============================================================================
// Dispatcher: resolves compound materials (Mask/BumpMap/TwoSided) using
// bounded-depth loops (SPIR-V compute shaders do not support recursion).
// ============================================================================

static const int MAX_UNWRAP = 6;

float3 bsdf_eval(uint mat_id, SurfaceHit hit, float3 wo)
{
    float3 extra_weight = float3(1, 1, 1);
    SurfaceHit h = hit;
    float3 w = wo;

    for (int i = 0; i < MAX_UNWRAP; ++i) {
        if (mat_id == INVALID_INDEX) return float3(0, 0, 0);
        GPUMaterialGPU m = materials[mat_id];

        if (m.type == MAT_TWOSIDED) {
            if (h.wi.z > 0.0) {
                mat_id = m.front_bsdf;
            } else {
                h = flip_hit(h);
                w = flip_vec(w);
                mat_id = m.back_bsdf;
            }
            continue;
        } else if (m.type == MAT_MASK) {
            extra_weight *= tex_eval(m.tex_opacity, h.uv).x;
            mat_id = m.nested_bsdf;
            continue;
        } else if (m.type == MAT_BUMPMAP) {
            SurfaceHit ph = bump_perturbed_hit(h, m.tex_bump, m.scale);
            float3 pwo = frame_to_local(ph, frame_to_world(h, w));
            if (!(pwo.z * w.z > 0.0)) return float3(0, 0, 0);
            h = ph; w = pwo;
            mat_id = m.nested_bsdf;
            continue;
        }

        float3 result;
        if (m.type == MAT_DIFFUSE) result = eval_diffuse(m, h, w);
        else if (m.type == MAT_PLASTIC) result = eval_plastic(m, h, w);
        else if (m.type == MAT_ROUGHCONDUCTOR) result = eval_roughconductor(m, h, w);
        else if (m.type == MAT_ROUGHDIELECTRIC) result = eval_roughdielectric(m, h, w);
        else if (m.type == MAT_ROUGHPLASTIC) result = eval_roughplastic(m, h, w);
        else result = float3(0, 0, 0);

        return result * extra_weight;
    }
    return float3(0, 0, 0);
}

float bsdf_pdf(uint mat_id, SurfaceHit hit, float3 wo)
{
    float extra_weight = 1.0;
    SurfaceHit h = hit;
    float3 w = wo;

    for (int i = 0; i < MAX_UNWRAP; ++i) {
        if (mat_id == INVALID_INDEX) return 0.0;
        GPUMaterialGPU m = materials[mat_id];

        if (m.type == MAT_TWOSIDED) {
            if (h.wi.z > 0.0) {
                mat_id = m.front_bsdf;
            } else {
                h = flip_hit(h);
                w = flip_vec(w);
                mat_id = m.back_bsdf;
            }
            continue;
        } else if (m.type == MAT_MASK) {
            extra_weight *= tex_eval(m.tex_opacity, h.uv).x;
            mat_id = m.nested_bsdf;
            continue;
        } else if (m.type == MAT_BUMPMAP) {
            SurfaceHit ph = bump_perturbed_hit(h, m.tex_bump, m.scale);
            float3 pwo = frame_to_local(ph, frame_to_world(h, w));
            if (!(pwo.z * w.z > 0.0)) return 0.0;
            h = ph; w = pwo;
            mat_id = m.nested_bsdf;
            continue;
        }

        float result;
        if (m.type == MAT_DIFFUSE) result = pdf_diffuse(m, h, w);
        else if (m.type == MAT_PLASTIC) result = pdf_plastic(m, h, w);
        else if (m.type == MAT_ROUGHCONDUCTOR) result = pdf_roughconductor(m, h, w);
        else if (m.type == MAT_ROUGHDIELECTRIC) result = pdf_roughdielectric(m, h, w);
        else if (m.type == MAT_ROUGHPLASTIC) result = pdf_roughplastic(m, h, w);
        else result = 0.0;

        return result * extra_weight;
    }
    return 0.0;
}

BSDFSample bsdf_sample(uint mat_id, SurfaceHit hit, float sample1, float2 sample2, out float3 value)
{
    value = float3(0, 0, 0);
    float3 extra_weight = float3(1, 1, 1);
    float mask_pdf_scale = 1.0;
    SurfaceHit h = hit;
    int op_stack[MAX_UNWRAP];
    SurfaceHit frame_stack[MAX_UNWRAP + 1];
    int depth = 0;

    for (int i = 0; i < MAX_UNWRAP; ++i) {
        if (mat_id == INVALID_INDEX) { BSDFSample bs = make_bsdf_sample(float3(0,0,0)); return bs; }
        GPUMaterialGPU m = materials[mat_id];

        if (m.type == MAT_TWOSIDED) {
            if (h.wi.z > 0.0) {
                mat_id = m.front_bsdf;
            } else {
                op_stack[depth] = 0; frame_stack[depth] = h; depth++;
                h = flip_hit(h);
                mat_id = m.back_bsdf;
            }
            continue;
        } else if (m.type == MAT_MASK) {
            float opacity = tex_eval(m.tex_opacity, h.uv).x;
            bool enter = (sample1 < opacity) && (opacity > 0.0);

            if (!enter) {
                BSDFSample bs = make_bsdf_sample(-h.wi);
                bs.eta = 1.0;
                bs.pdf = (1.0 - opacity) * mask_pdf_scale;
                bs.delta = true;
                value = extra_weight;

                for (int k = depth - 1; k >= 0; --k) {
                    if (op_stack[k] == 0) bs.wo = flip_vec(bs.wo);
                    else bs.wo = frame_to_local(frame_stack[k], frame_to_world(frame_stack[k + 1], bs.wo));
                }
                return bs;
            }

            sample1 = sample1 / opacity;
            mask_pdf_scale *= opacity;
            mat_id = m.nested_bsdf;
            continue;
        } else if (m.type == MAT_BUMPMAP) {
            SurfaceHit ph = bump_perturbed_hit(h, m.tex_bump, m.scale);
            op_stack[depth] = 1; frame_stack[depth] = h; depth++;
            frame_stack[depth] = ph;
            h = ph;
            mat_id = m.nested_bsdf;
            continue;
        }

        BSDFSample bs;
        float3 leaf_value;
        if (m.type == MAT_DIFFUSE) bs = sample_diffuse(m, h, sample2, leaf_value);
        else if (m.type == MAT_CONDUCTOR) bs = sample_conductor(m, h, leaf_value);
        else if (m.type == MAT_DIELECTRIC) bs = sample_dielectric(m, h, sample1, leaf_value);
        else if (m.type == MAT_THINDIELECTRIC) bs = sample_thindielectric(m, h, sample1, leaf_value);
        else if (m.type == MAT_PLASTIC) bs = sample_plastic(m, h, sample1, sample2, leaf_value);
        else if (m.type == MAT_ROUGHCONDUCTOR) bs = sample_roughconductor(m, h, sample2, leaf_value);
        else if (m.type == MAT_ROUGHDIELECTRIC) bs = sample_roughdielectric(m, h, sample1, sample2, leaf_value);
        else if (m.type == MAT_ROUGHPLASTIC) bs = sample_roughplastic(m, h, sample1, sample2, leaf_value);
        else { bs = make_bsdf_sample(float3(0,0,0)); leaf_value = float3(0,0,0); }

        value = leaf_value * extra_weight;
        bs.pdf *= mask_pdf_scale;

        for (int k = depth - 1; k >= 0; --k) {
            if (op_stack[k] == 0) {
                bs.wo = flip_vec(bs.wo);
            } else {
                bs.wo = frame_to_local(frame_stack[k], frame_to_world(frame_stack[k + 1], bs.wo));
            }
        }
        return bs;
    }

    BSDFSample bs = make_bsdf_sample(float3(0, 0, 0));
    return bs;
}

bool bsdf_is_smooth(uint mat_id)
{
    const uint SMOOTH_MASK = 0x2u | 0x4u | 0x8u | 0x10u;
    if (mat_id == INVALID_INDEX) return false;
    return (materials[mat_id].flags & SMOOTH_MASK) != 0u;
}

// ============================================================================
// Light sampling (mirrors Scene::sample_emitter_direction / Mesh::
// sample_direction / Area::sample_direction / EnvironmentMap)
// ============================================================================

struct DirSample {
    float3 p, n, d;
    float2 uv;
    float t, pdf, dist;
    bool delta;
    uint light_id;
};

uint cdf_sample(uint offset, uint count, float u, out float rescaled_u)
{
    uint lo = 0, hi = count;
    while (lo < hi) {
        uint mid = (lo + hi) / 2u;
        if (light_triangle_cdf[offset + mid] < u) lo = mid + 1u; else hi = mid;
    }
    uint index = min(lo, count - 1u);
    float cdf_here = light_triangle_cdf[offset + index];
    float cdf_prev = index > 0u ? light_triangle_cdf[offset + index - 1u] : 0.0;
    float pmf = cdf_here - cdf_prev;
    rescaled_u = pmf > 0.0 ? (u - cdf_prev) / pmf : 0.0;
    return index;
}

DirSample sample_area_light(uint light_idx, float3 ref_p, float2 sample_)
{
    GPULightGPU L = lights[light_idx];
    float rescaled_v;
    uint local_tri = cdf_sample(L.cdf_offset, L.cdf_count, sample_.y, rescaled_v);
    uint global_tri = L.base_triangle + local_tri;

    uint i0 = indices[global_tri * 3 + 0];
    uint i1 = indices[global_tri * 3 + 1];
    uint i2 = indices[global_tri * 3 + 2];
    float3 p0 = positions[i0].xyz, p1 = positions[i1].xyz, p2 = positions[i2].xyz;
    float3 e0 = p1 - p0, e1 = p2 - p0;
    float2 b = square_to_uniform_triangle(float2(sample_.x, rescaled_v));

    DirSample ds;
    ds.p = e0 * b.x + e1 * b.y + p0;
    float3 geo_n = normalize(cross(e0, e1));
    float3 n0 = normals[i0].xyz, n1 = normals[i1].xyz, n2 = normals[i2].xyz;
    float3 bary = float3(1.0 - b.x - b.y, b.x, b.y);
    float3 shading_n = (dot(n0,n0) + dot(n1,n1) + dot(n2,n2)) > 0.0 ? normalize(n0 * bary.x + n1 * bary.y + n2 * bary.z) : geo_n;
    ds.n = shading_n;
    ds.uv = float2(0, 0);
    ds.t = 0.0;
    ds.pdf = L.inv_total_area;
    ds.delta = false;

    ds.d = ds.p - ref_p;
    float dist2 = dot(ds.d, ds.d);
    ds.dist = sqrt(dist2);
    ds.d /= ds.dist;

    float dp = abs(dot(ds.d, ds.n));
    ds.pdf *= dp != 0.0 ? dist2 / dp : 0.0;
    ds.light_id = light_idx;
    return ds;
}

float pdf_area_light(uint light_idx, float3 ds_d, float3 ds_n, float ds_dist)
{
    GPULightGPU L = lights[light_idx];
    float pdf = L.inv_total_area;
    float dp = abs(dot(ds_d, ds_n));
    pdf *= dp != 0.0 ? ds_dist * ds_dist / dp : 0.0;
    return pdf;
}

float3 eval_envmap(uint light_idx, float3 wi_world_neg)
{
    GPULightGPU L = lights[light_idx];
    if (L.spatial_varying == 0u) {
        return tex_eval(L.radiance_tex, float2(0, 0));
    }
    float4 r0 = L.to_world_inv0, r1 = L.to_world_inv1, r2 = L.to_world_inv2;
    float3 d = float3(dot(r0.xyz, wi_world_neg), dot(r1.xyz, wi_world_neg), dot(r2.xyz, wi_world_neg));
    float phi = atan2(d.x, -d.z);
    float theta = acos(clamp(d.y, -1.0, 1.0));
    float2 uv = float2(phi * M_INV_TWOPI_F, theta * M_INV_PI_F);
    uv.x -= 0.5 / float(max(int(textures[L.radiance_tex].width), 1) - 1 + 1e-6);
    uv -= floor(uv);
    return tex_eval(L.radiance_tex, uv);
}

DirSample sample_envmap(uint light_idx, float3 ref_p, float2 sample_)
{
    GPULightGPU L = lights[light_idx];
    float radius = max(L.bounding_sphere_radius, length(ref_p - float3(L.bsc_x, L.bsc_y, L.bsc_z)));
    float dist = 2.0 * radius;

    float3 d = square_to_uniform_sphere(sample_);
    DirSample ds;
    ds.p = ref_p + d * dist;
    ds.n = -d;
    ds.t = 0.0;
    ds.pdf = square_to_uniform_sphere_pdf();
    ds.delta = false;
    ds.d = d;
    ds.dist = dist;
    ds.light_id = light_idx;
    ds.uv = float2(0, 0);
    return ds;
}

float pdf_envmap() { return square_to_uniform_sphere_pdf(); }

DirSample sample_emitter_direction(float3 ref_p, inout Rng rng, out float3 weight)
{
    uint n = pc.num_lights;
    weight = float3(0, 0, 0);
    DirSample ds;
    ds.pdf = 0.0; ds.delta = false; ds.light_id = INVALID_INDEX;
    if (n == 0u) return ds;

    uint idx = min(uint(rng_next(rng) * float(n)), n - 1u);
    float2 s2 = rng_next2(rng);

    GPULightGPU L = lights[idx];
    if (L.type == LIGHT_AREA) {
        ds = sample_area_light(idx, ref_p, s2);
        if (ds.pdf > 0.0 && dot(ds.d, ds.n) < 0.0) {
            float3 rad = tex_eval(L.radiance_tex, float2(0, 0));
            weight = rad / ds.pdf;
        } else {
            ds.pdf = 0.0;
        }
    } else {
        ds = sample_envmap(idx, ref_p, s2);
        // eval_envmap() expects the direction POINTING AT the environment
        // (i.e. CPU-side `-si.wi`), which is exactly ds.d here - passing
        // -ds.d would look up the antipodal texel (phi+pi, pi-theta), i.e.
        // a horizontally rotated + vertically flipped sky.
        float3 rad = eval_envmap(idx, ds.d);
        weight = ds.pdf > 0.0 ? rad / ds.pdf : float3(0, 0, 0);
    }

    ds.pdf *= 1.0 / float(n);
    weight *= float(n);

    return ds;
}

float pdf_emitter_direction(uint light_id, float3 ds_d, float3 ds_n, float ds_dist)
{
    if (light_id == INVALID_INDEX || pc.num_lights == 0u) return 0.0;
    GPULightGPU L = lights[light_id];
    float pdf = (L.type == LIGHT_AREA) ? pdf_area_light(light_id, ds_d, ds_n, ds_dist) : pdf_envmap();
    return pdf / float(pc.num_lights);
}

// ============================================================================
// Camera ray generation (mirrors PerspectiveCamera::sample_ray)
// ============================================================================

float3 transform_point(float4 r0, float4 r1, float4 r2, float4 r3, float3 p)
{
    float4 hp = float4(p, 1.0);
    float w = dot(r3, hp);
    return float3(dot(r0, hp), dot(r1, hp), dot(r2, hp)) / w;
}
float3 transform_vector(float4 r0, float4 r1, float4 r2, float4 r3, float3 v)
{
    float4 hv = float4(v, 0.0);
    return float3(dot(r0, hv), dot(r1, hv), dot(r2, hv));
}

// ============================================================================
// BVH traversal (mirrors BVHAccel::ray_intersect / ray_test)
// ============================================================================

static const uint BVH_META_COUNT_MASK = 0x3FFFFFFFu;
static const uint BVH_META_AXIS_SHIFT = 30u;

bool ray_aabb(float3 ro, float3 rd_inv, float3 bmin, float3 bmax, float t_min, float t_max, out float near_t, out float far_t)
{
    float3 t0 = (bmin - ro) * rd_inv;
    float3 t1 = (bmax - ro) * rd_inv;
    float3 tsmall = min(t0, t1);
    float3 tbig = max(t0, t1);
    near_t = max(max(tsmall.x, tsmall.y), max(tsmall.z, t_min));
    far_t = min(min(tbig.x, tbig.y), min(tbig.z, t_max));
    return near_t <= far_t;
}

bool ray_triangle(float3 ro, float3 rd, float3 p0, float3 p1, float3 p2, out float t, out float u, out float v)
{
    float3 e1 = p1 - p0, e2 = p2 - p0;
    float3 pvec = cross(rd, e2);
    float det = dot(e1, pvec);
    t = 0.0; u = 0.0; v = 0.0;
    if (abs(det) < 1e-9) return false;
    float inv_det = 1.0 / det;
    float3 tvec = ro - p0;
    u = dot(tvec, pvec) * inv_det;
    if (u < 0.0 || u > 1.0) return false;
    float3 qvec = cross(tvec, e1);
    v = dot(rd, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0) return false;
    t = dot(e2, qvec) * inv_det;
    return t > 1e-6;
}

struct HitResult {
    bool hit;
    float t;
    uint tri;
    float u, v;
};

HitResult trace_ray(float3 ro, float3 rd, float t_min, float t_max)
{
    float3 rd_inv = float3(1.0, 1.0, 1.0) / rd;
    uint stack[64];
    int sp = 0;
    uint node_idx = 0;

    HitResult result;
    result.hit = false;
    result.t = t_max;

    for (;;) {
        GPUBVHNodeGPU node = bvh_nodes[node_idx];
        float near_t, far_t;
        bool box_hit = ray_aabb(ro, rd_inv, node.bmin, node.bmax, t_min, result.t, near_t, far_t);

        if (box_hit) {
            uint count = node.meta & BVH_META_COUNT_MASK;
            if (count > 0u) {
                uint base = node.offset;
                for (uint i = 0u; i < count; ++i) {
                    uint tri = bvh_primitives[base + i];
                    float3 p0 = positions[indices[tri * 3 + 0]].xyz;
                    float3 p1 = positions[indices[tri * 3 + 1]].xyz;
                    float3 p2 = positions[indices[tri * 3 + 2]].xyz;
                    float t, u, v;
                    if (ray_triangle(ro, rd, p0, p1, p2, t, u, v) && t < result.t && t >= t_min) {
                        result.t = t; result.tri = tri; result.u = u; result.v = v; result.hit = true;
                    }
                }
                if (sp == 0) break;
                node_idx = stack[--sp];
            } else {
                uint axis = node.meta >> BVH_META_AXIS_SHIFT;
                bool neg = rd[axis] < 0.0;
                uint near_child = neg ? (node_idx + 1u) : node.offset;
                uint far_child = neg ? node.offset : (node_idx + 1u);
                stack[sp++] = far_child;
                node_idx = near_child;
            }
        } else {
            if (sp == 0) break;
            node_idx = stack[--sp];
        }
    }
    return result;
}

bool trace_shadow_ray(float3 ro, float3 rd, float t_min, float t_max)
{
    float3 rd_inv = float3(1.0, 1.0, 1.0) / rd;
    uint stack[64];
    int sp = 0;
    uint node_idx = 0;

    for (;;) {
        GPUBVHNodeGPU node = bvh_nodes[node_idx];
        float near_t, far_t;
        bool box_hit = ray_aabb(ro, rd_inv, node.bmin, node.bmax, t_min, t_max, near_t, far_t);

        if (box_hit) {
            uint count = node.meta & BVH_META_COUNT_MASK;
            if (count > 0u) {
                uint base = node.offset;
                for (uint i = 0u; i < count; ++i) {
                    uint tri = bvh_primitives[base + i];
                    float3 p0 = positions[indices[tri * 3 + 0]].xyz;
                    float3 p1 = positions[indices[tri * 3 + 1]].xyz;
                    float3 p2 = positions[indices[tri * 3 + 2]].xyz;
                    float t, u, v;
                    if (ray_triangle(ro, rd, p0, p1, p2, t, u, v) && t >= t_min && t <= t_max) return true;
                }
                if (sp == 0) break;
                node_idx = stack[--sp];
            } else {
                uint axis = node.meta >> BVH_META_AXIS_SHIFT;
                bool neg = rd[axis] < 0.0;
                uint near_child = neg ? (node_idx + 1u) : node.offset;
                uint far_child = neg ? node.offset : (node_idx + 1u);
                stack[sp++] = far_child;
                node_idx = near_child;
            }
        } else {
            if (sp == 0) break;
            node_idx = stack[--sp];
        }
    }
    return false;
}

// offset_p() (intersection.h) - self-intersection avoidance
float3 offset_p(float3 p, float3 n, float3 d)
{
    float mag = (1.0 + max(max(abs(p.x), abs(p.y)), abs(p.z))) * 1e-4;
    mag = dot(n, d) >= 0.0 ? mag : -mag;
    return n * mag + p;
}
