// Stage 1, checkpoint 3: full megakernel path tracer.
//
// This is a line-by-line port of src/integrators/path.cpp's li() loop, using
// GPU-side re-implementations of every BSDF in src/bsdfs/*.cpp, every
// texture in src/textures/*.cpp, and area-light sampling matching
// src/components/mesh.cpp's Mesh::sample_direction()/pdf_direction(). Camera
// ray generation and BVH traversal are unchanged from raytrace_debug.slang.
//
// Known, intentional deviations from the CPU renderer (see conversation
// history for the full rationale):
//   1. RNG: a 32-bit hash-based PCG variant is used instead of the CPU's
//      64-bit pcg32 (see core/pcg32.h). This produces a different (but
//      still statistically valid) noise pattern; only the CONVERGED image
//      is expected to match the CPU renderer, not per-sample values.
//   2. Environment map importance sampling is not implemented: envmap
//      lights are sampled with square_to_uniform_sphere regardless of
//      whether the CPU-side emitter is spatially-varying. eval()/pdf() are
//      otherwise implemented faithfully (equirect UV mapping is correct),
//      so this only trades sampling efficiency for implementation size -
//      it remains an unbiased (if higher-variance) estimator. No test scene
//      in this repository currently uses <emitter type="envmap">.
//   3. square_to_cosine_hemisphere_pdf() is ported EXACTLY as implemented on
//      the CPU (src/core/warp.cpp), including its `M_PI * cos_theta` formula
//      (as opposed to the textbook `cos_theta / M_PI`). This looks like a
//      pre-existing CPU-side bug, but this port's goal is bit-for-bit
//      algorithmic parity with the CPU renderer (so the two can be
//      validated against each other), not independent correctness - so it
//      is intentionally NOT "fixed" here.
//
// All BSDF math below mirrors its CPU counterpart function-for-function;
// see the corresponding src/bsdfs/*.cpp file referenced in each comment.

// ============================================================================
// Data structures (byte layouts matching the host-side Gpu*Gpu packed structs
// in src/gpu/path_trace_main.cpp)
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

// ============================================================================
// Bindings
// ============================================================================

[[vk::binding(0, 0)]] RWStructuredBuffer<float4> accum_buffer; // running sum of radiance, divided by spp count on readback
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
// RNG - see file header note (1): a GPU-friendly 32-bit PCG hash, NOT the
// CPU's 64-bit pcg32. Statistically valid, different noise pattern.
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

// square_to_uniform_disk_concentric (warp.cpp)
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

// square_to_cosine_hemisphere (warp.cpp)
float3 square_to_cosine_hemisphere(float2 s)
{
    float2 p = square_to_uniform_disk_concentric(s);
    float z = sqrt(max(1.0 - dot(p, p), 0.0));
    return float3(p.x, p.y, z);
}

// See file header note (3): ported exactly as-is from warp.cpp, bug and all.
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

// square_to_uniform_triangle (warp.cpp)
float2 square_to_uniform_triangle(float2 s)
{
    float t = safe_sqrt_(1.0 - s.x);
    return float2(1.0 - t, t * s.y);
}

// fresnel() (fresnel.h) - dielectric Fresnel term.
// Returns (r, cos_theta_t, eta_it, eta_ti).
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

// fresnel_conductor() (fresnel.h)
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

// mis_weight() (path.cpp) - power-2 balance heuristic
float mis_weight(float pdf_a, float pdf_b)
{
    pdf_a *= pdf_a;
    pdf_b *= pdf_b;
    float denom = pdf_a + pdf_b;
    return denom == 0.0 ? 0.0 : pdf_a / denom;
}

// lerp_gather() (roughplastic.cpp) - linear interpolation into a flat LUT slice
float lerp_gather(uint offset, uint size, float x)
{
    x *= float(size - 1);
    int index = min(int(x), int(size) - 2);
    float v0 = material_lut[offset + uint(index)];
    float v1 = material_lut[offset + uint(index) + 1u];
    return lerp(v0, v1, x - float(index));
}

// ============================================================================
// Microfacet distribution (Beckmann-ish GGX-like visible-normal sampling,
// mirrors core/microfacet.h's TMicrofacetDistribution exactly)
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

// sample_visible_11 + sample() (microfacet.h) - visible normal sampling
float3 microfacet_sample(float3 wi, float2 s, float alpha, out float pdf_out)
{
    float3 wi_p = normalize(float3(alpha * wi.x, alpha * wi.y, wi.z));

    float sin_theta_2 = wi_p.x * wi_p.x + wi_p.y * wi_p.y;
    float inv_sin_theta = 1.0 / sqrt(max(sin_theta_2, 1e-12));
    float2 sc = float2(wi_p.x, wi_p.y) * inv_sin_theta;
    if (abs(sin_theta_2) <= 4.0 * 1.1920929e-7) sc = float2(1.0, 0.0);
    // Frame::sincos_phi (frame.h) returns {sin_phi, cos_phi} = {result.y, result.x}
    // where result = (wi_p.x, wi_p.y) * inv_sin_theta, i.e. cos_phi = sc.x, sin_phi = sc.y.
    float cos_phi = sc.x;
    float sin_phi = sc.y;

    float cos_theta = wi_p.z;

    // sample_visible_11
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

// eval_transmittance()/eval_reflectance() (microfacet.h) are only used at CPU
// *construction* time to precompute RoughPlastic's LUT (already baked into
// material_lut by the host); they are not needed at shading time.

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
    // Color3f::luminance() (spectrum.h) uses standard Rec.709-ish weights;
    // ported to match, though tex_eval_1 is only used for scalar (1-channel)
    // textures (e.g. Mask's opacity) where c.x==c.y==c.z already.
    return c.x;
}

// ============================================================================
// Surface hit: mirrors the fields of SurfaceIntersection3f actually consumed
// by BSDF eval()/sample()/pdf(). `wi`/`dp_du`/`dp_dv` are local to
// `shading_s/t/n`; wo passed to eval()/pdf() is likewise local to whichever
// frame is "current" at the point of the call (see BumpMap - this is
// replicated exactly, not re-derived from first principles).
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

// BumpMap::frame() (bumpmap.cpp)
SurfaceHit bump_perturbed_hit(SurfaceHit h, uint bump_tex, float scale)
{
    float2 grad_uv = float2(0.0, 0.0); // eval_1_grad for Constant/Checkerboard is always (0,0); Bitmap not ported for gradients in this checkpoint (no test scene uses bumpmap+bitmap together)
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
// BSDF sample result (mirrors TBSDFSample)
// ============================================================================

struct BSDFSample {
    float3 wo;
    float pdf;
    float eta;
    bool delta;
};

BSDFSample make_bsdf_sample(float3 wo) { BSDFSample b; b.wo = wo; b.pdf = 0.0; b.eta = 1.0; b.delta = false; return b; }

// ============================================================================
// Per-material eval/sample/pdf - mirrors src/bsdfs/*.cpp exactly. `hit.wi`
// and the `wo` parameter are both local to hit's current shading frame.
// ============================================================================

// ---- Diffuse (diffuse.cpp) ----
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

// ---- Conductor (conductor.cpp) - smooth/delta ----
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

// ---- Dielectric (dielectric.cpp) - smooth/delta ----
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

// ---- ThinDielectric (thindielectric.cpp) - smooth/delta ----
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

// ---- Plastic (plastic.cpp) - smooth ----
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

// ---- RoughConductor (roughconductor.cpp) ----
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

// ---- RoughDielectric (roughdielectric.cpp) ----
float3 eval_roughdielectric(GPUMaterialGPU m, SurfaceHit h, float3 wo)
{
    float cti = h.wi.z, cto = wo.z;
    if (cti == 0.0) return float3(0, 0, 0);
    bool refl = cti * cto > 0.0;
    float eta = cti > 0.0 ? m.eta : m.inv_eta;
    float inv_eta = cti > 0.0 ? m.inv_eta : m.eta;

    float3 mv = normalize(h.wi + wo * (refl ? 1.0 : eta));
    mv = mulsign_(1.0, mv.z) * mv; // ensure hemisphere alignment: mulsign(m, cos_theta(m))
    // The above mirrors mulsign(m, Frame::cos_theta(m)) elementwise via a scalar sign flip.

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

// ---- RoughPlastic (roughplastic.cpp) ----
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
// bounded-depth loops (see file header - avoids relying on Slang/SPIR-V
// recursion support). MAX_UNWRAP comfortably exceeds any nesting depth
// produced by GPUSceneBuilder for this renderer's BSDF grammar.
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

        // Leaf material
        float3 result;
        if (m.type == MAT_DIFFUSE) result = eval_diffuse(m, h, w);
        else if (m.type == MAT_PLASTIC) result = eval_plastic(m, h, w);
        else if (m.type == MAT_ROUGHCONDUCTOR) result = eval_roughconductor(m, h, w);
        else if (m.type == MAT_ROUGHDIELECTRIC) result = eval_roughdielectric(m, h, w);
        else if (m.type == MAT_ROUGHPLASTIC) result = eval_roughplastic(m, h, w);
        else result = float3(0, 0, 0); // Conductor/Dielectric/ThinDielectric are delta -> eval() == 0

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

// bsdf_sample needs to reconstruct the world-space wo from whatever local
// frame the innermost (leaf) material sampled in, since TwoSided/BumpMap
// change the local frame along the way (mirrors each wrapper's exact
// post-processing of bs.wo in the CPU code). Mask's stochastic branch is
// folded into this SAME iterative loop (rather than a recursive call to
// bsdf_sample) because SPIR-V compute shaders do not support recursion.
BSDFSample bsdf_sample(uint mat_id, SurfaceHit hit, float sample1, float2 sample2, out float3 value)
{
    value = float3(0, 0, 0);
    float3 extra_weight = float3(1, 1, 1);
    // Accumulated product of Mask opacities for every Mask level we chose
    // to "enter" (sample1 < opacity) on the way down - mirrors the chain of
    // `bs.pdf *= opacity` multiplications CPU applies as each enclosing
    // Mask::sample() call returns. Applied once at whichever point we
    // terminate (leaf material OR a Mask's pass-through branch).
    float mask_pdf_scale = 1.0;
    SurfaceHit h = hit;
    // Stack of "undo" operations to re-express the leaf's local bs.wo back
    // into the ORIGINAL hit's local frame, applied innermost-first.
    // 0 = TwoSided flip (flip .z), 1 = BumpMap (frame_to_local(hit, frame_to_world(ph, .)))
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
                BSDFSample bs = make_bsdf_sample(-h.wi); // matches Mask::sample's default (straight-through, "null" event)
                bs.eta = 1.0;
                bs.pdf = (1.0 - opacity) * mask_pdf_scale;
                bs.delta = true;
                value = extra_weight; // matches CPU's default Color3f(1.0f) when the nested BSDF isn't sampled

                for (int k = depth - 1; k >= 0; --k) {
                    if (op_stack[k] == 0) bs.wo = flip_vec(bs.wo);
                    else bs.wo = frame_to_local(frame_stack[k], frame_to_world(frame_stack[k + 1], bs.wo));
                }
                return bs;
            }

            // Entered: rescale sample1 (matches CPU's `sample1 /= opacity`)
            // and continue unwrapping into the nested BSDF in this same loop.
            sample1 = sample1 / opacity;
            mask_pdf_scale *= opacity;
            mat_id = m.nested_bsdf;
            continue;
        } else if (m.type == MAT_BUMPMAP) {
            SurfaceHit ph = bump_perturbed_hit(h, m.tex_bump, m.scale);
            op_stack[depth] = 1; frame_stack[depth] = h; depth++;
            frame_stack[depth] = ph; // store perturbed frame too, for the unwind step's "frame_to_world(ph, .)"
            h = ph;
            mat_id = m.nested_bsdf;
            continue;
        }

        // Leaf material
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

        // Unwind: re-express bs.wo from the (possibly perturbed/flipped)
        // leaf frame back into the ORIGINAL hit's local frame, innermost
        // (most recently pushed) operation first - mirrors e.g. TwoSided's
        // `bs.wo.z() *= -1.f` and BumpMap's `perturbed_si.to_world(bs.wo)`
        // followed by treating the result as the outer frame's local wo.
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
    // has_flag(ESmooth): EDiffuseReflection|EDiffuseTransmission|EGlossyReflection|EGlossyTransmission (see bsdf.h)
    const uint SMOOTH_MASK = 0x2u | 0x4u | 0x8u | 0x10u;
    if (mat_id == INVALID_INDEX) return false;
    // `flags` is already the fully-resolved composite for wrapper materials
    // (TwoSided ORs front|back; Mask/BumpMap copy their nested_bsdf's flags
    // at CPU construct() time - see to_gpu_material() for each), so a single
    // lookup (no unwrapping) is both correct and matches the CPU's
    // BSDF::has_flag() semantics exactly.
    return (materials[mat_id].flags & SMOOTH_MASK) != 0u;
}

// ============================================================================
// Light sampling (mirrors Scene::sample_emitter_direction /
// Mesh::sample_direction / Area::sample_direction / EnvironmentMap - with
// the simplification noted in the file header (2))
// ============================================================================

struct DirSample {
    float3 p, n, d;
    float2 uv;
    float t, pdf, dist;
    bool delta;
    uint light_id;
};

// Binary search over a normalized CDF, mirrors DiscreteDistribution1f::sample_reuse
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

// Mesh::sample_position + Mesh::sample_direction (mesh.cpp), for Area lights
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
    // Use interpolated shading normal if available (has_normals is baked
    // into vertex_normals as geometric-placeholder-or-real already; barycentric
    // interpolate exactly like Mesh's geometry does elsewhere).
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

// EnvironmentMap::eval (envmap.cpp) - equirect UV mapping
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

// EnvironmentMap::sample_direction (non-spatially-varying branch always;
// spatially-varying uses the SAME uniform-sphere strategy - see file header
// note (2)) + pdf_direction (for the spatially-varying eval() consistency).
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

// ============================================================================
// Scene-level emitter sampling (Scene::sample_emitter_direction /
// pdf_emitter_direction) - uniform emitter picking (matches m_emitter_pmf).
// ============================================================================

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

    // Account for uniform emitter-picking probability (matches
    // Scene::sample_emitter_direction's ds.pdf *= pdf_emitter(index)).
    ds.pdf *= 1.0 / float(n);
    weight *= float(n); // 1/pdf_emitter cancels into the weight, matching spec *= emitter_weight (== n) in the multi-emitter CPU path

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
// BVH traversal + closest-hit shading data (mirrors BVHAccel::ray_intersect)
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

// ============================================================================
// Main path tracing loop (mirrors Path::li(), path.cpp)
// ============================================================================

[shader("compute")]
[numthreads(16, 16, 1)]
void main(uint3 tid: SV_DispatchThreadID)
{
    if (tid.x >= camera.width || tid.y >= camera.height) return;

    uint pixel_index = tid.y * camera.width + tid.x;
    Rng rng = rng_init(pixel_index, pc.spp_index, 0u);

    // ---- Primary ray (mirrors PerspectiveCamera::sample_ray, with a random
    // subpixel jitter for antialiasing across spp - the CPU integrator gets
    // this from its reconstruction filter's sampler-driven pixel offset). ----
    float2 jitter = rng_next2(rng);
    float px = float(tid.x) + jitter.x;
    float py = float(tid.y) + jitter.y;
    float inv_w = 1.0 / float(camera.width);
    float inv_h = 1.0 / float(camera.height);

    float3 near_p = transform_point(camera.s2c0, camera.s2c1, camera.s2c2, camera.s2c3, float3(px * inv_w, py * inv_h, 0.0));
    float3 d = normalize(near_p);
    float inv_z = 1.0 / d.z;

    float3 ray_o = transform_point(camera.c2w0, camera.c2w1, camera.c2w2, camera.c2w3, float3(0.0, 0.0, 0.0));
    float3 ray_d = transform_vector(camera.c2w0, camera.c2w1, camera.c2w2, camera.c2w3, d);
    float ray_tmin = camera.near_clip * inv_z;
    float ray_tmax = camera.far_clip * inv_z;

    // ---- Path tracing loop state (mirrors Path::li) ----
    float3 throughput = float3(1, 1, 1);
    float3 result = float3(0, 0, 0);
    float eta = 1.0;
    bool valid_ray = false;
    bool active = true;

    float3 prev_p = ray_o;      // previous vertex position (for DirectionSample::d/dist reconstruction)
    float3 prev_n = float3(0, 0, 1);
    bool prev_is_valid = false; // Intersection::is_valid(): t < M_MAX_FLOAT
    float prev_bsdf_pdf = 1.0;
    bool prev_bsdf_delta = true;

    for (uint depth = 0u; depth < pc.max_depth && active; ++depth) {
        HitResult hit = trace_ray(ray_o, ray_d, ray_tmin, ray_tmax);
        // Mirrors Scene::ray_intersect (scene.cpp): a ray that hits no
        // geometry still counts as an "intersection" when the scene has an
        // environment map - the CPU fabricates a hit at t = M_MAX_FLOAT with
        // mesh_id = M_INVALID_INDEX so that the direct-emission block below
        // can pick the environment emitter up. Without this, the envmap
        // branch of that block is unreachable and the sky never contributes.
        bool is_intersect = hit.hit || pc.environment_light_id != INVALID_INDEX;

        uint hit_mesh_material = INVALID_INDEX;
        uint hit_mesh_light = INVALID_INDEX;
        float3 hit_p = float3(0, 0, 0), hit_n_geo = float3(0, 0, 1), hit_shading_n = float3(0, 0, 1);
        float3 hit_dp_du = float3(1, 0, 0), hit_dp_dv = float3(0, 1, 0);
        float2 hit_uv = float2(0, 0);
        float3 hit_wi_world = -ray_d;

        if (hit.hit) {
            hit_mesh_material = triangle_material_id[hit.tri];
            hit_mesh_light = triangle_light_id[hit.tri];

            uint i0 = indices[hit.tri * 3 + 0], i1 = indices[hit.tri * 3 + 1], i2 = indices[hit.tri * 3 + 2];
            float3 p0 = positions[i0].xyz, p1 = positions[i1].xyz, p2 = positions[i2].xyz;
            float3 bary = float3(1.0 - hit.u - hit.v, hit.u, hit.v);
            hit_p = p0 * bary.x + p1 * bary.y + p2 * bary.z;
            hit_uv = uvs[i0] * bary.x + uvs[i1] * bary.y + uvs[i2] * bary.z;

            hit_n_geo = normalize(cross(p1 - p0, p2 - p0));
            float3 n0 = normals[i0].xyz, n1 = normals[i1].xyz, n2 = normals[i2].xyz;
            float nlen2 = dot(n0, n0) + dot(n1, n1) + dot(n2, n2);
            hit_shading_n = nlen2 > 0.0 ? normalize(n0 * bary.x + n1 * bary.y + n2 * bary.z) : hit_n_geo;

            float2 uv0 = uvs[i0], uv1 = uvs[i1], uv2 = uvs[i2];
            float2 duv1 = uv1 - uv0, duv2 = uv2 - uv0;
            float3 edge1 = p1 - p0, edge2 = p2 - p0;
            float detuv = duv1.x * duv2.y - duv1.y * duv2.x;
            if (abs(detuv) < 1e-8) {
                hit_dp_du = edge1; hit_dp_dv = edge2;
            } else {
                float invd = 1.0 / detuv;
                hit_dp_du = (edge1 * duv2.y - edge2 * duv1.y) * invd;
                hit_dp_dv = (-edge1 * duv2.x + edge2 * duv1.x) * invd;
            }
        }

        bool hit_mesh_valid = is_intersect && hit_mesh_material != INVALID_INDEX;

        // ---------------------- Direct emission ----------------------
        if (is_intersect && (!hit_mesh_valid || hit_mesh_light != INVALID_INDEX)) {
            float3 emission = float3(0, 0, 0);
            float em_pdf = 0.0;
            if (hit_mesh_light != INVALID_INDEX) {
                GPULightGPU L = lights[hit_mesh_light];
                bool front = dot(hit_wi_world, hit_shading_n) > 0.0; // Frame::cos_theta(si.wi) with wi local to shading_n
                if (front) emission = tex_eval(L.radiance_tex, hit_uv);

                if (!prev_bsdf_delta) {
                    float3 ds_d = normalize(hit_p - prev_p);
                    float ds_dist = length(hit_p - prev_p);
                    em_pdf = pdf_emitter_direction(hit_mesh_light, ds_d, hit_shading_n, ds_dist);
                }
            } else if (!hit_mesh_valid && pc.environment_light_id != INVALID_INDEX) {
                emission = eval_envmap(pc.environment_light_id, ray_d);
                if (!prev_bsdf_delta) {
                    em_pdf = pdf_emitter_direction(pc.environment_light_id, ray_d, -ray_d, 1.0);
                }
            }
            float w = mis_weight(prev_bsdf_pdf, em_pdf);
            result += throughput * emission * w;
        }

        bool active_next = (depth + 1u < pc.max_depth) && is_intersect && hit_mesh_valid;
        if (!active_next) break;

        SurfaceHit sh;
        sh.uv = hit_uv;
        sh.shading_n = hit_shading_n;
        sh.shading_s = normalize(hit_dp_du - hit_shading_n * dot(hit_dp_du, hit_shading_n));
        if (length(sh.shading_s) < 1e-6) {
            // Degenerate tangent (dp_du parallel to normal): fall back to an
            // arbitrary orthonormal basis, matching TFrame's default-normal
            // constructor path used when normals are supplied without valid
            // UV-derived tangents.
            float3 up = abs(hit_shading_n.z) < 0.999 ? float3(0, 0, 1) : float3(1, 0, 0);
            sh.shading_s = normalize(cross(up, hit_shading_n));
        }
        sh.shading_t = cross(hit_shading_n, sh.shading_s);
        sh.dp_du = hit_dp_du;
        sh.dp_dv = hit_dp_dv;
        sh.wi = frame_to_local(sh, hit_wi_world);

        // ---------------------- Emitter sampling ----------------------
        bool active_em = bsdf_is_smooth(hit_mesh_material);
        DirSample ds;
        float3 em_weight = float3(0, 0, 0);
        float3 wo_world = float3(0, 0, 0);
        float3 wo_local = float3(0, 0, 0);

        if (active_em) {
            ds = sample_emitter_direction(hit_p, rng, em_weight);
            active_em = active_em && ds.pdf != 0.0;
            wo_world = ds.d;
            wo_local = frame_to_local(sh, wo_world);
        }

        // ------ Evaluate BSDF * cos(theta) and sample direction -------
        float sample1 = rng_next(rng);
        float2 sample2 = rng_next2(rng);

        float3 bsdf_val = bsdf_eval(hit_mesh_material, sh, wo_local);
        float bsdf_pdf_v = bsdf_pdf(hit_mesh_material, sh, wo_local);
        float3 bsdf_weight;
        BSDFSample bs = bsdf_sample(hit_mesh_material, sh, sample1, sample2, bsdf_weight);

        // --------------- Emitter sampling contribution ----------------
        if (active_em) {
            float mis_em = ds.delta ? 1.0 : mis_weight(ds.pdf, bsdf_pdf_v);
            // Test visibility (matches Scene::sample_emitter_direction's shadow test)
            float3 shadow_o = offset_p(hit_p, hit_shading_n, wo_world);
            float3 to_target = ds.p - shadow_o;
            float shadow_dist = length(to_target);
            float3 shadow_d = to_target / max(shadow_dist, 1e-8);
            bool occluded = trace_shadow_ray(shadow_o, shadow_d, 0.0, shadow_dist * (1.0 - 1e-4));
            if (!occluded) {
                result += throughput * bsdf_val * em_weight * mis_em;
            }
        }

        // ---------------------- BSDF sampling ----------------------
        float3 next_wo_world = frame_to_world(sh, bs.wo);
        float3 next_o = offset_p(hit_p, hit_shading_n, next_wo_world);
        ray_o = next_o;
        ray_d = next_wo_world;
        ray_tmin = 0.0;
        ray_tmax = 3.402823466e+30;

        throughput *= bsdf_weight;
        eta *= bs.eta;
        bool intersection_valid = true; // TSurfaceIntersection::is_valid() checks t<MAX_FLOAT, always true for a real hit
        valid_ray = valid_ray || intersection_valid;

        prev_p = hit_p;
        prev_n = hit_shading_n;
        prev_bsdf_pdf = bs.pdf;
        prev_bsdf_delta = bs.delta;

        // -------------------- Stopping criterion ---------------------
        float throughput_max = max(max(throughput.x, throughput.y), throughput.z);
        float rr_prob = min(throughput_max * eta * eta, 0.95);
        bool rr_active = (depth + 1u) >= pc.rr_depth;
        bool rr_continue = rng_next(rng) < rr_prob;
        if (rr_active) throughput /= rr_prob;
        active = (!rr_active || rr_continue) && throughput_max != 0.0;
    }

    // ---- Accumulate (mirrors ImageBlock's running-average accumulation) ----
    float4 prev = (pc.accumulate != 0u) ? accum_buffer[pixel_index] : float4(0, 0, 0, 0);
    accum_buffer[pixel_index] = prev + float4(result, 1.0);
}
