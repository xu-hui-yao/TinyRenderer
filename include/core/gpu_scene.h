#pragma once

// ============================================================================
// GPU-ready flattened scene representation (Stage 0 of the CPU -> GPU port).
//
// This header defines plain-data (POD-ish, no shared_ptr/virtual dispatch)
// structures that mirror the information currently held by the virtual-
// dispatch Texture / BSDF / Emitter / Mesh class hierarchies. They are
// produced by an *additive* export layer (Texture::to_gpu_texture(),
// BSDF::to_gpu_material(), Emitter::to_gpu_light(), Scene::build_gpu_scene())
// that does not alter any existing rendering code path: the CPU renderer
// keeps using its virtual-dispatch classes as the ground-truth reference,
// while this flattened copy is what will eventually be uploaded to the GPU
// (Stage 1). Tagged-union "type" fields replace virtual dispatch, and
// indices replace shared_ptr for cross-references (textures referenced by
// materials, nested materials referenced by compound materials, etc.),
// exactly like Accel::ray_intersect's mesh_id already does for geometry.
// ============================================================================

#include <core/common.h>
#include <unordered_map>

M_NAMESPACE_BEGIN

// Texture is not among the forward declarations in core/common.h; declare it
// here since GPUSceneBuilder only ever needs a pointer/reference to it.
class Texture;

// ----------------------------------------------------------------------
// Textures
// ----------------------------------------------------------------------

enum class GPUTextureType : uint32_t { Constant = 0, Checkerboard, Bitmap, Count };

struct GPUTexture {
    GPUTextureType type   = GPUTextureType::Constant;
    bool spatial_varying  = false;
    int channels          = 3; // 1 (scalar) or 3 (RGB)

    // Constant: color0 holds the color (or a broadcast scalar in .x() when channels == 1)
    // Checkerboard: color0 / color1 are the two checker colors
    Color3f color0{ 1.0f };
    Color3f color1{ 0.0f };
    float scale_u = 1.0f, scale_v = 1.0f; // Checkerboard UV tiling

    // Bitmap: raw pixel data, row-major, `channels` floats per pixel.
    // Kept as an owned copy for now; Stage 1 (GPU upload) will decide how to
    // pack this into a shared texture atlas / image array.
    int width = 0, height = 0;
    std::vector<float> pixels;
};

// ----------------------------------------------------------------------
// Materials (BSDFs)
// ----------------------------------------------------------------------

enum class GPUMaterialType : uint32_t {
    Diffuse = 0,
    Conductor,
    Dielectric,
    ThinDielectric,
    Plastic,
    RoughConductor,
    RoughDielectric,
    RoughPlastic,
    Mask,
    BumpMap,
    TwoSided,
    Count
};

// A tagged union describing every BSDF type currently implemented by the
// renderer. Compound materials (Mask / BumpMap / TwoSided) reference their
// nested material(s) by index into GPUScene::materials instead of holding a
// shared_ptr<BSDF>, which keeps the representation a flat, GPU-uploadable
// array (a bounded-depth "resolve child index" loop on the GPU, instead of
// unrepresentable recursive virtual calls).
struct GPUMaterial {
    GPUMaterialType type = GPUMaterialType::Diffuse;
    uint32_t flags       = 0; // BSDFFlags bitmask (see components/bsdf.h)

    // Texture references (index into GPUScene::textures, M_INVALID_INDEX if unused)
    uint32_t tex_reflectance            = M_INVALID_INDEX; // Diffuse.reflectance; Plastic/RoughPlastic.diffuse_reflectance
    uint32_t tex_specular_reflectance   = M_INVALID_INDEX; // Conductor/RoughConductor/Dielectric/ThinDielectric/Plastic/RoughDielectric/RoughPlastic
    uint32_t tex_specular_transmittance = M_INVALID_INDEX; // Dielectric/ThinDielectric/RoughDielectric
    uint32_t tex_eta                    = M_INVALID_INDEX; // Conductor/RoughConductor (complex IOR, real part)
    uint32_t tex_k                      = M_INVALID_INDEX; // Conductor/RoughConductor (complex IOR, imaginary part)
    uint32_t tex_alpha                  = M_INVALID_INDEX; // RoughConductor/RoughDielectric (roughness texture)
    uint32_t tex_opacity                = M_INVALID_INDEX; // Mask
    uint32_t tex_bump                   = M_INVALID_INDEX; // BumpMap (height-field, sampled via eval_1_grad)

    // Nested material references (index into GPUScene::materials, M_INVALID_INDEX if unused)
    uint32_t nested_bsdf = M_INVALID_INDEX; // Mask / BumpMap
    uint32_t front_bsdf  = M_INVALID_INDEX; // TwoSided
    uint32_t back_bsdf   = M_INVALID_INDEX; // TwoSided (equals front_bsdf's index when not explicitly provided)

    // Scalar parameters (only the ones relevant to `type` are meaningful)
    float eta                       = 1.5f; // relative IOR: Dielectric/ThinDielectric/Plastic/RoughDielectric/RoughPlastic
    float inv_eta                   = 1.0f / 1.5f; // RoughDielectric
    float alpha                     = 0.1f; // RoughPlastic's fixed scalar roughness (RoughConductor/RoughDielectric use tex_alpha instead)
    float scale                     = 1.0f; // BumpMap displacement scale
    bool has_reflection             = true; // Dielectric/ThinDielectric/RoughDielectric
    bool has_transmission           = true; // Dielectric/ThinDielectric/RoughDielectric
    bool nonlinear                  = false; // Plastic/RoughPlastic
    float fdr_int                   = 0.0f; // Plastic
    float fdr_ext                   = 0.0f; // Plastic
    float inv_eta_2                 = 1.0f; // Plastic/RoughPlastic
    float specular_sampling_weight  = 0.5f; // Plastic/RoughPlastic lobe-selection weight

    // RoughPlastic: precomputed rough Fresnel transmittance LUT (indexed by cos_theta in [0,1]).
    // Kept as an owned copy for now; Stage 1 will pack this into a shared LUT buffer.
    std::vector<float> external_transmittance;
    float internal_reflectance = 0.0f;
};

// ----------------------------------------------------------------------
// Lights (Emitters)
// ----------------------------------------------------------------------

enum class GPULightType : uint32_t { Area = 0, Envmap, Count };

struct GPULight {
    GPULightType type = GPULightType::Area;

    // Area: owning mesh (index into GPUScene::meshes); its geometry is
    // already available in GPUScene for area sampling.
    uint32_t mesh_id      = M_INVALID_INDEX;
    uint32_t radiance_tex = M_INVALID_INDEX; // index into GPUScene::textures

    // Area: normalized CDF over the owning mesh's LOCAL triangle indices
    // (0..triangle_count-1), for exactly replicating Mesh::sample_position()'s
    // area-weighted triangle picking (see DiscreteDistribution1f::sample) on
    // the GPU via binary search. `cdf_offset`/`cdf_count` index into
    // GPUScene::light_triangle_cdf; `inv_total_area` == Mesh's m_area_pmf
    // normalization() (the resulting uniform-over-surface-area density).
    uint32_t cdf_offset    = M_INVALID_INDEX;
    uint32_t cdf_count     = 0;
    float inv_total_area   = 0.0f;
    // Global triangle id of this mesh's local triangle 0 (i.e. GPUMeshInfo's
    // base_index/3), needed to convert cdf_sample()'s LOCAL triangle index
    // result into a global id usable to index GPUScene::indices/positions.
    uint32_t base_triangle = 0;

    // Envmap-specific fields
    Transform4f to_world;
    Point3f bounding_sphere_center{ 0.0f };
    float bounding_sphere_radius = 0.0f;
    bool spatial_varying         = false;
    int distribution_width       = 0;
    int distribution_height      = 0;
    // Row-major, un-normalized level-0 luminance grid of the importance
    // sampling distribution. Building an actual GPU alias table/CDF from
    // this is deferred to Stage 1 (it is an algorithm, not just a data
    // layout, and is best validated against its GPU consumer directly).
    std::vector<float> distribution_luminance;
};

// ----------------------------------------------------------------------
// BVH (flattened acceleration structure, for in-shader traversal)
// ----------------------------------------------------------------------

// Mirrors BVHAccel's internal LinearBVHNode (see src/accelerate/bvh.cpp),
// packed to exactly 32 bytes (half a cache line) for GPU upload: bmin/offset
// form one float4-sized block, bmax/meta the other. `meta` packs the leaf
// primitive count into its low 30 bits (0 => interior node, matching the
// CPU-side n_primitives==0 convention) and the split axis into its top 2
// bits, since a leaf's primitive count comfortably fits well under 2^30 and
// axis only ever needs values 0/1/2.
struct GPUBVHNode {
    float3 bmin;
    uint32_t offset = 0; // leaf: index into GPUScene::bvh_primitives / interior: index of the right child (left child is always this node's index + 1)
    float3 bmax;
    uint32_t meta = 0;
};
static_assert(sizeof(GPUBVHNode) == 32, "GPUBVHNode must stay tightly packed to 32 bytes to match the Slang-side "
                                        "StructuredBuffer<GPUBVHNode> layout assumed by raytrace_debug.slang");

constexpr uint32_t M_BVH_META_COUNT_MASK = 0x3FFFFFFFu;
constexpr uint32_t M_BVH_META_AXIS_SHIFT = 30u;

inline uint32_t gpu_bvh_leaf_count(uint32_t meta) { return meta & M_BVH_META_COUNT_MASK; }

inline uint32_t gpu_bvh_axis(uint32_t meta) { return meta >> M_BVH_META_AXIS_SHIFT; }

inline uint32_t gpu_bvh_pack_meta(uint32_t n_primitives, uint32_t axis) {
    return (n_primitives & M_BVH_META_COUNT_MASK) | (axis << M_BVH_META_AXIS_SHIFT);
}

// ----------------------------------------------------------------------
// Camera (mirrors Camera::sample_ray()'s ray generation math)
// ----------------------------------------------------------------------

struct GPUCamera {
    // Row-major 4x4 matrices, i.e. element (r, c) is at index [r * 4 + c] -
    // matches TMatrix::operator()(row, col)'s m_data[row * Cols + col]
    // layout exactly, so the GPU-side matrix-vector multiply must also be
    // written as an explicit row-major dot product (see the Slang shader)
    // rather than relying on any particular shading-language default.
    float camera_to_world[16]  = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    float sample_to_camera[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    uint32_t width = 0, height = 0;
    float near_clip = 1e-4f, far_clip = 1e4f;
};
// This layout is byte-for-byte identical to a Slang/HLSL constant buffer
// struct made of 8 consecutive float4's (2x float[16], i.e. the two
// matrices decomposed into row vectors) followed by 4 tightly-packed 4-byte
// scalars - see raytrace_debug.slang's GPUCameraGPU. Uploading this struct's
// raw bytes directly therefore requires no host-side repacking. If this
// struct's field order/types ever change, double check that equivalence
// still holds (or repack explicitly instead of relying on it).
static_assert(sizeof(GPUCamera) == 144, "GPUCamera's byte layout is relied upon to exactly match "
                                        "raytrace_debug.slang's GPUCameraGPU constant buffer struct");

// ----------------------------------------------------------------------
// Geometry (concatenated across all meshes into global buffers)
// ----------------------------------------------------------------------

struct GPUMeshInfo {
    uint32_t base_vertex     = 0; // offset into GPUScene::vertex_* arrays
    uint32_t vertex_count    = 0;
    uint32_t base_index      = 0; // offset into GPUScene::indices (3 per triangle); base_index/3 is this mesh's first global triangle id
    uint32_t triangle_count  = 0;
    uint32_t material_id     = M_INVALID_INDEX; // index into GPUScene::materials
    uint32_t light_id         = M_INVALID_INDEX; // index into GPUScene::lights, if this mesh is an emitter
    bool has_normals          = false; // false => vertex_normals for this mesh's range are geometric placeholders
    bool has_uvs              = false; // false => vertex_uvs for this mesh's range are (0,0) placeholders
};

// ----------------------------------------------------------------------
// The full flattened scene
// ----------------------------------------------------------------------

struct GPUScene {
    // Geometry, concatenated across all meshes. `indices` stores GLOBAL
    // vertex indices (i.e. GPUMeshInfo::base_vertex has already been added),
    // exactly like a typical GPU vertex/index buffer pair. Global triangle
    // id `g` (used by bvh_primitives and triangle_material_id/
    // triangle_light_id below) refers to indices[g*3 + 0/1/2].
    std::vector<Point3f> vertex_positions;
    std::vector<Normal3f> vertex_normals;
    std::vector<Point2f> vertex_uvs;
    std::vector<uint32_t> indices; // 3 per triangle

    std::vector<GPUMeshInfo> meshes;
    std::vector<GPUMaterial> materials;
    std::vector<GPUTexture> textures;
    std::vector<GPULight> lights;
    GPUCamera camera;

    // Index into `lights`, or M_INVALID_INDEX if the scene has no environment emitter.
    uint32_t environment_light_id = M_INVALID_INDEX;

    // Flattened BVH (see Accel::export_gpu_bvh); empty if the scene's Accel
    // implementation does not expose a flat BVH (e.g. NaiveAccel/KDTreeAccel/
    // OctreeAccel - only BVHAccel currently does).
    std::vector<GPUBVHNode> bvh_nodes;
    std::vector<uint32_t> bvh_primitives; // global triangle ids, indexed by GPUBVHNode leaf offset/count

    // Denormalized per-global-triangle lookups (one entry per triangle,
    // parallel to the implicit global triangle id space), so a GPU shading
    // kernel can resolve a hit triangle's material/light in O(1) without
    // searching for which mesh it belongs to.
    std::vector<uint32_t> triangle_material_id;
    std::vector<uint32_t> triangle_light_id;

    // Shared backing storage for GPULight::cdf_offset/cdf_count (Area
    // lights' per-triangle area-weighted-sampling CDF; see GPULight).
    std::vector<float> light_triangle_cdf;
};


// ----------------------------------------------------------------------
// Builder: de-duplicates shared_ptr identity -> flat index while recursively
// invoking each object's to_gpu_*() export method. Implemented out-of-line
// in src/core/gpu_scene.cpp to avoid a header dependency cycle with
// components/texture.h, components/bsdf.h and components/emitter.h (which
// each need to reference GPUTexture / GPUMaterial / GPULight from this
// header for their own to_gpu_*() declarations).
// ----------------------------------------------------------------------

class GPUSceneBuilder {
public:
    // Resolves (and appends on first use) `tex` to an index into result.textures.
    // Returns M_INVALID_INDEX if `tex` is null.
    uint32_t add_texture(const std::shared_ptr<Texture> &tex);

    // Resolves (and appends on first use) `bsdf` to an index into result.materials.
    // May recurse for compound materials (Mask / BumpMap / TwoSided).
    // Returns M_INVALID_INDEX if `bsdf` is null.
    uint32_t add_material(const std::shared_ptr<BSDF> &bsdf);

    // Resolves (and appends on first use) `emitter` to an index into result.lights.
    // `mesh_id` is forwarded to Emitter::to_gpu_light() (meaningful for Area lights).
    // Returns M_INVALID_INDEX if `emitter` is null.
    uint32_t add_light(const std::shared_ptr<Emitter> &emitter, uint32_t mesh_id);

    GPUScene result;

private:
    std::unordered_map<const Texture *, uint32_t> texture_cache;
    std::unordered_map<const BSDF *, uint32_t> material_cache;
    std::unordered_map<const Emitter *, uint32_t> light_cache;
};

M_NAMESPACE_END
