#pragma once

// ============================================================================
// Host-side packing + Vulkan upload of a GPUScene into the shared
// "bindings 1..15" scene layout that every compute backend in this renderer
// consumes.
//
// This is the host-side counterpart of src/gpu/shaders/common/scene_common.slang:
// the packed structs below are byte-for-byte identical to that file's
// GPU*GPU structs, and upload_scene_buffers() fills the descriptor-set slots
// in exactly the binding order it declares. The two files must be kept in
// sync - if you add a field to one, add it to the other.
//
// WHY THIS IS ITS OWN TRANSLATION UNIT: this code used to live in an
// anonymous namespace inside src/gpu/gpu_renderer.cpp, where it was reachable
// only by the offline megakernel/wavefront render paths. The real-time
// interactive backend (include/gpu/gpu_session.h) needs the very same packing
// and upload, and duplicating ~250 lines of struct-field copying is exactly
// how the host and shader layouts silently drift apart. Both now call into
// this single implementation.
//
// NOTE: the standalone diagnostic executables (src/gpu/path_trace_main.cpp,
// src/gpu/wavefront_main.cpp) deliberately keep their own private copies of
// these structs - they recompile the whole core-sources glob into separate
// binaries on purpose (see src/gpu/CMakeLists.txt) and are not part of the
// tiny-renderer link, so they are intentionally left untouched here.
// ============================================================================

// core/gpu_scene.h only pulls in core/common.h, which merely FORWARD-DECLARES
// TArray/TRGBSpectrum/TTransform - it relies on its includer having already
// brought in the definitions. gpu_renderer.cpp got them transitively via
// gpu/gpu_renderer.h -> components/scene.h; this header has no such
// dependency, so it names what it actually needs explicitly.
#include <core/array.h>
#include <core/spectrum.h>
#include <core/transform.h>

#include <core/gpu_scene.h>
#include <gpu/gpu_buffer.h>
#include <gpu/gpu_triangles.h>
#include <gpu/vk_context.h>
#include <vector>

namespace tiny_renderer::gpu {

// ---- Packed structs, byte-for-byte matching common/scene_common.slang's
// GPU*GPU structs (shared by every backend's binding-1..15 scene layout). ----
#pragma pack(push, 1)
struct GPUTextureGPU {
    uint32_t type;
    uint32_t channels;
    uint32_t width;
    uint32_t height;
    uint32_t pixel_offset;
    float color0_r, color0_g, color0_b;
    float color1_r, color1_g, color1_b;
    float scale_u, scale_v;
};

struct GPUMaterialGPU {
    uint32_t type;
    uint32_t flags;
    uint32_t tex_reflectance;
    uint32_t tex_specular_reflectance;
    uint32_t tex_specular_transmittance;
    uint32_t tex_eta;
    uint32_t tex_k;
    uint32_t tex_alpha;
    uint32_t tex_opacity;
    uint32_t tex_bump;
    uint32_t nested_bsdf;
    uint32_t front_bsdf;
    uint32_t back_bsdf;
    float eta;
    float inv_eta;
    float alpha;
    float scale;
    uint32_t has_reflection;
    uint32_t has_transmission;
    uint32_t nonlinear;
    float fdr_int;
    float fdr_ext;
    float inv_eta_2;
    float specular_sampling_weight;
    uint32_t transmittance_lut_offset;
    uint32_t transmittance_lut_count;
    float internal_reflectance;
};

struct GPULightGPU {
    uint32_t type;
    uint32_t mesh_id;
    uint32_t radiance_tex;
    uint32_t cdf_offset;
    uint32_t cdf_count;
    float inv_total_area;
    uint32_t base_triangle;
    uint32_t spatial_varying;
    float to_world[16];
    float to_world_inv[16];
    float bsc_x, bsc_y, bsc_z;
    float bounding_sphere_radius;
};
#pragma pack(pop)

// Everything needed to upload a GPUScene's shared scene bindings - the
// variable-length arrays are flattened/float4-padded here so the upload step
// is a straight memcpy per buffer.
struct PackedScene {
    std::vector<float> packed_positions; // float4-padded
    std::vector<float> packed_normals;   // float4-padded
    // Triangle vertices denormalised into a flat array: 3 x float4 per
    // triangle, in the order the BVH's leaf primitives reference them.
    //
    // BVH traversal is the single hottest loop in this renderer, and doing
    // positions[indices[tri*3+k]] there costs three dependent gather chains (an
    // index load, then a vertex load that depends on it) for every triangle
    // test. Pre-expanding them turns that into one contiguous read the hardware
    // can coalesce. Duplicated vertices cost memory but no correctness: it is
    // the same positions, just laid out for the access pattern.
    std::vector<float> packed_triangles;
    std::vector<float> packed_uvs;
    std::vector<GPUTextureGPU> gpu_textures;
    std::vector<float> texture_pixels;
    std::vector<GPUMaterialGPU> gpu_materials;
    std::vector<float> material_lut;
    std::vector<GPULightGPU> gpu_lights;
};

PackedScene pack_gpu_scene(const GPUScene &gs);

// Deepest level of the BVH, to check against the shaders' fixed-size traversal
// stack (M_BVH_STACK_DEPTH in common/scene_common.slang, currently 32). An
// overflow there is silent, so it is measured rather than assumed.
int compute_bvh_depth(const GPUScene &gs);

// The uploaded scene buffers, in binding order. `camera_buf` is a uniform
// buffer (binding 1); everything else is a storage buffer (bindings 2..15).
struct SceneBufferSet {
    GpuBuffer positions_buf, normals_buf, uvs_buf, indices_buf;
    GpuBuffer bvh_nodes_buf, bvh_prims_buf, tri_mat_buf, tri_light_buf;
    GpuBuffer materials_buf, textures_buf, texture_pixels_buf, material_lut_buf;
    GpuBuffer lights_buf, light_cdf_buf, camera_buf;
    GpuBuffer packed_triangles_buf;

    void destroy(VulkanContext &ctx);
};

SceneBufferSet upload_scene_buffers(VulkanContext &ctx, const GPUScene &gs, const PackedScene &p);

} // namespace tiny_renderer::gpu
