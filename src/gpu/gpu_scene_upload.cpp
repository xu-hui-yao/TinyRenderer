// Host-side GPUScene packing + upload. See include/gpu/gpu_scene_upload.h for
// why this lives in its own translation unit rather than in an anonymous
// namespace inside gpu_renderer.cpp.
//
// The field-by-field copying below is verbatim from that former anonymous
// namespace - it is a pure move, so the offline megakernel/wavefront backends
// produce bit-identical buffers to before.

#include <gpu/gpu_scene_upload.h>

namespace tiny_renderer::gpu {

// The shaders' traversal stack is a fixed-size array (M_BVH_STACK_DEPTH in
// common/scene_common.slang). Overflowing it is not caught by anything: the
// write just goes out of bounds and traversal silently returns wrong hits,
// which looks like holes in the geometry rather than a crash. So the depth is
// measured here, where it can actually be reported.
int compute_bvh_depth(const GPUScene &gs) {
    if (gs.bvh_nodes.empty()) return 0;

    // Iterative walk with an explicit stack: a recursive version would risk
    // overflowing the HOST stack on a deep tree, which is the exact failure this
    // function exists to detect.
    std::vector<std::pair<uint32_t, int>> todo; // (node index, depth)
    todo.emplace_back(0, 1);
    int max_depth = 0;

    while (!todo.empty()) {
        const auto [idx, depth] = todo.back();
        todo.pop_back();
        if (idx >= gs.bvh_nodes.size()) continue;
        if (depth > max_depth) max_depth = depth;

        const GPUBVHNode &n = gs.bvh_nodes[idx];
        // Leaf nodes have a primitive count; interior nodes instead point at a
        // first child and leave the other implicit (see the traversal code).
        if ((n.meta & 0x3FFFFFFFu) == 0u) {
            todo.emplace_back(n.offset, depth + 1);
            todo.emplace_back(idx + 1, depth + 1);
        }
    }
    return max_depth;
}

PackedScene pack_gpu_scene(const GPUScene &gs) {
    PackedScene p;

    p.packed_positions.resize(gs.vertex_positions.size() * 4);
    p.packed_normals.resize(gs.vertex_normals.size() * 4);
    for (size_t i = 0; i < gs.vertex_positions.size(); ++i) {
        const auto &pos               = gs.vertex_positions[i];
        p.packed_positions[i * 4 + 0] = pos.x();
        p.packed_positions[i * 4 + 1] = pos.y();
        p.packed_positions[i * 4 + 2] = pos.z();
        p.packed_positions[i * 4 + 3] = 0.0f;
        const auto &n                 = gs.vertex_normals[i];
        p.packed_normals[i * 4 + 0]   = n.x();
        p.packed_normals[i * 4 + 1]   = n.y();
        p.packed_normals[i * 4 + 2]   = n.z();
        p.packed_normals[i * 4 + 3]   = 0.0f;
    }
    p.packed_uvs.resize(gs.vertex_uvs.size() * 2);
    for (size_t i = 0; i < gs.vertex_uvs.size(); ++i) {
        p.packed_uvs[i * 2 + 0] = gs.vertex_uvs[i].x();
        p.packed_uvs[i * 2 + 1] = gs.vertex_uvs[i].y();
    }

    p.packed_triangles = build_packed_triangles(gs);

    p.gpu_textures.resize(gs.textures.size());
    for (size_t i = 0; i < gs.textures.size(); ++i) {
        const GPUTexture &t = gs.textures[i];
        GPUTextureGPU &g    = p.gpu_textures[i];
        g.type              = static_cast<uint32_t>(t.type);
        g.channels          = static_cast<uint32_t>(t.channels);
        g.width             = static_cast<uint32_t>(t.width);
        g.height            = static_cast<uint32_t>(t.height);
        g.color0_r = t.color0(0); g.color0_g = t.color0(1); g.color0_b = t.color0(2);
        g.color1_r = t.color1(0); g.color1_g = t.color1(1); g.color1_b = t.color1(2);
        g.scale_u = t.scale_u; g.scale_v = t.scale_v;
        if (!t.pixels.empty()) {
            g.pixel_offset = static_cast<uint32_t>(p.texture_pixels.size());
            p.texture_pixels.insert(p.texture_pixels.end(), t.pixels.begin(), t.pixels.end());
        } else {
            g.pixel_offset = 0;
        }
    }

    p.gpu_materials.resize(gs.materials.size());
    for (size_t i = 0; i < gs.materials.size(); ++i) {
        const GPUMaterial &m = gs.materials[i];
        GPUMaterialGPU &g    = p.gpu_materials[i];
        g.type                       = static_cast<uint32_t>(m.type);
        g.flags                      = m.flags;
        g.tex_reflectance            = m.tex_reflectance;
        g.tex_specular_reflectance   = m.tex_specular_reflectance;
        g.tex_specular_transmittance = m.tex_specular_transmittance;
        g.tex_eta                    = m.tex_eta;
        g.tex_k                      = m.tex_k;
        g.tex_alpha                  = m.tex_alpha;
        g.tex_opacity                = m.tex_opacity;
        g.tex_bump                   = m.tex_bump;
        g.nested_bsdf                = m.nested_bsdf;
        g.front_bsdf                 = m.front_bsdf;
        g.back_bsdf                  = m.back_bsdf;
        g.eta                        = m.eta;
        g.inv_eta                    = m.inv_eta;
        g.alpha                      = m.alpha;
        g.scale                      = m.scale;
        g.has_reflection             = m.has_reflection ? 1u : 0u;
        g.has_transmission           = m.has_transmission ? 1u : 0u;
        g.nonlinear                  = m.nonlinear ? 1u : 0u;
        g.fdr_int                    = m.fdr_int;
        g.fdr_ext                    = m.fdr_ext;
        g.inv_eta_2                  = m.inv_eta_2;
        g.specular_sampling_weight   = m.specular_sampling_weight;
        g.internal_reflectance       = m.internal_reflectance;
        if (!m.external_transmittance.empty()) {
            g.transmittance_lut_offset = static_cast<uint32_t>(p.material_lut.size());
            g.transmittance_lut_count  = static_cast<uint32_t>(m.external_transmittance.size());
            p.material_lut.insert(p.material_lut.end(), m.external_transmittance.begin(),
                                   m.external_transmittance.end());
        } else {
            g.transmittance_lut_offset = 0;
            g.transmittance_lut_count  = 0;
        }
    }
    if (p.material_lut.empty()) p.material_lut.push_back(0.0f);
    if (p.texture_pixels.empty()) p.texture_pixels.push_back(0.0f);

    p.gpu_lights.resize(gs.lights.size());
    for (size_t i = 0; i < gs.lights.size(); ++i) {
        const GPULight &l = gs.lights[i];
        GPULightGPU &g    = p.gpu_lights[i];
        g.type            = static_cast<uint32_t>(l.type);
        g.mesh_id         = l.mesh_id;
        g.radiance_tex    = l.radiance_tex;
        g.cdf_offset      = l.cdf_offset;
        g.cdf_count       = l.cdf_count;
        g.inv_total_area  = l.inv_total_area;
        g.base_triangle   = l.base_triangle;
        g.spatial_varying = l.spatial_varying ? 1u : 0u;
        for (int r = 0; r < 4; ++r) {
            for (int c = 0; c < 4; ++c) {
                g.to_world[r * 4 + c]     = l.to_world.get_transform()(r, c);
                g.to_world_inv[r * 4 + c] = l.to_world.get_inverse()(r, c);
            }
        }
        g.bsc_x = l.bounding_sphere_center.x();
        g.bsc_y = l.bounding_sphere_center.y();
        g.bsc_z = l.bounding_sphere_center.z();
        g.bounding_sphere_radius = l.bounding_sphere_radius;
    }
    if (p.gpu_lights.empty()) p.gpu_lights.emplace_back();

    return p;
}

void SceneBufferSet::destroy(VulkanContext &ctx) {
    for (GpuBuffer *b : { &positions_buf, &normals_buf, &uvs_buf, &indices_buf, &bvh_nodes_buf, &bvh_prims_buf,
                          &tri_mat_buf, &tri_light_buf, &materials_buf, &textures_buf, &texture_pixels_buf,
                          &material_lut_buf, &lights_buf, &light_cdf_buf, &camera_buf, &packed_triangles_buf }) {
        destroy_buffer(ctx, *b);
    }
}

SceneBufferSet upload_scene_buffers(VulkanContext &ctx, const GPUScene &gs, const PackedScene &p) {
    SceneBufferSet s;
    s.positions_buf = create_buffer_with_data(ctx, p.packed_positions.data(), p.packed_positions.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.normals_buf   = create_buffer_with_data(ctx, p.packed_normals.data(), p.packed_normals.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.uvs_buf       = create_buffer_with_data(ctx, p.packed_uvs.data(), p.packed_uvs.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.indices_buf   = create_buffer_with_data(ctx, gs.indices.data(), gs.indices.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.bvh_nodes_buf = create_buffer_with_data(ctx, gs.bvh_nodes.data(), gs.bvh_nodes.size() * sizeof(GPUBVHNode), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.bvh_prims_buf = create_buffer_with_data(ctx, gs.bvh_primitives.data(), gs.bvh_primitives.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.tri_mat_buf   = create_buffer_with_data(ctx, gs.triangle_material_id.data(), gs.triangle_material_id.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.tri_light_buf = create_buffer_with_data(ctx, gs.triangle_light_id.data(), gs.triangle_light_id.size() * sizeof(uint32_t), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.materials_buf = create_buffer_with_data(ctx, p.gpu_materials.data(), p.gpu_materials.size() * sizeof(GPUMaterialGPU), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.textures_buf  = create_buffer_with_data(ctx, p.gpu_textures.data(), p.gpu_textures.size() * sizeof(GPUTextureGPU), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.texture_pixels_buf = create_buffer_with_data(ctx, p.texture_pixels.data(), p.texture_pixels.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.material_lut_buf   = create_buffer_with_data(ctx, p.material_lut.data(), p.material_lut.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.lights_buf    = create_buffer_with_data(ctx, p.gpu_lights.data(), p.gpu_lights.size() * sizeof(GPULightGPU), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    // Scene::build_gpu_scene() may return an empty CDF array for light-less scenes;
    // create_buffer_with_data already handles size_bytes==0 with a placeholder,
    // so no const_cast/push_back workaround is needed here (unlike the standalone tools).
    s.light_cdf_buf = create_buffer_with_data(ctx, gs.light_triangle_cdf.data(), gs.light_triangle_cdf.size() * sizeof(float), VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    s.camera_buf    = create_buffer_with_data(ctx, &gs.camera, sizeof(GPUCamera), VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT);
    s.packed_triangles_buf = create_buffer_with_data(ctx, p.packed_triangles.data(),
                                                     p.packed_triangles.size() * sizeof(float),
                                                     VK_BUFFER_USAGE_STORAGE_BUFFER_BIT);
    return s;
}

} // namespace tiny_renderer::gpu
