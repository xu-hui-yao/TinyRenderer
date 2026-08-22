#include <components/bsdf.h>
#include <components/emitter.h>
#include <components/texture.h>
#include <core/gpu_scene.h>

M_NAMESPACE_BEGIN

uint32_t GPUSceneBuilder::add_texture(const std::shared_ptr<Texture> &tex) {
    if (!tex) {
        return M_INVALID_INDEX;
    }

    auto it = texture_cache.find(tex.get());
    if (it != texture_cache.end()) {
        return it->second;
    }

    // Textures never nest other textures in this renderer, so there is no
    // risk of the recursive-append issue handled in add_material() below.
    GPUTexture gpu_tex = tex->to_gpu_texture(*this);

    auto index = static_cast<uint32_t>(result.textures.size());
    result.textures.push_back(std::move(gpu_tex));
    texture_cache[tex.get()] = index;
    return index;
}

uint32_t GPUSceneBuilder::add_material(const std::shared_ptr<BSDF> &bsdf) {
    if (!bsdf) {
        return M_INVALID_INDEX;
    }

    auto it = material_cache.find(bsdf.get());
    if (it != material_cache.end()) {
        return it->second;
    }

    // Reserve the slot and register it in the cache BEFORE recursing into
    // to_gpu_material(): compound materials (Mask / BumpMap / TwoSided) call
    // add_material() again for their nested BSDF(s), which may append
    // further entries to result.materials (and could in principle reallocate
    // the vector). We only ever refer back to our slot by integer index
    // (never by reference/pointer), so this ordering is safe either way.
    auto index = static_cast<uint32_t>(result.materials.size());
    result.materials.emplace_back();
    material_cache[bsdf.get()] = index;

    result.materials[index] = bsdf->to_gpu_material(*this);
    return index;
}

uint32_t GPUSceneBuilder::add_light(const std::shared_ptr<Emitter> &emitter, uint32_t mesh_id) {
    if (!emitter) {
        return M_INVALID_INDEX;
    }

    auto it = light_cache.find(emitter.get());
    if (it != light_cache.end()) {
        return it->second;
    }

    auto index = static_cast<uint32_t>(result.lights.size());
    result.lights.push_back(emitter->to_gpu_light(*this, mesh_id));
    light_cache[emitter.get()] = index;
    return index;
}

M_NAMESPACE_END
