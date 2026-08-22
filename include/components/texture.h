#pragma once

#include <components/object.h>
#include <core/gpu_scene.h>

M_NAMESPACE_BEGIN
class Texture : public Object {
public:
    explicit Texture(bool spatial_varying) : m_spatial_varying(spatial_varying) {}

    ~Texture() override = default;

    void construct() override = 0;

    void add_child(const std::shared_ptr<Object> &child) override = 0;

    virtual Color3f eval(const SurfaceIntersection3f &si, bool active) = 0;

    virtual float eval_1(const SurfaceIntersection3f &si, bool active) = 0;

    virtual Vector2f eval_1_grad(const SurfaceIntersection3f &si, bool active) = 0;

    virtual Color3f mean() = 0;

    [[nodiscard]] EClassType get_class_type() const override { return ETexture; }

    [[nodiscard]] std::string to_string() const override = 0;

    [[nodiscard]] bool is_spatial_varying() const { return m_spatial_varying; }

    // Export a flattened, GPU-uploadable description of this texture (Stage 0
    // of the CPU -> GPU port). See core/gpu_scene.h for details. This is
    // purely additive: it does not affect eval()/eval_1()/mean() or any
    // other existing rendering behavior.
    [[nodiscard]] virtual GPUTexture to_gpu_texture(GPUSceneBuilder &builder) const = 0;

protected:
    bool m_spatial_varying;
};

class Constant;

M_NAMESPACE_END
