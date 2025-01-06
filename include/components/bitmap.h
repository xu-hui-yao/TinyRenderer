#pragma once

#include <components/texture.h>
#include <core/tensor.h>
#include <fstream>

M_NAMESPACE_BEGIN
class Bitmap : public Texture {
public:
    explicit Bitmap(const PropertyList &properties);

    Bitmap(int height, int width, int channels);

    ~Bitmap() override;

    void add_child(const std::shared_ptr<Object> &child) override;

    void construct() override;

    [[nodiscard]] std::shared_ptr<TensorXf> get_data();

    float operator()(int row, int col, int channel) const;

    float &operator()(int row, int col, int channel);

    void save_exr(const std::string &filename) const;

    void save_png(const std::string &filename) const;

    void load_exr(const std::string &filename);

    void load_image(const std::string &filename);

    void load_hdr(const std::string &filename);

    [[nodiscard]] int get_rows() const;

    [[nodiscard]] int get_cols() const;

    Color3f eval(const SurfaceIntersection3f &si, bool &active) override;

    float eval_1(const SurfaceIntersection3f &si, bool &active) override;

    Vector2f eval_1_grad(const SurfaceIntersection3f &si, bool &active) override;

    Color3f mean() override;

    [[nodiscard]] std::string to_string() const override;

private:
    std::shared_ptr<TensorXf> m_data = nullptr;
};

M_NAMESPACE_END
