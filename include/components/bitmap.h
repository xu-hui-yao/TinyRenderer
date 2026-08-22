#pragma once

#include <components/texture.h>
#include <core/tensor.h>
#include <fstream>

M_NAMESPACE_BEGIN

// Output tonemapping applied to LINEAR radiance before the sRGB encoding
// step in Bitmap::save_png(). `None` (the default, preserving all prior
// behavior byte-for-byte) hard-clips values above 1.0 - correct/necessary
// for numerical ground-truth comparisons (e.g. CPU vs GPU renderer
// verification), but any physically-plausible HDR scene with strong light
// sources will show blown-out highlights. `ACES` instead compresses
// highlights with a smooth filmic shoulder (see
// TRGBSpectrum::aces_filmic()), trading exact linear-radiance recoverability
// for a more film-like, non-clipped look - a purely aesthetic choice, opt-in
// via `--tonemap=aces` on the command line (see src/main/main.cpp).
enum class ToneMapMode { None, ACES };

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

    void save_png(const std::string &filename, ToneMapMode tonemap = ToneMapMode::None) const;

    void load_exr(const std::string &filename);

    // Loads an 8-bit LDR image (PNG/JPEG/BMP/TGA). If `m_raw` is set (see
    // the `raw` XML property below, which defaults to true - this
    // project's XML scenes/assets already store LDR bitmap textures as
    // linear data, not sRGB-encoded), the normalized [0,1] values are
    // stored verbatim. Otherwise (raw explicitly set to false, for a
    // texture that IS actually sRGB-encoded) the loaded values are decoded
    // to linear light via Color3f::from_srgb() before being stored -
    // matching what the sRGB encoding step on the OUTPUT side
    // (Bitmap::save_png(), via to_srgb()) does in reverse.
    void load_image(const std::string &filename);

    void load_hdr(const std::string &filename);

    [[nodiscard]] int get_rows() const;

    [[nodiscard]] int get_cols() const;

    Color3f eval(const SurfaceIntersection3f &si, bool active) override;

    float eval_1(const SurfaceIntersection3f &si, bool active) override;

    Vector2f eval_1_grad(const SurfaceIntersection3f &si, bool active) override;

    Color3f mean() override;

    [[nodiscard]] GPUTexture to_gpu_texture(GPUSceneBuilder &builder) const override;

    [[nodiscard]] std::string to_string() const override;

private:
    std::shared_ptr<TensorXf> m_data = nullptr;

    // When true, load_image() skips the sRGB->linear decode: the image's
    // raw normalized [0,1] values ARE the data, verbatim. Mirrors Mitsuba's
    // <boolean name="raw"> texture property, but defaults to true here
    // since this project's XML scenes/assets already store LDR bitmap
    // textures (color or data, e.g. roughness/alpha/opacity/bump) as
    // linear data, not sRGB-encoded. Set <boolean name="raw"
    // value="false"/> to opt into the sRGB->linear decode for a texture
    // that IS actually sRGB-encoded. Irrelevant for load_exr()/load_hdr():
    // those formats are always already linear and are never decoded
    // regardless of this flag.
    bool m_raw = true;
};


M_NAMESPACE_END
