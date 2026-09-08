#pragma once

// ============================================================================
// InteractiveSession - a PERSISTENT Vulkan session for real-time rendering.
//
// The offline renderers (gpu_renderer.cpp's render_megakernel / render_wavefront)
// create a VulkanContext, re-upload the whole scene, build pipelines, render
// every requested sample, and tear it all down again inside a single call. That
// is fine for a batch render but makes per-frame invocation impossible.
//
// This class hoists all of that one-time setup into the constructor so that
// render_frame() is cheap enough to call 30+ times a second:
//
//     InteractiveSession session(scene->build_gpu_scene(), 1280, 720);
//     while (!window.should_close()) {
//         camera_controller.update();
//         const uint8_t* rgba = session.render_frame(render, display);
//         upload_and_present(rgba);
//     }
//
// Per-frame GPU work in P1 is two dispatches in one command buffer:
//   1. rt_path_trace  - spp_per_frame samples, accumulated into accum_buffer
//   2. svgf_modulate  - divide by sample count, tonemap, pack to RGBA8
//
// Only the RGBA8 display buffer is ever read by the host, and it is read as a
// raw persistently-mapped pointer already in glTexSubImage2D's layout, so the
// host does zero per-pixel work per frame (the offline preview path's CPU
// convert_to_ldr() is what would otherwise blow the frame budget).
//
// P2/P3 add the denoiser passes between these two dispatches.
// ============================================================================

// core/gpu_scene.h itself only pulls in core/common.h, which FORWARD-DECLARES
// TArray/TRGBSpectrum/TTransform and relies on the includer having brought in
// the definitions. Same fix as in gpu_scene_upload.h: name what we need.
#include <core/array.h>
#include <core/spectrum.h>
#include <core/transform.h>

#include <array>
#include <core/gpu_scene.h>
#include <cstdint>
#include <gpu/gpu_buffer.h>
#include <gpu/gpu_image.h>
#include <gpu/gpu_timings.h>
#include <gpu/gpu_scene_upload.h>
#include <gpu/vk_context.h>
#include <vector>

namespace tiny_renderer::gpu {

class InteractiveSession {
public:
    // Everything that may change from frame to frame without forcing the
    // scene buffers to be rebuilt.
    struct RenderParams {
        // Row-major 4x4s (16 floats each), matching GPUCamera's layout, which
        // is byte-identical to the shader's GPUCameraGPU constant buffer.
        const float *camera_to_world  = nullptr;
        const float *sample_to_camera = nullptr;

        uint32_t max_depth = 3; // 1 + indirect bounce + NEE; the interactive default
        uint32_t rr_depth  = 2;
        uint32_t spp_per_frame = 1;
        float radiance_clamp   = 10.0f; // 0 disables clamping

        // Divide out the surface albedo before filtering and multiply it back
        // at display time, so texture detail survives denoising. Turning this
        // off writes raw radiance and forces albedo to 1.0, which makes the
        // whole pipeline algebraically transparent - that is what lets
        // interactive_verify_main.cpp assert bit-exact agreement with the
        // offline megakernel.
        bool demodulate = true;

        // Set when the camera or any render setting changed: the running
        // average is stale for the new view and must be restarted.
        bool reset_history = false;
    };

    // Spatiotemporal variance-guided filtering (SVGF) parameters. These do not
    // invalidate history when changed - they only affect how the next frame
    // blends in - EXCEPT that a big change will look wrong for a frame or two
    // until the EMA catches up.
    struct DenoiseParams {
        bool enabled        = true;
        bool clamp_history  = true; // rectify history into the current 3x3 range (kills ghosting)
        // How much wider the rectification radius gets per accumulated sample.
        // 0 keeps it fixed at 1.5 sigma forever, which biases a static view and
        // stalls convergence; larger values let a long history converge at the
        // cost of weaker ghost removal. See svgf_temporal.slang.
        float clamp_growth  = 0.25f;
        float alpha_color   = 0.05f;   // minimum blend weight for new colour samples
        float alpha_moments = 0.05f;   // ...and for the luminance moments
        float phi_depth     = 1.0f;    // depth-rejection tolerance, in units of the depth gradient
        float phi_normal    = 0.5f;    // normal-rejection threshold, in RADIANS
        float max_hist_len  = 0.0f;    // cap on history length; 0 = uncapped (keeps converging)

        // ---- Spatial (a-trous) pass ----
        // 0 disables it and leaves the temporal result unsmoothed - useful both
        // as a baseline and as a debugging step, since it isolates how much of
        // the cleanup came from each pass.
        //
        // Each iteration keeps a fixed 5x5 stencil and doubles the stride, so N
        // iterations reach a (4*2^(N-1)+1)^2 neighbourhood for 25*N taps.
        // Beyond ~5 iterations the gathers are so wide they start over-blurring
        // for little noise reduction.
        uint32_t atrous_iterations = 4;
        float phi_color    = 4.0f;   // luminance edge-stopping scale (sigma_l)
        float phi_normal_a = 128.0f; // normal edge-stopping EXPONENT (sigma_n)
        float phi_depth_a  = 1.0f;   // depth edge-stopping scale (sigma_z)

        // Writes a status code instead of illumination, for diagnosing why
        // history is (not) being reused. See svgf_temporal.slang.
        bool debug_mode = false;
    };

    // Display-side settings. Note these live in the modulate pass, so changing
    // them NEVER invalidates the accumulation (no reset_history needed) -
    // dragging exposure at 1000 accumulated spp costs one cheap dispatch.
    struct DisplayParams {
        float exposure   = 1.0f;
        uint32_t tonemap = 2;   // 0 = clamp, 1 = Reinhard, 2 = ACES
        bool use_srgb    = true; // matches Bitmap::save_png()'s transfer function
        float gamma      = 2.2f; // only used when use_srgb == false

        // 0 = final image; otherwise one of the DEBUG_* views in
        // svgf_modulate.slang. Essential for tuning the filter.
        uint32_t debug_view = 0;
        float variance_scale = 1.0f; // display gain for the variance view
    };

    // `gs` is taken by value because its camera's width/height are overwritten
    // with the render resolution before upload (the shaders bounds-check
    // against camera.width/height).
    InteractiveSession(GPUScene gs, uint32_t width, uint32_t height);
    ~InteractiveSession();

    InteractiveSession(const InteractiveSession &)            = delete;
    InteractiveSession &operator=(const InteractiveSession &) = delete;

    // Dispatches one frame and returns a pointer to `width*height*4` bytes of
    // RGBA8 pixel data. The pointer stays valid until the next render_frame()
    // call (it is the display buffer's persistent mapping).
    const uint8_t *render_frame(const RenderParams &render, const DenoiseParams &denoise,
                                const DisplayParams &display);

    [[nodiscard]] uint32_t width() const { return m_width; }
    [[nodiscard]] uint32_t height() const { return m_height; }

    // Samples accumulated so far for the current view. Resets to
    // spp_per_frame whenever reset_history is passed.
    [[nodiscard]] uint32_t accumulated_spp() const { return m_accumulated_spp; }

    // Reads back the DENOISED (temporally accumulated, still demodulated)
    // illumination as one RGB triple per pixel. This is what the display pass
    // shows before remodulation, so with demodulate disabled it is directly
    // comparable to an offline render of the same sample count.
    //
    // Non-const: it issues a transfer to read the image back to the host.
    [[nodiscard]] std::vector<float> read_filtered_radiance();

    // Reads back the current accumulation as LINEAR radiance, one RGB triple
    // per pixel (row-major), i.e. the sum divided by the sample count.
    //
    // This is the correct source for anything that must match the offline
    // renderer's output (notably Bitmap::save_png(), which applies its own
    // tonemap + sRGB transfer). The RGBA8 pointer returned by render_frame()
    // has ALREADY been tonemapped by the modulate pass, so saving that instead
    // would double-tonemap - a mistake this codebase has hit before (see the
    // warning in include/gpu/gpu_renderer.h).
    [[nodiscard]] std::vector<float> read_linear_radiance() const;

    // The (possibly resolution-overridden) scene, for callers that need e.g.
    // the light count or the camera's original matrices.
    [[nodiscard]] const GPUScene &gpu_scene() const { return m_gs; }

    // Per-pass GPU timings for the most recent frame. See PassTimings.
    [[nodiscard]] const PassTimings &timings() const { return m_timings; }

private:
    struct Pipeline {
        VkDescriptorSetLayout set_layout = VK_NULL_HANDLE;
        VkPipelineLayout pipeline_layout = VK_NULL_HANDLE;
        VkPipeline pipeline              = VK_NULL_HANDLE;
        VkShaderModule module            = VK_NULL_HANDLE;
        // More than one set when a pass ping-pongs between two buffer sets and
        // needs the alternation baked into the bindings (temporal history,
        // a-trous iterations) rather than branched on in the shader.
        std::vector<VkDescriptorSet> sets;
    };

    // A descriptor slot that is either a storage/uniform buffer or a storage
    // image. The denoiser's screen-sized buffers are images (see
    // include/gpu/gpu_image.h for why), so pipeline setup has to handle both.
    struct DescriptorSlot {
        VkDescriptorType type       = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        VkBuffer buffer             = VK_NULL_HANDLE; // set when type is a buffer type
        VkImageView image_view      = VK_NULL_HANDLE; // set when type is STORAGE_IMAGE
    };

    // `slot_sets` is one entry per descriptor set to allocate.
    void create_pipeline(Pipeline &out, const std::string &spv_name, size_t push_constants_size,
                         const std::vector<VkDescriptorSetLayoutBinding> &bindings,
                         const std::vector<std::vector<DescriptorSlot>> &slot_sets);
    void destroy_pipeline(Pipeline &p);
    // Reads back the timestamp queries written during the last render_frame().
    void read_pass_timings();

    GPUScene m_gs;
    uint32_t m_width = 0, m_height = 0;

    VulkanContext m_ctx;
    PackedScene m_packed;
    SceneBufferSet m_scene_buffers;

    GpuBuffer m_accum_buf;   // float4 per pixel: linear radiance sum + sample count
    GpuBuffer m_display_buf; // uint per pixel: packed RGBA8

    // ---- G-Buffer (written by rt_path_trace) ----
    GpuImage m_gbuffer_albedo;       // RGBA16F
    GpuImage m_gbuffer_normal_depth; // RGBA32F: (oct normal.xy, linear_z, depth_grad)
    GpuImage m_gbuffer_pos_id;       // RGBA32F: (world pos.xyz, asfloat(material_id))

    // ---- Illumination + moments ----
    GpuImage m_illum_cur;   // RGBA16F, this frame's demodulated illumination
    GpuImage m_moments_cur; // RG32F, this frame's (lum, lum^2)

    // Ping-pong history. Index [0] and [1] alternate: the temporal pass reads
    // one and writes the other, so no copy is ever needed between frames.
    GpuImage m_illum_hist[2];    // RGBA16F
    GpuImage m_moments_hist[2];  // RG32F
    GpuImage m_hist_len[2];      // R32F

    // Ping-pong scratch for the a-trous iterations. Separate from the history
    // pair because the temporal output must be preserved as next frame's input;
    // overwriting it with filtered data would feed blurred results back into the
    // accumulator and compound the blur every frame.
    GpuImage m_atrous_illum[2];  // RGBA16F

    // The PREVIOUS frame's camera, for the temporal pass's reprojection:
    // world_to_camera then camera_to_sample (32 floats) + width/height.
    GpuBuffer m_prev_camera_buf;

    Pipeline m_path_trace;
    Pipeline m_temporal;
    Pipeline m_atrous;
    Pipeline m_modulate;

    // Which image the display pass shows this frame, so read_filtered_radiance()
    // returns what the user is actually looking at rather than an intermediate.
    const GpuImage *m_display_source = nullptr;

    VkDescriptorPool m_descriptor_pool = VK_NULL_HANDLE;

    // ---- GPU timestamp queries ----
    VkQueryPool m_query_pool      = VK_NULL_HANDLE;
    float m_timestamp_period_ns   = 0.0f; // 0 => unavailable
    PassTimings m_timings;
    // Number of timestamps written per frame: one before the first pass, then
    // one after each of the four.
    static constexpr uint32_t kTimestampCount = 5;

    uint32_t m_frame_index     = 0;
    uint32_t m_accumulated_spp = 0;
    // Which ping-pong pair currently holds the VALID history (0 or 1).
    uint32_t m_history_index = 0;
    // Whether the history buffers hold anything usable yet.
    bool m_history_valid = false;

    // Frames since the last history reset, i.e. how long the view has been
    // stable. Used to relax the blend floor so a static view keeps converging
    // instead of freezing into a fixed EMA - see render_frame().
    uint32_t m_still_frames = 0;

    // The previous frame's camera, in the form the temporal pass reprojects
    // with: world_to_camera (inverse of the camera_to_world we were handed) and
    // camera_to_sample. Keeping them here avoids re-inverting every frame for
    // the (common) case where the camera has not moved.
    std::array<float, 16> m_prev_w2c{};
    std::array<float, 16> m_prev_c2s{};
};

} // namespace tiny_renderer::gpu
