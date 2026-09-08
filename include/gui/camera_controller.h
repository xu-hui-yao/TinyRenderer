#pragma once

// ============================================================================
// CameraController - runtime camera manipulation for the interactive renderer.
//
// Deliberately does NOT touch the Camera class (include/components/camera.h):
// that class's camera_to_world is protected with no setter, and adding one
// would drag GLFW/input state into the offline renderer. Instead this class
// keeps its OWN eye/target/orientation and hands the caller two ready-to-upload
// row-major 4x4 matrices. The initial view is read once from the scene's
// Camera::to_gpu_camera(), so an interactive session opens exactly where the
// XML said the camera should be.
//
// Input-agnostic by design: it exposes orbit()/pan()/dolly()/move_local() and
// knows nothing about GLFW, so the window layer owns event handling and this
// stays trivially testable.
//
// ============================================================================
// MATRIX LAYOUT - subtle, and the reason this file defers to TTransform::look_at.
//
// GPUCamera stores camera_to_world row-major, i.e. index [r*4+c] is
// TMatrix::operator()(r, c), and the shaders multiply it as M * column_vector
// (see transform_point/transform_vector in common/scene_common.slang). So the
// camera's basis vectors live in the matrix COLUMNS:
//
//     col 0 = screen +X basis     col 2 = forward
//     col 1 = up                  col 3 = eye position
//
// Working through PerspectiveCamera's sample_to_camera gives
//
//     u = 0.5 - 0.5 * cot(fov/2) * (x_cam / z_cam)
//
// so u=0 (screen LEFT) implies x_cam > 0. Screen-left being +X while forward is
// +Z makes this a LEFT-handed frame, which is why TTransform::look_at names its
// first basis vector `left` and computes it as cross(up, forward) - the mirror
// of the right-handed cross(forward, up). Reusing that function is far safer
// than re-deriving the basis here: getting it backwards yields a mirrored
// image, which is a baffling bug to chase.
// ============================================================================

#include <core/array.h>
#include <core/transform.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>

namespace tiny_renderer::gui {

class CameraController {
public:
    enum class Mode { Orbit, Fly };

    struct Settings {
        float move_speed   = 1.0f;   // world units per second, scaled by scene radius
        float rotate_speed = 0.005f; // radians per pixel of mouse drag
        float zoom_speed   = 0.1f;   // fraction of distance per wheel notch
        float sensitivity  = 1.0f;   // global mouse multiplier
    };

    // A minimal float3. TArray only overloads arithmetic against scalars (plus
    // operator-= against another TArray), so doing the vector math with TArray
    // directly would need a pile of .x()/.y()/.z() noise; this keeps the
    // kinematics readable and dependency-free.
    struct Vec3 {
        float x = 0.0f, y = 0.0f, z = 0.0f;
        Vec3 operator+(const Vec3 &o) const { return { x + o.x, y + o.y, z + o.z }; }
        Vec3 operator-(const Vec3 &o) const { return { x - o.x, y - o.y, z - o.z }; }
        Vec3 operator*(float s) const { return { x * s, y * s, z * s }; }
        Vec3 &operator+=(const Vec3 &o) { x += o.x; y += o.y; z += o.z; return *this; }
    };
    static Vec3 cross(const Vec3 &a, const Vec3 &b) {
        return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
    }
    static Vec3 normalize(const Vec3 &v) {
        float l = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        return l > 1e-8f ? v * (1.0f / l) : Vec3{ 0.0f, 0.0f, 0.0f };
    }

    // Recovers the horizontal FOV (degrees) from a sample_to_camera matrix.
    //
    // PerspectiveCamera keeps `fov` private with no getter, and GPUCamera only
    // carries the baked matrix, so this is the only way for an interactive
    // session to start at the scene's actual FOV (otherwise the first touch of
    // the FOV slider would jump).
    //
    // Derivation: for the projection chain in PerspectiveCamera::construct(),
    // the clip-space x row works out to
    //     x_screen = -0.5*cot*x_cam + 0.5*z_cam,   w = z_cam
    // so inverting for x_cam gives
    //     x_cam = (w - 2*x_screen) / cot
    // hence sample_to_camera(0,0) = -2/cot, i.e. cot = -2 / s2c[0].
    static float fov_from_sample_to_camera(const float s2c[16]) {
        float cot = -2.0f / s2c[0];
        if (!(cot > 1e-6f)) return 30.0f; // degenerate/unknown: fall back to the renderer default
        return 2.0f * std::atan(1.0f / cot) * 180.0f / static_cast<float>(M_PI);
    }

    // `initial_camera_to_world`: the scene camera's row-major 4x4 (from
    // GPUCamera::camera_to_world). `fov` is HORIZONTAL, in degrees.
    // `scene_center`/`scene_radius` seed the orbit target, the default orbit
    // distance, and the pan/zoom step size; `radius <= 0` falls back to 1.
    CameraController(const float initial_camera_to_world[16], float fov, float near_clip, float far_clip,
                     const float scene_center[3], float scene_radius)
        : m_fov(fov), m_near_clip(near_clip), m_far_clip(far_clip) {

        std::copy_n(initial_camera_to_world, 16, m_initial_c2w);

        // Basis vectors are in the COLUMNS: eye = col 3, forward = col 2.
        Vec3 eye{ initial_camera_to_world[3], initial_camera_to_world[7], initial_camera_to_world[11] };
        Vec3 fwd = normalize(Vec3{ initial_camera_to_world[2], initial_camera_to_world[6],
                                   initial_camera_to_world[10] });
        if (fwd.x == 0.0f && fwd.y == 0.0f && fwd.z == 0.0f)
            fwd = { 0.0f, 0.0f, 1.0f };

        m_eye     = eye;
        m_yaw     = std::atan2(fwd.x, fwd.z);
        m_pitch   = std::asin(std::clamp(fwd.y, -1.0f, 1.0f));

        m_radius = scene_radius > 0.0f ? scene_radius : 1.0f;
        m_center = Vec3{ scene_center[0], scene_center[1], scene_center[2] };

        // Orbit distance defaults to the initial eye->center distance, which is
        // a sensible scale for this scene regardless of where the XML put the
        // camera.
        Vec3 d = m_eye - m_center;
        m_distance = std::max(std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z), m_radius * 0.1f);

        // The orbit target MUST lie on the initial view ray at exactly
        // `m_distance`, because Orbit mode recomputes eye = target -
        // forward*distance every rebuild. Anchoring it at the scene center
        // instead would silently move the camera on the very first frame, so
        // the interactive view would not match the scene XML (and resetting
        // would not restore it either).
        m_target = m_eye + fwd * m_distance;

        m_initial_eye    = m_eye;
        m_initial_yaw    = m_yaw;
        m_initial_pitch  = m_pitch;
        m_initial_target = m_target;

        rebuild();
    }

    // ---- Input handlers (called by the window layer) ----

    void orbit(float dx, float dy) {
        m_yaw -= dx * m_settings.rotate_speed * m_settings.sensitivity;
        m_pitch -= dy * m_settings.rotate_speed * m_settings.sensitivity;
        // Clamp just shy of the poles: at exactly +/-90 deg, forward is parallel
        // to world-up and cross(up, forward) degenerates, flipping the view.
        const float limit = static_cast<float>(M_PI_2) - 0.01f;
        m_pitch = std::clamp(m_pitch, -limit, limit);
        m_dirty = true;
    }

    void pan(float dx, float dy) {
        // Scaled by orbit distance so a drag feels the same at any zoom level.
        float scale = m_distance * 0.0015f * m_settings.sensitivity;
        // Dragging right should slide the scene right, i.e. move the camera
        // toward screen-LEFT, which is the +X basis (see the header note).
        Vec3 delta = left() * (dx * scale) + up() * (dy * scale);
        if (m_mode == Mode::Orbit)
            m_target += delta;
        else
            m_eye += delta;
        m_dirty = true;
    }

    void dolly(float notches) {
        if (m_mode == Mode::Orbit) {
            m_distance = std::clamp(m_distance * std::exp(-notches * m_settings.zoom_speed),
                                    m_radius * 0.01f, m_radius * 100.0f);
        } else {
            m_eye += forward() * (m_radius * m_settings.zoom_speed * notches * 4.0f);
        }
        m_dirty = true;
    }

    // Fly-mode WASD/QE translation along the camera's own axes.
    void move_local(float right_amount, float up_amount, float forward_amount, float dt) {
        float speed = m_settings.move_speed * m_radius * dt;
        m_eye += right() * (right_amount * speed) + up() * (up_amount * speed) +
                 forward() * (forward_amount * speed);
        m_dirty = true;
    }

    void set_mode(Mode mode) {
        if (mode == m_mode) return;
        // Entering Orbit adopts the point currently under the crosshair, so the
        // view doesn't jump when switching modes.
        if (mode == Mode::Orbit)
            m_target = m_eye + forward() * m_distance;
        m_mode  = mode;
        m_dirty = true;
    }
    [[nodiscard]] Mode mode() const { return m_mode; }

    // Marking dirty unconditionally would be a silent disaster: the render loop
    // calls set_fov()/set_aspect() EVERY frame with the panel's current value,
    // and treats `dirty` as "the view moved, drop the accumulation". Flagging a
    // no-op write reset the temporal history on every single frame, so the
    // interactive image was stuck at 1 spp forever.
    void set_fov(float fov_degrees) {
        const float f = std::clamp(fov_degrees, 5.0f, 170.0f);
        if (f == m_fov) return;
        m_fov   = f;
        m_dirty = true;
    }
    [[nodiscard]] float fov() const { return m_fov; }

    // Render resolution - feeds sample_to_camera's aspect handling only.
    void set_aspect(uint32_t width, uint32_t height) {
        const float a = height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;
        if (a == m_aspect) return;
        m_aspect = a;
        m_dirty  = true;
    }

    void reset() {
        m_eye     = m_initial_eye;
        m_yaw     = m_initial_yaw;
        m_pitch   = m_initial_pitch;
        m_target  = m_initial_target; // restores the initial view ray, see the constructor note
        m_mode    = Mode::Orbit;
        m_dirty   = true;
    }

    Settings &settings() { return m_settings; }

    // ---- Output ----

    void update() {
        if (m_dirty) rebuild();
    }

    [[nodiscard]] const float *camera_to_world() const { return m_c2w; }
    [[nodiscard]] const float *sample_to_camera() const { return m_s2c; }

    // True if the view changed since the last call. The caller uses this to
    // decide whether to reset temporal accumulation history.
    [[nodiscard]] bool consume_changed() {
        bool c = m_dirty;
        m_dirty = false;
        return c;
    }

    // Copy-pasteable <transform name="to_world"> matrix for pasting back into a
    // scene XML, which is how an interactively-found view gets saved.
    [[nodiscard]] std::string to_string() const {
        char buf[512];
        std::snprintf(buf, sizeof(buf),
                      "eye=(%.4f, %.4f, %.4f) yaw=%.1f pitch=%.1f dist=%.4f\nmatrix=%g %g %g %g %g %g %g %g %g %g %g "
                      "%g %g %g %g %g",
                      m_eye.x, m_eye.y, m_eye.z, m_yaw * 180.0f / static_cast<float>(M_PI),
                      m_pitch * 180.0f / static_cast<float>(M_PI), m_distance,
                      m_c2w[0], m_c2w[1], m_c2w[2], m_c2w[3], m_c2w[4], m_c2w[5], m_c2w[6], m_c2w[7],
                      m_c2w[8], m_c2w[9], m_c2w[10], m_c2w[11], m_c2w[12], m_c2w[13], m_c2w[14], m_c2w[15]);
        return std::string(buf);
    }

private:
    [[nodiscard]] Vec3 forward() const {
        float cp = std::cos(m_pitch);
        return { cp * std::sin(m_yaw), std::sin(m_pitch), cp * std::cos(m_yaw) };
    }
    // Screen +X basis. TTransform::look_at computes it as cross(up, forward).
    [[nodiscard]] Vec3 left() const {
        Vec3 l = cross(Vec3{ 0.0f, 1.0f, 0.0f }, forward());
        float len = std::sqrt(l.x * l.x + l.y * l.y + l.z * l.z);
        return len > 1e-5f ? l * (1.0f / len) : Vec3{ 1.0f, 0.0f, 0.0f };
    }
    [[nodiscard]] Vec3 right() const { return left() * -1.0f; }
    [[nodiscard]] Vec3 up() const { return cross(forward(), left()); }

    void rebuild() {
        Vec3 f = forward();

        if (m_mode == Mode::Orbit)
            m_eye = m_target - f * m_distance;

        // Defer to TTransform::look_at for the basis construction - it already
        // encodes this renderer's left-handed convention (see header note).
        // `valid` is a STICKY flag: TArray/TMatrix's check_zero() only ever does
        // `valid &= false`, it never sets it true. Every helper that takes it
        // (norm(), look_at(), ...) therefore expects the caller to seed it with
        // `true` and reads it as "nothing degenerate was encountered".
        //
        // Seeding it with `false` made look_at() bail on its very first
        // `if (!valid) return TTransform();` and hand back the IDENTITY
        // transform - a camera sitting at the world origin looking down +Z,
        // with the scene's actual view discarded. The offline path was
        // unaffected (it uses the camera baked by the XML), which is exactly
        // why only the interactive view rendered the wrong thing.
        bool valid = true;
        Transform4f c2w = Transform4f::look_at(Vector3f(m_eye.x, m_eye.y, m_eye.z),
                                               Vector3f(m_eye.x + f.x, m_eye.y + f.y, m_eye.z + f.z),
                                               Vector3f(0.0f, 1.0f, 0.0f), valid);
        if (!valid) {
            // Degenerate basis (forward parallel to world up, or a zero-length
            // view vector). Keep the last good view matrix rather than silently
            // teleporting the camera to the origin - the projection below is
            // independent of it and is still refreshed.
        } else {
            const Matrix4f M = c2w.get_transform();
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    m_c2w[r * 4 + c] = M(r, c);
        }

        // sample_to_camera, ported verbatim from PerspectiveCamera::construct()
        // (src/sensors/perspective.cpp) so the interactive projection matches
        // the offline renderer's exactly. `fov` is horizontal; the aspect
        // correction lives in the scale/translate pair.
        float recip  = 1.0f / (m_far_clip - m_near_clip);
        float cot    = 1.0f / std::tan(m_fov * 0.5f * static_cast<float>(M_PI) / 180.0f);
        float aspect = m_aspect;

        Matrix4f perspective  = Matrix4f::zero();
        perspective(0, 0)     = cot;
        perspective(1, 1)     = cot;
        perspective(2, 2)     = m_far_clip * recip;
        perspective(2, 3)     = -m_near_clip * m_far_clip * recip;
        perspective(3, 2)     = 1.0f;

        Transform4f s2c = Transform4f(Matrix4f::scale(Vector4f(-0.5f, -0.5f * aspect, 1.0f, 1.0f)) *
                                      Matrix4f::translate(Vector4f(-1.0f, -1.0f / aspect, 0.0f, 0.0f)) *
                                      perspective)
                              .inverse();
        const Matrix4f S = s2c.get_transform();
        for (int r = 0; r < 4; ++r)
            for (int c = 0; c < 4; ++c)
                m_s2c[r * 4 + c] = S(r, c);
    }

    Settings m_settings;
    Mode m_mode = Mode::Orbit;

    float m_fov = 30.0f, m_near_clip = 1e-4f, m_far_clip = 1e4f;
    float m_aspect = 16.0f / 9.0f;

    float m_yaw = 0.0f, m_pitch = 0.0f, m_distance = 1.0f;
    Vec3 m_eye{ 0.0f, 0.0f, 0.0f };
    Vec3 m_target{ 0.0f, 0.0f, 0.0f };
    Vec3 m_center{ 0.0f, 0.0f, 0.0f };
    float m_radius = 1.0f;

    float m_initial_c2w[16] = {};
    Vec3 m_initial_eye{ 0.0f, 0.0f, 0.0f };
    Vec3 m_initial_target{ 0.0f, 0.0f, 0.0f };
    float m_initial_yaw = 0.0f, m_initial_pitch = 0.0f;

    float m_c2w[16] = {};
    float m_s2c[16] = {};
    bool m_dirty = true;
};

} // namespace tiny_renderer::gui
