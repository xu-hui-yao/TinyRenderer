#pragma once

M_NAMESPACE_BEGIN
template <typename Scalar_, DeviceType Device_> class TFrame {
public:
    static constexpr DeviceType Device = Device_;

    typedef Scalar_ Scalar;
    typedef TArray<Scalar, 3, ArrayType::Vector, Device> VectorType;
    typedef TArray<Scalar, 2, ArrayType::Vector, Device> VectorType2;
    typedef TArray<Scalar, 3, ArrayType::Point, Device> PointType;
    typedef TArray<Scalar, 3, ArrayType::Normal, Device> NormalType;

    VectorType s, t;
    NormalType n;

    // ============================== Constructor ==============================

    // Default constructor -- performs no initialization!
    M_HOST_DEVICE explicit TFrame()
        : n(NormalType(static_cast<Scalar>(0), static_cast<Scalar>(0),
                       static_cast<Scalar>(1))) {
        bool valid = true;
        this->n.self_norm(valid);
        coordinate_system(n, s, t, valid);
        this->s.self_norm(valid);
        this->t.self_norm(valid);
    }

    // Given a normal and tangent vectors, construct a new coordinate frame
    M_HOST_DEVICE TFrame(const VectorType &s, const VectorType &t,
                         const NormalType &n)
        : s(s), t(t), n(n) {
        bool valid = true;
        this->n.self_norm(valid);
        this->s.self_norm(valid);
        this->t.self_norm(valid);
    }

    // Construct a new coordinate frame from a single vector
    M_HOST_DEVICE explicit TFrame(const NormalType &n) : n(n) {
        bool valid = true;
        this->n.self_norm(valid);
        coordinate_system(this->n, s, t, valid);
        this->s.self_norm(valid);
        this->t.self_norm(valid);
    }

    // ============================== Check ==============================

    M_HOST_DEVICE static void check_zero(Scalar &scalar, bool &valid) {
        if (scalar == static_cast<Scalar>(0)) {
            scalar += static_cast<Scalar>(M_EPSILON);
            valid &= false;
        }
    }

    // ============================== Function ==============================
    M_HOST_DEVICE void coordinate_system(const NormalType &n, VectorType &b,
                                         VectorType &c, bool &valid) {
        Scalar sign  = n.z() >= 0.0f ? 1.0f : -1.0f;
        Scalar a_val = -1.0f / (sign + n.z());
        Scalar b_val = n.x() * n.y() * a_val;
        b            = VectorType({ mulsign(n.x() * n.x() * a_val, n.z()) + 1,
                                    mulsign(b_val, n.z()), mulsign(-n.x(), n.z()) });
        c = VectorType({ b_val, n.y() * n.y() * a_val + sign, -n.y() });
    }

    // Convert from world coordinates to local coordinates
    M_HOST_DEVICE [[nodiscard]] VectorType to_local(const VectorType &v) const {
        return { v.dot(s), v.dot(t), v.dot(n) };
    }

    // Convert from local coordinates to world coordinates
    M_HOST_DEVICE [[nodiscard]] VectorType to_world(const VectorType &v) const {
        return s * v.x() + t * v.y() + VectorType(n * v.z());
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the cosine of the angle between the normal and v */
    M_HOST_DEVICE static Scalar cos_theta(const VectorType &v, bool &valid) {
        return v.z();
    }

    M_HOST_DEVICE static Scalar cos_theta(const NormalType &v, bool &valid) {
        return v.z();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the angle between the normal and v */
    M_HOST_DEVICE static Scalar sin_theta(const VectorType &v, bool &valid) {
        return safe_sqrt(sin_theta2(v, valid));
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared sine of the angle between the normal and v */
    M_HOST_DEVICE static Scalar sin_theta2(const VectorType &v, bool &valid) {
        return v.x() * v.x() + v.y() * v.y();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the tangent of the angle between the normal and v */
    M_HOST_DEVICE static Scalar tan_theta(const VectorType &v, bool &valid) {
        Scalar temp = static_cast<Scalar>(1) - v.z() * v.z();
        return safe_sqrt(temp) / v.z();
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the tangent of the angle between the normal and v */
    M_HOST_DEVICE static Scalar tan_theta2(const VectorType &v, bool &valid) {
        Scalar temp = static_cast<Scalar>(1) - v.z() * v.z();
        return M_MAX(Scalar(0), temp) / (v.z() * v.z());
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the sine of the phi parameter in spherical coordinates */
    M_HOST_DEVICE static Scalar sin_phi(const VectorType &v, bool &valid) {
        Scalar sin_theta = TFrame::sin_theta(v, valid);
        check_zero(sin_theta, valid);
        return clamp(v.y() / sin_theta, static_cast<Scalar>(-1),
                     static_cast<Scalar>(1));
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the cosine of the phi parameter in spherical coordinates
     */
    M_HOST_DEVICE static Scalar cos_phi(const VectorType &v, bool &valid) {
        Scalar sin_theta = TFrame::sin_theta(v, valid);
        check_zero(sin_theta, valid);
        return clamp(v.x() / sin_theta, static_cast<Scalar>(-1),
                     static_cast<Scalar>(1));
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared sine of the phi parameter in  spherical
     * coordinates */
    M_HOST_DEVICE static Scalar sin_phi2(const VectorType &v, bool &valid) {
        auto norm_v        = v.norm(valid);
        auto sin_theta_2_v = sin_theta2(norm_v, valid);
        check_zero(sin_theta_2_v, valid);
        return clamp(norm_v.y() * norm_v.y() / sin_theta_2_v,
                     static_cast<Scalar>(0), static_cast<Scalar>(1));
    }

    /** \brief Assuming that the given direction is in the local coordinate
     * system, return the squared cosine of the phi parameter in  spherical
     * coordinates */
    M_HOST_DEVICE static Scalar cos_phi2(const VectorType &v, bool &valid) {
        auto norm_v        = v.norm(valid);
        auto sin_theta_2_v = sin_theta2(norm_v, valid);
        check_zero(sin_theta_2_v, valid);
        return clamp(norm_v.x() * norm_v.x() / sin_theta_2_v,
                     static_cast<Scalar>(0), static_cast<Scalar>(1));
    }

    /** \brief Give a unit direction, this function returns the sine and cosine
     * of the azimuth in a reference spherical coordinate system (see the \ref
     * Frame description)
     */
    static std::pair<Scalar, Scalar> sincos_phi(const VectorType &v) {
        bool valid           = true;
        Scalar sin_theta_2   = sin_theta2(v, valid);
        Scalar inv_sin_theta = Scalar(1) / sqrt(sin_theta_2);

        VectorType2 result = VectorType2({ v.x(), v.y() }) * inv_sin_theta;

        result = abs(sin_theta_2) <= Scalar(4) * M_EPSILON
                     ? VectorType2({ 1.f, 0.f })
                     : result.clamp(Scalar(-1), Scalar(1));

        return { result.y(), result.x() };
    }

    // Equality test
    M_HOST_DEVICE bool operator==(const TFrame &frame) const {
        return (frame.s == s).all() && (frame.t == t).all() &&
               (frame.n == n).all();
    }

    // Inequality test
    M_HOST_DEVICE bool operator!=(const TFrame &frame) const {
        return !operator==(frame);
    }

    // Return a human-readable string summary of this frame
    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "Frame[\n"
               "  s = "
            << s.to_string()
            << ",\n"
               "  t = "
            << t.to_string()
            << ",\n"
               "  n = "
            << n.to_string()
            << "\n"
               "]\n";
        return oss.str();
    }
};

M_NAMESPACE_END
