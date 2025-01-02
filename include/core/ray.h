#pragma once

#include <core/array.h>

M_NAMESPACE_BEGIN
template <typename Scalar, int Dimension_> class TRay {
public:
    static constexpr int Dimension = Dimension_;
    typedef TArray<Scalar, Dimension, ArrayType::Point> Point;
    typedef TArray<Scalar, Dimension, ArrayType::Vector> Vector;

    // ============================== Constructor ==============================

    TRay() : m_min_t(M_EPSILON), m_max_t(M_MAX_FLOAT) {
        bool valid  = true;
        m_direction = Vector(static_cast<Scalar>(1), static_cast<Scalar>(1), static_cast<Scalar>(1)).norm(valid);
    }

    TRay(const Point &origin, const Vector &direction, bool &valid)
        : m_origin(origin), m_min_t(M_EPSILON), m_max_t(M_MAX_FLOAT) {
        this->m_direction = direction.norm(valid);
        update(valid);
    }

    TRay(const Point &origin, const Vector &direction, Scalar min_t, Scalar max_t, bool &valid)
        : m_origin(origin), m_min_t(min_t), m_max_t(max_t) {
        this->m_direction = direction.norm(valid);
        update(valid);
    }

    TRay(const TRay &ray)
        : m_origin(ray.m_origin), m_direction(ray.m_direction), m_direction_reciprocal(ray.m_direction_reciprocal),
          m_min_t(ray.m_min_t), m_max_t(ray.m_max_t) {}

    TRay(const TRay &ray, Scalar min_t, Scalar max_t)
        : m_origin(ray.m_origin), m_direction(ray.m_direction), m_direction_reciprocal(ray.m_direction_reciprocal),
          m_min_t(min_t), m_max_t(max_t) {}

    void update(bool &valid) { m_direction_reciprocal = m_direction.wise_inverse(valid); }

    // ============================== Getter and setter
    // ==============================

    Point o() const { return m_origin; }

    Point &o() { return m_origin; }

    Vector d() const { return m_direction; }

    Vector &d() { return m_direction; }

    Vector d_reciprocal() const { return m_direction_reciprocal; }

    Vector &d_reciprocal() { return m_direction_reciprocal; }

    Scalar min_t() const { return m_min_t; }

    Scalar &min_t() { return m_min_t; }

    Scalar max_t() const { return m_max_t; }

    Scalar &max_t() { return m_max_t; }

    // ============================== Function ==============================

    Point operator()(Scalar t) const { return m_origin + m_direction * t; }

    TRay reverse() const {
        TRay result;
        result.m_origin               = m_origin;
        result.m_direction            = -m_direction;
        result.m_direction_reciprocal = -m_direction_reciprocal;
        result.m_min_t                = m_min_t;
        result.m_max_t                = m_max_t;
        return result;
    }

    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "Ray[\n";
        oss << "  Origin: " << m_origin.to_string() << "\n";
        oss << "  Direction: " << m_direction.to_string() << "\n";
        oss << "  Min: " << m_min_t << "\n";
        oss << "  Max: " << m_max_t << "\n";
        oss << "]\n";
        return oss.str();
    }

private:
    Point m_origin;
    Vector m_direction;
    Vector m_direction_reciprocal;
    Scalar m_min_t;
    Scalar m_max_t;
};

M_NAMESPACE_END
