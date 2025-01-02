#pragma once

M_NAMESPACE_BEGIN
template <typename Scalar, int Dimension_>
class TBoundingBox {
    static constexpr int Dimension     = Dimension_;
    typedef TArray<Scalar, Dimension, ArrayType::Point> Point;
    typedef TArray<Scalar, Dimension, ArrayType::Vector> Vector;

public:
    // ============================== Constructor ==============================

    TBoundingBox() {
        m_min_p.set_constant(M_MAX_FLOAT);
        m_max_p.set_constant(-M_MAX_FLOAT);
    }

    explicit TBoundingBox(const Point &p)
        : m_min_p(p), m_max_p(p) {}

    TBoundingBox(const Point &min, const Point &max)
        : m_min_p(min), m_max_p(max) {}

    // ============================== Getter and setter
    // ==============================

    Point &get_min() { return m_min_p; }

    Point get_min() const { return m_min_p; }

    Point &get_max() { return m_max_p; }

    Point get_max() const { return m_max_p; }

    // ============================== Function ==============================

    bool operator==(const TBoundingBox &bounding_box) const {
        return m_min_p == bounding_box.m_min_p &&
               m_max_p == bounding_box.m_max_p;
    }

    bool operator!=(const TBoundingBox &bbox) const {
        return m_min_p != bbox.m_min_p || m_max_p != bbox.m_max_p;
    }

    Scalar get_volume() const {
        return (m_max_p - m_min_p).prod();
    }

    [[nodiscard]] Scalar get_surface_area() const {
        Vector d    = m_max_p - m_min_p;
        auto result = static_cast<Scalar>(0);
        for (int i = 0; i < Dimension; ++i) {
            auto term = static_cast<Scalar>(1);
            for (int j = 0; j < Dimension; ++j) {
                if (i == j)
                    continue;
                term *= d(j);
            }
            result += term;
        }
        return static_cast<Scalar>(2) * result;
    }

    Point get_center() const {
        return (m_max_p + m_min_p) * static_cast<Scalar>(0.5);
    }

    template <bool Strict = false>
    bool contains(const Point &p) const {
        if constexpr (Strict) {
            return (p > m_min_p).all() && (p < m_max_p).all();
        } else {
            return (p >= m_min_p).all() && (p <= m_max_p).all();
        }
    }

    template <bool Strict = false>
    bool contains(const TBoundingBox &bbox) const {
        if constexpr (Strict) {
            return (bbox.m_min_p > m_min_p).all() &&
                   (bbox.m_max_p < m_max_p).all();
        } else {
            return (bbox.m_min_p >= m_min_p).all() &&
                   (bbox.m_max_p <= m_max_p).all();
        }
    }

    template <bool Strict = false>
    bool overlaps(const TBoundingBox &bbox) const {
        if constexpr (Strict) {
            return (bbox.m_min_p < m_max_p).all() &&
                   (bbox.m_max_p > m_min_p).all();
        } else {
            return (bbox.m_min_p <= m_max_p).all() &&
                   (bbox.m_max_p >= m_min_p).all();
        }
    }

    Scalar squared_distance_to(const Point &p) const {
        Scalar result = 0;

        for (int i = 0; i < Dimension; ++i) {
            Scalar value = 0;
            if (p(i) < m_min_p(i))
                value = m_min_p(i) - p(i);
            else if (p(i) > m_max_p(i))
                value = p(i) - m_max_p(i);
            result += value * value;
        }

        return result;
    }

    Scalar distance_to(const Point &p) const {
        return sqrt(squared_distance_to(p));
    }

    Scalar squared_distance_to(const TBoundingBox &bbox) const {
        Scalar result = 0;

        for (int i = 0; i < Dimension; ++i) {
            Scalar value = 0;
            if (bbox.m_max_p(i) < m_min_p(i))
                value = m_min_p(i) - bbox.m_max_p(i);
            else if (bbox.m_min_p(i) > m_max_p(i))
                value = bbox.m_min_p(i) - m_max_p(i);
            result += value * value;
        }

        return result;
    }

    Scalar distance_to(const TBoundingBox &bbox) const {
        return sqrt(squared_distance_to(bbox));
    }

    [[nodiscard]] bool is_valid() const {
        return (m_max_p >= m_min_p).all();
    }

    [[nodiscard]] bool is_point() const {
        return (m_max_p == m_min_p).all();
    }

    [[nodiscard]] bool has_volume() const {
        return (m_max_p > m_min_p).all();
    }

    [[nodiscard]] int get_major_axis() const {
        Vector d    = m_max_p - m_min_p;
        int largest = 0;
        for (int i = 1; i < Dimension; ++i)
            if (d(i) > d(largest))
                largest = i;
        return largest;
    }

    [[nodiscard]] int get_minor_axis() const {
        Vector d     = m_max_p - m_min_p;
        int shortest = 0;
        for (int i = 1; i < Dimension; ++i)
            if (d(i) < d(shortest))
                shortest = i;
        return shortest;
    }

    Vector get_extents() const { return m_max_p - m_min_p; }

    void clip(const TBoundingBox &bounding_box) {
        m_min_p = m_min_p.wise_max(bounding_box.m_min_p);
        m_max_p = m_max_p.wise_min(bounding_box.m_max_p);
    }

    void expand_by(const Point &p) {
        m_min_p = m_min_p.wise_min(p);
        m_max_p = m_max_p.wise_max(p);
    }

    void expand_by(const TBoundingBox &bbox) {
        m_min_p = m_min_p.wise_min(bbox.m_min_p);
        m_max_p = m_max_p.wise_max(bbox.m_max_p);
    }

    static TBoundingBox merge(const TBoundingBox &bbox1,
                                            const TBoundingBox &bbox2) {
        return TBoundingBox(bbox1.m_min_p.wise_min(bbox2.m_min_p),
                            bbox1.m_max_p.wise_max(bbox2.m_max_p));
    }

    Point get_corner(int index) const {
        Point result;
        for (int i = 0; i < Dimension; ++i)
            result(i) = index & 1 << i ? m_max_p(i) : m_min_p(i);
        return result;
    }

    [[nodiscard]] bool ray_intersect(const Ray3f &ray) const {
        Scalar near_t = -M_MAX_FLOAT;
        Scalar far_t  = M_MAX_FLOAT;

        for (int i = 0; i < 3; i++) {
            Scalar origin  = ray.o()(i);
            Scalar min_val = m_min_p(i), max_val = m_max_p(i);

            Scalar dir = ray.d()(i);
            if (std::abs(dir) < M_EPSILON) {
                if (origin < min_val || origin > max_val)
                    return false;
            } else {
                Scalar inv_dir = ray.d_reciprocal()(i);
                Scalar t1      = (min_val - origin) * inv_dir;
                Scalar t2      = (max_val - origin) * inv_dir;

                near_t = M_MAX(M_MIN(t1, t2), near_t);
                far_t  = M_MIN(M_MAX(t1, t2), far_t);

                if (near_t > far_t + M_EPSILON)
                    return false;
            }
        }

        return ray.min_t() <= far_t && near_t <= ray.max_t();
    }

    bool ray_intersect(const Ray3f &ray, Scalar &near_t,
                                     Scalar &far_t) const {
        near_t = -M_MAX_FLOAT;
        far_t  = M_MAX_FLOAT;

        for (int i = 0; i < 3; i++) {
            Scalar origin  = ray.o()(i);
            Scalar min_val = m_min_p(i), max_val = m_max_p(i);

            if (ray.d()(i) == 0) {
                if (origin < min_val || origin > max_val)
                    return false;
            } else {
                Scalar t1 = (min_val - origin) * ray.d_reciprocal()(i);
                Scalar t2 = (max_val - origin) * ray.d_reciprocal()(i);

                if (t1 > t2) {
                    auto temp = t2;
                    t2        = t1;
                    t1        = temp;
                }

                near_t = M_MAX(t1, near_t);
                far_t  = M_MIN(t2, far_t);

                if (!(near_t <= far_t))
                    return false;
            }
        }

        return true;
    }

    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "Bounding box[\n";
        oss << "  min point: " << m_min_p.to_string() << "\n";
        oss << "  max point: " << m_max_p.to_string() << "\n";
        oss << "]";
        return oss.str();
    }

private:
    Point m_min_p;
    Point m_max_p;
};

M_NAMESPACE_END
