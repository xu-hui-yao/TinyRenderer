#pragma once

#include <cassert>
#include <core/common.h>
#include <string>

M_NAMESPACE_BEGIN
template <typename Scalar_, int Dimension_, ArrayType ArrayType_> class TArray {
public:
    template <typename, int, ArrayType> friend class TArray;

    static constexpr int Dimension       = Dimension_;
    static constexpr ArrayType ArrayType = ArrayType_;
    typedef Scalar_ Scalar;

    // ============================== Constructor ==============================
    // Create new vector with constant component values
    explicit TArray(Scalar value = 0) { std::fill(m_data, m_data + Dimension, value); }

    template <enum ArrayType OtherArrayType> explicit TArray(const TArray<Scalar, Dimension, OtherArrayType> &other) {
        std::copy(other.m_data, other.m_data + Dimension, m_data);
    }

    // Deep copy operator
    template <enum ArrayType OtherArrayType> TArray &operator=(const TArray<Scalar, Dimension, OtherArrayType> &other) {
        if ((*this == other).all()) {
            return *this;
        }
        std::copy(other.m_data, other.m_data + Dimension, m_data);
        return *this;
    }

    // Create a new 2D vector (type error if Dimension != 2)
    TArray(Scalar x, Scalar y) {
        check_dimension(2);
        m_data[0] = x;
        m_data[1] = y;
    }

    // Create a new 3D vector (type error if Dimension != 3)
    TArray(Scalar x, Scalar y, Scalar z) {
        check_dimension(3);
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
    }

    // Create a new 4D vector (type error if Dimension != 4)
    TArray(Scalar x, Scalar y, Scalar z, Scalar w) {
        check_dimension(4);
        m_data[0] = x;
        m_data[1] = y;
        m_data[2] = z;
        m_data[3] = w;
    }

    // ============================== Getter and setter
    // ==============================

    Scalar operator()(int index) const {
        check_range(index);
        return m_data[index];
    }

    Scalar &operator()(int index) {
        check_range(index);
        return m_data[index];
    }

    Scalar x() const { return this->operator()(0); }
    Scalar &x() { return this->operator()(0); }
    Scalar y() const { return this->operator()(1); }
    Scalar &y() { return this->operator()(1); }
    Scalar z() const { return this->operator()(2); }
    Scalar &z() { return this->operator()(2); }
    Scalar w() const { return this->operator()(3); }
    Scalar &w() { return this->operator()(3); }

    // ============================== Function ==============================
    void set_constant(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] = value;
        }
    }

    template <enum ArrayType OtherArrayType>
    TArray wise_min(const TArray<Scalar, Dimension, OtherArrayType> &other) const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = M_MIN(m_data[i], other(i));
        }
        return result;
    }

    template <enum ArrayType OtherArrayType>
    TArray wise_max(const TArray<Scalar, Dimension, OtherArrayType> &other) const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = M_MAX(m_data[i], other(i));
        }
        return result;
    }

    template <int ViewDimension> TArray<Scalar, ViewDimension, ArrayType> view(int start_index) const {
        check_range(start_index + ViewDimension - 1);

        TArray<Scalar, ViewDimension, ArrayType> sub_vector;

        for (int i = 0; i < ViewDimension; ++i) {
            sub_vector(i) = m_data[start_index + i];
        }
        return sub_vector;
    }

    TArray wise_inverse(bool &valid) const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            Scalar temp = m_data[i];
            check_zero(temp, valid);
            result.m_data[i] = static_cast<Scalar>(1) / temp;
        }

        return result;
    }

    template <enum ArrayType OtherArrayType>
    TArray cross(const TArray<Scalar, Dimension, OtherArrayType> &other) const {
        check_dimension(3);

        TArray result(0, 0, 0);
        result.x() = y() * other.z() - z() * other.y();
        result.y() = z() * other.x() - x() * other.z();
        result.z() = x() * other.y() - y() * other.x();
        return result;
    }

    template <enum ArrayType OtherArrayType> Scalar dot(const TArray<Scalar, Dimension, OtherArrayType> &other) const {
        Scalar result = 0;

        for (int i = 0; i < Dimension; ++i) {
            result += m_data[i] * other(i);
        }
        return result;
    }

    TArray norm(bool &valid) const {
        Scalar mag = this->magnitude();
        check_zero(mag, valid);
        return *this / mag;
    }

    TArray abs() const {
        TArray result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = m_data[i] < 0 ? -m_data[i] : m_data[i];
        }
        return result;
    }

    TArray clamp(Scalar min, Scalar max) const {
        TArray result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = M_MAX(M_MIN(m_data[i], max), min);
        }
        return result;
    }

    Scalar max_value() const {
        Scalar max = -M_MAX_FLOAT;
        for (int i = 0; i < Dimension; ++i) {
            if (m_data[i] > max) {
                max = m_data[i];
            }
        }
        return max;
    }

    Scalar min_value() const {
        Scalar min = M_MAX_FLOAT;
        for (int i = 0; i < Dimension; ++i) {
            if (m_data[i] < min) {
                min = m_data[i];
            }
        }
        return min;
    }

    void self_norm(bool &valid) {
        Scalar mag = this->magnitude();
        check_zero(mag, valid);

        for (int i = 0; i < Dimension; ++i) {
            m_data[i] /= mag;
        }
    }

    Scalar magnitude() const {
        Scalar sum = 0;

        for (int i = 0; i < Dimension; ++i) {
            auto temp = this->operator()(i);
            sum += temp * temp;
        }
        return sqrt(sum);
    }

    Scalar square_magnitude() const {
        Scalar sum = 0;

        for (int i = 0; i < Dimension; ++i) {
            auto temp = this->operator()(i);
            sum += temp * temp;
        }
        return sum;
    }

    // ============================== Operator ==============================
    TArray operator+(Scalar scalar) const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = m_data[i] + scalar;
        }
        return result;
    }

    TArray &operator+=(Scalar scalar) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] += scalar;
        }
        return *this;
    }

    TArray operator-() const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = -m_data[i];
        }
        return result;
    }

    TArray operator-(Scalar scalar) const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = m_data[i] - scalar;
        }
        return result;
    }

    TArray &operator-=(Scalar scalar) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] -= scalar;
        }
        return *this;
    }

    TArray operator*(Scalar scalar) const {
        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = m_data[i] * scalar;
        }
        return result;
    }

    TArray &operator*=(Scalar scalar) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] *= scalar;
        }
        return *this;
    }

    TArray operator/(Scalar scalar) const {
        bool valid = true;
        check_zero(scalar, valid);

        TArray result;

        for (int i = 0; i < Dimension; ++i) {
            result(i) = m_data[i] / scalar;
        }
        return result;
    }

    TArray &operator/=(Scalar scalar) {
        bool valid = true;
        check_zero(scalar, valid);

        for (int i = 0; i < Dimension; ++i) {
            m_data[i] /= scalar;
        }
        return *this;
    }

    TArray get_floor() {
        TArray result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = floor(m_data[i]);
        }
        return result;
    }

    TArray get_ceil() {
        TArray result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = ceil(m_data[i]);
        }
        return result;
    }

    // Element reduce
    Scalar prod() const {
        Scalar product = 1;

        for (int i = 0; i < Dimension; ++i) {
            product *= m_data[i];
        }
        return product;
    }

    Scalar add() const {
        Scalar sum = 0;

        for (int i = 0; i < Dimension; ++i) {
            sum += m_data[i];
        }
        return sum;
    }

    [[nodiscard]] bool all() const {
        bool result = true;
        for (int i = 0; i < Dimension; ++i) {
            result &= m_data[i] != static_cast<Scalar>(0);
        }
        return result;
    }

    [[nodiscard]] bool any() const {
        bool result = false;
        for (int i = 0; i < Dimension; ++i) {
            result |= m_data[i] != static_cast<Scalar>(0);
        }
        return result;
    }

    template <enum ArrayType OtherArrayType>
    TArray lerp(const TArray<Scalar, Dimension, OtherArrayType> &other, Scalar t) const {
        return *this * (1 - t) + other * t;
    }

    template <enum ArrayType OtherArrayType>
    TArray projection(const TArray<Scalar, Dimension, OtherArrayType> &other) const {
        return other * (dot(other) / other.dot(other));
    }

    [[nodiscard]] std::string to_string() const {

        std::ostringstream oss;
        switch (ArrayType) {
            case ArrayType::Vector:
                oss << "Vector";
                break;
            case ArrayType::Normal:
                oss << "Normal";
                break;
            case ArrayType::Point:
                oss << "Point";
                break;
            default:
                break;
        }
        oss << "[";
        for (int i = 0; i < Dimension - 1; ++i) {
            oss << m_data[i] << ", ";
        }
        oss << m_data[Dimension - 1] << "]";
        return oss.str();
    }

private:
    Scalar m_data[Dimension];

    // ============================== Check ==============================

    static void check_range(int index) {
#ifdef M_DEBUG
        assert(index < Dimension && "Index out of range");
#endif
    }

    static void check_dimension(int expected) {
#ifdef M_DEBUG
        assert(Dimension == expected && "Vector dimensions do not match");
#endif
    }

    static void check_zero(Scalar &scalar, bool &valid) {
        if (scalar == static_cast<Scalar>(0)) {
            scalar += static_cast<Scalar>(M_EPSILON);
            valid &= false;
        }
    }
};

#define M_ARRAY_OP(op, return_type, op1_type, op2_type)                                                                \
    template <typename Scalar, int Dimension> return_type operator op(const op1_type &op1, const op2_type &op2) {      \
        return_type result;                                                                                            \
        for (int i = 0; i < Dimension; ++i) {                                                                          \
            result(i) = op1(i) op op2(i);                                                                              \
        }                                                                                                              \
        return result;                                                                                                 \
    }
#define M_ARRAY_SCALAR_OP(op, return_type, op1_type)                                                                   \
    template <typename Scalar, int Dimension> return_type operator op(const op1_type &op1, Scalar op2) {               \
        return_type result;                                                                                            \
        for (int i = 0; i < Dimension; ++i) {                                                                          \
            result(i) = op1(i) op op2;                                                                                 \
        }                                                                                                              \
        return result;                                                                                                 \
    }
#define M_VECTOR_TYPE TArray<Scalar, Dimension, ArrayType::Vector>
#define M_POINT_TYPE TArray<Scalar, Dimension, ArrayType::Point>
#define M_NORMAL_TYPE TArray<Scalar, Dimension, ArrayType::Normal>
#define M_BOOL_TYPE TArray<int, Dimension, ArrayType::Vector>

M_ARRAY_OP(+, M_POINT_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(+, M_POINT_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(+, M_POINT_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(+, M_VECTOR_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(+, M_VECTOR_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(+, M_VECTOR_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(+, M_POINT_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(+, M_VECTOR_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(+, M_VECTOR_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(-, M_POINT_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(-, M_POINT_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(-, M_VECTOR_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(*, M_POINT_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(*, M_POINT_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(*, M_POINT_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(*, M_VECTOR_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(*, M_VECTOR_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(*, M_VECTOR_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(*, M_VECTOR_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(*, M_VECTOR_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(*, M_VECTOR_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(/, M_POINT_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(/, M_POINT_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(/, M_POINT_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(/, M_VECTOR_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(/, M_VECTOR_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(/, M_VECTOR_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(/, M_VECTOR_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(/, M_VECTOR_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(/, M_VECTOR_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(<, M_BOOL_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(<=, M_BOOL_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(>, M_BOOL_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(>=, M_BOOL_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(==, M_BOOL_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_VECTOR_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_VECTOR_TYPE, M_POINT_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_VECTOR_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_POINT_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_POINT_TYPE, M_POINT_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_POINT_TYPE, M_NORMAL_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_NORMAL_TYPE, M_VECTOR_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_NORMAL_TYPE, M_POINT_TYPE)
M_ARRAY_OP(!=, M_BOOL_TYPE, M_NORMAL_TYPE, M_NORMAL_TYPE)
M_ARRAY_SCALAR_OP(<, M_BOOL_TYPE, M_VECTOR_TYPE)
M_ARRAY_SCALAR_OP(<, M_BOOL_TYPE, M_POINT_TYPE)
M_ARRAY_SCALAR_OP(<, M_BOOL_TYPE, M_NORMAL_TYPE)
M_ARRAY_SCALAR_OP(<=, M_BOOL_TYPE, M_VECTOR_TYPE)
M_ARRAY_SCALAR_OP(<=, M_BOOL_TYPE, M_POINT_TYPE)
M_ARRAY_SCALAR_OP(<=, M_BOOL_TYPE, M_NORMAL_TYPE)
M_ARRAY_SCALAR_OP(>, M_BOOL_TYPE, M_VECTOR_TYPE)
M_ARRAY_SCALAR_OP(>, M_BOOL_TYPE, M_POINT_TYPE)
M_ARRAY_SCALAR_OP(>, M_BOOL_TYPE, M_NORMAL_TYPE)
M_ARRAY_SCALAR_OP(>=, M_BOOL_TYPE, M_VECTOR_TYPE)
M_ARRAY_SCALAR_OP(>=, M_BOOL_TYPE, M_POINT_TYPE)
M_ARRAY_SCALAR_OP(>=, M_BOOL_TYPE, M_NORMAL_TYPE)
M_ARRAY_SCALAR_OP(==, M_BOOL_TYPE, M_VECTOR_TYPE)
M_ARRAY_SCALAR_OP(==, M_BOOL_TYPE, M_POINT_TYPE)
M_ARRAY_SCALAR_OP(==, M_BOOL_TYPE, M_NORMAL_TYPE)
M_ARRAY_SCALAR_OP(!=, M_BOOL_TYPE, M_VECTOR_TYPE)
M_ARRAY_SCALAR_OP(!=, M_BOOL_TYPE, M_POINT_TYPE)
M_ARRAY_SCALAR_OP(!=, M_BOOL_TYPE, M_NORMAL_TYPE)

M_NAMESPACE_END
