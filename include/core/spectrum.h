#pragma once

#include <cassert>
#include <core/common.h>
#include <cuda_runtime.h>
#include <string>

M_NAMESPACE_BEGIN
template <typename Scalar_, int Dimension_> class TRGBSpectrum {
public:
    static constexpr int Dimension = Dimension_;
    typedef Scalar_ Scalar;
    template <typename, int> friend class TRGBSpectrum;

    // ============================== Constructor ==============================

    explicit TRGBSpectrum(Scalar value = 0) { std::fill(m_data, m_data + Dimension, value); }

    template <int OtherDimension>
    explicit TRGBSpectrum(TRGBSpectrum<Scalar, OtherDimension> other_spectrum, Scalar padding = 1) {
        static_assert(OtherDimension <= Dimension, "OtherDimension must be less than Dimension");
        constexpr int CopyDimension = OtherDimension < Dimension ? OtherDimension : Dimension;
        for (int i = 0; i < CopyDimension; i++) {
            m_data[i] = other_spectrum.m_data[i];
        }
        for (int i = CopyDimension; i < Dimension; ++i) {
            m_data[i] = padding;
        }
    }

    TRGBSpectrum(std::initializer_list<Scalar> values) {
        int i = 0;
        for (auto it = values.begin(); it != values.end() && i < Dimension; ++it, ++i) {
            m_data[i] = *it;
        }
        for (; i < Dimension; ++i) {
            m_data[i] = static_cast<Scalar>(0); // Fill remaining elements with 0
        }
    }

    TRGBSpectrum(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; i++) {
            m_data[i] = other.m_data[i];
        }
    }

    TRGBSpectrum &operator=(const TRGBSpectrum &other) {
        if (this != &other) {
            for (int i = 0; i < Dimension; i++) {
                m_data[i] = other.m_data[i];
            }
        }
        return *this;
    }

    // ============================== Function ==============================
    TRGBSpectrum operator*(Scalar value) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] * value;
        }
        return result;
    }

    TRGBSpectrum operator*(const TRGBSpectrum &other) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] * other.m_data[i];
        }
        return result;
    }

    TRGBSpectrum &operator*=(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] *= value;
        }
        return *this;
    }

    TRGBSpectrum &operator*=(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] *= other.m_data[i];
        }
        return *this;
    }

    TRGBSpectrum operator/(Scalar value) const {
        bool valid = true;
        check_zero(value, valid);
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] / value;
        }
        return result;
    }

    TRGBSpectrum operator/(const TRGBSpectrum &other) const {
        bool valid = true;
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            check_zero(other.m_data[i], valid);
            result.m_data[i] = this->m_data[i] / other.m_data[i];
        }
        return result;
    }

    TRGBSpectrum &operator/=(Scalar value) {
        bool valid = true;
        check_zero(value, valid);
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] /= value;
        }
        return *this;
    }

    TRGBSpectrum &operator/=(const TRGBSpectrum &other) {
        bool valid = true;
        for (int i = 0; i < Dimension; ++i) {
            check_zero(other.m_data[i], valid);
            m_data[i] /= other.m_data[i];
        }
        return *this;
    }

    TRGBSpectrum operator+(Scalar value) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] + value;
        }
        return result;
    }

    TRGBSpectrum operator+(const TRGBSpectrum &other) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] + other.m_data[i];
        }
        return result;
    }

    TRGBSpectrum &operator+=(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] += value;
        }
        return *this;
    }

    TRGBSpectrum &operator+=(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] += other.m_data[i];
        }
        return *this;
    }

    TRGBSpectrum operator-(Scalar value) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] - value;
        }
        return result;
    }

    TRGBSpectrum operator-(const TRGBSpectrum &other) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] - other.m_data[i];
        }
        return result;
    }

    TRGBSpectrum &operator-=(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] -= value;
        }
        return *this;
    }

    TRGBSpectrum &operator-=(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] -= other.m_data[i];
        }
        return *this;
    }

    TRGBSpectrum<bool, Dimension> operator==(const TRGBSpectrum &other) const {
        TRGBSpectrum<bool, Dimension> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] == other.m_data[i];
        }
        return result;
    }

    TRGBSpectrum<bool, Dimension> operator!=(const TRGBSpectrum &other) const {
        TRGBSpectrum<bool, Dimension> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] != other.m_data[i];
        }
        return result;
    }

    TRGBSpectrum<bool, Dimension> operator==(Scalar data) const {
        TRGBSpectrum<bool, Dimension> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] == data;
        }
        return result;
    }

    TRGBSpectrum<bool, Dimension> operator!=(Scalar data) const {
        TRGBSpectrum<bool, Dimension> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] != data;
        }
        return result;
    }

    bool all() {
        bool result = true;
        for (int i = 0; i < Dimension; ++i) {
            result = result && m_data[i];
        }
        return result;
    }

    bool any() {
        bool result = false;
        for (int i = 0; i < Dimension; ++i) {
            result = result || m_data[i];
        }
        return result;
    }

    Scalar add() const {
        Scalar sum = 0;

        for (int i = 0; i < Dimension; ++i) {
            sum += m_data[i];
        }
        return sum;
    }

    [[nodiscard]] bool is_black() const {
        for (int i = 0; i < Dimension; ++i) {
            if (m_data[i] != static_cast<Scalar>(0)) {
                return false;
            }
        }
        return true;
    }

    [[nodiscard]] bool is_valid() const {
        for (int i = 0; i < Dimension; ++i) {
            auto value = m_data[i];
            if (value < 0)
                return false;
        }
        return true;
    }

    [[nodiscard]] TRGBSpectrum clamp(Scalar low, Scalar high) {
        TRGBSpectrum result;

        for (int i = 0; i < Dimension; ++i) {
            m_data[i] = clamp(m_data[i], low, high);
        }
        return result;
    }

    [[nodiscard]] Scalar max_value() {
        Scalar max_value = 0;
        for (int i = 0; i < Dimension; ++i) {
            if (m_data[i] > max_value) {
                max_value = m_data[i];
            }
        }
        return max_value;
    }

    TRGBSpectrum<Scalar, 3> to_srgb() const {
        TRGBSpectrum<Scalar, 3> result;

        if constexpr (Dimension == 3) {
            for (int i = 0; i < 3; ++i) {
                Scalar value = m_data[i];
                if (value <= 0.0031308) {
                    result(i) = 12.92 * value;
                } else {
                    result(i) = (1.0 + 0.055) * pow(value, static_cast<Scalar>(1.0 / 2.4)) - 0.055;
                }
            }
        } else {
            printf("error, unimplemented yet.");
        }

        return result;
    }

    Scalar luminance() const {
        static_assert(Dimension == 3, "luminance");
        return m_data[0] * Scalar(0.212671) + m_data[1] * Scalar(0.715160) + m_data[2] * Scalar(0.072169);
    }

    [[nodiscard]] std::string to_string() const {
        std::ostringstream oss;
        oss << "Spectrum[";
        for (int i = 0; i < Dimension - 1; ++i) {
            oss << m_data[i] << ", ";
        }
        oss << m_data[Dimension - 1];
        oss << "]";
        return oss.str();
    }

    Scalar operator()(int index) const { return this->m_data[index]; }

    Scalar &operator()(int index) { return this->m_data[index]; }

    TRGBSpectrum<Scalar, Dimension - 1> divide_by_weight() {
        TRGBSpectrum<Scalar, Dimension - 1> result;
        for (int i = 0; i < Dimension - 1; ++i) {
            if (m_data[Dimension - 1] > M_EPSILON) {
                result.m_data[i] = m_data[i] / m_data[Dimension - 1];
            } else {
                result.m_data[i] = 0;
            }
        }
        return result;
    }

private:
    Scalar m_data[Dimension];

    static void check_zero(Scalar scalar, bool &valid) {
        if (scalar == static_cast<Scalar>(0)) {
            valid &= false;
        }
    }
};

M_NAMESPACE_END
