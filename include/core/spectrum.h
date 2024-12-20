#pragma once

#include <cassert>
#include <core/common.h>
#include <cuda_runtime.h>
#include <string>

M_NAMESPACE_BEGIN
template <typename Scalar_, int Dimension_, DeviceType Device_>
class TRGBSpectrum {
public:
    static constexpr DeviceType Device = Device_;
    static constexpr int Dimension     = Dimension_;
    typedef Scalar_ Scalar;
    template <typename, int, DeviceType> friend class TRGBSpectrum;

    // ============================== Constructor ==============================

    M_HOST_DEVICE explicit TRGBSpectrum(Scalar value = 0) {
#ifndef __CUDA_ARCH__
        if constexpr (Device == DeviceType::CPU) {
            std::fill(m_data, m_data + Dimension, value);
        } else {
            auto *temp = new Scalar[Dimension];
            for (auto &i : temp) {
                i = value;
            }
            check_cuda_error(cudaMemcpyAsync(m_data, temp,
                                             Dimension * sizeof(Scalar),
                                             cudaMemcpyHostToDevice),
                             "Memory copy failed");
        }
#else
        if constexpr (Device == DeviceType::CPU) {
            assert(
                false &&
                "Error: Can not construct Spectrum on CPU by device function");
        } else {
            for (int i = 0; i < Dimension; i++) {
                m_data[i] = value;
            }
        }
#endif
    }

    M_HOST_DEVICE explicit TRGBSpectrum(Scalar *input_data,
                                        DeviceType input_device) {
#ifndef __CUDA_ARCH__
        if constexpr (Device == DeviceType::CPU) {
            if (input_device == DeviceType::CPU) {
                std::copy(input_data, input_data + Dimension, m_data);
            } else {
                check_cuda_error(cudaMemcpyAsync(m_data, input_data,
                                                 Dimension * sizeof(Scalar),
                                                 cudaMemcpyDeviceToHost),
                                 "Memory copy failed");
            }
        } else {
            if (input_device == DeviceType::CPU) {
                check_cuda_error(cudaMemcpyAsync(m_data, input_data,
                                                 Dimension * sizeof(Scalar),
                                                 cudaMemcpyHostToDevice),
                                 "Memory copy failed");
            } else {
                check_cuda_error(cudaMemcpyAsync(m_data, input_data,
                                                 Dimension * sizeof(Scalar),
                                                 cudaMemcpyDeviceToDevice),
                                 "Memory copy failed");
            }
        }
#else
        if constexpr (Device == DeviceType::CPU) {
            assert(
                false &&
                "Error: Can not construct Spectrum on CPU by device function");
        } else {
            if (input_device == DeviceType::CPU) {
                assert(false && "Error: Can not construct Spectrum on GPU by "
                                "device function through a CPU m_data");
            } else {
                for (int i = 0; i < Dimension; i++) {
                    m_data[i] = input_data[i];
                }
            }
        }
#endif
    }

    template <int OtherDimension>
    M_HOST_DEVICE explicit TRGBSpectrum(
        TRGBSpectrum<Scalar, OtherDimension, Device> other_spectrum,
        Scalar padding = 1) {
#ifndef __CUDA_ARCH__
        static_assert(OtherDimension <= Dimension,
                      "OtherDimension must be less than Dimension");
        constexpr int CopyDimension =
            OtherDimension < Dimension ? OtherDimension : Dimension;
        if constexpr (Device == DeviceType::CPU) {
            for (int i = 0; i < CopyDimension; i++) {
                m_data[i] = other_spectrum.m_data[i];
            }
            for (int i = CopyDimension; i < Dimension; ++i) {
                m_data[i] = padding;
            }
        } else {
            check_cuda_error(cudaMemcpyAsync(m_data, other_spectrum.m_data,
                                             CopyDimension * sizeof(Scalar),
                                             cudaMemcpyDeviceToDevice),
                             "Memory copy failed");
            check_cuda_error(
                cudaMemset(m_data + CopyDimension, padding,
                           (Dimension - CopyDimension) * sizeof(Scalar)),
                "Memory set failed");
        }
#else
        if constexpr (Device == DeviceType::CPU) {
            assert(
                false &&
                "Error: Can not construct Spectrum on CPU by device function");
        } else {
            for (int i = 0; i < CopyDimension; i++) {
                m_data[i] = other_spectrum.m_data[i];
            }
            for (int i = CopyDimension; i < Dimension; i++) {
                m_data[i] = padding;
            }
        }
#endif
    }

    M_HOST_DEVICE TRGBSpectrum(std::initializer_list<Scalar> values) {
        int i = 0;
#ifndef __CUDA_ARCH__
        if constexpr (Device == DeviceType::CPU) {
            for (auto it = values.begin(); it != values.end() && i < Dimension;
                 ++it, ++i) {
                m_data[i] = *it;
            }
        } else {
            Scalar temp[Dimension] = { 0 };
            for (auto it = values.begin(); it != values.end() && i < Dimension;
                 ++it, ++i) {
                temp[i] = *it;
            }
            check_cuda_error(
                cudaMemcpyAsync(m_data, temp, Dimension * sizeof(Scalar),
                                cudaMemcpyHostToDevice),
                "Memory copy failed in initializer_list constructor");
        }
#else
        if constexpr (Device == DeviceType::CPU) {
            assert(
                false &&
                "Error: Cannot construct Spectrum on CPU in device function");
        } else {
            for (auto it = values.begin(); it != values.end() && i < Dimension;
                 ++it, ++i) {
                m_data[i] = *it;
            }
        }
#endif
        for (; i < Dimension; ++i) {
            m_data[i] =
                static_cast<Scalar>(0); // Fill remaining elements with 0
        }
    }

    M_HOST_DEVICE TRGBSpectrum(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; i++) {
            m_data[i] = other.m_data[i];
        }
    }

    M_HOST_DEVICE TRGBSpectrum &operator=(const TRGBSpectrum &other) {
        if (this != &other) {
            for (int i = 0; i < Dimension; i++) {
                m_data[i] = other.m_data[i];
            }
        }
        return *this;
    }

    // ============================== Function ==============================
    M_HOST_DEVICE TRGBSpectrum operator*(Scalar value) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] * value;
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum operator*(const TRGBSpectrum &other) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] * other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum &operator*=(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] *= value;
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum &operator*=(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] *= other.m_data[i];
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum operator/(Scalar value) const {
        bool valid = true;
        check_zero(value, valid);
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] / value;
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum operator/(const TRGBSpectrum &other) const {
        bool valid = true;
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            check_zero(other.m_data[i], valid);
            result.m_data[i] = this->m_data[i] / other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum &operator/=(Scalar value) {
        bool valid = true;
        check_zero(value, valid);
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] /= value;
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum &operator/=(const TRGBSpectrum &other) {
        bool valid = true;
        for (int i = 0; i < Dimension; ++i) {
            check_zero(other.m_data[i], valid);
            m_data[i] /= other.m_data[i];
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum operator+(Scalar value) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] + value;
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum operator+(const TRGBSpectrum &other) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] + other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum &operator+=(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] += value;
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum &operator+=(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] += other.m_data[i];
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum operator-(Scalar value) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] - value;
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum operator-(const TRGBSpectrum &other) const {
        TRGBSpectrum result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] - other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum &operator-=(Scalar value) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] -= value;
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum &operator-=(const TRGBSpectrum &other) {
        for (int i = 0; i < Dimension; ++i) {
            m_data[i] -= other.m_data[i];
        }
        return *this;
    }

    M_HOST_DEVICE TRGBSpectrum<bool, Dimension, Device>
    operator==(const TRGBSpectrum &other) const {
        TRGBSpectrum<bool, Dimension, Device> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] == other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum<bool, Dimension, Device>
    operator!=(const TRGBSpectrum &other) const {
        TRGBSpectrum<bool, Dimension, Device> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] != other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum<bool, Dimension, Device>
    operator==(Scalar data) const {
        TRGBSpectrum<bool, Dimension, Device> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] == data;
        }
        return result;
    }

    M_HOST_DEVICE TRGBSpectrum<bool, Dimension, Device>
    operator!=(Scalar data) const {
        TRGBSpectrum<bool, Dimension, Device> result;
        for (int i = 0; i < Dimension; ++i) {
            result.m_data[i] = this->m_data[i] != data;
        }
        return result;
    }

    M_HOST_DEVICE bool all() {
        bool result = true;
        for (int i = 0; i < Dimension; ++i) {
            result = result && m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE bool any() {
        bool result = false;
        for (int i = 0; i < Dimension; ++i) {
            result = result || m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE Scalar add() const {
        Scalar sum = 0;

        for (int i = 0; i < Dimension; ++i) {
            sum += m_data[i];
        }
        return sum;
    }

    M_HOST_DEVICE [[nodiscard]] bool is_black() const {
        for (int i = 0; i < Dimension; ++i) {
            if (m_data[i] != static_cast<Scalar>(0)) {
                return false;
            }
        }
        return true;
    }

    M_HOST_DEVICE [[nodiscard]] bool is_valid() const {
        for (int i = 0; i < Dimension; ++i) {
            auto value = m_data[i];
            if (value < 0)
                return false;
        }
        return true;
    }

    M_HOST_DEVICE [[nodiscard]] TRGBSpectrum clamp(Scalar low, Scalar high) {
        TRGBSpectrum result;

        for (int i = 0; i < Dimension; ++i) {
            m_data[i] = clamp(m_data[i], low, high);
        }
        return result;
    }

    M_HOST_DEVICE [[nodiscard]] Scalar max_value() {
        Scalar max_value = 0;
        for (int i = 0; i < Dimension; ++i) {
            if (m_data[i] > max_value) {
                max_value = m_data[i];
            }
        }
        return max_value;
    }

    M_HOST_DEVICE TRGBSpectrum<Scalar, 3, Device> to_srgb() const {
        TRGBSpectrum<Scalar, 3, Device> result;

        if constexpr (Dimension == 3) {
            for (int i = 0; i < 3; ++i) {
                Scalar value = m_data[i];
                if (value <= 0.0031308) {
                    result(i) = 12.92 * value;
                } else {
                    result(i) = (1.0 + 0.055) *
                                    pow(value, static_cast<Scalar>(1.0 / 2.4)) -
                                0.055;
                }
            }
        } else {
            printf("error, unimplemented yet.");
        }

        return result;
    }

    M_HOST_DEVICE Scalar luminance() const {
        static_assert(Dimension == 3, "luminance");
        return m_data[0] * Scalar(0.212671) + m_data[1] * Scalar(0.715160) +
               m_data[2] * Scalar(0.072169);
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

    M_HOST_DEVICE Scalar operator()(int index) const {
        return this->m_data[index];
    }

    M_HOST_DEVICE Scalar &operator()(int index) { return this->m_data[index]; }

    M_HOST_DEVICE TRGBSpectrum<Scalar, Dimension - 1, Device>
    divide_by_weight() {
        TRGBSpectrum<Scalar, Dimension - 1, Device> result;
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

    M_HOST_DEVICE static void check_cuda_error(cudaError_t err,
                                               const char *message) {
#ifdef M_DEBUG
        if (err != cudaSuccess) {
            assert(false && message);
        }
#endif
    }

    M_HOST_DEVICE static void check_zero(Scalar scalar, bool &valid) {
        if (scalar == static_cast<Scalar>(0)) {
            // scalar += static_cast<Scalar>(M_EPSILON);
            valid &= false;
        }
    }
};

M_NAMESPACE_END
