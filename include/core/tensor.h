#pragma once

#include <cassert>
#include <core/common.h>
#include <string>

M_NAMESPACE_BEGIN
template <typename Scalar> class TTensor {
public:
    // ============================== Constructor ==============================

    // Set matrix to value
    explicit TTensor(int rows, int cols, int channels, Scalar value = 0)
        : m_rows(rows), m_cols(cols), m_channels(channels) {
        assert(rows > 0 && cols > 0 && channels > 0 && "Rows and columns and channels must be positive integers");
        m_data = new Scalar[m_rows * m_cols * m_channels];
        std::fill(m_data, m_data + m_rows * m_cols * m_channels, value);
    }

    // Initialize matrix by pointer
    explicit TTensor(int rows, int cols, int channels, Scalar *input_m_data)
        : m_rows(rows), m_cols(cols), m_channels(channels) {
        assert(input_m_data != nullptr && "Input data pointer must not be null");
        assert(rows > 0 && cols > 0 && "Rows and columns must be positive integers");
        m_data = new Scalar[m_rows * m_cols * m_channels];
        std::copy(input_m_data, input_m_data + m_rows * m_cols * m_channels, m_data);
    }

    // Deep copy constructor
    TTensor(const TTensor &other) {
        m_rows     = other.m_rows;
        m_cols     = other.m_cols;
        m_channels = other.m_channels;
        delete[] m_data;
        m_data = new Scalar[m_rows * m_cols * m_channels];
        std::copy(other.m_data, other.m_data + m_rows * m_cols * m_channels, m_data);
    }

    // Move constructor
    TTensor(TTensor &&other) noexcept
        : m_rows(other.m_rows), m_cols(other.m_cols), m_channels(other.m_channels), m_data(other.m_data) {
        other.m_data     = nullptr;
        other.m_rows     = 0;
        other.m_cols     = 0;
        other.m_channels = 0;
    }

    // Move assignment
    TTensor &operator=(TTensor &&other) noexcept {
        if (this == &other)
            return *this;

        delete[] m_data;

        m_data     = other.m_data;
        m_rows     = other.m_rows;
        m_cols     = other.m_cols;
        m_channels = other.m_channels;

        other.m_data     = nullptr;
        other.m_rows     = 0;
        other.m_cols     = 0;
        other.m_channels = 0;

        return *this;
    }

    ~TTensor() {
        delete[] m_data;
        m_data = nullptr;
    }

    // Deep copy operator
    TTensor &operator=(const TTensor &other) {
        if (this == &other) {
            return *this;
        }

        delete[] m_data;
        m_data = new Scalar[m_rows * m_cols * m_channels];
        std::copy(other.m_data, other.m_data + m_rows * m_cols * m_channels, m_data);

        return *this;
    }

    // ============================== Getter and setter
    // ==============================

    [[nodiscard]] int get_rows() const { return m_rows; }

    [[nodiscard]] int get_cols() const { return m_cols; }

    [[nodiscard]] int get_channels() const { return m_channels; }

    [[nodiscard]] Scalar *get_data() const { return m_data; }

    [[nodiscard]] Scalar *&get_data() { return m_data; }

    Scalar operator()(int row, int col, int channel) const {
        check_range(row, col, channel);

        return m_data[(row * m_cols + col) * m_channels + channel];
    }

    Scalar &operator()(int row, int col, int channel) {
        check_range(row, col, channel);

        return m_data[(row * m_cols + col) * m_channels + channel];
    }

    void set_constant(Scalar value) {
        for (int i = 0; i < m_rows; ++i) {
            for (int j = 0; j < m_cols; ++j) {
                for (int k = 0; k < m_channels; ++k) {
                    m_data[i * m_cols * m_channels + j * m_channels + k] = value;
                }
            }
        }
    }

    // ============================== Function ==============================

    // Matrix-Scalar multiply
    TTensor operator*(Scalar scalar) const {
        TTensor result(m_rows, m_cols, m_channels, 0);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = m_data[i] * scalar;
        }
        return result;
    }

    TTensor &operator*=(Scalar scalar) {
        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            m_data[i] *= scalar;
        }
        return *this;
    }

    // Matrix-Scalar division
    TTensor operator/(Scalar scalar) const {
        bool valid = true;
        check_zero(scalar, valid);

        TTensor result(m_rows, m_cols, m_channels, 0);
        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = m_data[i] / scalar;
        }
        return result;
    }

    TTensor &operator/=(Scalar scalar) {
        bool valid = true;
        check_zero(scalar, valid);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            m_data[i] /= scalar;
        }
        return *this;
    }

    // Matrix-Matrix add
    TTensor operator+(const TTensor &other) const {
        assert(m_cols == other.get_cols() && m_rows == other.get_rows() && m_channels == other.get_channels() &&
               "Rows and columns must be equal");
        TTensor result(m_rows, m_cols, m_channels, 0);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = m_data[i] + other.m_data[i];
        }
        return result;
    }

    TTensor &operator+=(const TTensor &other) {
        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            m_data[i] += other.m_data[i];
        }
        return *this;
    }

    // Matrix-Matrix sub
    TTensor operator-(const TTensor &other) const {
        assert(m_cols == other.get_cols() && m_rows == other.get_rows() && m_channels == other.get_channels() &&
               "Rows and columns and channels must be equal");
        TTensor result(m_rows, m_cols, m_channels, 0);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = m_data[i] - other.m_data[i];
        }
        return result;
    }

    TTensor &operator-=(const TTensor &other) {
        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            m_data[i] -= other.m_data[i];
        }
        return *this;
    }

    // Matrix-Scalar add
    TTensor operator+(Scalar scalar) const {
        TTensor result(m_rows, m_cols, m_channels, 0);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = m_data[i] + scalar;
        }
        return result;
    }

    TTensor &operator+=(Scalar scalar) {
        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            m_data[i] += scalar;
        }
        return *this;
    }

    // Matrix-Scalar sub
    TTensor operator-(Scalar scalar) const {
        TTensor result(m_rows, m_cols, m_channels, 0);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = m_data[i] - scalar;
        }
        return result;
    }

    TTensor &operator-=(Scalar scalar) {
        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            m_data[i] -= scalar;
        }
        return *this;
    }

    // Operator -
    TTensor operator-() const {
        TTensor result(m_rows, m_cols, m_channels, 0);

        for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
            result.m_data[i] = -m_data[i];
        }
        return result;
    }

    [[nodiscard]] std::string to_string() const {
        if (m_rows == 0 || m_cols == 0 || m_channels == 0 || m_data == nullptr) {
            return "[]";
        }

        std::ostringstream oss;
        oss << "[\n";

        for (int i = 0; i < m_rows; ++i) {
            oss << "  [";
            for (int j = 0; j < m_cols; ++j) {
                oss << "[";
                for (int k = 0; k < m_channels; ++k) {
                    oss << m_data[i * m_cols * m_channels + j * m_channels + k];
                    if (k < m_channels - 1) {
                        oss << ", ";
                    }
                }
                oss << "]";
                if (j < m_cols - 1) {
                    oss << ", ";
                }
            }
            oss << "]";
            if (i < m_rows - 1) {
                oss << ",\n";
            }
        }
        oss << "\n]";
        return oss.str();
    }

private:
    int m_rows = 0, m_cols = 0, m_channels = 0;
    Scalar *m_data = nullptr;

    // ============================== Check ==============================

    static void check_range(int row, int col, int channel) {
#ifdef M_DEBUG
        assert(row < m_rows && col < m_cols && channel < m_channels && "Index out of range");
#endif
    }

    static void check_zero(Scalar &scalar, bool &valid) {
        if (scalar == static_cast<Scalar>(0)) {
            scalar += static_cast<Scalar>(M_EPSILON);
            valid &= false;
        }
    }
};

M_NAMESPACE_END
