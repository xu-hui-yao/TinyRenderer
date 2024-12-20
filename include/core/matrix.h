#pragma once

#include <cassert>
#include <core/common.h>
#include <cuda_runtime.h>
#include <string>

M_NAMESPACE_BEGIN
template <typename Scalar, int Rows_, int Cols_, DeviceType Device_>
class TMatrix {
public:
    static constexpr DeviceType Device = Device_;
    static constexpr int Rows          = Rows_;
    static constexpr int Cols          = Cols_;

    // ============================== Constructor ==============================

    // Set matrix to value
    M_HOST_DEVICE explicit TMatrix(Scalar value = 0) {
#ifndef __CUDA_ARCH__
        if constexpr (Device == DeviceType::CPU) {
            std::fill(m_data, m_data + Rows * Cols, value);
        } else {
            auto *temp = new Scalar[Rows * Cols];
            for (auto &i : temp) {
                i = value;
            }
            check_cuda_error(cudaMemcpyAsync(m_data, temp,
                                             Rows * Cols * sizeof(Scalar),
                                             cudaMemcpyHostToDevice),
                             "Memory copy failed");
            delete[] temp;
        }
#else
        if constexpr (Device == DeviceType::CPU) {
            printf("Error: Can not construct Matrix on CPU by device function");
        } else {
            for (int i = 0; i < Rows * Cols; i++) {
                m_data[i] = value;
            }
        }
#endif
    }

    // Initialize matrix by pointer
    M_HOST_DEVICE explicit TMatrix(Scalar *input_m_data,
                                   DeviceType input_device) {
#ifndef __CUDA_ARCH__
        if constexpr (Device == DeviceType::CPU) {
            if (input_device == DeviceType::CPU) {
                std::copy(input_m_data, input_m_data + Rows * Cols, m_data);
            } else {
                check_cuda_error(cudaMemcpyAsync(m_data, input_m_data,
                                                 Rows * Cols * sizeof(Scalar),
                                                 cudaMemcpyDeviceToHost),
                                 "Memory copy failed");
            }
        } else {
            if (input_device == DeviceType::CPU) {
                check_cuda_error(cudaMemcpyAsync(m_data, input_m_data,
                                                 Rows * Cols * sizeof(Scalar),
                                                 cudaMemcpyHostToDevice),
                                 "Memory copy failed");
            } else {
                check_cuda_error(cudaMemcpyAsync(m_data, input_m_data,
                                                 Rows * Cols * sizeof(Scalar),
                                                 cudaMemcpyDeviceToDevice),
                                 "Memory copy failed");
            }
        }
#else
        if constexpr (Device == DeviceType::CPU) {
            printf("Error: Can not construct Matrix on CPU by device function");
        } else {
            if (input_device == DeviceType::CPU) {
                printf("Error: Can not construct Matrix on GPU by device "
                       "function through a CPU m_data");
            } else {
                for (int i = 0; i < Rows * Cols; i++) {
                    m_data[i] = input_m_data[i];
                }
            }
        }
#endif
    }

    // Deep copy constructor
    M_HOST_DEVICE TMatrix(const TMatrix &other) {
#ifndef __CUDA_ARCH__
        if constexpr (Device == DeviceType::CPU) {
            std::copy(other.m_data, other.m_data + Rows * Cols, m_data);
        } else {
            assert(false &&
                   "Error: Can not construct Matrix on GPU by host function");
        }
#else
        if (Device == DeviceType::CPU) {
            printf("Error: Can not construct Matrix on CPU by device function");
        } else {
            for (int i = 0; i < Rows * Cols; i++) {
                m_data[i] = other.m_data[i];
            }
        }
#endif
    }

    // Deep copy operator
    M_HOST_DEVICE TMatrix &operator=(const TMatrix &other) {
#ifndef __CUDA_ARCH__
        if (this == &other) {
            return *this;
        }

        if constexpr (Device == DeviceType::CPU) {
            std::copy(other.m_data, other.m_data + Rows * Cols, m_data);
        } else {
            check_cuda_error(cudaMemcpyAsync(m_data, other.m_data,
                                             Rows * Cols * sizeof(Scalar),
                                             cudaMemcpyDeviceToDevice),
                             "CUDA memory copy failed");
        }

        return *this;
#else
        if (this == &other) {
            return *this;
        }

        if (Device == DeviceType::CPU) {
            printf("Error: Can not construct Matrix on CPU by device function");
        } else {
            for (int i = 0; i < Rows * Cols; i++) {
                m_data[i] = other.m_data[i];
            }
        }

        return *this;
#endif
    }

    // ============================== Getter and setter
    // ==============================

    M_HOST_DEVICE Scalar operator()(int row, int col) const {
        check_range(row, col);

        return m_data[row * Cols + col];
    }

    M_HOST_DEVICE Scalar &operator()(int row, int col) {
        check_range(row, col);

        return m_data[row * Cols + col];
    }

    M_HOST_DEVICE void set_constant(Scalar value) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m_data[i * Cols + j] = value;
            }
        }
    }

    M_HOST_DEVICE void set_identity() {
        check_square();

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m_data[i * Cols + j] =
                    i == j ? static_cast<Scalar>(1) : static_cast<Scalar>(0);
            }
        }
    }

    M_HOST_DEVICE TMatrix<Scalar, Rows, 1, Device> get_diagonal() {
        check_square();

        TMatrix<Scalar, Rows, 1, Device> result;
        for (int i = 0; i < Rows; ++i) {
            result(i) = m_data[i * Cols + i];
        }
        return result;
    }

    // ============================== Function ==============================

    M_HOST_DEVICE TMatrix<Scalar, Cols, Rows, Device> transpose() const {
        TMatrix<Scalar, Cols, Rows, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    M_HOST_DEVICE TMatrix inverse(bool &valid) const {
        check_square();

        TMatrix result;
        // Step 1: LU decomposition
        TMatrix L;
        TMatrix U;

        // Initialize L as an identity matrix
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                L(i, j) =
                    i == j ? static_cast<Scalar>(1) : static_cast<Scalar>(0);
                U(i, j) = 0;
            }
        }

        // Compute LU decomposition (no parallelism due to m_data dependency)
        for (int i = 0; i < Rows; ++i) {
            // Step 1: Find the pivot row
            int pivot_row = i;
            for (int j = i + 1; j < Rows; ++j) {
                if (abs(U(j, i)) > abs(U(pivot_row, i))) {
                    pivot_row = j;
                }
            }

            // Step 2: Swap rows in U and L if necessary
            if (pivot_row != i) {
                for (int k = 0; k < Cols; ++k) {
                    auto temp       = U(i, k);
                    U(i, k)         = U(pivot_row, k);
                    U(pivot_row, k) = temp;
                }
                for (int k = 0; k < i; ++k) {
                    auto temp       = L(i, k);
                    L(i, k)         = L(pivot_row, k);
                    L(pivot_row, k) = temp;
                }
            }

            // Step 3: Compute U (upper triangular matrix)
            for (int j = i; j < Cols; ++j) {
                Scalar sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += L(i, k) * U(k, j);
                }
                U(i, j) = (*this)(i, j) - sum;
            }

            // Step 4: Compute L (lower triangular matrix)
            for (int j = i + 1; j < Rows; ++j) {
                Scalar sum = 0;
                for (int k = 0; k < i; ++k) {
                    sum += L(j, k) * U(k, i);
                }
                check_zero(U(i, i), valid);
                Scalar denominator = U(i, i);
                L(j, i)            = ((*this)(j, i) - sum) / denominator;
            }
        }

        // Step 2: Compute the inverse of L and U
        TMatrix L_inv;
        TMatrix U_inv;

        // Inverse of L (no parallelism due to dependency)
        for (int i = 0; i < Rows; ++i) {
            check_zero(L(i, i), valid);
            Scalar denominator = L(i, i);
            L_inv(i, i)        = static_cast<Scalar>(1) / denominator;
            for (int j = i + 1; j < Rows; ++j) {
                Scalar sum = 0;
                for (int k = i; k < j; ++k) {
                    sum += L(j, k) * L_inv(k, i);
                }
                check_zero(L(j, j), valid);
                denominator = L(j, j);
                L_inv(j, i) = -sum / denominator;
            }
        }

        // Inverse of U (no parallelism due to dependency)
        for (int i = Rows - 1; i >= 0; --i) {
            check_zero(U(i, i), valid);
            Scalar denominator = U(i, i);
            U_inv(i, i)        = static_cast<Scalar>(1) / denominator;
            for (int j = i - 1; j >= 0; --j) {
                Scalar sum = 0;
                for (int k = j + 1; k <= i; ++k) {
                    sum += U(j, k) * U_inv(k, i);
                }
                check_zero(U(j, j), valid);
                denominator = U(j, j);
                U_inv(j, i) = -sum / denominator;
            }
        }

        // Step 3: Multiply U_inv and L_inv to get the inverse
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                Scalar sum = 0;
                for (int k = 0; k < Rows; ++k) {
                    sum += U_inv(i, k) * L_inv(k, j);
                }
                result(i, j) = sum;
            }
        }

        return result;
    }

    template <int ViewRow, int ViewCol>
    M_HOST_DEVICE TMatrix<Scalar, ViewRow, ViewCol, Device>
    view(int start_row, int start_col) const {
        check_range(start_row + ViewRow - 1, start_col + ViewCol - 1);

        TMatrix<Scalar, ViewRow, ViewCol, Device> sub_matrix;

        for (int i = 0; i < ViewRow; ++i) {
            for (int j = 0; j < ViewCol; ++j) {
                sub_matrix(i, j) = (*this)(start_row + i, start_col + j);
            }
        }
        return sub_matrix;
    }

    M_HOST_DEVICE TMatrix wise_inverse(bool &valid) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            check_zero(m_data[i], valid);
            result.m_data[i] = static_cast<Scalar>(1) / m_data[i];
        }

        return result;
    }

    M_HOST_DEVICE TMatrix wise_min(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = M_MIN(m_data[i], other.m_data[i]);
        }
        return result;
    }

    M_HOST_DEVICE TMatrix wise_max(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = M_MAX(m_data[i], other.m_data[i]);
        }
        return result;
    }

    // Matrix-Matrix multiply
    template <int OtherCols>
    M_HOST_DEVICE TMatrix<Scalar, Rows, OtherCols, Device>
    operator*(const TMatrix<Scalar, Cols, OtherCols, Device> &other) const {
        TMatrix<Scalar, Rows, OtherCols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < OtherCols; ++j) {
                Scalar sum = 0;
                for (int k = 0; k < Cols; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    template <int OtherCols>
    M_HOST_DEVICE TMatrix<Scalar, Rows, OtherCols, Device> &
    operator*=(const TMatrix<Scalar, Cols, OtherCols, Device> &other) {
        *this = *this * other;
        return *this;
    }

    // Matrix-Vector multiply
    template <ArrayType OtherArrayType>
    M_HOST_DEVICE TArray<Scalar, Rows, OtherArrayType, Device> operator*(
        const TArray<Scalar, Cols, OtherArrayType, Device> &vector) const {
        TArray<Scalar, Rows, OtherArrayType, Device> result;

        for (int i = 0; i < Rows; ++i) {
            Scalar sum = 0;
            for (int k = 0; k < Cols; ++k) {
                sum += (*this)(i, k) * vector(k);
            }
            result(i) = sum;
        }
        return result;
    }

    template <ArrayType OtherArrayType>
    M_HOST_DEVICE TArray<Scalar, Rows, OtherArrayType, Device> &
    operator*=(const TArray<Scalar, Cols, OtherArrayType, Device> &vector) {
        *this = *this * vector;
        return *this;
    }

    // Matrix-Scalar multiply
    M_HOST_DEVICE TMatrix operator*(Scalar scalar) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] * scalar;
        }
        return result;
    }

    M_HOST_DEVICE TMatrix &operator*=(Scalar scalar) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] *= scalar;
        }
        return *this;
    }

    // Matrix-Scalar division
    M_HOST_DEVICE TMatrix operator/(Scalar scalar) const {
        bool valid = true;
        check_zero(scalar, valid);

        TMatrix result;
        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] / scalar;
        }
        return result;
    }

    M_HOST_DEVICE TMatrix &operator/=(Scalar scalar) {
        bool valid = true;
        check_zero(scalar, valid);

        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] /= scalar;
        }
        return *this;
    }

    // Matrix-Matrix add
    M_HOST_DEVICE TMatrix operator+(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] + other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TMatrix &operator+=(const TMatrix &other) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] += other.m_data[i];
        }
        return *this;
    }

    // Matrix-Matrix sub
    M_HOST_DEVICE TMatrix operator-(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] - other.m_data[i];
        }
        return result;
    }

    M_HOST_DEVICE TMatrix &operator-=(const TMatrix &other) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] -= other.m_data[i];
        }
        return *this;
    }

    // Matrix-Scalar add
    M_HOST_DEVICE TMatrix operator+(Scalar scalar) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] + scalar;
        }
        return result;
    }

    M_HOST_DEVICE TMatrix &operator+=(Scalar scalar) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] += scalar;
        }
        return *this;
    }

    // Matrix-Scalar sub
    M_HOST_DEVICE TMatrix operator-(Scalar scalar) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] - scalar;
        }
        return result;
    }

    M_HOST_DEVICE TMatrix &operator-=(Scalar scalar) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] -= scalar;
        }
        return *this;
    }

    // Operator -
    M_HOST_DEVICE TMatrix operator-() const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = -m_data[i];
        }
        return result;
    }

    // Comparator
    M_HOST_DEVICE TMatrix<int, Rows, Cols, Device>
    operator<(const TMatrix &other) const {
        TMatrix<int, Rows, Cols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] < other(i, j) - M_EPSILON;
            }
        }
        return result;
    }

    M_HOST_DEVICE TMatrix<int, Rows, Cols, Device>
    operator<=(const TMatrix &other) const {
        TMatrix<int, Rows, Cols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] <= other(i, j) - M_EPSILON;
            }
        }
        return result;
    }

    M_HOST_DEVICE TMatrix<int, Rows, Cols, Device>
    operator>(const TMatrix &other) const {
        TMatrix<int, Rows, Cols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] > other(i, j) + M_EPSILON;
            }
        }
        return result;
    }

    M_HOST_DEVICE TMatrix<int, Rows, Cols, Device>
    operator>=(const TMatrix &other) const {
        TMatrix<int, Rows, Cols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] >= other(i, j) + M_EPSILON;
            }
        }
        return result;
    }

    M_HOST_DEVICE TMatrix<int, Rows, Cols, Device>
    operator==(const TMatrix &other) const {
        TMatrix<int, Rows, Cols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) =
                    abs(m_data[i * Cols + j] - other(i, j)) < M_EPSILON;
            }
        }
        return result;
    }

    M_HOST_DEVICE TMatrix<int, Rows, Cols, Device>
    operator!=(const TMatrix &other) const {
        TMatrix<int, Rows, Cols, Device> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = abs(m_data[i * Cols + j] -
                                   other.m_data[i * Cols + j]) >= M_EPSILON;
            }
        }
        return result;
    }

    // Element reduce
    M_HOST_DEVICE Scalar prod() const {
        Scalar product = 1;

        for (int i = 0; i < Rows * Cols; ++i) {
            product *= m_data[i];
        }
        return product;
    }

    M_HOST_DEVICE Scalar add() const {
        Scalar sum = 0;

        for (int i = 0; i < Rows * Cols; ++i) {
            sum += m_data[i];
        }
        return sum;
    }

    M_HOST_DEVICE [[nodiscard]] bool all() const {
        bool result = true;
        for (int i = 0; i < Rows * Cols; ++i) {
            result &= m_data[i] != static_cast<Scalar>(0);
        }
        return result;
    }

    M_HOST_DEVICE [[nodiscard]] bool any() const {
        bool result = false;
        for (int i = 0; i < Rows * Cols; ++i) {
            result |= m_data[i] != static_cast<Scalar>(0);
        }
        return result;
    }

    M_HOST_DEVICE Scalar square_distance(const TMatrix &other) const {
        Scalar distance = 0;

        for (int i = 0; i < Rows * Cols; ++i) {
            Scalar diff = m_data[i] - other.m_data[i];
            distance += diff * diff;
        }
        return distance;
    }

    [[nodiscard]] std::string to_string() const {
        Scalar temp[Rows * Cols];
        if constexpr (Device == DeviceType::GPU) {
            check_cuda_error(cudaMemcpyAsync(temp, m_data,
                                             Rows * Cols * sizeof(Scalar),
                                             cudaMemcpyDeviceToHost));
            cudaDeviceSynchronize();
        } else {
            std::copy(m_data, m_data + Rows * Cols, temp);
        }

        std::string matrix;
        matrix += "[\n";
        for (int i = 0; i < Rows; ++i) {
            matrix += "  [";
            for (int j = 0; j < Cols - 1; ++j) {
                matrix += std::to_string(temp[i * Cols + j]) + ", ";
            }
            matrix += std::to_string(temp[i * Cols + Cols - 1]) + "]\n";
        }
        matrix += "]";
        return matrix;
    }

    static constexpr TMatrix identity() {
        static_assert(Rows == Cols,
                      "Matrix identity requires a square matrix.");

        TMatrix result(0);
        for (int i = 0; i < Rows; ++i) {
            result(i, i) = static_cast<Scalar>(1);
        }
        return result;
    }

    static constexpr TMatrix zero() {
        TMatrix result(0);
        return result;
    }

    template <ArrayType ArrayType_>
    static TMatrix translate(TArray<Scalar, Rows, ArrayType_, Device> vector) {
        static_assert(Rows == Cols,
                      "Matrix translate requires a square matrix.");

        TMatrix result = identity();
        for (int i = 0; i < Rows; ++i) {
            result(i, Rows - 1) = vector(i);
        }
        return result;
    }

    template <ArrayType ArrayType_>
    static TMatrix scale(TArray<Scalar, Rows, ArrayType_, Device> vector) {
        static_assert(Rows == Cols, "Matrix scale requires a square matrix.");

        TMatrix result = identity();
        for (int i = 0; i < Rows; ++i) {
            result(i, i) = vector(i);
        }
        return result;
    }

private:
    Scalar m_data[Rows * Cols];

    // ============================== Check ==============================

    M_HOST_DEVICE static void check_cuda_error(cudaError_t err,
                                               const char *message) {
#ifdef M_DEBUG
        if (err != cudaSuccess) {
            assert(false && message);
        }
#endif
    }

    M_HOST_DEVICE static void check_range(int row, int col) {
#ifdef M_DEBUG
        assert(row < Rows && col < Cols && "Index out of range");
#endif
    }

    M_HOST_DEVICE static void check_square() {
#ifdef M_DEBUG
#ifdef __CUDA_ARCH__
        if (Rows != Cols) {
            printf("Matrix must be square");
        }
#else
        static_assert(Rows == Cols && "Matrix must be square");
#endif
#endif
    }

    M_HOST_DEVICE static void check_zero(Scalar &scalar, bool &valid) {
        if (scalar == static_cast<Scalar>(0)) {
            scalar += static_cast<Scalar>(M_EPSILON);
            valid &= false;
        }
    }
};

M_NAMESPACE_END
