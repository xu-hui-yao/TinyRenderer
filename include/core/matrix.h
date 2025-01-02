#pragma once

#include <cassert>
#include <core/common.h>
#include <string>

M_NAMESPACE_BEGIN
template <typename Scalar, int Rows_, int Cols_> class TMatrix {
public:
    static constexpr int Rows          = Rows_;
    static constexpr int Cols          = Cols_;

    // ============================== Constructor ==============================

    // Set matrix to value
    explicit TMatrix(Scalar value = 0) { std::fill(m_data, m_data + Rows * Cols, value); }

    // Deep copy constructor
    TMatrix(const TMatrix &other) { std::copy(other.m_data, other.m_data + Rows * Cols, m_data); }

    // Deep copy operator
    TMatrix &operator=(const TMatrix &other) {
        if (this == &other) {
            return *this;
        }

        std::copy(other.m_data, other.m_data + Rows * Cols, m_data);

        return *this;
    }

    // ============================== Getter and setter ==============================

    Scalar operator()(int row, int col) const {
        check_range(row, col);

        return m_data[row * Cols + col];
    }

    Scalar &operator()(int row, int col) {
        check_range(row, col);

        return m_data[row * Cols + col];
    }

    void set_constant(Scalar value) {
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m_data[i * Cols + j] = value;
            }
        }
    }

    void set_identity() {
        check_square();

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                m_data[i * Cols + j] = i == j ? static_cast<Scalar>(1) : static_cast<Scalar>(0);
            }
        }
    }

    TMatrix<Scalar, Rows, 1> get_diagonal() {
        check_square();

        TMatrix<Scalar, Rows, 1> result;
        for (int i = 0; i < Rows; ++i) {
            result(i) = m_data[i * Cols + i];
        }
        return result;
    }

    // ============================== Function ==============================

    TMatrix<Scalar, Cols, Rows> transpose() const {
        TMatrix<Scalar, Cols, Rows> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(j, i) = (*this)(i, j);
            }
        }
        return result;
    }

    TMatrix inverse(bool &valid) const {
        check_square();

        TMatrix result;
        // Step 1: LU decomposition
        TMatrix L;
        TMatrix U;

        // Initialize L as an identity matrix
        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                L(i, j) = i == j ? static_cast<Scalar>(1) : static_cast<Scalar>(0);
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
    TMatrix<Scalar, ViewRow, ViewCol> view(int start_row, int start_col) const {
        check_range(start_row + ViewRow - 1, start_col + ViewCol - 1);

        TMatrix<Scalar, ViewRow, ViewCol> sub_matrix;

        for (int i = 0; i < ViewRow; ++i) {
            for (int j = 0; j < ViewCol; ++j) {
                sub_matrix(i, j) = (*this)(start_row + i, start_col + j);
            }
        }
        return sub_matrix;
    }

    TMatrix wise_inverse(bool &valid) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            check_zero(m_data[i], valid);
            result.m_data[i] = static_cast<Scalar>(1) / m_data[i];
        }

        return result;
    }

    TMatrix wise_min(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = M_MIN(m_data[i], other.m_data[i]);
        }
        return result;
    }

    TMatrix wise_max(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = M_MAX(m_data[i], other.m_data[i]);
        }
        return result;
    }

    // Matrix-Matrix multiply
    template <int OtherCols>
    TMatrix<Scalar, Rows, OtherCols> operator*(const TMatrix<Scalar, Cols, OtherCols> &other) const {
        TMatrix<Scalar, Rows, OtherCols> result;

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
    TMatrix<Scalar, Rows, OtherCols> &operator*=(const TMatrix<Scalar, Cols, OtherCols> &other) {
        *this = *this * other;
        return *this;
    }

    // Matrix-Vector multiply
    template <ArrayType OtherArrayType>
    TArray<Scalar, Rows, OtherArrayType>
    operator*(const TArray<Scalar, Cols, OtherArrayType> &vector) const {
        TArray<Scalar, Rows, OtherArrayType> result;

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
    TArray<Scalar, Rows, OtherArrayType> &
    operator*=(const TArray<Scalar, Cols, OtherArrayType> &vector) {
        *this = *this * vector;
        return *this;
    }

    // Matrix-Scalar multiply
    TMatrix operator*(Scalar scalar) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] * scalar;
        }
        return result;
    }

    TMatrix &operator*=(Scalar scalar) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] *= scalar;
        }
        return *this;
    }

    // Matrix-Scalar division
    TMatrix operator/(Scalar scalar) const {
        bool valid = true;
        check_zero(scalar, valid);

        TMatrix result;
        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] / scalar;
        }
        return result;
    }

    TMatrix &operator/=(Scalar scalar) {
        bool valid = true;
        check_zero(scalar, valid);

        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] /= scalar;
        }
        return *this;
    }

    // Matrix-Matrix add
    TMatrix operator+(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] + other.m_data[i];
        }
        return result;
    }

    TMatrix &operator+=(const TMatrix &other) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] += other.m_data[i];
        }
        return *this;
    }

    // Matrix-Matrix sub
    TMatrix operator-(const TMatrix &other) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] - other.m_data[i];
        }
        return result;
    }

    TMatrix &operator-=(const TMatrix &other) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] -= other.m_data[i];
        }
        return *this;
    }

    // Matrix-Scalar add
    TMatrix operator+(Scalar scalar) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] + scalar;
        }
        return result;
    }

    TMatrix &operator+=(Scalar scalar) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] += scalar;
        }
        return *this;
    }

    // Matrix-Scalar sub
    TMatrix operator-(Scalar scalar) const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = m_data[i] - scalar;
        }
        return result;
    }

    TMatrix &operator-=(Scalar scalar) {
        for (int i = 0; i < Rows * Cols; ++i) {
            m_data[i] -= scalar;
        }
        return *this;
    }

    // Operator -
    TMatrix operator-() const {
        TMatrix result;

        for (int i = 0; i < Rows * Cols; ++i) {
            result.m_data[i] = -m_data[i];
        }
        return result;
    }

    // Comparator
    TMatrix<int, Rows, Cols> operator<(const TMatrix &other) const {
        TMatrix<int, Rows, Cols> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] < other(i, j) - M_EPSILON;
            }
        }
        return result;
    }

    TMatrix<int, Rows, Cols> operator<=(const TMatrix &other) const {
        TMatrix<int, Rows, Cols> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] <= other(i, j) - M_EPSILON;
            }
        }
        return result;
    }

    TMatrix<int, Rows, Cols> operator>(const TMatrix &other) const {
        TMatrix<int, Rows, Cols> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] > other(i, j) + M_EPSILON;
            }
        }
        return result;
    }

    TMatrix<int, Rows, Cols> operator>=(const TMatrix &other) const {
        TMatrix<int, Rows, Cols> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = m_data[i * Cols + j] >= other(i, j) + M_EPSILON;
            }
        }
        return result;
    }

    TMatrix<int, Rows, Cols> operator==(const TMatrix &other) const {
        TMatrix<int, Rows, Cols> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = abs(m_data[i * Cols + j] - other(i, j)) < M_EPSILON;
            }
        }
        return result;
    }

    TMatrix<int, Rows, Cols> operator!=(const TMatrix &other) const {
        TMatrix<int, Rows, Cols> result;

        for (int i = 0; i < Rows; ++i) {
            for (int j = 0; j < Cols; ++j) {
                result(i, j) = abs(m_data[i * Cols + j] - other.m_data[i * Cols + j]) >= M_EPSILON;
            }
        }
        return result;
    }

    // Element reduce
    Scalar prod() const {
        Scalar product = 1;

        for (int i = 0; i < Rows * Cols; ++i) {
            product *= m_data[i];
        }
        return product;
    }

    Scalar add() const {
        Scalar sum = 0;

        for (int i = 0; i < Rows * Cols; ++i) {
            sum += m_data[i];
        }
        return sum;
    }

    [[nodiscard]] bool all() const {
        bool result = true;
        for (int i = 0; i < Rows * Cols; ++i) {
            result &= m_data[i] != static_cast<Scalar>(0);
        }
        return result;
    }

    [[nodiscard]] bool any() const {
        bool result = false;
        for (int i = 0; i < Rows * Cols; ++i) {
            result |= m_data[i] != static_cast<Scalar>(0);
        }
        return result;
    }

    Scalar square_distance(const TMatrix &other) const {
        Scalar distance = 0;

        for (int i = 0; i < Rows * Cols; ++i) {
            Scalar diff = m_data[i] - other.m_data[i];
            distance += diff * diff;
        }
        return distance;
    }

    [[nodiscard]] std::string to_string() const {
        std::string matrix;
        matrix += "[\n";
        for (int i = 0; i < Rows; ++i) {
            matrix += "  [";
            for (int j = 0; j < Cols - 1; ++j) {
                matrix += std::to_string(m_data[i * Cols + j]) + ", ";
            }
            matrix += std::to_string(m_data[i * Cols + Cols - 1]) + "]\n";
        }
        matrix += "]";
        return matrix;
    }

    static constexpr TMatrix identity() {
        static_assert(Rows == Cols, "Matrix identity requires a square matrix.");

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

    template <ArrayType ArrayType_> static TMatrix translate(TArray<Scalar, Rows, ArrayType_> vector) {
        static_assert(Rows == Cols, "Matrix translate requires a square matrix.");

        TMatrix result = identity();
        for (int i = 0; i < Rows; ++i) {
            result(i, Rows - 1) = vector(i);
        }
        return result;
    }

    template <ArrayType ArrayType_> static TMatrix scale(TArray<Scalar, Rows, ArrayType_> vector) {
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
    static void check_range(int row, int col) {
#ifdef M_DEBUG
        assert(row < Rows && col < Cols && "Index out of range");
#endif
    }

    static void check_square() {
#ifdef M_DEBUG
        static_assert(Rows == Cols && "Matrix must be square");
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
