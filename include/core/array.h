#pragma once

#include <cassert>
#include <string>
#include <cuda_runtime.h>
#include <core/common.h>

M_NAMESPACE_BEGIN
	template <typename Scalar_, int Dimension_, ArrayType ArrayType_, DeviceType Device_>
	class TArray {
	public:
		template <typename, int, ArrayType, DeviceType>
		friend class TArray;

		static constexpr DeviceType Device = Device_;
		static constexpr int Dimension = Dimension_;
		static constexpr ArrayType ArrayType = ArrayType_;
		typedef Scalar_ Scalar;

		// ============================== Constructor ==============================
		// Create new vector with constant component values
		M_HOST_DEVICE explicit TArray(Scalar value = 0) {
#ifndef __CUDA_ARCH__
			if constexpr (Device == DeviceType::CPU) {
				std::fill(m_data, m_data + Dimension, value);
			} else {
				auto* temp = new Scalar[Dimension];
				for (auto& i : temp) {
					i = value;
				}
				check_cuda_error(
					cudaMemcpyAsync(m_data, temp, Dimension * sizeof(Scalar), cudaMemcpyHostToDevice),
					"Memory copy failed");
			}
#else
			if constexpr (Device == DeviceType::CPU) {
				assert(false && "Error: Can not construct Vector on CPU by device function");
			} else {
				for (int i = 0; i < Dimension; i++) {
					m_data[i] = value;
				}
			}
#endif
		}

		M_HOST_DEVICE explicit TArray(Scalar* input_data, DeviceType input_device) {
#ifndef __CUDA_ARCH__
			if constexpr (Device == DeviceType::CPU) {
				if (input_device == DeviceType::CPU) {
					std::copy(input_data, input_data + Dimension, m_data);
				} else {
					check_cuda_error(
						cudaMemcpyAsync(m_data, input_data, Dimension * sizeof(Scalar), cudaMemcpyDeviceToHost),
						"Memory copy failed");
				}
			} else {
				if (input_device == DeviceType::CPU) {
					check_cuda_error(
						cudaMemcpyAsync(m_data, input_data, Dimension * sizeof(Scalar), cudaMemcpyHostToDevice),
						"Memory copy failed");
				} else {
					check_cuda_error(
						cudaMemcpyAsync(m_data, input_data, Dimension * sizeof(Scalar), cudaMemcpyDeviceToDevice),
						"Memory copy failed");
				}
			}
#else
			if constexpr (Device == DeviceType::CPU) {
				assert(false && "Error: Can not construct Vector on CPU by device function");
			} else {
				if (input_device == DeviceType::CPU) {
					assert(false && "Error: Can not construct Vector on GPU by device function through a CPU m_data");
				} else {
					for (int i = 0; i < Dimension; i++) {
						m_data[i] = input_data[i];
					}
				}
			}
#endif
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE explicit TArray(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) {
#ifndef __CUDA_ARCH__
			if constexpr (Device == DeviceType::CPU) {
				std::copy(other.m_data, other.m_data + Dimension, m_data);
			} else {
				assert(false && "Error: Can not construct Vector on GPU by host function");
			}
#else
			if (Device == DeviceType::CPU) {
				printf("Error: Can not construct Vector on CPU by device function");
			} else {
				for (int i = 0; i < Dimension; i++) {
					m_data[i] = other.m_data[i];
				}
			}
#endif
		}

		// Deep copy operator
		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE TArray& operator=(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) {
#ifndef __CUDA_ARCH__
			if ((*this == other).all()) {
				return *this;
			}

			if constexpr (Device == DeviceType::CPU) {
				std::copy(other.m_data, other.m_data + Dimension, m_data);
			} else {
				check_cuda_error(
					cudaMemcpyAsync(m_data, other.m_data, Dimension * sizeof(Scalar), cudaMemcpyDeviceToDevice),
					"CUDA memory copy failed");
			}

			return *this;
#else
			if ((*this == other),all()) {
				return *this;
			}

			if (Device == DeviceType::CPU) {
				printf("Error: Can not construct Matrix on CPU by device function");
			} else {
				for (int i = 0; i < Dimension; i++) {
					m_data[i] = other.m_data[i];
				}
			}

			return *this;
#endif
		}

		// Create a new 2D vector (type error if Dimension != 2)
		M_HOST_DEVICE TArray(Scalar x, Scalar y) {
			check_dimension(2);
			m_data[0] = x;
			m_data[1] = y;
		}

		// Create a new 3D vector (type error if Dimension != 3)
		M_HOST_DEVICE TArray(Scalar x, Scalar y, Scalar z) {
			check_dimension(3);
			m_data[0] = x;
			m_data[1] = y;
			m_data[2] = z;
		}

		// Create a new 4D vector (type error if Dimension != 4)
		M_HOST_DEVICE TArray(Scalar x, Scalar y, Scalar z, Scalar w) {
			check_dimension(4);
			m_data[0] = x;
			m_data[1] = y;
			m_data[2] = z;
			m_data[3] = w;
		}

		// ============================== Getter and setter ==============================

		M_HOST_DEVICE Scalar operator()(int index) const {
			check_range(index);
			return m_data[index];
		}

		M_HOST_DEVICE Scalar& operator()(int index) {
			check_range(index);
			return m_data[index];
		}

		M_HOST_DEVICE Scalar x() const { return this->operator()(0); }
		M_HOST_DEVICE Scalar& x() { return this->operator()(0); }
		M_HOST_DEVICE Scalar y() const { return this->operator()(1); }
		M_HOST_DEVICE Scalar& y() { return this->operator()(1); }
		M_HOST_DEVICE Scalar z() const { return this->operator()(2); }
		M_HOST_DEVICE Scalar& z() { return this->operator()(2); }
		M_HOST_DEVICE Scalar w() const { return this->operator()(3); }
		M_HOST_DEVICE Scalar& w() { return this->operator()(3); }

		// ============================== Function ==============================
		M_HOST_DEVICE void set_constant(Scalar value) {
			for (int i = 0; i < Dimension; ++i) {
				m_data[i] = value;
			}
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE TArray wise_min(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = M_MIN(m_data[i], other(i));
			}
			return result;
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE TArray wise_max(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = M_MAX(m_data[i], other(i));
			}
			return result;
		}

		template <int ViewDimension>
		M_HOST_DEVICE TArray<Scalar, ViewDimension, ArrayType, Device> view(int start_index) const {
			check_range(start_index + ViewDimension - 1);

			TArray<Scalar, ViewDimension, ArrayType, Device> sub_vector;

			for (int i = 0; i < ViewDimension; ++i) {
				sub_vector(i) = m_data[start_index + i];
			}
			return sub_vector;
		}

		M_HOST_DEVICE TArray wise_inverse(bool& valid) const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				Scalar temp = m_data[i];
				check_zero(temp, valid);
				result.m_data[i] = static_cast<Scalar>(1) / temp;
			}

			return result;
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE TArray cross(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) const {
			check_dimension(3);

			TArray result(0, 0, 0);
			result.x() = y() * other.z() - z() * other.y();
			result.y() = z() * other.x() - x() * other.z();
			result.z() = x() * other.y() - y() * other.x();
			return result;
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE Scalar dot(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) const {
			Scalar result = 0;

			for (int i = 0; i < Dimension; ++i) {
				result += m_data[i] * other(i);
			}
			return result;
		}

		M_HOST_DEVICE TArray norm(bool& valid) const {
			Scalar mag = this->magnitude();
			check_zero(mag, valid);
			return *this / mag;
		}

		M_HOST_DEVICE TArray abs() const {
			TArray result;
			for (int i = 0; i < Dimension; ++i) {
				result.m_data[i] = m_data[i] < 0 ? -m_data[i] : m_data[i];
			}
			return result;
		}

		M_HOST_DEVICE TArray clamp(Scalar min, Scalar max) const {
			TArray result;
			for (int i = 0; i < Dimension; ++i) {
				result.m_data[i] = M_MAX(M_MIN(m_data[i], max), min);
			}
			return result;
		}

		M_HOST_DEVICE Scalar max_value() const {
			Scalar max = -M_MAX_FLOAT;
			for (int i = 0; i < Dimension; ++i) {
				if (m_data[i] > max) {
					max = m_data[i];
				}
			}
			return max;
		}

		M_HOST_DEVICE Scalar min_value() const {
			Scalar min = M_MAX_FLOAT;
			for (int i = 0; i < Dimension; ++i) {
				if (m_data[i] < min) {
					min = m_data[i];
				}
			}
			return min;
		}

		M_HOST_DEVICE void self_norm(bool& valid) {
			Scalar mag = this->magnitude();
			check_zero(mag, valid);

			for (int i = 0; i < Dimension; ++i) {
				m_data[i] /= mag;
			}
		}

		M_HOST_DEVICE Scalar magnitude() const {
			Scalar sum = 0;

			for (int i = 0; i < Dimension; ++i) {
				auto temp = this->operator()(i);
				sum += temp * temp;
			}
			return sqrt(sum);
		}

		M_HOST_DEVICE Scalar square_magnitude() const {
			Scalar sum = 0;

			for (int i = 0; i < Dimension; ++i) {
				auto temp = this->operator()(i);
				sum += temp * temp;
			}
			return sum;
		}

		// ============================== Operator ==============================
		M_HOST_DEVICE TArray operator+(Scalar scalar) const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = m_data[i] + scalar;
			}
			return result;
		}

		M_HOST_DEVICE TArray& operator+=(Scalar scalar) {
			for (int i = 0; i < Dimension; ++i) {
				m_data[i] += scalar;
			}
			return *this;
		}

		M_HOST_DEVICE TArray operator-() const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = -m_data[i];
			}
			return result;
		}

		M_HOST_DEVICE TArray operator-(Scalar scalar) const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = m_data[i] - scalar;
			}
			return result;
		}

		M_HOST_DEVICE TArray& operator-=(Scalar scalar) {
			for (int i = 0; i < Dimension; ++i) {
				m_data[i] -= scalar;
			}
			return *this;
		}

		M_HOST_DEVICE TArray operator*(Scalar scalar) const {
			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = m_data[i] * scalar;
			}
			return result;
		}

		M_HOST_DEVICE TArray& operator*=(Scalar scalar) {
			for (int i = 0; i < Dimension; ++i) {
				m_data[i] *= scalar;
			}
			return *this;
		}

		M_HOST_DEVICE TArray operator/(Scalar scalar) const {
			bool valid = true;
			check_zero(scalar, valid);

			TArray result;

			for (int i = 0; i < Dimension; ++i) {
				result(i) = m_data[i] / scalar;
			}
			return result;
		}

		M_HOST_DEVICE TArray& operator/=(Scalar scalar) {
			bool valid = true;
			check_zero(scalar, valid);

			for (int i = 0; i < Dimension; ++i) {
				m_data[i] /= scalar;
			}
			return *this;
		}

		// Element reduce
		M_HOST_DEVICE Scalar prod() const {
			Scalar product = 1;

			for (int i = 0; i < Dimension; ++i) {
				product *= m_data[i];
			}
			return product;
		}

		M_HOST_DEVICE Scalar add() const {
			Scalar sum = 0;

			for (int i = 0; i < Dimension; ++i) {
				sum += m_data[i];
			}
			return sum;
		}

		M_HOST_DEVICE [[nodiscard]] bool all() const {
			bool result = true;
			for (int i = 0; i < Dimension; ++i) {
				result &= m_data[i] != static_cast<Scalar>(0);
			}
			return result;
		}

		M_HOST_DEVICE [[nodiscard]] bool any() const {
			bool result = false;
			for (int i = 0; i < Dimension; ++i) {
				result |= m_data[i] != static_cast<Scalar>(0);
			}
			return result;
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE TArray lerp(const TArray<Scalar, Dimension, OtherArrayType, Device>& other, Scalar t) const {
			return *this * (1 - t) + other * t;
		}

		template <enum ArrayType OtherArrayType>
		M_HOST_DEVICE TArray projection(const TArray<Scalar, Dimension, OtherArrayType, Device>& other) const {
			return other * (dot(other) / other.dot(other));
		}

		[[nodiscard]] std::string to_string() const {
			Scalar temp[Dimension];
			if constexpr (Device == DeviceType::GPU) {
				check_cuda_error(cudaMemcpyAsync(temp, m_data, Dimension * sizeof(Scalar), cudaMemcpyDeviceToHost));
				cudaDeviceSynchronize();
			} else {
				std::copy(m_data, m_data + Dimension, temp);
			}

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
				oss << temp[i] << ", ";
			}
			oss << temp[Dimension - 1] << "]";
			return oss.str();
		}

	private:
		Scalar m_data[Dimension];

		// ============================== Check ==============================

		M_HOST_DEVICE static void check_range(int index) {
#ifdef M_DEBUG
			assert(index < Dimension && "Index out of range");
#endif
		}

		M_HOST_DEVICE static void check_dimension(int expected) {
#ifdef M_DEBUG
			assert(Dimension == expected && "Vector dimensions do not match");
#endif
		}

		M_HOST_DEVICE static void check_zero(Scalar& scalar, bool& valid) {
			if (scalar == static_cast<Scalar>(0)) {
				scalar += static_cast<Scalar>(M_EPSILON);
				valid &= false;
			}
		}

		M_HOST_DEVICE static void check_cuda_error(cudaError_t err, const char* message) {
#ifdef M_DEBUG
			if (err != cudaSuccess) {
				assert(false && message);
			}
#endif
		}
	};

#define M_ARRAY_OP(op, return_type, op1_type, op2_type) \
template <typename Scalar, int Dimension, DeviceType Device> \
M_HOST_DEVICE return_type operator op (const op1_type &op1, const op2_type &op2) { \
	return_type result; \
	for (int i = 0; i < Dimension; ++i) { \
		result(i) = op1(i) op op2(i); \
	} \
	return result; \
}
#define M_ARRAY_SCALAR_OP(op, return_type, op1_type) \
template <typename Scalar, int Dimension, DeviceType Device> \
M_HOST_DEVICE return_type operator op (const op1_type &op1, Scalar op2) { \
return_type result; \
for (int i = 0; i < Dimension; ++i) { \
result(i) = op1(i) op op2; \
} \
return result; \
}
#define M_VECTOR_TYPE TArray<Scalar, Dimension, ArrayType::Vector, Device>
#define M_POINT_TYPE TArray<Scalar, Dimension, ArrayType::Point, Device>
#define M_NORMAL_TYPE TArray<Scalar, Dimension, ArrayType::Normal, Device>
#define M_BOOL_TYPE TArray<int, Dimension, ArrayType::Vector, Device>

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
