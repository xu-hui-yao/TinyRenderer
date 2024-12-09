#pragma once

#include <iostream>
#include <core/ray.h>
#include <core/matrix.h>

M_NAMESPACE_BEGIN
	template <typename Scalar_, DeviceType Device_>
	class TTransform {
		static constexpr DeviceType Device = Device_;

		typedef Scalar_ Scalar;
		typedef TMatrix<Scalar, 4, 4, Device> Matrix4;
		typedef TMatrix<Scalar, 3, 3, Device> Matrix3;
		typedef TArray<Scalar, 3, ArrayType::Vector, Device> VectorType3;
		typedef TArray<Scalar, 4, ArrayType::Vector, Device> VectorType4;
		typedef TArray<Scalar, 3, ArrayType::Point, Device> PointType3;
		typedef TArray<Scalar, 3, ArrayType::Normal, Device> NormalType3;

	public:
		// ============================== Constructor ==============================

		// Create the identity transform
		M_HOST_DEVICE TTransform() : transform(Matrix4::identity()), inv(Matrix4::identity()) {}

		// Create a new transform instance for the given matrix
		M_HOST_DEVICE explicit TTransform(const Matrix4& transform) : transform(transform) {
			bool valid = true;
			inv = transform.inverse(valid);
		}

		// Create a new transform instance for the given matrix and its inverse
		M_HOST_DEVICE TTransform(Matrix4 transform, Matrix4 inv): transform(transform), inv(inv) {}

		// ============================== Function ==============================

		// Return the underlying matrix
		M_HOST_DEVICE [[nodiscard]] Matrix4 get_transform() const {
			return transform;
		}

		M_HOST_DEVICE [[nodiscard]] Matrix4& get_transform() {
			return transform;
		}

		// Return the inverse of the underlying matrix
		M_HOST_DEVICE [[nodiscard]] Matrix4 get_inverse() const {
			return inv;
		}

		M_HOST_DEVICE [[nodiscard]] Matrix4& get_inverse() {
			return inv;
		}

		// Return the inverse transformation
		M_HOST_DEVICE [[nodiscard]] TTransform inverse() const {
			return {inv, transform};
		}

		// Concatenate with another transform
		M_HOST_DEVICE TTransform operator*(const TTransform& t) const {
			return {transform * t.transform, t.inv * inv};
		}

		// Apply the homogeneous transformation to a 3D vector
		M_HOST_DEVICE VectorType3 operator*(const VectorType3& v) const {
			TMatrix<Scalar, 3, 3, Device> sub_matrix = transform.template view<3, 3>(0, 0);
			return sub_matrix * v;
		}

		// Apply the homogeneous transformation to a 3D normal
		M_HOST_DEVICE NormalType3 operator*(const NormalType3& n) const {
			TMatrix<Scalar, 3, 3, Device> sub_matrix = inv.template view<3, 3>(0, 0);
			return sub_matrix.transpose() * n;
		}

		// Transform a point by an arbitrary matrix in homogeneous coordinates
		M_HOST_DEVICE PointType3 operator*(const PointType3& p) const {
			VectorType4 homo_vector = VectorType4(p(0), p(1), p(2), static_cast<Scalar>(1));
			VectorType4 result = transform * homo_vector;
			VectorType3 sub_vector = result.template view<3>(0) / result(3);
			return PointType3(sub_vector(0), sub_vector(1), sub_vector(2));
		}

		// Apply the homogeneous transformation to a ray
		M_HOST_DEVICE TRay<Scalar, 3, Device> operator*(const TRay<Scalar, 3, Device>& r) const {
			bool valid = true;
			return {operator*(r.o()), operator*(r.d()), r.min_t(), r.max_t(), valid};
		}

		// Return a human-readable string summary of this frame
		[[nodiscard]] std::string to_string() {
			static_assert(Device == DeviceType::CPU, "Device is not on CPU");
			std::ostringstream oss;
			oss << "Transform[\n" << transform.to_string() << "]\n";
			return oss.str();
		}

		M_HOST_DEVICE static TTransform translate(Scalar x, Scalar y, Scalar z) {
			Matrix4 transform_matrix = Matrix4::identity();
			transform_matrix(0, 3) = x;
			transform_matrix(1, 3) = y;
			transform_matrix(2, 3) = z;
			Matrix4 inverse_matrix = Matrix4::identity();
			inverse_matrix(0, 3) = -x;
			inverse_matrix(1, 3) = -y;
			inverse_matrix(2, 3) = -z;
			return {transform_matrix, inverse_matrix};
		}

		M_HOST_DEVICE static TTransform translate(VectorType3 vector) {
			Matrix4 transform_matrix = Matrix4::identity();
			transform_matrix(0, 3) = vector(0);
			transform_matrix(1, 3) = vector(1);
			transform_matrix(2, 3) = vector(2);
			Matrix4 inverse_matrix = Matrix4::identity();
			inverse_matrix(0, 3) = -vector(0);
			inverse_matrix(1, 3) = -vector(1);
			inverse_matrix(2, 3) = -vector(2);
			return {transform_matrix, inverse_matrix};
		}

		M_HOST_DEVICE static TTransform scale(Scalar x, Scalar y, Scalar z, bool& valid) {
			check_zero(x, valid);
			check_zero(y, valid);
			check_zero(z, valid);
			Matrix4 transform_matrix = Matrix4::identity();
			transform_matrix(0, 0) = x;
			transform_matrix(1, 1) = y;
			transform_matrix(2, 2) = z;
			Matrix4 inverse_matrix = Matrix4::identity();
			inverse_matrix(0, 0) = 1 / x;
			inverse_matrix(1, 1) = 1 / y;
			inverse_matrix(2, 2) = 1 / z;
			return {transform_matrix, inverse_matrix};
		}

		M_HOST_DEVICE static TTransform scale(VectorType3 vector, bool& valid) {
			check_zero(vector(0), valid);
			check_zero(vector(1), valid);
			check_zero(vector(2), valid);
			Matrix4 transform_matrix = Matrix4::identity();
			transform_matrix(0, 0) = vector(0);
			transform_matrix(1, 1) = vector(1);
			transform_matrix(2, 2) = vector(2);
			Matrix4 inverse_matrix = Matrix4::identity();
			inverse_matrix(0, 0) = 1 / vector(0);
			inverse_matrix(1, 1) = 1 / vector(1);
			inverse_matrix(2, 2) = 1 / vector(2);
			return {transform_matrix, inverse_matrix};
		}

		M_HOST_DEVICE static TTransform rotate(Scalar radiance, VectorType3 axis, bool& valid) {
			axis = axis.norm(valid);
			Scalar cos_theta = cos(radiance);
			Scalar sin_theta = sin(radiance);
			Scalar one_minus_cos_theta = static_cast<Scalar>(1) - cos_theta;

			// Components of the axis
			Scalar x = axis(0);
			Scalar y = axis(1);
			Scalar z = axis(2);

			Matrix3 k = Matrix3::zero();
			k(0, 0) = 0;
			k(0, 1) = -z;
			k(0, 2) = y;
			k(1, 0) = z;
			k(1, 1) = 0;
			k(1, 2) = -x;
			k(2, 0) = -y;
			k(2, 1) = x;
			k(2, 2) = 0;

			Matrix3 rotate_matrix = Matrix3::identity() + k * sin_theta + k * k * one_minus_cos_theta;

			// Construct rotation matrix using Rodrigues' formula
			Matrix4 transform_matrix = Matrix4::identity();
			transform_matrix(0, 0) = rotate_matrix(0, 0);
			transform_matrix(0, 1) = rotate_matrix(0, 1);
			transform_matrix(0, 2) = rotate_matrix(0, 2);
			transform_matrix(1, 0) = rotate_matrix(1, 0);
			transform_matrix(1, 1) = rotate_matrix(1, 1);
			transform_matrix(1, 2) = rotate_matrix(1, 2);
			transform_matrix(2, 0) = rotate_matrix(2, 0);
			transform_matrix(2, 1) = rotate_matrix(2, 1);
			transform_matrix(2, 2) = rotate_matrix(2, 2);

			// The inverse of a rotation matrix is its transpose
			Matrix4 inverse_matrix = transform_matrix.transpose();

			return {transform_matrix, inverse_matrix};
		}

		M_HOST_DEVICE static TTransform look_at(VectorType3 origin, VectorType3 target, VectorType3 up,
		                                        bool& valid) {
			// Calculate the forward direction (camera looks towards -Z by default in right-handed coordinate system)
			VectorType3 forward = (target - origin).norm(valid);
			if (!valid) return TTransform();

			// Calculate the right direction
			VectorType3 left = up.norm(valid).cross(forward).norm(valid);
			if (!valid) return TTransform();

			// Recalculate up to ensure orthogonality
			VectorType3 adjusted_up = forward.cross(left).norm(valid);
			if (!valid) return TTransform();

			// Construct the look-at matrix
			Matrix4 rotate_matrix = Matrix4::identity();
			rotate_matrix(0, 0) = left(0);
			rotate_matrix(1, 0) = left(1);
			rotate_matrix(2, 0) = left(2);
			rotate_matrix(0, 1) = adjusted_up(0);
			rotate_matrix(1, 1) = adjusted_up(1);
			rotate_matrix(2, 1) = adjusted_up(2);
			rotate_matrix(0, 2) = forward(0);
			rotate_matrix(1, 2) = forward(1);
			rotate_matrix(2, 2) = forward(2);
			rotate_matrix(0, 3) = origin(0);
			rotate_matrix(1, 3) = origin(1);
			rotate_matrix(2, 3) = origin(2);

			Matrix4 transform_matrix = rotate_matrix;

			// The inverse of a look-at transformation is itself a look-at in reverse
			Matrix4 inverse_matrix = transform_matrix.inverse(valid);

			return {transform_matrix, inverse_matrix};
		}

		[[nodiscard]] std::string to_string() const {
			std::ostringstream oss;
			oss << "Transform[\n  matrix = " << indent(transform.to_string()) << ",\n]";
			return oss.str();
		}

	protected:
		Matrix4 transform;
		Matrix4 inv;

	private:
		M_HOST_DEVICE static void check_zero(Scalar& scalar, bool& valid) {
			if (scalar == static_cast<Scalar>(0)) {
				scalar += static_cast<Scalar>(M_EPSILON);
				valid &= false;
			}
		}
	};

M_NAMESPACE_END
