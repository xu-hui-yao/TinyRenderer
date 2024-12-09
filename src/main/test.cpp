#include <iostream>
#include <core/matrix.h>
#include <core/array.h>
#include <core/common.h>
#include <core/transform.h>
#include <core/frame.h>
#include <core/bounding_box.h>
#include <core/tensor.h>
#include <parse/parser.h>

using namespace tiny_renderer;

void core_test() {
	bool is_valid = true;

	// 测试TMatrix的基本操作
	std::cout << "is valid: " << is_valid << std::endl;
	float data[16]{6, 2, -1, 2, 5, 0, 4, -10, -9, -4, 2, 0, 6, 8, 2, -10};
	TMatrix<float, 4, 4, DeviceType::CPU> matrix_a(data, DeviceType::CPU);
	std::cout << "Matrix A:\n" << matrix_a.to_string() << "\n\n";

	TMatrix<float, 4, 4, DeviceType::CPU> matrix_b = matrix_a.inverse(is_valid);
	std::cout << "Matrix A Inverse (B):\n" << matrix_b.to_string() << "\n\n";

	std::cout << "Matrix A * Matrix B (should be identity):\n" << (matrix_a * matrix_b).to_string() << "\n\n";

	std::cout << "Matrix A Transpose:\n" << matrix_a.transpose().to_string() << "\n\n";

	// TTensor
	float tensor_data[24]{
		1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f,
		7.0f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f,
		13.0f, 14.0f, 15.0f, 16.0f, 17.0f, 18.0f,
		19.0f, 20.0f, 21.0f, 22.0f, 23.0f, 24.0f
	};

	TTensor<float, DeviceType::CPU> tensor_a(2, 2, 6, tensor_data, DeviceType::CPU);
	std::cout << "Tensor A:\n" << tensor_a.to_string() << "\n\n";
	TTensor<float, DeviceType::CPU> tensor_b = tensor_a * 2.0f;
	std::cout << "Tensor B (A * 2):\n" << tensor_b.to_string() << "\n\n";
	TTensor<float, DeviceType::CPU> tensor_c = tensor_a + 3.0f;
	std::cout << "Tensor C (A + 3):\n" << tensor_c.to_string() << "\n\n";
	TTensor<float, DeviceType::CPU> tensor_d = tensor_c - tensor_b;
	std::cout << "Tensor D (C - B):\n" << tensor_d.to_string() << "\n\n";
	TTensor tensor_copy(tensor_a);
	std::cout << "Tensor Copy (Deep Copy of A):\n" << tensor_copy.to_string() << "\n\n";
	tensor_a += 1.0f;
	std::cout << "Tensor A:\n" << tensor_a.to_string() << "\n\n";

	// 测试TVector的基本操作
	std::cout << "is valid: " << is_valid << std::endl;
	Vector3f vec_a(3.0f, -2.0f, 5.0f);
	Vector3f vec_b(1.0f, 4.0f, -2.0f);

	std::cout << "Vector A: (" << vec_a.x() << ", " << vec_a.y() << ", " << vec_a.z() << ")\n";
	std::cout << "Vector B: (" << vec_b.x() << ", " << vec_b.y() << ", " << vec_b.z() << ")\n";

	auto vec_cross = vec_a.cross(vec_b);
	std::cout << "Vector A cross Vector B: (" << vec_cross.x() << ", " << vec_cross.y() << ", " << vec_cross.z() <<
		")\n";

	float dot_product = vec_a.dot(vec_b);
	std::cout << "Vector A dot Vector B: " << dot_product << "\n";

	auto vec_norm = vec_a.norm(is_valid);
	std::cout << "Vector A normalized: (" << vec_norm.x() << ", " << vec_norm.y() << ", " << vec_norm.z() << ")\n";

	// 测试TPoint的基本操作
	std::cout << "is valid: " << is_valid << std::endl;
	Point3f point_a(2.0f, 3.0f, 5.0f);
	Point3f point_b(1.0f, -4.0f, 2.0f);

	std::cout << "Point A: (" << point_a.x() << ", " << point_a.y() << ", " << point_a.z() << ")\n";
	std::cout << "Point B: (" << point_b.x() << ", " << point_b.y() << ", " << point_b.z() << ")\n";

	auto point_diff = point_a - point_b;
	std::cout << "Vector from Point A to Point B: (" << point_diff.x() << ", " << point_diff.y() << ", " << point_diff.
		z() << ")\n";

	// 测试TNormal的基本操作
	std::cout << "is valid: " << is_valid << std::endl;
	Normal3f normal_a(1.0f, 1.0f, 1.0f);
	std::cout << "Normal A: (" << normal_a.x() << ", " << normal_a.y() << ", " << normal_a.z() << ")\n";

	Normal3f normal_b(0.2f);
	std::cout << "Normal A: (" << normal_b.x() << ", " << normal_b.y() << ", " << normal_b.z() << ")\n";

	// 测试TRay的基本操作
	std::cout << "is valid: " << is_valid << std::endl;
	Point3f origin(1.0f, 2.0f, 3.0f);
	Vector3f direction(0.1f, 1.0f, 0.1f);

	Ray3f ray_default;
	std::cout << "Default Ray: Origin = (" << ray_default.o().x() << ", " << ray_default.o().y() << ", " << ray_default.
		o().z() << "), "
		<< "Direction = (" << ray_default.d().x() << ", " << ray_default.d().y() << ", " << ray_default.d().z() <<
		")\n";

	Ray3f ray(origin, direction, is_valid);
	std::cout << "Ray with Origin and Direction:\n"
		<< "Origin = (" << ray.o().x() << ", " << ray.o().y() << ", " << ray.o().z() << "), "
		<< "Direction = (" << ray.d().x() << ", " << ray.d().y() << ", " << ray.d().z() << ")\n";

	std::cout << "Direction Reciprocal = (" << ray.d_reciprocal().x() << ", " << ray.d_reciprocal().y() << ", " << ray.
		d_reciprocal().z() << ")\n";

	float t = 2.0f;
	auto point_on_ray = ray(t);
	std::cout << "Point on Ray at t = " << t << ": (" << point_on_ray.x() << ", " << point_on_ray.y() << ", " <<
		point_on_ray.z() << ")\n";

	ray.min_t() = 1.0f;
	ray.max_t() = 5.0f;
	std::cout << "Ray min_t: " << ray.min_t() << ", max_t: " << ray.max_t() << "\n";

	auto reversed_ray = ray.reverse();
	std::cout << "Reversed Ray:\n"
		<< "Origin = (" << reversed_ray.o().x() << ", " << reversed_ray.o().y() << ", " << reversed_ray.o().z() << "), "
		<< "Direction = (" << reversed_ray.d().x() << ", " << reversed_ray.d().y() << ", " << reversed_ray.d().z() <<
		")\n";

	// 测试Transform
	std::cout << "is valid: " << is_valid << std::endl;
	float data1[16]{1.0f, 0.0f, 0.0f, 2.0f, 0.0f, 1.0f, 0.0f, 3.0f, 0.0f, 0.0f, 1.0f, 4.0f, 0.0f, 0.0f, 0.0f, 1.0f};
	TMatrix<float, 4, 4, DeviceType::CPU> matrix(data1, DeviceType::CPU);

	TTransform<float, DeviceType::CPU> transform(matrix);

	std::cout << "Transform Matrix:\n" << transform.get_transform().to_string() << "\n\n";
	std::cout << "Inverse Transform Matrix:\n" << transform.get_inverse().to_string() << "\n\n";

	Vector3f vector(1.0f, 0.1f, 0.1f);
	auto transformed_vector = transform * vector;
	std::cout << "Original Vector: (" << vector.x() << ", " << vector.y() << ", " << vector.z() << ")\n";
	std::cout << "Transformed Vector: (" << transformed_vector.x() << ", "
		<< transformed_vector.y() << ", " << transformed_vector.z() << ")\n\n";

	Normal3f normal(0.0f, 1.0f, 0.0f);
	auto transformed_normal = transform * normal;
	std::cout << "Original Normal: (" << normal.x() << ", " << normal.y() << ", " << normal.z() << ")\n";
	std::cout << "Transformed Normal: (" << transformed_normal.x() << ", "
		<< transformed_normal.y() << ", " << transformed_normal.z() << ")\n\n";

	Point3f point(1.0f, 2.0f, 3.0f);
	auto transformed_point = transform * point;
	std::cout << "Original Point: (" << point.x() << ", " << point.y() << ", " << point.z() << ")\n";
	std::cout << "Transformed Point: (" << transformed_point.x() << ", "
		<< transformed_point.y() << ", " << transformed_point.z() << ")\n\n";

	Ray3f ray1(point, vector, is_valid);
	auto transformed_ray = transform * ray1;
	std::cout << "Original Ray Origin: (" << ray1.o().x() << ", " << ray1.o().y() << ", " << ray1.o().z() << ")\n";
	std::cout << "Original Ray Direction: (" << ray1.d().x() << ", " << ray1.d().y() << ", " << ray1.d().z() << ")\n";
	std::cout << "Transformed Ray Origin: (" << transformed_ray.o().x() << ", "
		<< transformed_ray.o().y() << ", " << transformed_ray.o().z() << ")\n";
	std::cout << "Transformed Ray Direction: (" << transformed_ray.d().x() << ", "
		<< transformed_ray.d().y() << ", " << transformed_ray.d().z() << ")\n\n";

	auto inverse_transform = transform.inverse();
	std::cout << "Inverse Transform:\n" << inverse_transform.get_transform().to_string() << "\n";

	// 测试Frame
	std::cout << "is valid: " << is_valid << std::endl;
	Normal3f normal_vector(0.2f, 1.0f, 1.0f);
	TFrame<float, DeviceType::CPU> frame(normal_vector);
	std::cout << "Frame:\n" << frame.to_string() << "\n";
	Vector3f world_vector(1.0f, 1.0f, 0.3f);
	auto local_vector = frame.to_local(world_vector);
	std::cout << "World Vector: (" << world_vector.x() << ", " << world_vector.y() << ", " << world_vector.z() << ")\n";
	std::cout << "Local Vector: (" << local_vector.x() << ", " << local_vector.y() << ", " << local_vector.z() <<
		")\n\n";
	auto reconstructed_world_vector = frame.to_world(local_vector);
	std::cout << "Reconstructed World Vector: (" << reconstructed_world_vector.x() << ", "
		<< reconstructed_world_vector.y() << ", " << reconstructed_world_vector.z() << ")\n\n";
	std::cout << "Cosine Theta: " << TFrame<float, DeviceType::CPU>::cos_theta(local_vector, is_valid) << "\n";
	std::cout << "Sine Theta: " << TFrame<float, DeviceType::CPU>::sin_theta(local_vector, is_valid) << "\n";
	std::cout << "Tangent Theta: " << TFrame<float, DeviceType::CPU>::tan_theta(local_vector, is_valid) << "\n";
	std::cout << "Cosine Phi: " << TFrame<float, DeviceType::CPU>::cos_phi(local_vector, is_valid) << "\n";
	std::cout << "Sine Phi: " << TFrame<float, DeviceType::CPU>::sin_phi(local_vector, is_valid) << "\n";
	std::cout << "Cosine Phi Squared: " << TFrame<float, DeviceType::CPU>::cos_phi2(local_vector, is_valid) << "\n";
	std::cout << "Sine Phi Squared: " << TFrame<float, DeviceType::CPU>::sin_phi2(local_vector, is_valid) << "\n\n";
	TFrame<float, DeviceType::CPU> frame2(normal_vector);
	std::cout << "Is Frame equal to Frame2? " << (frame == frame2 ? "Yes" : "No") << "\n";
	frame2 = TFrame<float, DeviceType::CPU>(Normal3f(1.0f, 0.01f, 0.01f));
	std::cout << "Is Frame still equal to Frame2 after modification? " << (frame != frame2 ? "No" : "Yes") << "\n";

	// 测试BoundingBox
	std::cout << "is valid: " << is_valid << std::endl;
	Point3f point1(0.0f, 0.0f, 0.0f);
	Point3f point2(2.0f, 2.0f, 2.0f);
	BoundingBox3f bbox1(point1, point2);
	std::cout << "BoundingBox Min: (" << bbox1.get_min().x() << ", " << bbox1.get_min().y() << ", " << bbox1.get_min().
		z() << ")\n";
	std::cout << "BoundingBox Max: (" << bbox1.get_max().x() << ", " << bbox1.get_max().y() << ", " << bbox1.get_max().
		z() << ")\n";
	Point3f point_inside(1.0f, 1.0f, 1.0f);
	Point3f point_outside(3.0f, 3.0f, 3.0f);
	std::cout << "Contains Point Inside: " << bbox1.contains(point_inside) << "\n";
	std::cout << "Contains Point Outside: " << bbox1.contains(point_outside) << "\n";
	BoundingBox3f bbox2(Point3f(0.5f, 0.5f, 0.5f),
	                                                      Point3f(1.5f, 1.5f, 1.5f));
	std::cout << "Contains BBox2: " << bbox1.contains(bbox2) << "\n";
	BoundingBox3f bbox3(Point3f(1.5f, 1.5f, 1.5f),
	                                                      Point3f(3.0f, 3.0f, 3.0f));
	std::cout << "Overlaps BBox3: " << bbox1.overlaps(bbox3) << "\n";

	bbox1.expand_by(point_outside);
	std::cout << "BoundingBox Max after expansion: (" << bbox1.get_max().x() << ", " << bbox1.get_max().y() << ", " <<
		bbox1.get_max().z() << ")\n";

	std::cout << "Surface Area: " << bbox1.get_surface_area() << "\n";
	std::cout << "Volume: " << bbox1.get_volume() << "\n";

	std::cout << "Distance to Point Outside: " << bbox1.distance_to(point_outside) << "\n";

	std::cout << "Major Axis: " << bbox1.get_major_axis() << "\n";
	std::cout << "Minor Axis: " << bbox1.get_minor_axis() << "\n";

	std::cout << "Is Valid: " << bbox1.is_valid() << "\n";
	std::cout << "Is Point: " << bbox1.is_point() << "\n";

	std::cout << "is valid: " << is_valid << std::endl;
}

void read_test() {
	std::string scene_name = "assets/bunny/bunny.xml";

	try {
		auto root(load_from_xml(scene_name));
		if (root->get_class_type() == Object::EScene) {
			std::cout << "initialize success" << std::endl;
		}
	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
	}
}

int main() {
	core_test();
	// read_test();
	return 0;
}
