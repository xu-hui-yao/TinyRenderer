#include <components/bitmap.h>
#ifndef STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#endif
#include <stb_image.h>
#ifndef STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#endif
#include <stb_image_write.h>
#ifndef TINYEXR_IMPLEMENTATION
#define TINYEXR_IMPLEMENTATION
#define TINYEXR_USE_MINIZ 0
#define TINYEXR_USE_STB_ZLIB 1
#undef max
#undef min
#endif
#include <iostream>
#include <tinyexr.h>
#include <core/intersection.h>

M_NAMESPACE_BEGIN
	Bitmap::Bitmap(const PropertyList& properties) : Texture(true) {
		filesystem::path filename = get_file_resolver()->resolve(
			filesystem::path(properties.get_string("filename")));

		std::string ext = filename.extension();
		std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

		if (ext == "exr") {
			load_exr(filename.str());
		} else if (ext == "png" || ext == "jpg" || ext == "jpeg" || ext == "bmp" || ext == "tga") {
			load_image(filename.str());
		} else {
			throw std::runtime_error("Unsupported image format: " + ext);
		}

		m_name = properties.get_string("name", "bitmap");
	}

	Bitmap::Bitmap(int height, int width, int channels) : Texture(true) {
		m_data = std::make_shared<TensorXf>(height, width, channels);
	}

	void Bitmap::add_child(const std::shared_ptr<Object>& child) {}

	void Bitmap::construct() {
#ifdef M_DEBUG
		std::cout << "Construct" << class_type_name(get_class_type()) << std::endl;
#endif
		if (m_data == nullptr) {
			throw std::runtime_error("Bitmap data is not loaded.");
		}
		if (m_data->get_channels() != 1 && m_data->get_channels() != 3) {
			throw std::runtime_error(
				"Bitmap data does not support channels: " + std::to_string(m_data->get_channels()));
		}
	}

	Bitmap::~Bitmap() = default;

	[[nodiscard]] std::shared_ptr<TensorXf> Bitmap::get_data() {
		return m_data;
	}

	float Bitmap::operator()(int row, int col, int channel) const {
		return m_data->operator()(row, col, channel);
	}

	float& Bitmap::operator()(int row, int col, int channel) {
		return m_data->operator()(row, col, channel);
	}

	void Bitmap::save_exr(const std::string& filename) const {
		// Convert the TensorXf data to a float array suitable for EXR
		int width = m_data->get_cols();
		int height = m_data->get_rows();

		EXRHeader header;
		InitEXRHeader(&header);

		EXRImage image;
		InitEXRImage(&image);

		image.num_channels = 3;

		std::vector<float> images[3];
		images[0].resize(width * height);
		images[1].resize(width * height);
		images[2].resize(width * height);

		for (int i = 0; i < height; i++) {
			for (int j = 0; j < width; j++) {
				images[0][i * width + j] = m_data->operator()(i, j, 0);
				images[1][i * width + j] = m_data->operator()(i, j, 1);
				images[2][i * width + j] = m_data->operator()(i, j, 2);
			}
		}

		float* image_ptr[3];
		image_ptr[0] = &images[2].at(0); // B
		image_ptr[1] = &images[1].at(0); // G
		image_ptr[2] = &images[0].at(0); // R

		image.images = reinterpret_cast<unsigned char**>(image_ptr);
		image.width = width;
		image.height = height;
		header.num_channels = 3;
		header.channels = static_cast<EXRChannelInfo*>(malloc(sizeof(EXRChannelInfo) * header.num_channels));
		// Must be (A)BGR order, since most of EXR viewers expect this channel order.
		strncpy_s(header.channels[0].name, "B", 255);
		header.channels[0].name[strlen("B")] = '\0';
		strncpy_s(header.channels[1].name, "G", 255);
		header.channels[1].name[strlen("G")] = '\0';
		strncpy_s(header.channels[2].name, "R", 255);
		header.channels[2].name[strlen("R")] = '\0';

		header.pixel_types = static_cast<int*>(malloc(sizeof(int) * header.num_channels));
		header.requested_pixel_types = static_cast<int*>(malloc(sizeof(int) * header.num_channels));
		for (int i = 0; i < header.num_channels; i++) {
			header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
			header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF;
			// pixel type of output image to be stored in .EXR
		}

		const char* err = nullptr;
		int ret = SaveEXRImageToFile(&image, &header, filename.c_str(), &err);
		if (ret != TINYEXR_SUCCESS) {
			std::cerr << "Save EXR err: " << err << std::endl;
			FreeEXRErrorMessage(err); // free's buffer for an error message
			return;
		}
		free(header.channels);
		free(header.pixel_types);
		free(header.requested_pixel_types);
	}

	void Bitmap::save_png(const std::string& filename) const {
		int width = m_data->get_cols();
		int height = m_data->get_rows();

		// Convert the float data to 8-bit (0-255) data for saving as PNG
		std::vector<unsigned char> img_data(width * height * 3);
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				Color3f color;
				color(0) = m_data->operator()(y, x, 0);
				color(1) = m_data->operator()(y, x, 1);
				color(2) = m_data->operator()(y, x, 2);
				auto srgb = color.to_srgb();
				for (int c = 0; c < 3; ++c) {
					// Convert to 8-bit by clamping to 0-255 range
					int index = (y * width + x) * 3 + c;
					auto value = static_cast<unsigned char>(clamp(srgb(c) * 255.0f, 0.0f, 255.0f));
					img_data[index] = value;
				}
			}
		}

		// Save as PNG
		int ret = stbi_write_png(filename.c_str(), width, height, 3, img_data.data(), width * 3);
		if (ret == 0) {
			throw std::runtime_error("Failed to save PNG file: " + filename);
		}
	}

	void Bitmap::load_exr(const std::string& filename) {
		EXRVersion exr_version;
		EXRHeader exr_header;
		const char* err = nullptr;
		InitEXRHeader(&exr_header);
		if (ParseEXRHeaderFromFile(&exr_header, &exr_version, filename.c_str(), &err) != TINYEXR_SUCCESS) {
			std::string error_message = err ? std::string(err) : "Unknown error";
			if (err) FreeEXRErrorMessage(err);
			throw std::runtime_error("Failed to parse EXR header: " + error_message);
		}

		float* exr_data;
		int width, height;
		int channels = exr_header.num_channels;

		FreeEXRHeader(&exr_header);

		int ret = LoadEXR(&exr_data, &width, &height, filename.c_str(), &err);
		if (ret != TINYEXR_SUCCESS) {
			std::string error_message = err ? std::string(err) : "Unknown error";
			if (err) FreeEXRErrorMessage(err);
			throw std::runtime_error("Failed to load EXR file: " + error_message);
		}

		m_data = std::make_shared<TensorXf>(height, width, channels, 0.0f);
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				for (int c = 0; c < channels; ++c) {
					m_data->operator()(y, x, c) = exr_data[(y * width + x) * 3 + c];
				}
			}
		}
		free(exr_data);
	}

	void Bitmap::load_image(const std::string& filename) {
		int width, height, channels;
		unsigned char* img_data = stbi_load(filename.c_str(), &width, &height, &channels, 0);
		if (!img_data) {
			throw std::runtime_error("Failed to load image file: " + std::string(stbi_failure_reason()));
		}

		m_data = std::make_shared<TensorXf>(height, width, channels, 0.0f);
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				for (int c = 0; c < channels; ++c) {
					m_data->operator()(y, x, c) = static_cast<float>(img_data[(y * width + x) * channels + c]) / 255.0f;
				}
			}
		}
		stbi_image_free(img_data);
	}

	int Bitmap::get_cols() const {
		return m_data->get_cols();
	}

	int Bitmap::get_rows() const {
		return m_data->get_rows();
	}

	Color3f Bitmap::eval(const SurfaceIntersection3f& si, bool& active) {
		if (m_data->get_channels() == 3) {
			Point2f uv = si.uv;

			// Ensure uv coordinates are within the [0, 1] range
			uv.x() = clamp(uv.x(), 0.0f, 1.0f);
			uv.y() = clamp(uv.y(), 0.0f, 1.0f);

			// Scale uv to the texture dimensions
			float x = uv.x() * (static_cast<float>(m_data->get_cols()) - 1);
			float y = uv.y() * (static_cast<float>(m_data->get_rows()) - 1);

			// Calculate the integer and fractional parts
			int x0 = static_cast<int>(floor(x));
			int y0 = static_cast<int>(floor(y));
			int x1 = clamp(x0 + 1, 0, m_data->get_cols() - 1);
			int y1 = clamp(y0 + 1, 0, m_data->get_rows() - 1);

			float fx = x - static_cast<float>(x0);
			float fy = y - static_cast<float>(y0);

			// Perform bi-linear interpolation
			Color3f c00({m_data->operator()(y0, x0, 0), m_data->operator()(y0, x0, 1), m_data->operator()(y0, x0, 2)});
			Color3f c10({m_data->operator()(y0, x1, 0), m_data->operator()(y0, x1, 1), m_data->operator()(y0, x1, 2)});
			Color3f c01({m_data->operator()(y1, x0, 0), m_data->operator()(y1, x0, 1), m_data->operator()(y1, x0, 2)});
			Color3f c11({m_data->operator()(y1, x1, 0), m_data->operator()(y1, x1, 1), m_data->operator()(y1, x1, 2)});

			// Interpolate horizontally and vertically
			Color3f c0 = c00 * (1 - fx) + c10 * fx;
			Color3f c1 = c01 * (1 - fx) + c11 * fx;
			Color3f result = c0 * (1 - fy) + c1 * fy;

			return result;
		} else {
			return Color3f(eval_1(si, active));
		}
	}

	float Bitmap::eval_1(const SurfaceIntersection3f& si, bool& active) {
		if (m_data->get_channels() == 1) {
			Point2f uv = si.uv;

			// Ensure uv coordinates are within the [0, 1] range
			uv.x() = clamp(uv.x(), 0.0f, 1.0f);
			uv.y() = clamp(uv.y(), 0.0f, 1.0f);

			// Scale uv to the texture dimensions
			float x = uv.x() * (static_cast<float>(m_data->get_cols()) - 1);
			float y = uv.y() * (static_cast<float>(m_data->get_rows()) - 1);

			// Calculate the integer and fractional parts
			int x0 = static_cast<int>(floor(x));
			int y0 = static_cast<int>(floor(y));
			int x1 = clamp(x0 + 1, 0, m_data->get_cols() - 1);
			int y1 = clamp(y0 + 1, 0, m_data->get_rows() - 1);

			float fx = x - static_cast<float>(x0);
			float fy = y - static_cast<float>(y0);

			// Perform bi-linear interpolation
			float c00 = m_data->operator()(y0, x0, 0);
			float c10 = m_data->operator()(y0, x1, 0);
			float c01 = m_data->operator()(y1, x0, 0);
			float c11 = m_data->operator()(y1, x1, 0);

			// Interpolate horizontally and vertically
			float c0 = c00 * (1 - fx) + c10 * fx;
			float c1 = c01 * (1 - fx) + c11 * fx;
			float result = c0 * (1 - fy) + c1 * fy;

			return result;
		} else {
			return eval(si, active).luminance();
		}
	}

	Vector2f Bitmap::eval_1_grad(const SurfaceIntersection3f& si, bool& active) {
		// Ensure UV coordinates are within [0, 1] range
		Point2f uv = si.uv;
		uv.x() = clamp(uv.x(), 0.0f, 1.0f);
		uv.y() = clamp(uv.y(), 0.0f, 1.0f);

		// Scale UV to texture dimensions
		float x = uv.x() * (static_cast<float>(m_data->get_cols()) - 1);
		float y = uv.y() * (static_cast<float>(m_data->get_rows()) - 1);

		// Integer and fractional parts
		int x0 = static_cast<int>(floor(x));
		int y0 = static_cast<int>(floor(y));
		int x1 = clamp(x0 + 1, 0, m_data->get_cols() - 1);
		int y1 = clamp(y0 + 1, 0, m_data->get_rows() - 1);

		float fx = x - static_cast<float>(x0);
		float fy = y - static_cast<float>(y0);

		// Load values from the bitmap
		float f00 = m_data->operator()(y0, x0, 0);
		float f10 = m_data->operator()(y0, x1, 0);
		float f01 = m_data->operator()(y1, x0, 0);
		float f11 = m_data->operator()(y1, x1, 0);

		// Compute gradients w.r.t. pixel coordinates (x, y)
		Vector2f df_xy;
		df_xy.x() = (1 - fy) * (f10 - f00) + fy * (f11 - f01); // Partial derivative w.r.t. x
		df_xy.y() = (1 - fx) * (f01 - f00) + fx * (f11 - f10); // Partial derivative w.r.t. y

		// Scale gradients to UV space
		auto scale_u = static_cast<float>(m_data->get_cols() - 1);
		auto scale_v = static_cast<float>(m_data->get_rows() - 1);
		Vector2f df_uv(df_xy.x() * scale_u, df_xy.y() * scale_v);

		return df_uv;
	}

	Color3f Bitmap::mean() {
		if (m_data->get_channels() == 1) {
			float sum = 0;
			for (int i = 0; i < m_data->get_cols(); i++) {
				for (int j = 0; j < m_data->get_rows(); j++) {
					sum += m_data->operator()(j, i, 0);
				}
			}
			sum /= static_cast<float>(m_data->get_cols() * m_data->get_rows());
			return Color3f(sum);
		} else {
			Color3f sum(0.0f);
			for (int i = 0; i < m_data->get_cols(); i++) {
				for (int j = 0; j < m_data->get_rows(); j++) {
					sum(0) += m_data->operator()(i, j, 0);
					sum(1) += m_data->operator()(i, j, 1);
					sum(2) += m_data->operator()(i, j, 2);
				}
			}
			sum /= static_cast<float>(m_data->get_cols() * m_data->get_rows());
			return sum;
		}
	}

	std::string Bitmap::to_string() const {
		return "Bitmap";
	}


	REGISTER_CLASS(Bitmap, "bitmap")

M_NAMESPACE_END
