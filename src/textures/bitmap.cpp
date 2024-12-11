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
		float* exr_data;
		int width, height;
		const char* err = nullptr;

		int ret = LoadEXR(&exr_data, &width, &height, filename.c_str(), &err);
		if (ret != TINYEXR_SUCCESS) {
			std::string error_message = err ? std::string(err) : "Unknown error";
			if (err) FreeEXRErrorMessage(err);
			throw std::runtime_error("Failed to load EXR file: " + error_message);
		}

		m_data = std::make_shared<TensorXf>(height, width, 3, 0.0f);
		for (int y = 0; y < height; ++y) {
			for (int x = 0; x < width; ++x) {
				for (int c = 0; c < 3; ++c) {
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

		m_data = std::make_shared<TensorXf>(height, width, 3, 0.0f);
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

	Color3f Bitmap::eval(const SurfaceIntersection3f& si, bool active) {
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
	}

	std::string Bitmap::to_string() const {
		return "Bitmap\n";
	}


	REGISTER_CLASS(Bitmap, "bitmap")

M_NAMESPACE_END
