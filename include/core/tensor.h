#pragma once

#include <cassert>
#include <string>
#include <cuda_runtime.h>
#include <core/common.h>

M_NAMESPACE_BEGIN
	template <typename Scalar, DeviceType Device_>
	class TTensor {
	public:
		static constexpr DeviceType Device = Device_;

		// ============================== Constructor ==============================

		// Set matrix to value
		M_HOST_DEVICE explicit TTensor(int rows, int cols, int channels, Scalar value = 0)
			: m_rows(rows), m_cols(cols), m_channels(channels) {
#ifndef __CUDA_ARCH__
			assert(rows > 0 && cols > 0 && channels > 0 && "Rows and columns and channels must be positive integers");
			if constexpr (Device == DeviceType::CPU) {
				m_data = new Scalar[m_rows * m_cols * m_channels];
				std::fill(m_data, m_data + m_rows * m_cols * m_channels, value);
			} else {
				auto* temp = new Scalar[m_rows * m_cols * m_channels];
				std::fill(temp, temp + m_rows * m_cols * m_channels, value);
				check_cuda_error(
					cudaMemcpy(m_data, temp, m_rows * m_cols * m_channels * sizeof(Scalar), cudaMemcpyHostToDevice),
					"Memory copy failed");
				delete[] temp;
			}
#else
			printf("Error: Can not construct Tensor on GPU by device function");
#endif
		}

		// Initialize matrix by pointer
		M_HOST_DEVICE explicit TTensor(int rows, int cols, int channels, Scalar* input_m_data, DeviceType input_device)
			: m_rows(rows), m_cols(cols), m_channels(channels) {
#ifndef __CUDA_ARCH__
			assert(input_m_data != nullptr && "Input data pointer must not be null");
			assert(rows > 0 && cols > 0 && "Rows and columns must be positive integers");
			if constexpr (Device == DeviceType::CPU) {
				m_data = new Scalar[m_rows * m_cols * m_channels];
				if (input_device == DeviceType::CPU) {
					std::copy(input_m_data, input_m_data + m_rows * m_cols * m_channels, m_data);
				} else {
					check_cuda_error(
						cudaMemcpy(m_data, input_m_data, m_rows * m_cols * m_channels * sizeof(Scalar),
						           cudaMemcpyDeviceToHost),
						"Memory copy failed");
				}
			} else {
				assert(input_m_data != nullptr && "Input data pointer must not be null");\
				assert(rows > 0 && cols > 0 && "Rows and columns must be positive integers");
				check_cuda_error(cudaMalloc(&m_data, m_rows * m_cols * m_channels * sizeof(Scalar)));
				if (input_device == DeviceType::CPU) {
					check_cuda_error(
						cudaMemcpy(m_data, input_m_data, m_rows * m_cols * m_channels * sizeof(Scalar),
						           cudaMemcpyHostToDevice),
						"Memory copy failed");
				} else {
					check_cuda_error(
						cudaMemcpy(m_data, input_m_data, m_rows * m_cols * m_channels * sizeof(Scalar),
						           cudaMemcpyDeviceToDevice),
						"Memory copy failed");
				}
			}
#else
			printf("Error: Can not construct Tensor on GPU by device function");
#endif
		}

		// Deep copy constructor
		M_HOST_DEVICE TTensor(const TTensor& other) {
#ifndef __CUDA_ARCH__
			m_rows = other.m_rows;
			m_cols = other.m_cols;
			m_channels = other.m_channels;
			if constexpr (Device == DeviceType::CPU) {
				delete[] m_data;
				m_data = new Scalar[m_rows * m_cols * m_channels];
				std::copy(other.m_data, other.m_data + m_rows * m_cols * m_channels, m_data);
			} else {
				check_cuda_error(cudaFree(m_data), "Memory copy failed");
				check_cuda_error(cudaMalloc(&m_data, m_rows * m_cols * m_channels * sizeof(Scalar)));
				check_cuda_error(
					cudaMemcpy(m_data, other.m_data, m_rows * m_cols * m_channels * sizeof(Scalar),
					           cudaMemcpyDeviceToDevice),
					"Memory copy failed");
			}
#else
			printf("Error: Can not construct Tensor on GPU by device function");
#endif
		}

		// Move constructor
		M_HOST_DEVICE TTensor(TTensor&& other) noexcept
			: m_rows(other.m_rows), m_cols(other.m_cols), m_data(other.m_data), m_channels(other.m_channels) {
			other.m_data = nullptr;
			other.m_rows = 0;
			other.m_cols = 0;
			other.m_channels = 0;
		}

		// Move assignment
		M_HOST_DEVICE TTensor& operator=(TTensor&& other) noexcept {
			if (this == &other) return *this;

			if constexpr (Device == DeviceType::CPU) {
				delete[] m_data;
			} else {
				check_cuda_error(cudaFree(m_data), "CUDA free failed");
			}

			m_data = other.m_data;
			m_rows = other.m_rows;
			m_cols = other.m_cols;
			m_channels = other.m_channels;

			other.m_data = nullptr;
			other.m_rows = 0;
			other.m_cols = 0;
			other.m_channels = 0;

			return *this;
		}

		~TTensor() {
			if constexpr (Device == DeviceType::GPU) {
				check_cuda_error(cudaFree(m_data), "CUDA free failed");
			} else {
				delete[] m_data;
			}
			m_data = nullptr;
		}

		// Deep copy operator
		M_HOST_DEVICE TTensor& operator=(const TTensor& other) {
#ifndef __CUDA_ARCH__
			if (this == &other) {
				return *this;
			}

			if constexpr (Device == DeviceType::CPU) {
				delete[] m_data;
				m_data = new Scalar[m_rows * m_cols * m_channels];
				std::copy(other.m_data, other.m_data + m_rows * m_cols * m_channels, m_data);
			} else {
				check_cuda_error(cudaFree(m_data), "Memory copy failed");
				check_cuda_error(cudaMalloc(&m_data, m_rows * m_cols * m_channels * sizeof(Scalar)));
				check_cuda_error(
					cudaMemcpy(m_data, other.m_data, m_rows * m_cols * m_channels * sizeof(Scalar),
					           cudaMemcpyDeviceToDevice),
					"CUDA memory copy failed");
			}

			return *this;
#else
			printf("Error: Can not construct Tensor on GPU by device function");
#endif
		}

		// ============================== Getter and setter ==============================

		M_HOST_DEVICE [[nodiscard]] int get_rows() const {
			return m_rows;
		}

		M_HOST_DEVICE [[nodiscard]] int get_cols() const {
			return m_cols;
		}

		M_HOST_DEVICE [[nodiscard]] int get_channels() const {
			return m_channels;
		}

		M_HOST_DEVICE [[nodiscard]] Scalar* get_data() const {
			return m_data;
		}

		M_HOST_DEVICE [[nodiscard]] Scalar* & get_data() {
			return m_data;
		}

		M_HOST_DEVICE Scalar operator()(int row, int col, int channel) const {
			check_range(row, col, channel);

			return m_data[(row * m_cols + col) * m_channels + channel];
		}

		M_HOST_DEVICE Scalar& operator()(int row, int col, int channel) {
			check_range(row, col, channel);

			return m_data[(row * m_cols + col) * m_channels + channel];
		}

		M_HOST_DEVICE void set_constant(Scalar value) {
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
		M_HOST_DEVICE TTensor operator*(Scalar scalar) const {
			TTensor result(m_rows, m_cols, m_channels, 0);

			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				result.m_data[i] = m_data[i] * scalar;
			}
			return result;
		}

		M_HOST_DEVICE TTensor& operator*=(Scalar scalar) {
			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				m_data[i] *= scalar;
			}
			return *this;
		}

		// Matrix-Scalar division
		M_HOST_DEVICE TTensor operator/(Scalar scalar) const {
			bool valid = true;
			check_zero(scalar, valid);

			TTensor result(m_rows, m_cols, m_channels, 0);
			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				result.m_data[i] = m_data[i] / scalar;
			}
			return result;
		}

		M_HOST_DEVICE TTensor& operator/=(Scalar scalar) {
			bool valid = true;
			check_zero(scalar, valid);

			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				m_data[i] /= scalar;
			}
			return *this;
		}

		// Matrix-Matrix add
		M_HOST_DEVICE TTensor operator+(const TTensor& other) const {
			assert(
				m_cols == other.get_cols() && m_rows == other.get_rows() && m_channels == other.get_channels() &&
				"Rows and columns must be equal");
			TTensor result(m_rows, m_cols, m_channels, 0);

			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				result.m_data[i] = m_data[i] + other.m_data[i];
			}
			return result;
		}

		M_HOST_DEVICE TTensor& operator+=(const TTensor& other) {
			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				m_data[i] += other.m_data[i];
			}
			return *this;
		}

		// Matrix-Matrix sub
		M_HOST_DEVICE TTensor operator-(const TTensor& other) const {
			assert(
				m_cols == other.get_cols() && m_rows == other.get_rows() && m_channels == other.get_channels() &&
				"Rows and columns and channels must be equal");
			TTensor result(m_rows, m_cols, m_channels, 0);

			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				result.m_data[i] = m_data[i] - other.m_data[i];
			}
			return result;
		}

		M_HOST_DEVICE TTensor& operator-=(const TTensor& other) {
			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				m_data[i] -= other.m_data[i];
			}
			return *this;
		}

		//Matrix-Scalar add
		M_HOST_DEVICE TTensor operator+(Scalar scalar) const {
			TTensor result(m_rows, m_cols, m_channels, 0);

			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				result.m_data[i] = m_data[i] + scalar;
			}
			return result;
		}

		M_HOST_DEVICE TTensor& operator+=(Scalar scalar) {
			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				m_data[i] += scalar;
			}
			return *this;
		}

		// Matrix-Scalar sub
		M_HOST_DEVICE TTensor operator-(Scalar scalar) const {
			TTensor result(m_rows, m_cols, m_channels, 0);

			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				result.m_data[i] = m_data[i] - scalar;
			}
			return result;
		}

		M_HOST_DEVICE TTensor& operator-=(Scalar scalar) {
			for (int i = 0; i < m_rows * m_cols * m_channels; ++i) {
				m_data[i] -= scalar;
			}
			return *this;
		}

		// Operator -
		M_HOST_DEVICE TTensor operator-() const {
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

			std::unique_ptr<Scalar[]> temp(new Scalar[m_rows * m_cols * m_channels]);
			if constexpr (Device == DeviceType::GPU) {
				check_cuda_error(
					cudaMemcpy(temp.get(), m_data, m_rows * m_cols * m_channels * sizeof(Scalar),
					           cudaMemcpyDeviceToHost),
					"Failed to copy tensor data from GPU to host");
				cudaDeviceSynchronize();
			} else {
				std::copy(m_data, m_data + m_rows * m_cols * m_channels, temp.get());
			}

			std::ostringstream oss;
			oss << "[\n";

			for (int i = 0; i < m_rows; ++i) {
				oss << "  [";
				for (int j = 0; j < m_cols; ++j) {
					oss << "[";
					for (int k = 0; k < m_channels; ++k) {
						oss << temp[i * m_cols * m_channels + j * m_channels + k];
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
		Scalar* m_data = nullptr;

		// ============================== Check ==============================

		M_HOST_DEVICE static void check_cuda_error(cudaError_t err, const char* message) {
#ifdef M_DEBUG
			if (err != cudaSuccess) {
				assert(false && message);
			}
#endif
		}

		M_HOST_DEVICE void check_range(int row, int col, int channel) const {
#ifdef M_DEBUG
			assert(row < m_rows && col < m_cols && channel < m_channels && "Index out of range");
#endif
		}

		M_HOST_DEVICE static void check_zero(Scalar& scalar, bool& valid) {
			if (scalar == static_cast<Scalar>(0)) {
				scalar += static_cast<Scalar>(M_EPSILON);
				valid &= false;
			}
		}
	};

M_NAMESPACE_END
