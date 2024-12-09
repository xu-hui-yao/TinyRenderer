#pragma once

#include <core/spectrum.h>
#include <core/common.h>
#include <core/array.h>
#include <components/bitmap.h>
#include <mutex>
#include <memory>
#include <vector>

#define M_BLOCK_SIZE 32 /* Block size used for parallelization */

M_NAMESPACE_BEGIN
	class ImageBlock {
	public:
		ImageBlock(const Vector2i& size, const std::shared_ptr<ReconstructionFilter>& filter);
		~ImageBlock();

		void set_offset(const Point2i& offset);
		[[nodiscard]] const Point2i& get_offset() const;

		void set_size(const Point2i& size);
		[[nodiscard]] const Vector2i& get_size() const;

		[[nodiscard]] int get_border_size() const;

		void clear(const Color3f& color = Color3f(0)) const;
		void put(const Point2f& pos, const Color3f& value);
		void put(const ImageBlock& block);

		void lock() const;
		void unlock() const;

		std::shared_ptr<Bitmap> to_bitmap() const;
		void from_bitmap(const Bitmap& bitmap);

		[[nodiscard]] std::string to_string() const;

	private:
		Point2i m_offset;
		Vector2i m_size;
		int m_border_size = 0;

		std::unique_ptr<float[]> m_filter;
		std::vector<float> m_weights_x, m_weights_y;
		float m_filter_radius = 0;
		float m_lookup_factor = 0;

		mutable std::mutex m_mutex;
		std::unique_ptr<Color4f[]> m_color;
		int m_rows, m_cols;
	};

	class BlockGenerator {
	public:
		BlockGenerator(const Vector2i& size, int block_size);

		bool next(ImageBlock& block);

		bool has_next() const;

		[[nodiscard]] int get_block_count() const;

	private:
		enum EDirection { ERight = 0, EDown, ELeft, EUp };

		Point2i m_block;
		Vector2i m_num_blocks, m_size;
		int m_block_size, m_num_steps, m_blocks_left, m_steps_left, m_direction;
		mutable std::mutex m_mutex;
	};

M_NAMESPACE_END
