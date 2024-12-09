#pragma once

#include <core/common.h>
#include <cinttypes>
#include <cassert>
#include <algorithm>

M_NAMESPACE_BEGIN
#define PCG32_DEFAULT_STATE  0x853c49e6748fea9bULL
#define PCG32_DEFAULT_STREAM 0xda3e39cb94b95bdbULL
#define PCG32_MULT           0x5851f42d4c957f2dULL

	struct pcg32 {
		M_HOST_DEVICE pcg32() : state(PCG32_DEFAULT_STATE), inc(PCG32_DEFAULT_STREAM) {}

		M_HOST_DEVICE explicit pcg32(uint64_t init_state, uint64_t init_seq = 1u) { seed(init_state, init_seq); }

		/**
		 * \brief Seed the pseudorandom number generator
		 *
		 * Specified in two parts: a state initializer and a sequence selection
		 * constant (a.k.a. stream id)
		 */
		M_HOST_DEVICE void seed(uint64_t init_state, uint64_t init_seq = 1) {
			state = 0U;
			inc = (init_seq << 1u) | 1u;
			next_uint();
			state += init_state;
			next_uint();
		}

		// Generate a uniformly distributed unsigned 32-bit random number
		M_HOST_DEVICE uint32_t next_uint() {
			uint64_t old_state = state;
			state = old_state * PCG32_MULT + inc;
			auto xor_shifted = static_cast<uint32_t>((old_state >> 18u ^ old_state) >> 27u);
			auto rot = static_cast<uint32_t>(old_state >> 59u);
			return xor_shifted >> rot | xor_shifted << (~rot + 1u & 31);
		}

		// Generate a uniformly distributed number, r, where 0 <= r < bound
		M_HOST_DEVICE uint32_t next_uint(uint32_t bound) {
			// To avoid bias, we need to make the range of the RNG a multiple of
			// bound, which we do by dropping output less than a threshold.
			// A naive scheme to calculate the threshold would be to do
			//
			//     uint32_t threshold = 0x100000000ull % bound;
			//
			// but 64-bit div/mod is slower than 32-bit div/mod (especially on
			// 32-bit platforms).  In essence, we do
			//
			//     uint32_t threshold = (0x100000000ull-bound) % bound;
			//
			// because this version will calculate the same modulus, but the LHS
			// value is less than 2^32.

			uint32_t threshold = (~bound + 1u) % bound;

			// Uniformity guarantees that this loop will terminate.  In practice, it
			// should usually terminate quickly; on average (assuming all bounds are
			// equally likely), 82.25% of the time, we can expect it to require just
			// one iteration.  In the worst case, someone passes a bound of 2^31 + 1
			// (i.e., 2147483649), which invalidates almost 50% of the range.  In
			// practice, bounds are typically small and only a tiny amount of the range
			// is eliminated.
			for (;;) {
				uint32_t r = next_uint();
				if (r >= threshold)
					return r % bound;
			}
		}

		// Generate a single precision floating point value on the interval [0, 1)
		M_HOST_DEVICE float next_float() {
			/* Trick from MTGP: generate a uniformly distributed
			   single precision number in [1,2) and subtract 1. */
			union {
				uint32_t u;
				float f;
			} x{};
			x.u = next_uint() >> 9 | 0x3f800000u;
			return x.f - 1.0f;
		}

		/**
		 * \brief Generate a double precision floating point value on the interval [0, 1)
		 *
		 * \remark Since the underlying random number generator produces 32 bit output,
		 * only the first 32 mantissa bits will be filled (however, the resolution is still
		 * finer than in \ref nextFloat(), which only uses 23 mantissa bits)
		 */
		double next_double() {
			/* Trick from MTGP: generate a uniformly distributed
			   double precision number in [1,2) and subtract 1. */
			union {
				uint64_t u;
				double d;
			} x{};
			x.u = static_cast<uint64_t>(next_uint() << 20) | 0x3ff0000000000000ULL;
			return x.d - 1.0;
		}

		/**
		 * \brief Multistep advance function (jump-ahead, jump-back)
		 *
		 * The method used here is based on Brown, "Random Number Generation
		 * with Arbitrary Stride", Transactions of the American Nuclear
		 * Society (Nov. 1994). The algorithm is very similar to fast
		 * exponentiation.
		 */
		M_HOST_DEVICE void advance(int64_t delta_) {
			uint64_t
				cur_mult = PCG32_MULT,
				cur_plus = inc,
				acc_mult = 1u,
				acc_plus = 0u;

			/* Even though delta is an unsigned integer, we can pass a signed
			   integer to go backwards, it just goes "the long way round". */
			auto delta = static_cast<uint64_t>(delta_);

			while (delta > 0) {
				if (delta & 1) {
					acc_mult *= cur_mult;
					acc_plus = acc_plus * cur_mult + cur_plus;
				}
				cur_plus = (cur_mult + 1) * cur_plus;
				cur_mult *= cur_mult;
				delta /= 2;
			}
			state = acc_mult * state + acc_plus;
		}

		/**
		 * \brief Draw uniformly distributed permutation and permute the
		 * given STL container
		 *
		 * From: Knuth, TAoCP Vol. 2 (3rd 3d), Section 3.4.2
		 */
		template <typename Iterator>
		M_HOST_DEVICE void shuffle(Iterator begin, Iterator end) {
			for (Iterator it = end - 1; it > begin; --it)
				std::iter_swap(it, begin + next_uint(static_cast<uint32_t>(it - begin + 1)));
		}

		/// Compute the distance between two PCG32 pseudorandom number generators
		M_HOST_DEVICE int64_t operator-(const pcg32& other) const {
			assert(inc == other.inc);

			uint64_t
				cur_mult = PCG32_MULT,
				cur_plus = inc,
				cur_state = other.state,
				the_bit = 1u,
				distance = 0u;

			while (state != cur_state) {
				if ((state & the_bit) != (cur_state & the_bit)) {
					cur_state = cur_state * cur_mult + cur_plus;
					distance |= the_bit;
				}
				assert((state & the_bit) == (cur_state & the_bit));
				the_bit <<= 1;
				cur_plus = (cur_mult + 1ULL) * cur_plus;
				cur_mult *= cur_mult;
			}

			return static_cast<int64_t>(distance);
		}

		// Equality operator
		M_HOST_DEVICE bool operator==(const pcg32& other) const { return state == other.state && inc == other.inc; }

		// Inequality operator
		M_HOST_DEVICE bool operator!=(const pcg32& other) const { return state != other.state || inc != other.inc; }

		uint64_t state{}; // RNG state.  All values are possible.
		uint64_t inc{}; // Controls which RNG sequence (stream) is selected. Must *always* be odd.
	};

M_NAMESPACE_END
