#pragma once

#include "ExtendedBch256x128Eq3.h"
#include "RiffleRm2SubS19.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

#include <immintrin.h>
#include <cryptoTools/Common/Defines.h>
#if __has_include(<cryptoTools/Crypto/PRNG.h>)
#include <cryptoTools/Crypto/PRNG.h>
#define LIBOTE_RIFFLE_HAS_PRNG 1
#endif

namespace osuCrypto
{
	namespace detail
	{
		struct RiffleFieldProduct128
		{
			u64 low;
			u64 high;
		};

		OC_FORCEINLINE RiffleFieldProduct128 riffleCarrylessMultiply(
			u64 left, u64 right) noexcept
		{
			const auto product = _mm_clmulepi64_si128(
				_mm_cvtsi64_si128(static_cast<long long>(left)),
				_mm_cvtsi64_si128(static_cast<long long>(right)), 0x00);
			return {
				static_cast<u64>(_mm_cvtsi128_si64(product)),
				static_cast<u64>(_mm_extract_epi64(product, 1))
			};
		}

		OC_FORCEINLINE __m256i riffleFieldMultiplyVector(
			__m256i value, u64 coefficient) noexcept
		{
			const auto laneMask = _mm256_set_epi64x(0, -1, 0, -1);
			const auto coefficientVector = _mm256_set_epi64x(
				0, static_cast<long long>(coefficient),
				0, static_cast<long long>(coefficient));
			const auto modulusVector = _mm256_set_epi64x(0, 0x807, 0, 0x807);
			auto product = _mm256_clmulepi64_epi128(value, coefficientVector, 0x00);
			const auto high = _mm256_srli_si256(product, 8);
			auto reduced = _mm256_xor_si256(
				_mm256_and_si256(product, laneMask),
				_mm256_clmulepi64_epi128(high, modulusVector, 0x00));
			const auto secondHigh = _mm256_srli_si256(reduced, 8);
			reduced = _mm256_xor_si256(
				_mm256_and_si256(reduced, laneMask),
				_mm256_clmulepi64_epi128(secondHigh, modulusVector, 0x00));
			return _mm256_and_si256(reduced, laneMask);
		}

		template<int Shift, u64 Mask>
		OC_FORCEINLINE void riffleTranspose64LargeStage(u64* rows) noexcept
		{
			const auto mask = _mm256_set1_epi64x(static_cast<long long>(Mask));
			for (u64 base = 0; base < 64; base += 2 * Shift)
			{
				for (u64 offset = 0; offset < Shift; offset += 4)
				{
					auto low = _mm256_load_si256(
						reinterpret_cast<const __m256i*>(rows + base + offset));
					auto high = _mm256_load_si256(
						reinterpret_cast<const __m256i*>(rows + base + Shift + offset));
					const auto difference = _mm256_and_si256(
						_mm256_xor_si256(_mm256_srli_epi64(low, Shift), high), mask);
					low = _mm256_xor_si256(low, _mm256_slli_epi64(difference, Shift));
					high = _mm256_xor_si256(high, difference);
					_mm256_store_si256(
						reinterpret_cast<__m256i*>(rows + base + offset), low);
					_mm256_store_si256(
						reinterpret_cast<__m256i*>(rows + base + Shift + offset), high);
				}
			}
		}

		OC_FORCEINLINE void riffleTranspose64Avx2(u64* rows) noexcept
		{
			riffleTranspose64LargeStage<32, 0x00000000ffffffffULL>(rows);
			riffleTranspose64LargeStage<16, 0x0000ffff0000ffffULL>(rows);
			riffleTranspose64LargeStage<8, 0x00ff00ff00ff00ffULL>(rows);
			riffleTranspose64LargeStage<4, 0x0f0f0f0f0f0f0f0fULL>(rows);

			const auto mask2 = _mm256_set1_epi64x(0x3333333333333333ULL);
			for (u64 base = 0; base < 64; base += 4)
			{
				auto rows4 = _mm256_load_si256(reinterpret_cast<const __m256i*>(rows + base));
				const auto partner = _mm256_permute4x64_epi64(rows4, 0x4e);
				auto difference = _mm256_and_si256(
					_mm256_xor_si256(_mm256_srli_epi64(rows4, 2), partner), mask2);
				difference = _mm256_permute4x64_epi64(difference, 0x44);
				const auto adjustment = _mm256_blend_epi32(
					_mm256_slli_epi64(difference, 2), difference, 0xf0);
				rows4 = _mm256_xor_si256(rows4, adjustment);
				_mm256_store_si256(reinterpret_cast<__m256i*>(rows + base), rows4);
			}

			const auto mask1 = _mm256_set1_epi64x(0x5555555555555555ULL);
			for (u64 base = 0; base < 64; base += 4)
			{
				auto rows4 = _mm256_load_si256(reinterpret_cast<const __m256i*>(rows + base));
				const auto partner = _mm256_permute4x64_epi64(rows4, 0xb1);
				auto difference = _mm256_and_si256(
					_mm256_xor_si256(_mm256_srli_epi64(rows4, 1), partner), mask1);
				difference = _mm256_permute4x64_epi64(difference, 0xa0);
				const auto adjustment = _mm256_blend_epi32(
					_mm256_slli_epi64(difference, 1), difference, 0xcc);
				rows4 = _mm256_xor_si256(rows4, adjustment);
				_mm256_store_si256(reinterpret_cast<__m256i*>(rows + base), rows4);
			}
		}

		inline void riffleTransformField(block* state, u64 coefficient) noexcept
		{
			alignas(64) std::array<u64, 128> lanes;
			for (u64 row = 0; row < 64; ++row)
			{
				lanes[row] = static_cast<u64>(_mm_cvtsi128_si64(state[row].mData));
				lanes[64 + row] = static_cast<u64>(_mm_extract_epi64(state[row].mData, 1));
			}
			riffleTranspose64Avx2(lanes.data());
			riffleTranspose64Avx2(lanes.data() + 64);
			for (u64 bit = 0; bit < 64; ++bit)
			{
				const auto packed = _mm256_set_epi64x(
					0, static_cast<long long>(lanes[64 + bit]),
					0, static_cast<long long>(lanes[bit]));
				const auto product = riffleFieldMultiplyVector(packed, coefficient);
				lanes[bit] = static_cast<u64>(
					_mm_cvtsi128_si64(_mm256_castsi256_si128(product)));
				lanes[64 + bit] = static_cast<u64>(
					_mm_cvtsi128_si64(_mm256_extracti128_si256(product, 1)));
			}
			riffleTranspose64Avx2(lanes.data());
			riffleTranspose64Avx2(lanes.data() + 64);
			for (u64 row = 0; row < 64; ++row)
			{
				state[row].mData = _mm_set_epi64x(
					static_cast<long long>(lanes[64 + row]),
					static_cast<long long>(lanes[row]));
			}
		}
	}

	// Fixed first production candidate for the exact-permutation construction.
	// One object owns immutable setup data; Workspace owns all per-call scratch.
	class RiffleExactPermFieldCheckpoint
	{
	public:
		static constexpr u64 messageBlocks = u64{ 1 } << 20;
		static constexpr u64 codeBlocks = u64{ 1 } << 21;
		static constexpr u64 outerDimension = 128;
		static constexpr u64 outerLength = 256;
		static constexpr u64 outerBlocks = messageBlocks / outerDimension;
		static constexpr u64 regions = outerLength;
		static constexpr u64 regionSize = outerBlocks;
		static constexpr u64 epochBlocks = 256;
		static constexpr u64 stateBlocks = 64;
		static constexpr u64 epochs = codeBlocks / epochBlocks;
		static constexpr u64 tileOuterBlocks = 2048;
		static constexpr u64 tileBlocks = tileOuterBlocks * outerLength;
		static constexpr u64 buckets = outerBlocks / tileOuterBlocks;
		static constexpr u64 parityFanoutLayers = 56;

		static_assert(outerBlocks * outerLength == codeBlocks);
		static_assert(buckets == 4);

		struct Workspace
		{
			std::vector<block> bucketValues;
			std::vector<block> tile;

			Workspace()
				: bucketValues(codeBlocks)
				, tile(tileBlocks)
			{
			}
		};

		struct ExperimentalWorkspace : Workspace
		{
			std::vector<u32> bucketOffsets;

			ExperimentalWorkspace()
				: bucketOffsets(codeBlocks)
			{
			}
		};

		struct ParityFanout31x33Schedule
		{
			std::array<u8, 33> targets;
			std::array<u8, 31> sources;
		};

		void init(
			u64 permutationSeed,
			u64 coefficientSeed,
			bool retainExperimentalSchedules = false)
		{
			std::mt19937_64 random(permutationSeed);
			initPermutation(random);
			mCoefficients.resize(epochs - 1);
			for (auto& coefficient : mCoefficients)
			{
				coefficient = splitmix64(coefficientSeed);
				coefficient |= static_cast<u64>(coefficient == 0);
			}
			if (!retainExperimentalSchedules)
				releaseExperimentalSchedules();
		}

		// Production initialization. The caller controls the cryptographic PRNG
		// seed and can retain it when schedule reproducibility is required.
#ifdef LIBOTE_RIFFLE_HAS_PRNG
		void init(PRNG& random, bool retainExperimentalSchedules = false)
		{
			initPermutation(random);
			mCoefficients.resize(epochs - 1);
			for (auto& coefficient : mCoefficients)
			{
				coefficient = random.get<u64>();
				coefficient |= static_cast<u64>(coefficient == 0);
			}
			if (!retainExperimentalSchedules)
				releaseExperimentalSchedules();
		}
#endif

		void init(
			u64 permutationSeed,
			const u64* coefficients,
			u64 count,
			bool retainExperimentalSchedules = false)
		{
			if (count != epochs - 1)
				throw std::invalid_argument("FieldCheckpoint coefficient count mismatch");
			std::mt19937_64 random(permutationSeed);
			initPermutation(random);
			mCoefficients.assign(coefficients, coefficients + count);
			for (const auto coefficient : mCoefficients)
				if (coefficient == 0)
					throw std::invalid_argument("FieldCheckpoint coefficients must be nonzero");
			if (!retainExperimentalSchedules)
				releaseExperimentalSchedules();
		}

		void dualEncodeTo(
			const block* __restrict input,
			u64 inputSize,
			block* __restrict output,
			u64 outputSize,
			Workspace& workspace) const
		{
			if (mInnerToSlot24.size() < 3 * codeBlocks ||
				mBucketOffsets24.size() < 3 * codeBlocks ||
				mCoefficients.size() != epochs - 1)
				throw std::logic_error("RiffleExactPermFieldCheckpoint is not initialized");
			if (inputSize != codeBlocks || outputSize != messageBlocks)
				throw std::invalid_argument("RiffleExactPermFieldCheckpoint span size mismatch");
			dualEncodeUnchecked(input, output, workspace);
		}

		void dualEncodeUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			dualEncodePackedSchedulesUnchecked<0>(input, output, workspace);
		}

		// Benchmark-only copy of the promoted path. The callback is invoked after
		// the state/value pass (phase 0), then after each tile scatter (phase 1)
		// and its BCH transpose (phase 2). Keeping this separate avoids putting
		// timers or callback branches in the production kernel.
		template<typename Boundary>
		void dualEncodeProfiledUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace,
			Boundary& boundary) const noexcept
		{
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			const u8* slotBytes = mInnerToSlot24.data() + 3 * codeBlocks;

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					slotBytes -= 3;
					u32 slot;
					std::memcpy(&slot, slotBytes, sizeof(slot));
					state[lane] ^= input[inner];
					values[slot & static_cast<u32>(codeBlocks - 1)] = state[lane];
				}
			}
			boundary(0);

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				const u8* offsetBytes =
					mBucketOffsets24.data() + 3 * bucketBase;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					u32 local;
					std::memcpy(&local, offsetBytes, sizeof(local));
					offsetBytes += 3;
					tile[local & static_cast<u32>(tileBlocks - 1)] =
						values[bucketBase + offset];
				}
				boundary(1);

				const u64 outerBase = bucket * tileOuterBlocks;
#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
				{
					ExtendedBch256x128Eq3::transposeBlock4(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						tile + (offset + 2) * outerLength,
						tile + (offset + 3) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension,
						output + (outerBase + offset + 2) * outerDimension,
						output + (outerBase + offset + 3) * outerDimension);
				}
#else
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
#endif
				boundary(2);
			}
		}

		template<u64 PrefetchDistance>
		void dualEncodePrefetchUnchecked(
			const block* __restrict input,
			block* __restrict output,
			ExperimentalWorkspace& workspace) const noexcept
		{
			static_assert(PrefetchDistance != 0);
			dualEncodeBucketUnchecked<PrefetchDistance>(input, output, workspace);
		}

		template<u64 PrefetchDistance>
		void dualEncodeBucketUnchecked(
			const block* __restrict input,
			block* __restrict output,
			ExperimentalWorkspace& workspace) const noexcept
		{
			std::array<u32, buckets> counts;
			counts.fill(static_cast<u32>(tileBlocks));
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			u32* __restrict offsets = workspace.bucketOffsets.data();
			const u32* __restrict inverse = mInnerToOuter.data();

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					state[lane] ^= input[inner];
					const u32 outer = inverse[inner];
					const u32 bucket = outer / tileBlocks;
					const u32 local = outer & (tileBlocks - 1);
					const u64 slot = static_cast<u64>(bucket) * tileBlocks + --counts[bucket];
					values[slot] = state[lane];
					offsets[slot] = local;
				}
			}

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					if constexpr (PrefetchDistance != 0)
					{
						if (offset + PrefetchDistance < tileBlocks)
							_mm_prefetch(
								reinterpret_cast<const char*>(
									workspace.tile.data() +
									offsets[bucketBase + offset + PrefetchDistance]),
								_MM_HINT_T0);
					}
					tile[offsets[bucketBase + offset]] = values[bucketBase + offset];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
				{
					ExtendedBch256x128Eq3::transposeBlock4(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						tile + (offset + 2) * outerLength,
						tile + (offset + 3) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension,
						output + (outerBase + offset + 2) * outerDimension,
						output + (outerBase + offset + 3) * outerDimension);
				}
#else
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
#endif
			}
		}

		template<u64 PrefetchDistance>
		void dualEncodePrecomputedRouteUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			const u32* __restrict slots = mInnerToSlot.data();

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					state[lane] ^= input[inner];
					values[slots[inner]] = state[lane];
				}
			}

			block* __restrict tile = workspace.tile.data();
			const u32* __restrict offsets = mBucketOffsets.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					if constexpr (PrefetchDistance != 0)
					{
						if (offset + PrefetchDistance < tileBlocks)
							_mm_prefetch(
								reinterpret_cast<const char*>(
									tile + offsets[bucketBase + offset + PrefetchDistance]),
								_MM_HINT_T0);
					}
					tile[offsets[bucketBase + offset]] = values[bucketBase + offset];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
				{
					ExtendedBch256x128Eq3::transposeBlock4(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						tile + (offset + 2) * outerLength,
						tile + (offset + 3) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension,
						output + (outerBase + offset + 2) * outerDimension,
						output + (outerBase + offset + 3) * outerDimension);
				}
#else
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
#endif
			}
		}

		template<u64 PrefetchDistance>
		void dualEncodeBucketIdRouteUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			std::array<u32, buckets> counts;
			counts.fill(static_cast<u32>(tileBlocks));
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			const u8* __restrict bucketIds = mInnerBucket.data();

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					state[lane] ^= input[inner];
					const u32 bucket = bucketIds[inner];
					const u64 slot =
						static_cast<u64>(bucket) * tileBlocks + --counts[bucket];
					values[slot] = state[lane];
				}
			}

			finishPrecomputedRoute<PrefetchDistance>(output, workspace);
		}

		template<u64 PrefetchDistance>
		void dualEncodePackedOffsetRouteUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			const u32* __restrict slots = mInnerToSlot.data();

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					state[lane] ^= input[inner];
					values[slots[inner]] = state[lane];
				}
			}

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					if constexpr (PrefetchDistance != 0)
					{
						if (offset + PrefetchDistance < tileBlocks)
							_mm_prefetch(
								reinterpret_cast<const char*>(
									tile + packedOffset(
										bucketBase + offset + PrefetchDistance)),
								_MM_HINT_T0);
					}
					tile[packedOffset(bucketBase + offset)] = values[bucketBase + offset];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
			}
		}

		template<u64 PrefetchDistance>
		void dualEncodePackedSchedulesUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			const u8* slotBytes = mInnerToSlot24.data() + 3 * codeBlocks;

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					slotBytes -= 3;
					u32 slot;
					std::memcpy(&slot, slotBytes, sizeof(slot));
					state[lane] ^= input[inner];
					values[slot & static_cast<u32>(codeBlocks - 1)] = state[lane];
				}
			}

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				const u8* offsetBytes =
					mBucketOffsets24.data() + 3 * bucketBase;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					if constexpr (PrefetchDistance != 0)
					{
						if (offset + PrefetchDistance < tileBlocks)
						{
							u32 future;
							std::memcpy(
								&future,
								offsetBytes + 3 * PrefetchDistance,
								sizeof(future));
							_mm_prefetch(
								reinterpret_cast<const char*>(
									tile + (future & static_cast<u32>(tileBlocks - 1))),
								_MM_HINT_T0);
						}
					}
					u32 local;
					std::memcpy(&local, offsetBytes, sizeof(local));
					offsetBytes += 3;
					tile[local & static_cast<u32>(tileBlocks - 1)] =
						values[bucketBase + offset];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
				{
					ExtendedBch256x128Eq3::transposeBlock4(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						tile + (offset + 2) * outerLength,
						tile + (offset + 3) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension,
						output + (outerBase + offset + 2) * outerDimension,
						output + (outerBase + offset + 3) * outerDimension);
				}
#else
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
#endif
			}
		}

		// Exact same 24-bit permutation route and BCH transpose as the promoted
		// encoder, with a compile-time inner emitter.  The emitter finalizes
		// values in reverse inner order, so no materialized inner word or
		// indirect callback is required in the hot path.
		template<u64 PrefetchDistance, typename Inner>
		void dualEncodePackedInnerUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace,
			const Inner& inner) const noexcept
		{
			dualEncodePackedInnerFanoutUnchecked<PrefetchDistance, 1>(
				input, output, workspace, inner);
		}

		template<u64 PrefetchDistance, u64 AppliedFanoutLayers, typename Inner>
		void dualEncodePackedInnerFanoutUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace,
			const Inner& inner) const noexcept
		{
			static_assert(AppliedFanoutLayers <= parityFanoutLayers);
			block* __restrict values = workspace.bucketValues.data();
			const u8* slotBytes = mInnerToSlot24.data() + 3 * codeBlocks;
			inner.emitReverse(input, codeBlocks, [&](u64 innerIndex, block value) {
				(void)innerIndex;
				slotBytes -= 3;
				u32 slot;
				std::memcpy(&slot, slotBytes, sizeof(slot));
				values[slot & static_cast<u32>(codeBlocks - 1)] = value;
			});

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				const u8* offsetBytes = mBucketOffsets24.data() + 3 * bucketBase;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					if constexpr (PrefetchDistance != 0)
					{
						if (offset + PrefetchDistance < tileBlocks)
						{
							u32 future;
							std::memcpy(&future, offsetBytes + 3 * PrefetchDistance, sizeof(future));
							_mm_prefetch(
								reinterpret_cast<const char*>(
									tile + (future & static_cast<u32>(tileBlocks - 1))),
								_MM_HINT_T0);
						}
					}
					u32 local;
					std::memcpy(&local, offsetBytes, sizeof(local));
					offsetBytes += 3;
					tile[local & static_cast<u32>(tileBlocks - 1)] =
						values[bucketBase + offset];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
#if defined(LIBOTE_RIFFLE_EXPERIMENTAL_AVX512_OUTER) && defined(__AVX512F__)
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 4)
				{
					for (u64 lane = 0; lane < 4; ++lane)
					{
						block* word = tile + (offset + lane) * outerLength;
						const u64 fanoutBase =
							(outerBase + offset + lane) * parityFanoutLayers;
						for (u64 layer = AppliedFanoutLayers; layer-- > 0;)
						{
							const auto& fanout = mParityFanouts[fanoutBase + layer];
							block pivot = word[fanout.targets[0]];
							for (unsigned index = 1; index < fanout.targets.size(); ++index)
								pivot ^= word[fanout.targets[index]];
							for (const u8 source : fanout.sources) word[source] ^= pivot;
						}
					}
					ExtendedBch256x128Eq3::transposeBlock4(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						tile + (offset + 2) * outerLength,
						tile + (offset + 3) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension,
						output + (outerBase + offset + 2) * outerDimension,
						output + (outerBase + offset + 3) * outerDimension);
				}
#else
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					block* word0 = tile + offset * outerLength;
					block* word1 = word0 + outerLength;
					const u64 fanoutBase0 =
						(outerBase + offset) * parityFanoutLayers;
					const u64 fanoutBase1 = fanoutBase0 + parityFanoutLayers;
					for (u64 layer = AppliedFanoutLayers; layer-- > 0;)
					{
						const auto& fanout0 = mParityFanouts[fanoutBase0 + layer];
						const auto& fanout1 = mParityFanouts[fanoutBase1 + layer];
						block pivot0 = word0[fanout0.targets[0]];
						block pivot1 = word1[fanout1.targets[0]];
						for (unsigned index = 1; index < fanout0.targets.size(); ++index)
						{
							pivot0 ^= word0[fanout0.targets[index]];
							pivot1 ^= word1[fanout1.targets[index]];
						}
						for (unsigned index = 0; index < fanout0.sources.size(); ++index)
						{
							word0[fanout0.sources[index]] ^= pivot0;
							word1[fanout1.sources[index]] ^= pivot1;
						}
					}
					ExtendedBch256x128Eq3::transposeBlock2(
						word0,
						word1,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
#endif
			}
		}

		void dualEncodePairPackedUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			constexpr u64 slotMask = codeBlocks - 1;
			constexpr u64 offsetMask = tileBlocks - 1;
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				const u64 epochBase = epoch * epochBlocks;
				for (u64 offset = epochBlocks; offset != 0; offset -= 2)
				{
					const u64 inner0 = epochBase + offset - 2;
					const u64 inner1 = inner0 + 1;
					u64 packed;
					std::memcpy(
						&packed,
						mInnerToSlotPair6.data() + 6 * (inner0 >> 1),
						sizeof(packed));
					const u32 slot0 = static_cast<u32>(packed & slotMask);
					const u32 slot1 = static_cast<u32>((packed >> 21) & slotMask);
					state[(offset - 1) & (stateBlocks - 1)] ^= input[inner1];
					values[slot1] = state[(offset - 1) & (stateBlocks - 1)];
					state[(offset - 2) & (stateBlocks - 1)] ^= input[inner0];
					values[slot0] = state[(offset - 2) & (stateBlocks - 1)];
				}
			}

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				for (u64 offset = 0; offset < tileBlocks; offset += 2)
				{
					u64 packed;
					std::memcpy(
						&packed,
						mBucketOffsetsPair5.data() + 5 * ((bucketBase + offset) >> 1),
						sizeof(packed));
					const u32 local0 = static_cast<u32>(packed & offsetMask);
					const u32 local1 = static_cast<u32>((packed >> 19) & offsetMask);
					tile[local0] = values[bucketBase + offset];
					tile[local1] = values[bucketBase + offset + 1];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
			}
		}

		// Minimal-traffic alternative: finalized checkpoint values scatter
		// directly into exact outer order, followed by one sequential outer pass.
		// bucketValues is reused as the N-element outer word.
		void dualEncodeDirectScatterUnchecked(
			const block* __restrict input,
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			std::array<block, stateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict outerWord = workspace.bucketValues.data();
			const u32* __restrict inverse = mInnerToOuter.data();

			for (u64 epoch = epochs; epoch-- > 0;)
			{
				if (epoch + 1 < epochs)
					detail::riffleTransformField(state.data(), mCoefficients[epoch]);
				for (u64 offset = epochBlocks; offset-- > 0;)
				{
					const u64 inner = epoch * epochBlocks + offset;
					const u64 lane = offset & (stateBlocks - 1);
					state[lane] ^= input[inner];
					outerWord[inverse[inner]] = state[lane];
				}
			}

			for (u64 outer = 0; outer < outerBlocks; outer += 2)
			{
				ExtendedBch256x128Eq3::transposeBlock2(
					outerWord + outer * outerLength,
					outerWord + (outer + 1) * outerLength,
					output + outer * outerDimension,
					output + (outer + 1) * outerDimension);
			}
		}

		const std::vector<u32>& innerToOuter() const noexcept { return mInnerToOuter; }
		const std::vector<u64>& coefficients() const noexcept { return mCoefficients; }
		const std::vector<ParityFanout31x33Schedule>& parityFanouts() const noexcept
		{
			return mParityFanouts;
		}

	private:
		std::vector<u32> mInnerToOuter;
		std::vector<u32> mInnerToSlot;
		std::vector<u8> mInnerToSlot24;
		std::vector<u8> mInnerToSlotPair6;
		std::vector<u8> mInnerBucket;
		std::vector<u32> mBucketOffsets;
		std::vector<u8> mBucketOffsets24;
		std::vector<u8> mBucketOffsetsPair5;
		std::vector<u64> mCoefficients;
		std::vector<ParityFanout31x33Schedule> mParityFanouts;

		template<typename Random>
		void initPermutation(Random& random)
		{
			std::vector<u32> outerToInner(codeBlocks);
			std::vector<u8> regionToCoordinate(codeBlocks);
			std::array<u16, outerLength> coordinateToRegion;
			for (u64 outerBlock = 0; outerBlock < outerBlocks; ++outerBlock)
			{
				std::iota(coordinateToRegion.begin(), coordinateToRegion.end(), u16{ 0 });
				std::shuffle(coordinateToRegion.begin(), coordinateToRegion.end(), random);
				for (u64 coordinate = 0; coordinate < outerLength; ++coordinate)
				{
					const u64 region = coordinateToRegion[coordinate];
					regionToCoordinate[region * regionSize + outerBlock] =
						static_cast<u8>(coordinate);
				}
			}

			std::vector<u32> blockToPosition(regionSize);
			for (u64 region = 0; region < regions; ++region)
			{
				std::iota(blockToPosition.begin(), blockToPosition.end(), u32{ 0 });
				std::shuffle(blockToPosition.begin(), blockToPosition.end(), random);
				for (u64 outerBlock = 0; outerBlock < outerBlocks; ++outerBlock)
				{
					const u64 coordinate =
						regionToCoordinate[region * regionSize + outerBlock];
					const u32 position = blockToPosition[outerBlock];
					outerToInner[outerBlock * outerLength + coordinate] =
						static_cast<u32>(region * regionSize + position);
				}
			}

			mInnerToOuter.resize(codeBlocks);
			for (u64 outer = 0; outer < codeBlocks; ++outer)
				mInnerToOuter[outerToInner[outer]] = static_cast<u32>(outer);

			mInnerToSlot.resize(codeBlocks);
			mInnerToSlot24.resize(3 * codeBlocks + sizeof(u32));
			mInnerToSlotPair6.resize(3 * codeBlocks + sizeof(u64));
			mInnerBucket.resize(codeBlocks);
			mBucketOffsets.resize(codeBlocks);
			mBucketOffsets24.resize(3 * codeBlocks + sizeof(u32));
			mBucketOffsetsPair5.resize(5 * (codeBlocks / 2) + sizeof(u64));
			std::array<u32, buckets> counts;
			counts.fill(static_cast<u32>(tileBlocks));
			for (u64 inner = codeBlocks; inner-- > 0;)
			{
				const u32 outer = mInnerToOuter[inner];
				const u32 bucket = outer / tileBlocks;
				const u32 local = outer & (tileBlocks - 1);
				const u32 slot = static_cast<u32>(
					static_cast<u64>(bucket) * tileBlocks + --counts[bucket]);
				mInnerToSlot[inner] = slot;
				auto* packedSlotBytes =
					mInnerToSlot24.data() + 3 * static_cast<u64>(inner);
				packedSlotBytes[0] = static_cast<u8>(slot);
				packedSlotBytes[1] = static_cast<u8>(slot >> 8);
				packedSlotBytes[2] = static_cast<u8>(slot >> 16);
				mInnerBucket[inner] = static_cast<u8>(bucket);
				mBucketOffsets[slot] = local;
				auto* packed = mBucketOffsets24.data() + 3 * static_cast<u64>(slot);
				packed[0] = static_cast<u8>(local);
				packed[1] = static_cast<u8>(local >> 8);
				packed[2] = static_cast<u8>(local >> 16);
			}
			for (u64 pair = 0; pair < codeBlocks / 2; ++pair)
			{
				const u64 slots = static_cast<u64>(mInnerToSlot[2 * pair]) |
					(static_cast<u64>(mInnerToSlot[2 * pair + 1]) << 21);
				auto* slotBytes = mInnerToSlotPair6.data() + 6 * pair;
				for (u64 byte = 0; byte < 6; ++byte)
					slotBytes[byte] = static_cast<u8>(slots >> (8 * byte));

				const u64 offsets = static_cast<u64>(mBucketOffsets[2 * pair]) |
					(static_cast<u64>(mBucketOffsets[2 * pair + 1]) << 19);
				auto* offsetBytes = mBucketOffsetsPair5.data() + 5 * pair;
				for (u64 byte = 0; byte < 5; ++byte)
					offsetBytes[byte] = static_cast<u8>(offsets >> (8 * byte));
			}

			mParityFanouts.resize(outerBlocks * parityFanoutLayers);
			std::uniform_int_distribution<unsigned> coordinateDistribution(0, outerLength - 1);
			for (auto& fanout : mParityFanouts)
			{
				std::array<u8, 64> sampled;
				for (unsigned pick = 0; pick < sampled.size(); ++pick)
				{
					u8 candidate;
					do candidate = static_cast<u8>(coordinateDistribution(random));
					while (std::find(sampled.begin(), sampled.begin() + pick, candidate) != sampled.begin() + pick);
					sampled[pick] = candidate;
				}
				std::copy(sampled.begin(), sampled.begin() + 33, fanout.targets.begin());
				std::copy(sampled.begin() + 33, sampled.end(), fanout.sources.begin());
			}
		}

		OC_FORCEINLINE u32 packedOffset(u64 index) const noexcept
		{
			u32 value;
			std::memcpy(&value, mBucketOffsets24.data() + 3 * index, sizeof(value));
			return value & static_cast<u32>(tileBlocks - 1);
		}

		OC_FORCEINLINE u32 packedSlot(u64 index) const noexcept
		{
			u32 value;
			std::memcpy(&value, mInnerToSlot24.data() + 3 * index, sizeof(value));
			return value & static_cast<u32>(codeBlocks - 1);
		}

		template<u64 PrefetchDistance>
		void finishPrecomputedRoute(
			block* __restrict output,
			Workspace& workspace) const noexcept
		{
			block* __restrict values = workspace.bucketValues.data();
			block* __restrict tile = workspace.tile.data();
			const u32* __restrict offsets = mBucketOffsets.data();
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
				{
					if constexpr (PrefetchDistance != 0)
					{
						if (offset + PrefetchDistance < tileBlocks)
							_mm_prefetch(
								reinterpret_cast<const char*>(
									tile + offsets[bucketBase + offset + PrefetchDistance]),
								_MM_HINT_T0);
					}
					tile[offsets[bucketBase + offset]] = values[bucketBase + offset];
				}

				const u64 outerBase = bucket * tileOuterBlocks;
				for (u64 offset = 0; offset < tileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * outerLength,
						tile + (offset + 1) * outerLength,
						output + (outerBase + offset) * outerDimension,
						output + (outerBase + offset + 1) * outerDimension);
				}
			}
		}

		static u64 splitmix64(u64& state) noexcept
		{
			u64 value = (state += 0x9e3779b97f4a7c15ULL);
			value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
			value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
			return value ^ (value >> 31);
		}

		void releaseExperimentalSchedules()
		{
			mInnerToOuter.clear();
			mInnerToOuter.shrink_to_fit();
			mInnerToSlot.clear();
			mInnerToSlot.shrink_to_fit();
			mInnerToSlotPair6.clear();
			mInnerToSlotPair6.shrink_to_fit();
			mInnerBucket.clear();
			mInnerBucket.shrink_to_fit();
			mBucketOffsets.clear();
			mBucketOffsets.shrink_to_fit();
			mBucketOffsetsPair5.clear();
			mBucketOffsetsPair5.shrink_to_fit();
		}
	};
}

#ifdef LIBOTE_RIFFLE_HAS_PRNG
#undef LIBOTE_RIFFLE_HAS_PRNG
#endif
