#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <immintrin.h>
#ifdef _WIN32
#include <intrin.h>
#endif

#include "libOTe/Tools/RiffleCode/ExtendedBch256x128Eq3.h"
#if defined(LIBOTE_RIFFLE_BCH_CIRCUIT_PAIR_BENCH)
#include "libOTe/Tools/RiffleCode/ExtendedBch256x128Eq3TransposeCircuit2Seed57.h"
#endif
#include "libOTe/Tools/RiffleCode/RiffleExactPermFieldCheckpoint.h"

#ifdef _WIN32
#include <Windows.h>
#else
#include <pthread.h>
#include <sched.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

using namespace osuCrypto;

namespace
{
	constexpr u64 K = u64{ 1 } << 20;
	constexpr u64 N = u64{ 1 } << 21;
	constexpr u64 OuterDimension = 128;
	constexpr u64 OuterLength = 256;
	constexpr u64 OuterBlocks = K / OuterDimension;
	constexpr u64 Regions = OuterLength;
	constexpr u64 RegionSize = OuterBlocks;
	constexpr u64 EpochBlocks = 256;
	constexpr u64 StateBlocks = 64;
	constexpr u64 Epochs = N / EpochBlocks;

	static_assert(OuterBlocks * OuterLength == N);
	static_assert(Regions * RegionSize == N);
	static_assert((StateBlocks & (StateBlocks - 1)) == 0);

	volatile u64 benchmarkSink = 0;

	struct Product128
	{
		u64 low;
		u64 high;
	};

	OC_FORCEINLINE Product128 carrylessMultiply(u64 left, u64 right) noexcept
	{
		const auto product = _mm_clmulepi64_si128(
			_mm_cvtsi64_si128(static_cast<long long>(left)),
			_mm_cvtsi64_si128(static_cast<long long>(right)), 0x00);
		return {
			static_cast<u64>(_mm_cvtsi128_si64(product)),
			static_cast<u64>(_mm_extract_epi64(product, 1))
		};
	}

	OC_FORCEINLINE u64 fieldMultiplyScalar(u64 value, u64 coefficient) noexcept
	{
		// x^64 + x^11 + x^2 + x + 1.
		constexpr u64 modulusLow = 0x807;
		auto product = carrylessMultiply(value, coefficient);
		auto fold = carrylessMultiply(product.high, modulusLow);
		product.low ^= fold.low;
		fold = carrylessMultiply(fold.high, modulusLow);
		return product.low ^ fold.low;
	}

	OC_FORCEINLINE __m256i fieldMultiplyVector(__m256i value, u64 coefficient) noexcept
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
	OC_FORCEINLINE void transpose64LargeStage(u64* rows) noexcept
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

	OC_FORCEINLINE void transpose64Avx2(u64* rows) noexcept
	{
		transpose64LargeStage<32, 0x00000000ffffffffULL>(rows);
		transpose64LargeStage<16, 0x0000ffff0000ffffULL>(rows);
		transpose64LargeStage<8, 0x00ff00ff00ff00ffULL>(rows);
		transpose64LargeStage<4, 0x0f0f0f0f0f0f0f0fULL>(rows);

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

	void transpose64Scalar(u64* rows) noexcept
	{
		u64 mask = 0x00000000ffffffffULL;
		for (u64 shift = 32; shift != 0; shift >>= 1)
		{
			for (u64 base = 0; base < 64; base = (base + shift + 1) & ~shift)
			{
				const u64 difference = ((rows[base] >> shift) ^ rows[base + shift]) & mask;
				rows[base] ^= difference << shift;
				rows[base + shift] ^= difference;
			}
			mask ^= mask << (shift >> 1);
		}
	}

	template<void (*Transpose)(u64*)>
	void transformField(block* state, u64 coefficient) noexcept
	{
		alignas(64) std::array<u64, 128> lanes;
		for (u64 row = 0; row < StateBlocks; ++row)
		{
			lanes[row] = static_cast<u64>(_mm_cvtsi128_si64(state[row].mData));
			lanes[64 + row] = static_cast<u64>(_mm_extract_epi64(state[row].mData, 1));
		}
		Transpose(lanes.data());
		Transpose(lanes.data() + 64);
		if constexpr (Transpose == transpose64Avx2)
		{
			for (u64 bit = 0; bit < 64; ++bit)
			{
				const auto packed = _mm256_set_epi64x(
					0, static_cast<long long>(lanes[64 + bit]),
					0, static_cast<long long>(lanes[bit]));
				const auto product = fieldMultiplyVector(packed, coefficient);
				lanes[bit] = static_cast<u64>(
					_mm_cvtsi128_si64(_mm256_castsi256_si128(product)));
				lanes[64 + bit] = static_cast<u64>(
					_mm_cvtsi128_si64(_mm256_extracti128_si256(product, 1)));
			}
		}
		else
		{
			for (auto& lane : lanes)
				lane = fieldMultiplyScalar(lane, coefficient);
		}
		Transpose(lanes.data());
		Transpose(lanes.data() + 64);
		for (u64 row = 0; row < StateBlocks; ++row)
		{
			state[row].mData = _mm_set_epi64x(
				static_cast<long long>(lanes[64 + row]),
				static_cast<long long>(lanes[row]));
		}
	}

	template<void (*Transform)(block*, u64)>
	void checkpointTranspose(block* word, u64 wordBlocks, const u64* coefficients) noexcept
	{
		std::array<block, StateBlocks> state;
		state.fill(block(0, 0));
		const u64 epochs = wordBlocks / EpochBlocks;
		for (u64 epoch = epochs; epoch-- > 0;)
		{
			if (epoch + 1 < epochs)
				Transform(state.data(), coefficients[epoch]);
			block* position = word + epoch * EpochBlocks;
			for (u64 offset = EpochBlocks; offset-- > 0;)
			{
				const u64 lane = offset & (StateBlocks - 1);
				state[lane] ^= position[offset];
				position[offset] = state[lane];
			}
		}
	}

	void checkpointOptimized(block* word, u64 wordBlocks, const u64* coefficients) noexcept
	{
		checkpointTranspose<transformField<transpose64Avx2>>(word, wordBlocks, coefficients);
	}

	void checkpointReference(block* word, u64 wordBlocks, const u64* coefficients) noexcept
	{
		checkpointTranspose<transformField<transpose64Scalar>>(word, wordBlocks, coefficients);
	}

	class StructuredSchedule
	{
	public:
		static constexpr u64 optimizedTileOuterBlocks = 512;
		static constexpr u64 optimizedPrefetchDistance = 48;

		void init(u64 seed)
		{
			std::mt19937_64 random(seed);
			mOuterToInner.resize(N);
			mInnerToOuter.resize(N);
			mRegionPlan.resize(N);
			std::vector<u8> regionToCoordinate(N);
			std::array<u16, OuterLength> coordinateToRegion;
			for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
			{
				std::iota(coordinateToRegion.begin(), coordinateToRegion.end(), u16{ 0 });
				std::shuffle(coordinateToRegion.begin(), coordinateToRegion.end(), random);
				for (u64 coordinate = 0; coordinate < OuterLength; ++coordinate)
				{
					const u64 region = coordinateToRegion[coordinate];
					regionToCoordinate[region * RegionSize + outerBlock] =
						static_cast<u8>(coordinate);
				}
			}

			std::vector<u32> blockToPosition(RegionSize);
			for (u64 region = 0; region < Regions; ++region)
			{
				std::iota(blockToPosition.begin(), blockToPosition.end(), u32{ 0 });
				std::shuffle(blockToPosition.begin(), blockToPosition.end(), random);
				for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
				{
					const u64 coordinate =
						regionToCoordinate[region * RegionSize + outerBlock];
					const u32 position = blockToPosition[outerBlock];
					mOuterToInner[outerBlock * OuterLength + coordinate] =
						static_cast<u32>(region * RegionSize + position);
					mRegionPlan[region * RegionSize + outerBlock] =
						static_cast<u32>(coordinate) | (position << 8);
				}
			}
			for (u64 outer = 0; outer < N; ++outer)
				mInnerToOuter[mOuterToInner[outer]] = static_cast<u32>(outer);
		}

		void validate() const
		{
			std::vector<u8> seen(N, 0);
			for (const auto index : mOuterToInner)
			{
				if (index >= N || seen[index])
					throw std::runtime_error("structured schedule is not a permutation");
				seen[index] = 1;
			}
			for (u64 inner = 0; inner < N; ++inner)
			{
				const u64 outer = mInnerToOuter[inner];
				if (outer >= N || mOuterToInner[outer] != inner)
					throw std::runtime_error("inverse structured schedule is inconsistent");
			}
		}

		void releaseFlatSchedule()
		{
			mOuterToInner.clear();
			mOuterToInner.shrink_to_fit();
		}

		const std::vector<u32>& regionPlan() const noexcept
		{
			return mRegionPlan;
		}

		void gather(const block* __restrict inner, block* __restrict outer) const noexcept
		{
			for (u64 index = 0; index < N; ++index)
				outer[index] = inner[mOuterToInner[index]];
		}

		void gatherUnrolled(const block* __restrict inner, block* __restrict outer) const noexcept
		{
			const u32* __restrict indices = mOuterToInner.data();
			for (u64 index = 0; index < N; index += 8)
			{
				const u32 i0 = indices[index + 0];
				const u32 i1 = indices[index + 1];
				const u32 i2 = indices[index + 2];
				const u32 i3 = indices[index + 3];
				const u32 i4 = indices[index + 4];
				const u32 i5 = indices[index + 5];
				const u32 i6 = indices[index + 6];
				const u32 i7 = indices[index + 7];
				outer[index + 0] = inner[i0];
				outer[index + 1] = inner[i1];
				outer[index + 2] = inner[i2];
				outer[index + 3] = inner[i3];
				outer[index + 4] = inner[i4];
				outer[index + 5] = inner[i5];
				outer[index + 6] = inner[i6];
				outer[index + 7] = inner[i7];
			}
		}

		template<u64 Distance>
		void gatherPrefetch(const block* __restrict inner, block* __restrict outer) const noexcept
		{
			static_assert(Distance % 8 == 0);
			const u32* __restrict indices = mOuterToInner.data();
			for (u64 index = 0; index < N; index += 8)
			{
				if (index + Distance + 7 < N)
				{
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 0]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 1]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 2]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 3]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 4]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 5]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 6]), _MM_HINT_T0);
					_mm_prefetch(reinterpret_cast<const char*>(inner + indices[index + Distance + 7]), _MM_HINT_T0);
				}
				const u32 i0 = indices[index + 0];
				const u32 i1 = indices[index + 1];
				const u32 i2 = indices[index + 2];
				const u32 i3 = indices[index + 3];
				const u32 i4 = indices[index + 4];
				const u32 i5 = indices[index + 5];
				const u32 i6 = indices[index + 6];
				const u32 i7 = indices[index + 7];
				outer[index + 0] = inner[i0];
				outer[index + 1] = inner[i1];
				outer[index + 2] = inner[i2];
				outer[index + 3] = inner[i3];
				outer[index + 4] = inner[i4];
				outer[index + 5] = inner[i5];
				outer[index + 6] = inner[i6];
				outer[index + 7] = inner[i7];
			}
		}

		void gatherAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message) const noexcept
		{
			alignas(64) std::array<block, OuterLength> local;
			for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
			{
				const u64 base = outerBlock * OuterLength;
				for (u64 coordinate = 0; coordinate < OuterLength; ++coordinate)
					local[coordinate] = inner[mOuterToInner[base + coordinate]];
				ExtendedBch256x128Eq3::transposeBlock(
					local.data(), message + outerBlock * OuterDimension);
			}
		}

		void gatherAndOuterTranspose2(
			const block* __restrict inner,
			block* __restrict message) const noexcept
		{
			alignas(64) std::array<block, 2 * OuterLength> local;
			for (u64 outerBlock = 0; outerBlock < OuterBlocks; outerBlock += 2)
			{
				for (u64 pair = 0; pair < 2; ++pair)
				{
					const u64 base = (outerBlock + pair) * OuterLength;
					for (u64 coordinate = 0; coordinate < OuterLength; ++coordinate)
						local[pair * OuterLength + coordinate] =
							inner[mOuterToInner[base + coordinate]];
				}
				ExtendedBch256x128Eq3::transposeBlock2(
					local.data(), local.data() + OuterLength,
					message + outerBlock * OuterDimension,
					message + (outerBlock + 1) * OuterDimension);
			}
		}

		void tiledGatherAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace,
			u64 batchSize) const noexcept
		{
			tiledGatherAndOuterTransposeImpl<0, false>(
				inner, message, workspace, batchSize);
		}

		void tiledPrefetchGatherAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace,
			u64 batchSize) const noexcept
		{
			tiledGatherAndOuterTransposeImpl<12, false>(
				inner, message, workspace, batchSize);
		}

		void tiledPackedGatherAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace,
			u64 batchSize) const noexcept
		{
			tiledGatherAndOuterTransposeImpl<0, true>(
				inner, message, workspace, batchSize);
		}

		void tiledPrefetchPackedGatherAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace,
			u64 batchSize) const noexcept
		{
			tiledGatherAndOuterTransposeImpl<12, true>(
				inner, message, workspace, batchSize);
		}

		template<u64 PrefetchDistance>
		void tiledPackedPrefetchDistance(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace,
			u64 batchSize) const noexcept
		{
			tiledGatherAndOuterTransposeImpl<PrefetchDistance, true>(
				inner, message, workspace, batchSize);
		}

		void optimizedGatherAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace) const noexcept
		{
			tiledGatherAndOuterTransposeImpl<optimizedPrefetchDistance, true>(
				inner, message, workspace, optimizedTileOuterBlocks);
		}

		template<u64 TileOuterBlocks, typename Offset>
		void bucketRouteAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict bucketValues,
			Offset* __restrict bucketOffsets,
			block* __restrict tile) const noexcept
		{
			static_assert(TileOuterBlocks != 0 &&
				(TileOuterBlocks & (TileOuterBlocks - 1)) == 0);
			static_assert(OuterBlocks % TileOuterBlocks == 0);
			constexpr u64 tileBlocks = TileOuterBlocks * OuterLength;
			constexpr u64 buckets = OuterBlocks / TileOuterBlocks;
			static_assert(tileBlocks <= u64{ 1 } +
				static_cast<u64>(std::numeric_limits<Offset>::max()));
			std::array<u32, buckets> counts{};
			const u32* __restrict inverse = mInnerToOuter.data();

			// Read the inner word once in order. Each bucket is an independent
			// sequential stream; only the small destination offset travels with it.
			for (u64 innerIndex = 0; innerIndex < N; ++innerIndex)
			{
				const u32 outerIndex = inverse[innerIndex];
				const u32 bucket = outerIndex / tileBlocks;
				const u32 local = outerIndex & (tileBlocks - 1);
				const u64 slot = static_cast<u64>(bucket) * tileBlocks + counts[bucket]++;
				bucketValues[slot] = inner[innerIndex];
				bucketOffsets[slot] = static_cast<Offset>(local);
			}
			finishBucketRoute<TileOuterBlocks>(
				message, bucketValues, bucketOffsets, tile);
		}

		template<u64 TileOuterBlocks, typename Offset>
		void checkpointBucketRouteAndOuterTranspose(
			const block* __restrict input,
			block* __restrict message,
			const u64* __restrict coefficients,
			block* __restrict bucketValues,
			Offset* __restrict bucketOffsets,
			block* __restrict tile) const noexcept
		{
			static_assert(TileOuterBlocks != 0 &&
				(TileOuterBlocks & (TileOuterBlocks - 1)) == 0);
			static_assert(OuterBlocks % TileOuterBlocks == 0);
			constexpr u64 tileBlocks = TileOuterBlocks * OuterLength;
			constexpr u64 buckets = OuterBlocks / TileOuterBlocks;
			static_assert(tileBlocks <= u64{ 1 } +
				static_cast<u64>(std::numeric_limits<Offset>::max()));
			std::array<u32, buckets> counts;
			counts.fill(static_cast<u32>(tileBlocks));
			std::array<block, StateBlocks> state;
			state.fill(block(0, 0));
			const u32* __restrict inverse = mInnerToOuter.data();

			// The checkpoint accumulator finalizes outputs in reverse source order.
			// Route each finalized value immediately, avoiding a 32 MiB write and
			// subsequent read of the materialized inner word.
			for (u64 epoch = Epochs; epoch-- > 0;)
			{
				if (epoch + 1 < Epochs)
					transformField<transpose64Avx2>(state.data(), coefficients[epoch]);
				for (u64 offset = EpochBlocks; offset-- > 0;)
				{
					const u64 innerIndex = epoch * EpochBlocks + offset;
					const u64 lane = offset & (StateBlocks - 1);
					state[lane] ^= input[innerIndex];
					const u32 outerIndex = inverse[innerIndex];
					const u32 bucket = outerIndex / tileBlocks;
					const u32 local = outerIndex & (tileBlocks - 1);
					const u64 slot = static_cast<u64>(bucket) * tileBlocks + --counts[bucket];
					bucketValues[slot] = state[lane];
					bucketOffsets[slot] = static_cast<Offset>(local);
				}
			}
			finishBucketRoute<TileOuterBlocks>(
				message, bucketValues, bucketOffsets, tile);
		}

	private:
		template<u64 TileOuterBlocks, typename Offset>
		static void finishBucketRoute(
			block* __restrict message,
			const block* __restrict bucketValues,
			const Offset* __restrict bucketOffsets,
			block* __restrict tile) noexcept
		{
			constexpr u64 tileBlocks = TileOuterBlocks * OuterLength;
			constexpr u64 buckets = OuterBlocks / TileOuterBlocks;
			for (u64 bucket = 0; bucket < buckets; ++bucket)
			{
				const u64 bucketBase = bucket * tileBlocks;
				for (u64 offset = 0; offset < tileBlocks; ++offset)
					tile[bucketOffsets[bucketBase + offset]] =
						bucketValues[bucketBase + offset];

				const u64 outerBase = bucket * TileOuterBlocks;
				for (u64 offset = 0; offset < TileOuterBlocks; offset += 2)
				{
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * OuterLength,
						tile + (offset + 1) * OuterLength,
						message + (outerBase + offset) * OuterDimension,
						message + (outerBase + offset + 1) * OuterDimension);
				}
			}
		}

		template<u64 PrefetchDistance, bool Packed>
		void tiledGatherAndOuterTransposeImpl(
			const block* __restrict inner,
			block* __restrict message,
			block* __restrict workspace,
			u64 batchSize) const noexcept
		{
			for (u64 outerBase = 0; outerBase < OuterBlocks; outerBase += batchSize)
			{
				for (u64 region = 0; region < Regions; ++region)
				{
					const u64 regionBase = region * RegionSize;
					const u32* __restrict plan =
						mRegionPlan.data() + regionBase + outerBase;
					for (u64 offset = 0; offset < batchSize; ++offset)
					{
						if constexpr (PrefetchDistance != 0)
						{
							if (offset + PrefetchDistance < batchSize)
							{
								const u32 future = plan[offset + PrefetchDistance];
								_mm_prefetch(
									reinterpret_cast<const char*>(
										inner + regionBase + (future >> 8)),
									_MM_HINT_T0);
							}
						}
						const u32 operation = plan[offset];
						workspace[offset * OuterLength + (operation & 255)] =
							inner[regionBase + (operation >> 8)];
					}
				}
				if constexpr (Packed)
				{
					for (u64 offset = 0; offset < batchSize; offset += 2)
					{
						ExtendedBch256x128Eq3::transposeBlock2(
							workspace + offset * OuterLength,
							workspace + (offset + 1) * OuterLength,
							message + (outerBase + offset) * OuterDimension,
							message + (outerBase + offset + 1) * OuterDimension);
					}
				}
				else
				{
					for (u64 offset = 0; offset < batchSize; ++offset)
					{
						ExtendedBch256x128Eq3::transposeBlock(
							workspace + offset * OuterLength,
							message + (outerBase + offset) * OuterDimension);
					}
				}
			}
		}

		std::vector<u32> mOuterToInner;
		// Inverse schedule for source-major exact routing. Keeping both maps is
		// an evaluator choice; it does not change the sampled permutation.
		std::vector<u32> mInnerToOuter;
		// Region-major packed records: low 8 bits are the destination outer
		// coordinate and the remaining 13 bits are the source position in the
		// region. This preserves the two-dimensional permutation structure.
		std::vector<u32> mRegionPlan;
	};

	// Narrow performance evaluator for the distinct g=4 packet-permutation
	// construction. Each region permutes 2048 fixed four-block groups. The
	// transposed hot path routes one 64-byte packet at a time while retaining
	// the current FieldCheckpoint and BCH constituents.
	class PacketShuffle4Schedule
	{
	public:
		static constexpr u64 PacketWidth = 4;
		static constexpr u64 PacketsPerRegion = RegionSize / PacketWidth;
		static constexpr u64 PacketCount = N / PacketWidth;
		static constexpr u64 TileOuterBlocks = 2048;
		static constexpr u64 TileBlocks = TileOuterBlocks * OuterLength;
		static constexpr u64 Buckets = OuterBlocks / TileOuterBlocks;
		static constexpr u64 TileOuterGroups = TileOuterBlocks / PacketWidth;
		static constexpr u64 TilePackets = TileBlocks / PacketWidth;

		static_assert(RegionSize % PacketWidth == 0);
		static_assert(TileOuterBlocks % PacketWidth == 0);
		static_assert(Buckets == 4);

		struct Workspace
		{
			std::vector<block> bucketValues = std::vector<block>(N);
			std::vector<block> tile = std::vector<block>(TileBlocks);
		};

		void init(u64 seed)
		{
			std::mt19937_64 random(seed);
			mOuterToInner.resize(N);
			std::vector<u8> regionToCoordinate(N);
			std::array<u16, OuterLength> coordinateToRegion;
			for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
			{
				std::iota(coordinateToRegion.begin(), coordinateToRegion.end(), u16{ 0 });
				std::shuffle(coordinateToRegion.begin(), coordinateToRegion.end(), random);
				for (u64 coordinate = 0; coordinate < OuterLength; ++coordinate)
				{
					const u64 region = coordinateToRegion[coordinate];
					regionToCoordinate[region * RegionSize + outerBlock] =
						static_cast<u8>(coordinate);
				}
			}

			std::vector<u16> innerPacketToOuterGroup(PacketCount);
			std::vector<u16> groupToPacket(PacketsPerRegion);
			for (u64 region = 0; region < Regions; ++region)
			{
				std::iota(groupToPacket.begin(), groupToPacket.end(), u16{ 0 });
				std::shuffle(groupToPacket.begin(), groupToPacket.end(), random);
				for (u64 outerGroup = 0; outerGroup < PacketsPerRegion; ++outerGroup)
				{
					const u64 innerPacket =
						region * PacketsPerRegion + groupToPacket[outerGroup];
					innerPacketToOuterGroup[innerPacket] = static_cast<u16>(outerGroup);
					for (u64 lane = 0; lane < PacketWidth; ++lane)
					{
						const u64 outerBlock = outerGroup * PacketWidth + lane;
						const u64 coordinate =
							regionToCoordinate[region * RegionSize + outerBlock];
						const u64 outer = outerBlock * OuterLength + coordinate;
						const u64 inner = innerPacket * PacketWidth + lane;
						mOuterToInner[outer] = static_cast<u32>(inner);
					}
				}
			}

			mInnerPacketToSlot24.resize(3 * PacketCount + sizeof(u32));
			mSlotRoute.resize(PacketCount);
			std::array<u32, Buckets> counts;
			counts.fill(static_cast<u32>(TilePackets));
			for (u64 innerPacket = PacketCount; innerPacket-- > 0;)
			{
				const u64 region = innerPacket / PacketsPerRegion;
				const u64 outerGroup = innerPacketToOuterGroup[innerPacket];
				const u64 bucket = outerGroup / TileOuterGroups;
				const u64 localGroup = outerGroup & (TileOuterGroups - 1);
				const u64 slotPacket = bucket * TilePackets + --counts[bucket];
				auto* slot = mInnerPacketToSlot24.data() + 3 * innerPacket;
				slot[0] = static_cast<u8>(slotPacket);
				slot[1] = static_cast<u8>(slotPacket >> 8);
				slot[2] = static_cast<u8>(slotPacket >> 16);

				u64 route = localGroup;
				for (u64 lane = 0; lane < PacketWidth; ++lane)
				{
					const u64 outerBlock = outerGroup * PacketWidth + lane;
					const u64 coordinate =
						regionToCoordinate[region * RegionSize + outerBlock];
					route |= coordinate << (9 + 8 * lane);
				}
				mSlotRoute[slotPacket] = route;
			}
		}

		void validate() const
		{
			std::vector<u8> seen(N, 0);
			for (const u32 inner : mOuterToInner)
			{
				if (inner >= N || seen[inner])
					throw std::runtime_error("packet schedule is not a permutation");
				seen[inner] = 1;
			}
		}

		void gather(
			const block* __restrict inner,
			block* __restrict outer) const noexcept
		{
			for (u64 index = 0; index < N; ++index)
				outer[index] = inner[mOuterToInner[index]];
		}

		void releaseOracle()
		{
			mOuterToInner.clear();
			mOuterToInner.shrink_to_fit();
		}

		void dualEncodeUnchecked(
			const block* __restrict input,
			block* __restrict output,
			const u64* __restrict coefficients,
			Workspace& workspace) const noexcept
		{
			std::array<block, StateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict values = workspace.bucketValues.data();
			const u8* slotBytes = mInnerPacketToSlot24.data() + 3 * PacketCount;

			for (u64 epoch = Epochs; epoch-- > 0;)
			{
				if (epoch + 1 < Epochs)
					transformField<transpose64Avx2>(state.data(), coefficients[epoch]);
				const u64 epochBase = epoch * EpochBlocks;
				for (u64 offset = EpochBlocks; offset != 0; offset -= PacketWidth)
				{
					const u64 innerBase = epochBase + offset - PacketWidth;
					slotBytes -= 3;
					u32 slotPacket;
					std::memcpy(&slotPacket, slotBytes, sizeof(slotPacket));
					block* __restrict destination =
						values + (slotPacket & static_cast<u32>(PacketCount - 1)) * PacketWidth;
					for (u64 lane = PacketWidth; lane-- > 0;)
					{
						const u64 stateLane = (offset - PacketWidth + lane) & (StateBlocks - 1);
						state[stateLane] ^= input[innerBase + lane];
						destination[lane] = state[stateLane];
					}
				}
			}

			block* __restrict tile = workspace.tile.data();
			for (u64 bucket = 0; bucket < Buckets; ++bucket)
			{
				const u64 packetBase = bucket * TilePackets;
				for (u64 packet = 0; packet < TilePackets; ++packet)
				{
					const u64 slotPacket = packetBase + packet;
					const u64 route = mSlotRoute[slotPacket];
					const u64 localGroup = route & (TileOuterGroups - 1);
					const block* __restrict source = values + slotPacket * PacketWidth;
					for (u64 lane = 0; lane < PacketWidth; ++lane)
					{
						const u64 coordinate = (route >> (9 + 8 * lane)) & 255;
						tile[(localGroup * PacketWidth + lane) * OuterLength + coordinate] =
							source[lane];
					}
				}

				const u64 outerBase = bucket * TileOuterBlocks;
				for (u64 offset = 0; offset < TileOuterBlocks; offset += 2)
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * OuterLength,
						tile + (offset + 1) * OuterLength,
						output + (outerBase + offset) * OuterDimension,
						output + (outerBase + offset + 1) * OuterDimension);
			}
		}

	private:
		std::vector<u32> mOuterToInner;
		std::vector<u8> mInnerPacketToSlot24;
		// 9 low bits select one of 512 local outer groups. The next four
		// bytes are the destination BCH coordinates for the packet lanes.
		std::vector<u64> mSlotRoute;
	};

	// An exact factorization of each uniform 8192-position region permutation.
	// Every permutation is represented as 256 independent 32-way permutations,
	// followed by 32 independent 256-way permutations, followed by 256 more
	// 32-way permutations. The factorization changes only the evaluator, not the
	// sampled permutation or its distribution.
	class ExactPermClos256x32
	{
	public:
		static constexpr u64 Groups = 256;
		static constexpr u64 GroupSize = 32;
		static constexpr u64 OuterTileBlocks = 256;
		static constexpr u64 MaxOuterTileBlocks = 512;

		struct Workspace
		{
			std::vector<block> regionMajor = std::vector<block>(N);
			std::vector<block> stage0 = std::vector<block>(RegionSize);
			std::vector<block> stage1 = std::vector<block>(RegionSize);
			std::vector<block> outerTile =
				std::vector<block>(MaxOuterTileBlocks * OuterLength);
		};

		void init(const StructuredSchedule& schedule)
		{
			mFirstDestinationColor.resize(N);
			mMiddleSourceGroup.resize(N);
			mLastSourceColor.resize(N);
			mRegionToCoordinate.resize(N);

			const auto& plan = schedule.regionPlan();
			std::array<u16, RegionSize> destinationForSource;
			std::array<u8, RegionSize> used{};
			std::array<int, Groups> matchedLeft;
			std::array<int, Groups> matchedEdge;
			std::array<u8, Groups> seen;
			std::array<int, Groups> selectedEdge;

			for (u64 region = 0; region < Regions; ++region)
			{
				const u64 base = region * RegionSize;
				for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
				{
					const u32 operation = plan[base + outerBlock];
					const u64 sourcePosition = operation >> 8;
					destinationForSource[sourcePosition] = static_cast<u16>(outerBlock);
					mRegionToCoordinate[base + outerBlock] =
						static_cast<u8>(operation);
				}
				used.fill(0);

				for (u64 color = 0; color < GroupSize; ++color)
				{
					matchedLeft.fill(-1);
					matchedEdge.fill(-1);
					selectedEdge.fill(-1);

					auto augment = [&](auto&& self, int left) -> bool {
						for (u64 inputOffset = 0; inputOffset < GroupSize; ++inputOffset)
						{
							const int edge = left * static_cast<int>(GroupSize) +
								static_cast<int>(inputOffset);
							if (used[edge])
								continue;
							const int right = destinationForSource[edge] /
								static_cast<int>(GroupSize);
							if (seen[right])
								continue;
							seen[right] = 1;
							if (matchedLeft[right] == -1 || self(self, matchedLeft[right]))
							{
								matchedLeft[right] = left;
								matchedEdge[right] = edge;
								return true;
							}
						}
						return false;
					};

					for (u64 left = 0; left < Groups; ++left)
					{
						seen.fill(0);
						if (!augment(augment, static_cast<int>(left)))
							throw std::runtime_error("failed to factor exact region permutation");
					}
					for (u64 right = 0; right < Groups; ++right)
					{
						const int edge = matchedEdge[right];
						const int left = matchedLeft[right];
						if (edge < 0 || left < 0 || selectedEdge[left] != -1)
							throw std::runtime_error("invalid perfect matching in exact factorization");
						selectedEdge[left] = edge;
						used[edge] = 1;
						mMiddleSourceGroup[base + right * GroupSize + color] =
							static_cast<u8>(left);
					}
					for (u64 left = 0; left < Groups; ++left)
					{
						const int edge = selectedEdge[left];
						const u64 inputOffset = static_cast<u64>(edge) & (GroupSize - 1);
						const u64 destination = destinationForSource[edge];
						const u64 right = destination / GroupSize;
						const u64 outputOffset = destination & (GroupSize - 1);
						mFirstDestinationColor[base + left * GroupSize + inputOffset] =
							static_cast<u8>(color);
						mLastSourceColor[base + right * GroupSize + outputOffset] =
							static_cast<u8>(color);
					}
				}
			}
		}

		void validate(const StructuredSchedule& schedule) const
		{
			const auto& plan = schedule.regionPlan();
			std::array<u16, RegionSize> result;
			std::array<u16, RegionSize> first;
			std::array<u16, RegionSize> middle;
			for (u64 region = 0; region < Regions; ++region)
			{
				const u64 base = region * RegionSize;
				for (u64 left = 0; left < Groups; ++left)
					for (u64 inputOffset = 0; inputOffset < GroupSize; ++inputOffset)
						first[left * GroupSize +
							mFirstDestinationColor[base + left * GroupSize + inputOffset]] =
							static_cast<u16>(left * GroupSize + inputOffset);
				for (u64 right = 0; right < Groups; ++right)
					for (u64 color = 0; color < GroupSize; ++color)
						middle[right * GroupSize + color] = first[
							mMiddleSourceGroup[base + right * GroupSize + color] * GroupSize + color];
				for (u64 right = 0; right < Groups; ++right)
					for (u64 outputOffset = 0; outputOffset < GroupSize; ++outputOffset)
						result[right * GroupSize + outputOffset] = middle[
							right * GroupSize +
							mLastSourceColor[base + right * GroupSize + outputOffset]];

				for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
				{
					const u32 operation = plan[base + outerBlock];
					if (result[outerBlock] != (operation >> 8) ||
						mRegionToCoordinate[base + outerBlock] != static_cast<u8>(operation))
						throw std::runtime_error("Clos factorization disagrees with exact permutation");
				}
			}
		}

		void permutationAndOuterTranspose(
			const block* __restrict inner,
			block* __restrict message,
			Workspace& workspace) const noexcept
		{
			regionPermutation(inner, workspace);
			regionMajorOuterTranspose(message, workspace);
		}

		void fusedCheckpointAndPermutation(
			const block* __restrict input,
			const u64* __restrict coefficients,
			Workspace& workspace) const noexcept
		{
			constexpr u64 epochsPerRegion = RegionSize / EpochBlocks;
			static_assert(RegionSize % EpochBlocks == 0);
			std::array<block, StateBlocks> state;
			state.fill(block(0, 0));
			block* __restrict stage0 = workspace.stage0.data();
			block* __restrict stage1 = workspace.stage1.data();
			block* __restrict regionMajor = workspace.regionMajor.data();

			for (u64 region = Regions; region-- > 0;)
			{
				const u64 base = region * RegionSize;
				const u8* __restrict first = mFirstDestinationColor.data() + base;
				for (u64 localEpoch = epochsPerRegion; localEpoch-- > 0;)
				{
					const u64 epoch = region * epochsPerRegion + localEpoch;
					if (epoch + 1 < Epochs)
						transformField<transpose64Avx2>(state.data(), coefficients[epoch]);
					const u64 epochBase = localEpoch * EpochBlocks;
					for (u64 offset = EpochBlocks; offset-- > 0;)
					{
						const u64 sourcePosition = epochBase + offset;
						const u64 lane = offset & (StateBlocks - 1);
						state[lane] ^= input[base + sourcePosition];
						const u64 groupBase = sourcePosition & ~(GroupSize - 1);
						stage0[groupBase + first[sourcePosition]] = state[lane];
					}
				}

				const u8* __restrict middle = mMiddleSourceGroup.data() + base;
				const u8* __restrict last = mLastSourceColor.data() + base;
				for (u64 right = 0; right < Groups; ++right)
				{
					const u64 groupBase = right * GroupSize;
					for (u64 color = 0; color < GroupSize; ++color)
						stage1[groupBase + color] = stage0[
							middle[groupBase + color] * GroupSize + color];
				}
				block* __restrict destination = regionMajor + base;
				for (u64 right = 0; right < Groups; ++right)
				{
					const u64 outputBase = right * GroupSize;
					for (u64 outputOffset = 0; outputOffset < GroupSize; ++outputOffset)
						destination[outputBase + outputOffset] = stage1[
							outputBase + last[outputBase + outputOffset]];
				}
			}
		}

		void fusedCheckpointPermutationAndOuterTranspose(
			const block* __restrict input,
			block* __restrict message,
			const u64* __restrict coefficients,
			Workspace& workspace) const noexcept
		{
			fusedCheckpointAndPermutation(input, coefficients, workspace);
			regionMajorOuterTranspose(message, workspace);
		}

		void regionPermutation(
			const block* __restrict inner,
			Workspace& workspace) const noexcept
		{
			block* __restrict stage0 = workspace.stage0.data();
			block* __restrict stage1 = workspace.stage1.data();
			block* __restrict regionMajor = workspace.regionMajor.data();

			for (u64 region = 0; region < Regions; ++region)
			{
				const u64 base = region * RegionSize;
				const block* __restrict source = inner + base;
				const u8* __restrict first = mFirstDestinationColor.data() + base;
				const u8* __restrict middle = mMiddleSourceGroup.data() + base;
				const u8* __restrict last = mLastSourceColor.data() + base;

				for (u64 left = 0; left < Groups; ++left)
				{
					const u64 groupBase = left * GroupSize;
					for (u64 inputOffset = 0; inputOffset < GroupSize; ++inputOffset)
						stage0[groupBase + first[groupBase + inputOffset]] =
							source[groupBase + inputOffset];
				}
				for (u64 right = 0; right < Groups; ++right)
				{
					const u64 groupBase = right * GroupSize;
					for (u64 color = 0; color < GroupSize; ++color)
						stage1[groupBase + color] = stage0[
							middle[groupBase + color] * GroupSize + color];
				}
				block* __restrict destination = regionMajor + base;
				for (u64 right = 0; right < Groups; ++right)
				{
					const u64 outputBase = right * GroupSize;
					for (u64 outputOffset = 0; outputOffset < GroupSize; ++outputOffset)
						destination[outputBase + outputOffset] = stage1[
							outputBase + last[outputBase + outputOffset]];
				}
			}
		}

		void regionMajorOuterTranspose(
			block* __restrict message,
			Workspace& workspace) const noexcept
		{
			regionMajorOuterTransposeTiled<OuterTileBlocks>(message, workspace);
		}

		template<u64 TileOuterBlocks>
		void regionMajorOuterTransposeTiled(
			block* __restrict message,
			Workspace& workspace) const noexcept
		{
			static_assert(TileOuterBlocks >= 2 &&
				(TileOuterBlocks & (TileOuterBlocks - 1)) == 0);
			static_assert(TileOuterBlocks <= MaxOuterTileBlocks);
			static_assert(OuterBlocks % TileOuterBlocks == 0);
			const block* __restrict regionMajor = workspace.regionMajor.data();
			block* __restrict tile = workspace.outerTile.data();
			for (u64 outerBase = 0; outerBase < OuterBlocks;
				outerBase += TileOuterBlocks)
			{
				for (u64 region = 0; region < Regions; ++region)
				{
					const u64 base = region * RegionSize + outerBase;
					const block* __restrict source = regionMajor + base;
					const u8* __restrict coordinate = mRegionToCoordinate.data() + base;
					for (u64 offset = 0; offset < TileOuterBlocks; ++offset)
						tile[offset * OuterLength + coordinate[offset]] = source[offset];
				}
				for (u64 offset = 0; offset < TileOuterBlocks; offset += 2)
					ExtendedBch256x128Eq3::transposeBlock2(
						tile + offset * OuterLength,
						tile + (offset + 1) * OuterLength,
						message + (outerBase + offset) * OuterDimension,
						message + (outerBase + offset + 1) * OuterDimension);
			}
		}

	private:
		// Group-major inverse of the first local permutation: each source
		// offset maps directly to its matching color. This lets checkpoint
		// finalization feed the first stage without materializing its output.
		std::vector<u8> mFirstDestinationColor;
		std::vector<u8> mMiddleSourceGroup;
		std::vector<u8> mLastSourceColor;
		std::vector<u8> mRegionToCoordinate;
	};

	void outerTranspose(const block* __restrict outer, block* __restrict message) noexcept
	{
		for (u64 outerBlock = 0; outerBlock < OuterBlocks; ++outerBlock)
		{
			ExtendedBch256x128Eq3::transposeBlock(
				outer + outerBlock * OuterLength,
				message + outerBlock * OuterDimension);
		}
	}

	void outerTranspose2(const block* __restrict outer, block* __restrict message) noexcept
	{
		for (u64 outerBlock = 0; outerBlock < OuterBlocks; outerBlock += 2)
		{
			ExtendedBch256x128Eq3::transposeBlock2(
				outer + outerBlock * OuterLength,
				outer + (outerBlock + 1) * OuterLength,
				message + outerBlock * OuterDimension,
				message + (outerBlock + 1) * OuterDimension);
		}
	}

	u64 splitmix64(u64& state) noexcept
	{
		u64 value = (state += 0x9e3779b97f4a7c15ULL);
		value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
		value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
		return value ^ (value >> 31);
	}

	block randomBlock(u64& state) noexcept
	{
		const u64 low = splitmix64(state);
		const u64 high = splitmix64(state);
		return block(low, high);
	}

	bool equalBlocks(const block* left, const block* right, u64 size) noexcept
	{
		return std::memcmp(left, right, size * sizeof(block)) == 0;
	}

	u64 checksum(const block* values, u64 size) noexcept
	{
		block sum(0, 0);
		for (u64 index = 0; index < size; ++index)
			sum ^= values[index];
		return sum.get<u64>()[0] ^ sum.get<u64>()[1];
	}

	void selfTest(const StructuredSchedule& schedule, const std::vector<u64>& coefficients)
	{
		u64 state = 0x726966666c652d31ULL;
		for (u64 trial = 0; trial < 4; ++trial)
		{
			std::array<block, OuterLength> word;
			std::array<block, OuterDimension> optimized;
			std::array<block, OuterDimension> reference;
			for (auto& value : word)
				value = randomBlock(state);
			ExtendedBch256x128Eq3::transposeBlock(word.data(), optimized.data());
			ExtendedBch256x128Eq3::transposeReference(word.data(), reference.data());
			if (!equalBlocks(optimized.data(), reference.data(), reference.size()))
				throw std::runtime_error("generated BCH transpose disagrees with dense reference");
		}
		{
			std::array<block, OuterLength> word0;
			std::array<block, OuterLength> word1;
			std::array<block, OuterDimension> packed0;
			std::array<block, OuterDimension> packed1;
			std::array<block, OuterDimension> reference0;
			std::array<block, OuterDimension> reference1;
			for (auto& value : word0)
				value = randomBlock(state);
			for (auto& value : word1)
				value = randomBlock(state);
			ExtendedBch256x128Eq3::transposeBlock2(
				word0.data(), word1.data(), packed0.data(), packed1.data());
			ExtendedBch256x128Eq3::transposeBlock(word0.data(), reference0.data());
			ExtendedBch256x128Eq3::transposeBlock(word1.data(), reference1.data());
			if (!equalBlocks(packed0.data(), reference0.data(), OuterDimension) ||
				!equalBlocks(packed1.data(), reference1.data(), OuterDimension))
				throw std::runtime_error("packed BCH transpose disagrees with scalar circuit");
		}

		constexpr u64 testBlocks = 3 * EpochBlocks;
		std::array<block, testBlocks> optimizedInner;
		for (auto& value : optimizedInner)
			value = randomBlock(state);
		auto referenceInner = optimizedInner;
		checkpointOptimized(optimizedInner.data(), testBlocks, coefficients.data());
		checkpointReference(referenceInner.data(), testBlocks, coefficients.data());
		if (!equalBlocks(optimizedInner.data(), referenceInner.data(), testBlocks))
			throw std::runtime_error("optimized checkpoint inner disagrees with scalar reference");

		std::vector<block> source(N);
		std::vector<block> stagedWord(N);
		std::vector<block> stagedMessage(K);
		std::vector<block> fusedMessage(K);
		std::vector<block> tiledMessage(K);
		std::vector<block> tiledWorkspace(1024 * OuterLength);
		for (auto& value : source)
			value = randomBlock(state);
		schedule.gather(source.data(), stagedWord.data());
		outerTranspose(stagedWord.data(), stagedMessage.data());
		schedule.gatherAndOuterTranspose(source.data(), fusedMessage.data());
		if (!equalBlocks(stagedMessage.data(), fusedMessage.data(), K))
			throw std::runtime_error("fused structured gather disagrees with staged path");
		schedule.gatherAndOuterTranspose2(source.data(), fusedMessage.data());
		if (!equalBlocks(stagedMessage.data(), fusedMessage.data(), K))
			throw std::runtime_error("packed fused gather disagrees with staged path");
		for (const u64 batchSize : { u64{ 64 } })
		{
			schedule.tiledGatherAndOuterTranspose(
				source.data(), tiledMessage.data(), tiledWorkspace.data(), batchSize);
			if (!equalBlocks(stagedMessage.data(), tiledMessage.data(), K))
				throw std::runtime_error("region-tiled gather disagrees with staged path");
			schedule.tiledPrefetchGatherAndOuterTranspose(
				source.data(), tiledMessage.data(), tiledWorkspace.data(), batchSize);
			if (!equalBlocks(stagedMessage.data(), tiledMessage.data(), K))
				throw std::runtime_error("prefetched region-tiled gather disagrees with staged path");
			schedule.tiledPackedGatherAndOuterTranspose(
				source.data(), tiledMessage.data(), tiledWorkspace.data(), batchSize);
			if (!equalBlocks(stagedMessage.data(), tiledMessage.data(), K))
				throw std::runtime_error("packed region-tiled gather disagrees with staged path");
			schedule.tiledPrefetchPackedGatherAndOuterTranspose(
				source.data(), tiledMessage.data(), tiledWorkspace.data(), batchSize);
			if (!equalBlocks(stagedMessage.data(), tiledMessage.data(), K))
				throw std::runtime_error("prefetched packed region-tiled gather disagrees with staged path");
		}
		benchmarkSink ^= checksum(fusedMessage.data(), fusedMessage.size());
	}

	template<typename Function>
	double milliseconds(Function&& function)
	{
		const auto begin = std::chrono::steady_clock::now();
		function();
		const auto end = std::chrono::steady_clock::now();
		return std::chrono::duration<double, std::milli>(end - begin).count();
	}

	double median(std::vector<double> samples)
	{
		std::sort(samples.begin(), samples.end());
		return samples[samples.size() / 2];
	}

	template<typename Prepare, typename Function>
	std::vector<double> benchmark(u64 trials, Prepare&& prepare, Function&& function)
	{
		std::vector<double> samples;
		samples.reserve(trials);
		for (u64 trial = 0; trial < trials; ++trial)
		{
			prepare();
			samples.push_back(milliseconds(function));
		}
		return samples;
	}

	void printSamples(const char* name, const std::vector<double>& samples)
	{
		std::cout << name << "_median_ms=" << median(samples) << ' ' << name << "_samples_ms=";
		for (const auto sample : samples)
			std::cout << sample << ',';
		std::cout << '\n';
	}

	void verifyRm2SubS19Transpose()
	{
		constexpr u64 testEpochs = 5;
		constexpr u64 testBlocks = 128 * testEpochs;
		u64 randomState = 0x524d325355423136ULL;
		std::vector<block> input(testBlocks);
		for (auto& value : input)
			value = randomBlock(randomState);

		RiffleRm2SubS19Transpose inner;
		inner.init(testEpochs, 0x5452414e53504f53ULL);
		std::vector<block> optimized(testBlocks);
		inner.emitReverse(input.data(), testBlocks,
			[&](u64 index, block value) { optimized[index] = value; });

		std::vector<block> reference(testBlocks);
		std::array<block, 19> state;
		state.fill(block(0, 0));
		for (u64 epoch = testEpochs; epoch-- > 0;)
		{
			const block* node = input.data() + 128 * epoch;
			for (u64 point = 0; point < 128; ++point)
			{
				block addend(0, 0);
				u32 column = detail::Rm2Sub19Columns[point];
				while (column)
				{
					const unsigned bit = std::countr_zero(column);
					addend ^= state[bit];
					column &= column - 1;
				}
				reference[128 * epoch + point] = node[point] ^ addend;
			}
			if (epoch == 0)
				break;

			std::array<block, 19> syndrome;
			syndrome.fill(block(0, 0));
			for (u64 point = 0; point < 128; ++point)
			{
				u32 column = detail::Rm2Sub19Columns[point];
				while (column)
				{
					const unsigned bit = std::countr_zero(column);
					syndrome[bit] ^= node[point];
					column &= column - 1;
				}
			}
			if (epoch + 1 == testEpochs)
				state = syndrome;
			else
			{
				std::array<block, 19> multiplied;
				const auto& masks = inner.schedules()[epoch].rowMasks;
				for (unsigned row = 0; row < 19; ++row)
				{
					block value(0, 0);
					u32 mask = masks[row];
					while (mask)
					{
						const unsigned bit = std::countr_zero(mask);
						value ^= state[bit];
						mask &= mask - 1;
					}
					multiplied[row] = value ^ syndrome[row];
				}
				state = multiplied;
			}
		}
		if (!equalBlocks(reference.data(), optimized.data(), testBlocks))
		{
			for (u64 index = 0; index < testBlocks; ++index)
				if (reference[index] != optimized[index])
					throw std::runtime_error(
						"optimized RM2Sub s19 transpose disagrees with dense reference at block " +
						std::to_string(index));
			throw std::runtime_error("optimized RM2Sub s19 transpose disagrees with dense reference");
		}
	}

#ifndef _WIN32
	bool adviseHugePages(void* data, u64 bytes) noexcept
	{
		const u64 page = static_cast<u64>(sysconf(_SC_PAGESIZE));
		const auto address = reinterpret_cast<std::uintptr_t>(data);
		const auto begin = (address + page - 1) & ~(page - 1);
		const auto end = (address + bytes) & ~(page - 1);
		if (end <= begin || madvise(reinterpret_cast<void*>(begin), end - begin, MADV_HUGEPAGE))
			return false;
#ifdef MADV_COLLAPSE
		return madvise(reinterpret_cast<void*>(begin), end - begin, MADV_COLLAPSE) == 0;
#else
		return true;
#endif
	}
#endif
}

int main(int argc, char** argv)
{
	try
	{
		const char* benchmarkCpuText = std::getenv("RIFFLE_BENCH_CPU");
		const u64 benchmarkCpu = benchmarkCpuText
			? std::strtoull(benchmarkCpuText, nullptr, 0)
			: u64{ 2 };
#ifdef _WIN32
		SetThreadAffinityMask(GetCurrentThread(), DWORD_PTR{ 1 } << benchmarkCpu);
		SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);
#else
		cpu_set_t affinity;
		CPU_ZERO(&affinity);
		CPU_SET(benchmarkCpu, &affinity);
		if (pthread_setaffinity_np(pthread_self(), sizeof(affinity), &affinity) != 0)
			throw std::runtime_error("failed to pin benchmark thread");
#endif
		const u64 trials = argc > 1 ? std::strtoull(argv[1], nullptr, 0) : 9;
		const bool selectedOnly = argc > 2 && std::string(argv[2]) == "selected";
		const bool bucketOnly = argc > 2 && std::string(argv[2]) == "bucket";
		const bool bucketSelectedOnly =
			argc > 2 && std::string(argv[2]) == "bucket-selected";
		const bool packedTuningOnly =
			argc > 2 && std::string(argv[2]) == "packed-tuning";
		const bool productionOnly =
			argc > 2 && std::string(argv[2]) == "production";
		const bool productionCounterOnly =
			argc > 2 && std::string(argv[2]) == "production-counter";
		const bool phaseProfileOnly =
			argc > 2 && std::string(argv[2]) == "phase-profile";
		const bool clos256x32Only =
			argc > 2 && std::string(argv[2]) == "clos256x32";
		const bool packet4Only =
			argc > 2 && std::string(argv[2]) == "packet4";
		const bool packet4CounterOnly =
			argc > 2 && std::string(argv[2]) == "packet4-counter";
		const bool rm2SubS19Only =
			argc > 2 && std::string(argv[2]) == "rm2sub-s19";
#if defined(LIBOTE_RIFFLE_BCH_CIRCUIT_PAIR_BENCH)
		const bool outerCircuitPairOnly =
			argc > 2 && std::string(argv[2]) == "outer-circuit-pair";
#endif
		if (trials < 3 || !(trials & 1))
			throw std::invalid_argument("trials must be an odd integer at least three");

		StructuredSchedule schedule;
		const double setupMs = milliseconds([&] { schedule.init(0x5354525543543235ULL); });
		schedule.validate();

		u64 randomState = 0x4649454c4443484bULL;
		std::vector<u64> coefficients(Epochs - 1);
		for (auto& coefficient : coefficients)
		{
			coefficient = splitmix64(randomState);
			coefficient |= static_cast<u64>(coefficient == 0);
		}
		selfTest(schedule, coefficients);

		std::vector<block> source(N);
		std::vector<block> working(N);
		std::vector<block> outerWord(N);
		std::vector<block> message(K);
		std::vector<block> tiledWorkspace(1024 * OuterLength);
		for (auto& value : source)
			value = randomBlock(randomState);

		if (rm2SubS19Only)
		{
			constexpr u64 rm2Epochs = N / RiffleRm2SubS19Transpose::stepBlocks;
			verifyRm2SubS19Transpose();
			std::cerr << "rm2sub: dense inner PASS\n";
			RiffleRm2SubS19Transpose inner;
			inner.init(rm2Epochs, 0x524d32434f454646ULL);
			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL,
				coefficients.data(), coefficients.size(), true);
			RiffleExactPermFieldCheckpoint::Workspace workspace;
			std::vector<block> innerWord(N);
			std::vector<block> stagedMessage(K);
			std::vector<block> fusedMessage(K);

			inner.emitReverse(source.data(), N,
				[&](u64 index, block value) { innerWord[index] = value; });
			std::cerr << "rm2sub: full inner materialization PASS\n";
			for (u64 innerIndex = 0; innerIndex < N; ++innerIndex)
				outerWord[promoted.innerToOuter()[innerIndex]] = innerWord[innerIndex];
			for (u64 outer = 0; outer < OuterBlocks; outer += 2)
			{
				for (u64 lane = 0; lane < 2; ++lane)
				{
					block* word = outerWord.data() + (outer + lane) * OuterLength;
					const auto& fanout = promoted.parityFanouts()[outer + lane];
					block pivot = word[fanout.targets[0]];
					for (unsigned index = 1; index < fanout.targets.size(); ++index)
						pivot ^= word[fanout.targets[index]];
					for (const u8 sourceIndex : fanout.sources) word[sourceIndex] ^= pivot;
				}
				ExtendedBch256x128Eq3::transposeBlock2(
					outerWord.data() + outer * OuterLength,
					outerWord.data() + (outer + 1) * OuterLength,
					stagedMessage.data() + outer * OuterDimension,
					stagedMessage.data() + (outer + 1) * OuterDimension);
			}
			std::cerr << "rm2sub: staged fanout31x33 complete PASS\n";
			promoted.dualEncodePackedInnerUnchecked<32>(
				source.data(), fusedMessage.data(), workspace, inner);
			std::cerr << "rm2sub: fused complete PASS\n";
			if (!equalBlocks(stagedMessage.data(), fusedMessage.data(), K))
				throw std::runtime_error(
					"integrated RM2Sub s16 encoder disagrees with staged complete encoder");

			auto integrated = benchmark(trials, [] {}, [&] {
				promoted.dualEncodePackedInnerUnchecked<32>(
					source.data(), fusedMessage.data(), workspace, inner);
			});
			auto materializedInner = benchmark(trials, [] {}, [&] {
				inner.emitReverse(source.data(), N,
					[&](u64 index, block value) { innerWord[index] = value; });
			});
			benchmarkSink ^= checksum(fusedMessage.data(), fusedMessage.size());
			benchmarkSink ^= checksum(innerWord.data(), innerWord.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ParityFanout31x33-BCHPerm-BitShuffle-RM2SubS19"
				<< " mode=rm2sub-s19"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " t=128 s=19 epochs=" << rm2Epochs
				<< " trials=" << trials
				<< " correctness=dense-inner+staged-fanout31x33-complete-PASS\n";
			printSamples("integrated_full", integrated);
			printSamples("materialized_inner", materializedInner);
			std::cout << "checksum=0x" << std::hex << benchmarkSink << std::dec << '\n';
			return 0;
		}

#if defined(LIBOTE_RIFFLE_BCH_CIRCUIT_PAIR_BENCH)
		if (outerCircuitPairOnly)
		{
			std::vector<block> candidateMessage(K);
			auto runSeed114 = [&] {
				return milliseconds([&] {
					for (u64 outer = 0; outer < OuterBlocks; outer += 2)
						ExtendedBch256x128Eq3::transposeBlock2(
							source.data() + outer * OuterLength,
							source.data() + (outer + 1) * OuterLength,
							message.data() + outer * OuterDimension,
							message.data() + (outer + 1) * OuterDimension);
				});
			};
			auto runSeed57 = [&] {
				return milliseconds([&] {
					for (u64 outer = 0; outer < OuterBlocks; outer += 2)
						detail::seed57::extendedBch256x128Eq3TransposeCircuit2(
							source.data() + outer * OuterLength,
							source.data() + (outer + 1) * OuterLength,
							candidateMessage.data() + outer * OuterDimension,
							candidateMessage.data() + (outer + 1) * OuterDimension);
				});
			};
			runSeed114();
			runSeed57();
			if (!equalBlocks(message.data(), candidateMessage.data(), K))
				throw std::runtime_error("seed-57 outer circuit disagrees with seed 114");
			std::vector<double> seed114;
			std::vector<double> seed57;
			seed114.reserve(trials);
			seed57.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					seed57.push_back(runSeed57());
					seed114.push_back(runSeed114());
				}
				else
				{
					seed114.push_back(runSeed114());
					seed57.push_back(runSeed57());
				}
			}
			benchmarkSink ^= checksum(message.data(), message.size());
			benchmarkSink ^= checksum(candidateMessage.data(), candidateMessage.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=ExtendedBch256x128Eq3"
				<< " mode=outer-circuit-pair"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " outer_blocks=" << OuterBlocks
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("seed114_xors2934", seed114);
			printSamples("seed57_xors2949", seed57);
			std::cout << "seed57_over_seed114="
				<< median(seed57) / median(seed114)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}
#endif

		if (packet4Only || packet4CounterOnly)
		{
			PacketShuffle4Schedule packet;
			const double packetSetupMs = milliseconds([&] {
				packet.init(0x5354525543543235ULL);
			});
			packet.validate();
			PacketShuffle4Schedule::Workspace packetWorkspace;

			std::memcpy(working.data(), source.data(), N * sizeof(block));
			checkpointOptimized(working.data(), N, coefficients.data());
			packet.gather(working.data(), outerWord.data());
			outerTranspose2(outerWord.data(), message.data());
			const auto packetReference = message;
			packet.dualEncodeUnchecked(
				source.data(), message.data(), coefficients.data(), packetWorkspace);
			if (!equalBlocks(packetReference.data(), message.data(), K))
				throw std::runtime_error("packet evaluator disagrees with its flat oracle");
			packet.releaseOracle();
			if (packet4CounterOnly)
			{
				packet.dualEncodeUnchecked(
					source.data(), message.data(), coefficients.data(), packetWorkspace);
				const auto packetFull = benchmark(
					trials,
					[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
					[&] {
						packet.dualEncodeUnchecked(
							working.data(), message.data(), coefficients.data(), packetWorkspace);
					});
				if (!equalBlocks(packetReference.data(), message.data(), K))
					throw std::runtime_error("isolated packet evaluator disagrees with its oracle");
				benchmarkSink ^= checksum(message.data(), message.size());
				std::cout << std::fixed << std::setprecision(6)
					<< "construction=Riffle-BCHPerm-TransposePacketShuffle-FieldCheckpoint-g4"
					<< " mode=packet4-counter"
					<< " benchmark_cpu=" << benchmarkCpu
					<< " K=" << K << " N=" << N
					<< " packet_width=4 packets_per_region=2048"
					<< " setup_ms=" << packetSetupMs
					<< " hot_schedule_mib=5.500000"
					<< " workspace_mib=40.000000"
					<< " trials=" << trials << " correctness=PASS\n";
				printSamples("packet4_full", packetFull);
				std::cout << "sink=" << benchmarkSink << '\n';
				return 0;
			}

			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size());
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			promoted.dualEncodeUnchecked(
				source.data(), message.data(), promotedWorkspace);

			auto runProduction = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					promoted.dualEncodeUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
			};
			auto runPacket = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					packet.dualEncodeUnchecked(
						working.data(), message.data(), coefficients.data(), packetWorkspace);
				});
			};
			std::vector<double> productionFull;
			std::vector<double> packetFull;
			productionFull.reserve(trials);
			packetFull.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					runPacket();
					packetFull.push_back(runPacket());
					runProduction();
					productionFull.push_back(runProduction());
				}
				else
				{
					runProduction();
					productionFull.push_back(runProduction());
					runPacket();
					packetFull.push_back(runPacket());
				}
			}

			if (!equalBlocks(packetReference.data(), message.data(), K))
				throw std::runtime_error("timed packet evaluator disagrees with its flat oracle");
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-BCHPerm-TransposePacketShuffle-FieldCheckpoint-g4"
				<< " mode=packet4"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " K=" << K << " N=" << N
				<< " packet_width=4 packets_per_region=2048"
				<< " sample_policy=self-warmed-alternating-pairs"
				<< " setup_ms=" << packetSetupMs
				<< " hot_schedule_mib=5.500000"
				<< " workspace_mib=40.000000"
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("production_full", productionFull);
			printSamples("packet4_full", packetFull);
			std::cout << "packet4_over_production="
				<< median(packetFull) / median(productionFull)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (clos256x32Only)
		{
			ExactPermClos256x32 clos;
			const double closSetupMs = milliseconds([&] { clos.init(schedule); });
			clos.validate(schedule);
			ExactPermClos256x32::Workspace closWorkspace;

			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size());
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			auto runPromoted = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					promoted.dualEncodeUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
			};
			auto runClos = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					clos.permutationAndOuterTranspose(
						working.data(), message.data(), closWorkspace);
				});
			};
			auto runFusedClos = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					clos.fusedCheckpointPermutationAndOuterTranspose(
						working.data(), message.data(), coefficients.data(), closWorkspace);
				});
			};

			std::memcpy(working.data(), source.data(), N * sizeof(block));
			checkpointOptimized(working.data(), N, coefficients.data());
			schedule.gather(working.data(), outerWord.data());
			outerTranspose2(outerWord.data(), message.data());
			const auto referenceMessage = message;
			runPromoted();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("promoted evaluator disagrees with exact oracle");
			runClos();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("Clos evaluator disagrees with exact oracle");
			runFusedClos();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("fused Clos evaluator disagrees with exact oracle");

			std::memcpy(working.data(), source.data(), N * sizeof(block));
			checkpointOptimized(working.data(), N, coefficients.data());
			clos.regionPermutation(working.data(), closWorkspace);
			clos.regionMajorOuterTranspose(message.data(), closWorkspace);
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("split Clos evaluator disagrees with exact oracle");

			const auto checkpointOnly = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] { checkpointOptimized(working.data(), N, coefficients.data()); });
			const auto closRegionPermutation = benchmark(
				trials, [] {}, [&] {
					clos.regionPermutation(working.data(), closWorkspace);
				});
			const auto fusedCheckpointPermutation = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					clos.fusedCheckpointAndPermutation(
						working.data(), coefficients.data(), closWorkspace);
				});
			struct TileSamples
			{
				u64 tileOuterBlocks;
				std::vector<double> outer;
			};
			auto runTile = [&]<u64 TileOuterBlocks>(
				std::integral_constant<u64, TileOuterBlocks>) {
				clos.regionMajorOuterTransposeTiled<TileOuterBlocks>(
					message.data(), closWorkspace);
				if (!equalBlocks(referenceMessage.data(), message.data(), K))
					throw std::runtime_error("tiled Clos handoff disagrees with exact oracle");
				return TileSamples{
					TileOuterBlocks,
					benchmark(trials, [] {}, [&] {
						clos.regionMajorOuterTransposeTiled<TileOuterBlocks>(
							message.data(), closWorkspace);
					})
				};
			};
			std::vector<TileSamples> tileSamples;
			tileSamples.emplace_back(runTile(std::integral_constant<u64, 32>{}));
			tileSamples.emplace_back(runTile(std::integral_constant<u64, 64>{}));
			tileSamples.emplace_back(runTile(std::integral_constant<u64, 128>{}));
			tileSamples.emplace_back(runTile(std::integral_constant<u64, 256>{}));
			tileSamples.emplace_back(runTile(std::integral_constant<u64, 512>{}));
			const auto promotedFull = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					promoted.dualEncodeUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
			const auto closFull = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					clos.permutationAndOuterTranspose(
						working.data(), message.data(), closWorkspace);
				});
			const auto fusedClosFull = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					clos.fusedCheckpointPermutationAndOuterTranspose(
						working.data(), message.data(), coefficients.data(), closWorkspace);
				});

			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-Clos256x32-FieldCheckpoint"
				<< " mode=clos256x32"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " K=" << K << " N=" << N
				<< " setup_ms=" << closSetupMs
				<< " setup_schedule_mib=8.000000"
				<< " sweep_workspace_mib=34.250000"
				<< " local_stages=256x32,32x256,256x32"
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("checkpoint_only", checkpointOnly);
			printSamples("clos_region_permutation", closRegionPermutation);
			printSamples("fused_checkpoint_permutation", fusedCheckpointPermutation);
			for (const auto& sample : tileSamples)
				printSamples(
					("clos_outer_t" + std::to_string(sample.tileOuterBlocks)).c_str(),
					sample.outer);
			printSamples("promoted_full", promotedFull);
			printSamples("clos_full", closFull);
			printSamples("fused_clos_full", fusedClosFull);
			std::cout << "clos_over_promoted="
				<< median(closFull) / median(promotedFull)
				<< " fused_clos_over_promoted="
				<< median(fusedClosFull) / median(promotedFull)
				<< " component_sum_ms="
				<< median(checkpointOnly) + median(closRegionPermutation) +
					median(tileSamples[3].outer)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (productionCounterOnly)
		{
			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size());
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			promoted.dualEncodeUnchecked(
				source.data(), message.data(), promotedWorkspace);
			std::vector<double> productionFull;
			productionFull.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
				productionFull.push_back(milliseconds([&] {
					promoted.dualEncodeUnchecked(
						source.data(), message.data(), promotedWorkspace);
				}));
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-FieldCheckpoint-v1"
				<< " mode=production-counter"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " K=" << K << " N=" << N
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("production_full", productionFull);
			std::cout << "sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (phaseProfileOnly)
		{
			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size());
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			promoted.dualEncodeUnchecked(
				source.data(), message.data(), promotedWorkspace);
			const auto referenceMessage = message;
			std::vector<double> stateValue;
			std::vector<double> tileScatter;
			std::vector<double> outerTranspose;
			std::vector<double> profiledFull;
			stateValue.reserve(trials);
			tileScatter.reserve(trials);
			outerTranspose.reserve(trials);
			profiledFull.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				std::array<double, 3> phaseMs{};
				auto previous = std::chrono::steady_clock::now();
				const auto begin = previous;
				auto boundary = [&](u64 phase) {
					const auto now = std::chrono::steady_clock::now();
					phaseMs[phase] +=
						std::chrono::duration<double, std::milli>(now - previous).count();
					previous = now;
				};
				promoted.dualEncodeProfiledUnchecked(
					source.data(), message.data(), promotedWorkspace, boundary);
				const auto end = std::chrono::steady_clock::now();
				stateValue.push_back(phaseMs[0]);
				tileScatter.push_back(phaseMs[1]);
				outerTranspose.push_back(phaseMs[2]);
				profiledFull.push_back(
					std::chrono::duration<double, std::milli>(end - begin).count());
			}
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("profiled evaluator disagrees with production evaluator");
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-FieldCheckpoint-v1"
				<< " mode=phase-profile"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " K=" << K << " N=" << N
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("state_value", stateValue);
			printSamples("tile_scatter", tileScatter);
			printSamples("outer_transpose", outerTranspose);
			printSamples("profiled_full", profiledFull);
			std::cout << "phase_sum_median_ms="
				<< median(stateValue) + median(tileScatter) + median(outerTranspose)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (productionOnly)
		{
			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size());
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			auto runBaseline = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					schedule.optimizedGatherAndOuterTranspose(
						working.data(), message.data(), tiledWorkspace.data());
				});
			};
			auto runProduction = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					promoted.dualEncodeUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
			};
			runBaseline();
			const auto referenceMessage = message;
			runProduction();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("production evaluator disagrees with frozen evaluator");
			std::vector<double> baselineFull;
			std::vector<double> productionFull;
			baselineFull.reserve(trials);
			productionFull.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					productionFull.push_back(runProduction());
					baselineFull.push_back(runBaseline());
				}
				else
				{
					baselineFull.push_back(runBaseline());
					productionFull.push_back(runProduction());
				}
			}
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-FieldCheckpoint-v1"
				<< " mode=production"
				<< " benchmark_cpu=" << benchmarkCpu
				<< " K=" << K << " N=" << N
				<< " tile_mib=8.000000"
				<< " setup_schedule_mib=12.000000"
				<< " workspace_mib=40.000000"
				<< " outer_xors="
				<< ExtendedBch256x128Eq3::transposedEncoderXorCount
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("baseline_full", baselineFull);
			printSamples("production_full", productionFull);
			std::cout << "production_over_baseline="
				<< median(productionFull) / median(baselineFull)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (packedTuningOnly)
		{
			std::memcpy(working.data(), source.data(), N * sizeof(block));
			checkpointOptimized(working.data(), N, coefficients.data());
			schedule.gather(working.data(), outerWord.data());
			outerTranspose2(outerWord.data(), message.data());
			const auto referenceMessage = message;

			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size(), true);
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			struct PackedSamples
			{
				u64 distance;
				std::vector<double> full;
			};
			auto runDistance = [&]<u64 Distance>(
				std::integral_constant<u64, Distance>) {
				promoted.dualEncodePackedSchedulesUnchecked<Distance>(
					source.data(), message.data(), promotedWorkspace);
				if (!equalBlocks(referenceMessage.data(), message.data(), K))
					throw std::runtime_error("packed schedule tuning path disagrees with staged oracle");
				return PackedSamples{
					Distance,
					benchmark(
						trials,
						[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
						[&] {
							promoted.dualEncodePackedSchedulesUnchecked<Distance>(
								working.data(), message.data(), promotedWorkspace);
						})
				};
			};
			std::vector<PackedSamples> samples;
			samples.emplace_back(runDistance(std::integral_constant<u64, 0>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 8>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 16>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 24>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 32>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 40>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 48>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 64>{}));
			samples.emplace_back(runDistance(std::integral_constant<u64, 96>{}));
			promoted.dualEncodePairPackedUnchecked(
				source.data(), message.data(), promotedWorkspace);
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("pair-packed schedule disagrees with staged oracle");
			const auto pairPacked = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					promoted.dualEncodePairPackedUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
#ifndef _WIN32
			const bool hugeValue = adviseHugePages(
				promotedWorkspace.bucketValues.data(),
				promotedWorkspace.bucketValues.size() * sizeof(block));
			const bool hugeTile = adviseHugePages(
				promotedWorkspace.tile.data(),
				promotedWorkspace.tile.size() * sizeof(block));
			const auto hugePageFull = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					promoted.dualEncodeUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
			RiffleExactPermFieldCheckpoint::Workspace normalWorkspace;
			auto runWorkspace = [&](RiffleExactPermFieldCheckpoint::Workspace& workspace) {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					promoted.dualEncodeUnchecked(working.data(), message.data(), workspace);
				});
			};
			std::vector<double> normalPaired;
			std::vector<double> hugePaired;
			normalPaired.reserve(trials);
			hugePaired.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					hugePaired.push_back(runWorkspace(promotedWorkspace));
					normalPaired.push_back(runWorkspace(normalWorkspace));
				}
				else
				{
					normalPaired.push_back(runWorkspace(normalWorkspace));
					hugePaired.push_back(runWorkspace(promotedWorkspace));
				}
			}
#endif
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-FieldCheckpoint-v1"
				<< " mode=packed-tuning K=" << K << " N=" << N
				<< " trials=" << trials << " correctness=PASS\n";
			for (const auto& sample : samples)
				printSamples(
					("packed_prefetch" + std::to_string(sample.distance) + "_full").c_str(),
					sample.full);
			printSamples("pair_packed_full", pairPacked);
#ifndef _WIN32
			std::cout << "huge_value=" << hugeValue << " huge_tile=" << hugeTile << '\n';
			printSamples("huge_page_full", hugePageFull);
			printSamples("normal_workspace_paired_full", normalPaired);
			printSamples("huge_workspace_paired_full", hugePaired);
#endif
			std::cout << "sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (bucketSelectedOnly)
		{
			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size(), true);
			RiffleExactPermFieldCheckpoint::Workspace promotedWorkspace;
			auto runBaseline = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					schedule.optimizedGatherAndOuterTranspose(
						working.data(), message.data(), tiledWorkspace.data());
				});
			};
			auto runBucket = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					promoted.dualEncodePackedSchedulesUnchecked<32>(
						working.data(), message.data(), promotedWorkspace);
				});
			};
			auto runSlot = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					promoted.dualEncodePackedOffsetRouteUnchecked<32>(
						working.data(), message.data(), promotedWorkspace);
				});
			};

			runBaseline();
			const auto referenceMessage = message;
			runBucket();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("selected bucket route disagrees with frozen evaluator");
			runSlot();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("precomputed-slot route disagrees with frozen evaluator");

			std::vector<double> baselineFull;
			std::vector<double> slotBaselineFull;
			std::vector<double> bucketFull;
			std::vector<double> slotFull;
			baselineFull.reserve(trials);
			slotBaselineFull.reserve(trials);
			bucketFull.reserve(trials);
			slotFull.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					slotBaselineFull.push_back(runSlot());
					baselineFull.push_back(runBaseline());
				}
				else
				{
					baselineFull.push_back(runBaseline());
					slotBaselineFull.push_back(runSlot());
				}
			}
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					slotFull.push_back(runSlot());
					bucketFull.push_back(runBucket());
				}
				else
				{
					bucketFull.push_back(runBucket());
					slotFull.push_back(runSlot());
				}
			}
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-FieldCheckpoint-v1"
				<< " mode=bucket-selected"
				<< " K=" << K << " N=" << N
				<< " tile_outer_blocks="
				<< RiffleExactPermFieldCheckpoint::tileOuterBlocks
				<< " tile_mib="
				<< double(RiffleExactPermFieldCheckpoint::tileBlocks * sizeof(block)) /
					double(u64{ 1 } << 20)
				<< " bucket_value_mib="
				<< double(N * sizeof(block)) / double(u64{ 1 } << 20)
				<< " bucket_offset_mib="
				<< double(N * sizeof(u32)) / double(u64{ 1 } << 20)
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("baseline_full", baselineFull);
			printSamples("packed_offset24_vs_baseline_full", slotBaselineFull);
			printSamples("packed_both24_full", bucketFull);
			printSamples("packed_offset24_full", slotFull);
			std::cout << "bucket_over_baseline="
				<< median(slotBaselineFull) / median(baselineFull)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (bucketOnly)
		{
			schedule.gather(source.data(), outerWord.data());
			outerTranspose2(outerWord.data(), message.data());
			const auto rawReferenceMessage = message;
			std::memcpy(working.data(), source.data(), N * sizeof(block));
			checkpointOptimized(working.data(), N, coefficients.data());
			schedule.gather(working.data(), outerWord.data());
			outerTranspose2(outerWord.data(), message.data());
			const auto fullReferenceMessage = message;
			RiffleExactPermFieldCheckpoint promoted;
			promoted.init(
				0x5354525543543235ULL, coefficients.data(), coefficients.size(), true);
			RiffleExactPermFieldCheckpoint::ExperimentalWorkspace promotedWorkspace;
			promoted.dualEncodeDirectScatterUnchecked(
				source.data(), message.data(), promotedWorkspace);
			if (!equalBlocks(fullReferenceMessage.data(), message.data(), K))
				throw std::runtime_error("direct scatter disagrees with staged oracle");
			const auto directFull = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					promoted.dualEncodeDirectScatterUnchecked(
						working.data(), message.data(), promotedWorkspace);
				});
			auto benchmarkWritePrefetch = [&]<u64 Distance>(
				std::integral_constant<u64, Distance>) {
				promoted.dualEncodePrefetchUnchecked<Distance>(
					source.data(), message.data(), promotedWorkspace);
				if (!equalBlocks(fullReferenceMessage.data(), message.data(), K))
					throw std::runtime_error("write-prefetched bucket route disagrees with staged oracle");
				return benchmark(
					trials,
					[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
					[&] {
						promoted.dualEncodePrefetchUnchecked<Distance>(
							working.data(), message.data(), promotedWorkspace);
					});
			};
			const auto writePrefetch8 = benchmarkWritePrefetch(
				std::integral_constant<u64, 8>{});
			const auto writePrefetch16 = benchmarkWritePrefetch(
				std::integral_constant<u64, 16>{});
			const auto writePrefetch32 = benchmarkWritePrefetch(
				std::integral_constant<u64, 32>{});
			const auto writePrefetch64 = benchmarkWritePrefetch(
				std::integral_constant<u64, 64>{});

			struct BucketSamples
			{
				u64 tileOuterBlocks;
				u64 offsetBits;
				std::vector<double> permutationOuter;
				std::vector<double> full;
				std::vector<double> fusedFull;
			};
			auto runBucket = [&]<u64 TileOuterBlocks, typename Offset>(
				std::integral_constant<u64, TileOuterBlocks>, Offset) {
				std::vector<block> bucketValues(N);
				std::vector<Offset> bucketOffsets(N);
				std::vector<block> tile(TileOuterBlocks * OuterLength);
				schedule.bucketRouteAndOuterTranspose<TileOuterBlocks>(
					source.data(), message.data(), bucketValues.data(),
					bucketOffsets.data(), tile.data());
				if (!equalBlocks(rawReferenceMessage.data(), message.data(), K))
					throw std::runtime_error("bucket route disagrees with staged oracle");
				schedule.checkpointBucketRouteAndOuterTranspose<TileOuterBlocks>(
					source.data(), message.data(), coefficients.data(), bucketValues.data(),
					bucketOffsets.data(), tile.data());
				if (!equalBlocks(fullReferenceMessage.data(), message.data(), K))
					throw std::runtime_error("fused checkpoint bucket route disagrees with staged oracle");

				BucketSamples samples{
					TileOuterBlocks, sizeof(Offset) * 8,
					benchmark(trials, [] {}, [&] {
						schedule.bucketRouteAndOuterTranspose<TileOuterBlocks>(
							source.data(), message.data(), bucketValues.data(),
							bucketOffsets.data(), tile.data());
					}),
					benchmark(
						trials,
						[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
						[&] {
							checkpointOptimized(working.data(), N, coefficients.data());
							schedule.bucketRouteAndOuterTranspose<TileOuterBlocks>(
								working.data(), message.data(), bucketValues.data(),
								bucketOffsets.data(), tile.data());
						})
					,
					benchmark(trials, [] {}, [&] {
						schedule.checkpointBucketRouteAndOuterTranspose<TileOuterBlocks>(
							source.data(), message.data(), coefficients.data(),
							bucketValues.data(), bucketOffsets.data(), tile.data());
					})
				};
				benchmarkSink ^= checksum(message.data(), message.size());
				return samples;
			};

			std::vector<BucketSamples> buckets;
			buckets.emplace_back(runBucket(
				std::integral_constant<u64, 256>{}, u16{}));
			buckets.emplace_back(runBucket(
				std::integral_constant<u64, 512>{}, u32{}));
			buckets.emplace_back(runBucket(
				std::integral_constant<u64, 1024>{}, u32{}));
			buckets.emplace_back(runBucket(
				std::integral_constant<u64, 2048>{}, u32{}));
			buckets.emplace_back(runBucket(
				std::integral_constant<u64, 4096>{}, u32{}));

			const auto selectedPermutationOuter = benchmark(
				trials, [] {}, [&] {
					schedule.optimizedGatherAndOuterTranspose(
						source.data(), message.data(), tiledWorkspace.data());
				});
			const auto selectedFull = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					schedule.optimizedGatherAndOuterTranspose(
						working.data(), message.data(), tiledWorkspace.data());
				});
			benchmarkSink ^= checksum(message.data(), message.size());

			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-ExactPerm-FieldCheckpoint-v1"
				<< " mode=bucket-route"
				<< " K=" << K << " N=" << N
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("selected_permutation_outer", selectedPermutationOuter);
			printSamples("selected_full", selectedFull);
			printSamples("direct_scatter_full", directFull);
			printSamples("write_prefetch8_full", writePrefetch8);
			printSamples("write_prefetch16_full", writePrefetch16);
			printSamples("write_prefetch32_full", writePrefetch32);
			printSamples("write_prefetch64_full", writePrefetch64);
			for (const auto& samples : buckets)
			{
				const std::string label = "bucket_b" +
					std::to_string(samples.tileOuterBlocks) + "_u" +
					std::to_string(samples.offsetBits);
				printSamples((label + "_permutation_outer").c_str(), samples.permutationOuter);
				printSamples((label + "_full").c_str(), samples.full);
				printSamples((label + "_fused_full").c_str(), samples.fusedFull);
			}
			std::cout << "sink=" << benchmarkSink << '\n';
			return 0;
		}

		if (selectedOnly)
		{
			auto runBaseline = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					schedule.gather(working.data(), outerWord.data());
					outerTranspose(outerWord.data(), message.data());
				});
			};
			auto runSelected = [&] {
				std::memcpy(working.data(), source.data(), N * sizeof(block));
				return milliseconds([&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					schedule.optimizedGatherAndOuterTranspose(
						working.data(), message.data(), tiledWorkspace.data());
				});
			};
			runBaseline();
			const auto referenceMessage = message;
			runSelected();
			if (!equalBlocks(referenceMessage.data(), message.data(), K))
				throw std::runtime_error("selected full path disagrees with staged baseline");
			std::vector<double> baselineFull;
			std::vector<double> selectedFull;
			baselineFull.reserve(trials);
			selectedFull.reserve(trials);
			for (u64 trial = 0; trial < trials; ++trial)
			{
				if (trial & 1)
				{
					selectedFull.push_back(runSelected());
					baselineFull.push_back(runBaseline());
				}
				else
				{
					baselineFull.push_back(runBaseline());
					selectedFull.push_back(runSelected());
				}
			}

			const auto selectedInner = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] { checkpointOptimized(working.data(), N, coefficients.data()); });
			const auto selectedPermutationOuter = benchmark(
				trials, [] {}, [&] {
					schedule.optimizedGatherAndOuterTranspose(
						working.data(), message.data(), tiledWorkspace.data());
				});
			benchmarkSink ^= checksum(message.data(), message.size());
			std::cout << std::fixed << std::setprecision(6)
				<< "construction=Riffle-TransposeBitShuffle-FieldCheckpoint"
				<< " mode=selected"
				<< " K=" << K << " N=" << N
				<< " outer=ExtendedBch256x128-Eq3-AVX2x2"
				<< " tile_outer_blocks=" << StructuredSchedule::optimizedTileOuterBlocks
				<< " prefetch_distance=" << StructuredSchedule::optimizedPrefetchDistance
				<< " schedule_mib=8.000000"
				<< " trials=" << trials << " correctness=PASS\n";
			printSamples("selected_inner", selectedInner);
			printSamples("selected_permutation_outer", selectedPermutationOuter);
			printSamples("baseline_full", baselineFull);
			printSamples("selected_full", selectedFull);
			std::cout << "selected_over_baseline="
				<< median(selectedFull) / median(baselineFull)
				<< " sink=" << benchmarkSink << '\n';
			return 0;
		}

		const auto inner = benchmark(
			trials,
			[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
			[&] { checkpointOptimized(working.data(), N, coefficients.data()); });
		benchmarkSink ^= checksum(working.data(), working.size());

		const auto permutation = benchmark(
			trials, [] {}, [&] { schedule.gather(working.data(), outerWord.data()); });
		benchmarkSink ^= checksum(outerWord.data(), outerWord.size());
		const auto unrolledPermutation = benchmark(
			trials, [] {}, [&] { schedule.gatherUnrolled(working.data(), outerWord.data()); });
		const auto prefetch32Permutation = benchmark(
			trials, [] {}, [&] { schedule.gatherPrefetch<32>(working.data(), outerWord.data()); });
		const auto prefetch64Permutation = benchmark(
			trials, [] {}, [&] { schedule.gatherPrefetch<64>(working.data(), outerWord.data()); });
		const auto prefetch128Permutation = benchmark(
			trials, [] {}, [&] { schedule.gatherPrefetch<128>(working.data(), outerWord.data()); });
		benchmarkSink ^= checksum(outerWord.data(), outerWord.size());

		const auto outer = benchmark(
			trials, [] {}, [&] { outerTranspose(outerWord.data(), message.data()); });
		benchmarkSink ^= checksum(message.data(), message.size());
		const auto packedOuter = benchmark(
			trials, [] {}, [&] { outerTranspose2(outerWord.data(), message.data()); });
		benchmarkSink ^= checksum(message.data(), message.size());

		const auto fusedPermutationOuter = benchmark(
			trials, [] {},
			[&] { schedule.gatherAndOuterTranspose(working.data(), message.data()); });
		benchmarkSink ^= checksum(message.data(), message.size());
		const auto packedFusedPermutationOuter = benchmark(
			trials, [] {},
			[&] { schedule.gatherAndOuterTranspose2(working.data(), message.data()); });
		benchmarkSink ^= checksum(message.data(), message.size());

		const auto stagedFull = benchmark(
			trials,
			[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
			[&] {
				checkpointOptimized(working.data(), N, coefficients.data());
				schedule.gather(working.data(), outerWord.data());
				outerTranspose(outerWord.data(), message.data());
			});
		benchmarkSink ^= checksum(message.data(), message.size());
		const auto packedStagedFull = benchmark(
			trials,
			[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
			[&] {
				checkpointOptimized(working.data(), N, coefficients.data());
				schedule.gather(working.data(), outerWord.data());
				outerTranspose2(outerWord.data(), message.data());
			});
		benchmarkSink ^= checksum(message.data(), message.size());
		auto benchmarkPackedStaged = [&](auto&& gatherFunction) {
			return benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					gatherFunction();
					outerTranspose2(outerWord.data(), message.data());
				});
		};
		const auto packedUnrolledStagedFull = benchmarkPackedStaged(
			[&] { schedule.gatherUnrolled(working.data(), outerWord.data()); });
		const auto packedPrefetch32StagedFull = benchmarkPackedStaged(
			[&] { schedule.gatherPrefetch<32>(working.data(), outerWord.data()); });
		const auto packedPrefetch64StagedFull = benchmarkPackedStaged(
			[&] { schedule.gatherPrefetch<64>(working.data(), outerWord.data()); });
		const auto packedPrefetch128StagedFull = benchmarkPackedStaged(
			[&] { schedule.gatherPrefetch<128>(working.data(), outerWord.data()); });
		benchmarkSink ^= checksum(message.data(), message.size());

		const auto fusedFull = benchmark(
			trials,
			[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
			[&] {
				checkpointOptimized(working.data(), N, coefficients.data());
				schedule.gatherAndOuterTranspose(working.data(), message.data());
			});
		benchmarkSink ^= checksum(message.data(), message.size());

		struct TiledSamples
		{
			u64 batchSize;
			bool prefetch;
			bool packed;
			std::vector<double> permutationOuter;
			std::vector<double> full;
		};
		std::vector<TiledSamples> tiled;
		struct TiledConfig
		{
			u64 batchSize;
			bool packed;
			bool prefetch;
		};
		for (const auto [batchSize, packed, prefetch] : std::array{
			TiledConfig{ 64, false, false },
			TiledConfig{ 256, false, true },
			TiledConfig{ 128, true, false },
			TiledConfig{ 256, true, true } })
		{
				TiledSamples samples{ batchSize, prefetch, packed };
				samples.permutationOuter = benchmark(
					trials, [] {}, [&] {
						if (packed && prefetch)
							schedule.tiledPrefetchPackedGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
						else if (packed)
							schedule.tiledPackedGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
						else if (prefetch)
							schedule.tiledPrefetchGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
						else
							schedule.tiledGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
					});
				samples.full = benchmark(
					trials,
					[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
					[&] {
						checkpointOptimized(working.data(), N, coefficients.data());
						if (packed && prefetch)
							schedule.tiledPrefetchPackedGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
						else if (packed)
							schedule.tiledPackedGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
						else if (prefetch)
							schedule.tiledPrefetchGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
						else
							schedule.tiledGatherAndOuterTranspose(
								working.data(), message.data(), tiledWorkspace.data(), batchSize);
					});
				benchmarkSink ^= checksum(message.data(), message.size());
				tiled.emplace_back(std::move(samples));
		}

		struct PrefetchSamples
		{
			u64 batchSize;
			u64 distance;
			std::vector<double> permutationOuter;
			std::vector<double> full;
		};
		auto benchmarkPrefetchDistance = [&](auto distanceTag, u64 batchSize) {
			constexpr u64 distance = decltype(distanceTag)::value;
			PrefetchSamples samples{ batchSize, distance };
			samples.permutationOuter = benchmark(
				trials, [] {}, [&] {
					schedule.tiledPackedPrefetchDistance<distance>(
						working.data(), message.data(), tiledWorkspace.data(), batchSize);
				});
			samples.full = benchmark(
				trials,
				[&] { std::memcpy(working.data(), source.data(), N * sizeof(block)); },
				[&] {
					checkpointOptimized(working.data(), N, coefficients.data());
					schedule.tiledPackedPrefetchDistance<distance>(
						working.data(), message.data(), tiledWorkspace.data(), batchSize);
				});
			benchmarkSink ^= checksum(message.data(), message.size());
			return samples;
		};
		std::vector<PrefetchSamples> tunedPrefetch;
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 4>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 8>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 16>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 24>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 32>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 40>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 48>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 64>{}, 256));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 8>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 12>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 16>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 24>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 32>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 40>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 48>{}, 512));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 24>{}, 1024));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 32>{}, 1024));
		tunedPrefetch.emplace_back(benchmarkPrefetchDistance(std::integral_constant<u64, 40>{}, 1024));

		std::cout << std::fixed << std::setprecision(6)
			<< "construction=Riffle-TransposeBitShuffle-FieldCheckpoint"
			<< " K=" << K << " N=" << N
			<< " outer=ExtendedBch256x128-Eq3"
			<< " outer_xors=" << ExtendedBch256x128Eq3::transposedEncoderXorCount
			<< " state=" << StateBlocks << " epoch=" << EpochBlocks
			<< " trials=" << trials << '\n'
			<< "schedule_setup_ms=" << setupMs
			<< " flat_schedule_mib=" << double(N * sizeof(u32)) / double(u64{ 1 } << 20)
			<< " tiled_schedule_mib=" << double(N * sizeof(u32)) / double(u64{ 1 } << 20)
			<< " correctness=PASS\n";
		printSamples("inner", inner);
		printSamples("structured_permutation", permutation);
		printSamples("unrolled_permutation", unrolledPermutation);
		printSamples("prefetch32_permutation", prefetch32Permutation);
		printSamples("prefetch64_permutation", prefetch64Permutation);
		printSamples("prefetch128_permutation", prefetch128Permutation);
		printSamples("outer", outer);
		printSamples("packed_outer", packedOuter);
		printSamples("fused_permutation_outer", fusedPermutationOuter);
		printSamples("packed_fused_permutation_outer", packedFusedPermutationOuter);
		printSamples("staged_full", stagedFull);
		printSamples("packed_staged_full", packedStagedFull);
		printSamples("packed_unrolled_staged_full", packedUnrolledStagedFull);
		printSamples("packed_prefetch32_staged_full", packedPrefetch32StagedFull);
		printSamples("packed_prefetch64_staged_full", packedPrefetch64StagedFull);
		printSamples("packed_prefetch128_staged_full", packedPrefetch128StagedFull);
		printSamples("fused_full", fusedFull);
		for (const auto& samples : tiled)
		{
			const std::string label = "tiled_b" + std::to_string(samples.batchSize) +
				(samples.packed ? "_packed" : "") +
				(samples.prefetch ? "_prefetch" : "");
			printSamples((label + "_permutation_outer").c_str(), samples.permutationOuter);
			printSamples((label + "_full").c_str(), samples.full);
		}
		for (const auto& samples : tunedPrefetch)
		{
			const std::string label = "tuned_b" + std::to_string(samples.batchSize) +
				"_packed_prefetch" + std::to_string(samples.distance);
			printSamples((label + "_permutation_outer").c_str(), samples.permutationOuter);
			printSamples((label + "_full").c_str(), samples.full);
		}
		std::cout << "component_sum_ms="
			<< median(inner) + median(permutation) + median(outer)
			<< " fused_sum_ms=" << median(inner) + median(fusedPermutationOuter)
			<< " sink=" << benchmarkSink << '\n';
		return 0;
	}
	catch (const std::exception& error)
	{
		std::cerr << "error=" << error.what() << '\n';
		return 1;
	}
}
