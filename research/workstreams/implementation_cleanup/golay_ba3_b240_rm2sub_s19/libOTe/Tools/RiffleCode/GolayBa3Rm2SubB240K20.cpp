#include "GolayBa3Rm2SubB240K20.h"

#include "GolayBa3Rm2Sub_TestAccess.hpp"

#include <libOTe/Tools/RiffleCode/RiffleRm2SubS19.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <memory>
#include <new>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace osuCrypto
{
	namespace
	{
		using Code = GolayBa3Rm2SubB240K20;
		constexpr u64 RouteBuckets = 8;
		constexpr u64 TileRows = Code::outerRows / RouteBuckets;
		constexpr u64 TileBlocks = TileRows * Code::outerLength;
		static_assert(Code::outerRows % RouteBuckets == 0);
		static_assert(TileRows == 1'104);
		static_assert(TileBlocks < (u64{ 1 } << 19));
		static_assert(Code::codeBlocks < (u64{ 1 } << 22));

		constexpr std::array<u32, 12> GolayGeneratorRows{
			0x800ae3U, 0x8015c6U, 0x802b8cU, 0x805718U,
			0x80ae30U, 0x815c60U, 0x82b8c0U, 0x857180U,
			0x8ae300U, 0x95c600U, 0xab8c00U, 0xd71800U,
		};

		struct Implementation
		{
			std::vector<u8> innerToSlot24;
			std::vector<u8> bucketOffsets24;
			std::vector<u8> baPermutations;
			RiffleRm2SubS19Transpose inner;
			GolayBa3OuterMode mode = GolayBa3OuterMode::Reused;
			bool initialized = false;
		};

		struct EncodeWorkspace
		{
			EncodeWorkspace()
				: bucketValues(Code::codeBlocks)
				, tile(TileBlocks)
			{
			}

			std::vector<block> bucketValues;
			std::vector<block> tile;
		};

		static_assert(sizeof(Implementation) <= 256);
		static_assert(alignof(Implementation) <= 64);
		static_assert(sizeof(EncodeWorkspace) <= 64);
		static_assert(alignof(EncodeWorkspace) <= 64);

		Implementation& implementation(std::byte* storage) noexcept
		{
			return *std::launder(reinterpret_cast<Implementation*>(storage));
		}

		const Implementation& implementation(const std::byte* storage) noexcept
		{
			return *std::launder(reinterpret_cast<const Implementation*>(storage));
		}

		EncodeWorkspace& encodeWorkspace(std::byte* storage) noexcept
		{
			return *std::launder(reinterpret_cast<EncodeWorkspace*>(storage));
		}

		u64 splitmix64(u64& state) noexcept
		{
			u64 value = (state += 0x9e3779b97f4a7c15ULL);
			value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
			value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
			return value ^ (value >> 31);
		}

		u64 uniformBelow(u64& state, u64 bound) noexcept
		{
			const u64 threshold = (u64{ 0 } - bound) % bound;
			for (;;)
			{
				const u64 value = splitmix64(state);
				if (value >= threshold)
					return value % bound;
			}
		}

		template<typename Integer, std::size_t Size>
		void samplePermutation(std::array<Integer, Size>& permutation, u64& state)
		{
			std::iota(permutation.begin(), permutation.end(), Integer{ 0 });
			for (std::size_t last = Size; last > 1; --last)
				std::swap(permutation[last - 1], permutation[uniformBelow(state, last)]);
		}

		OC_FORCEINLINE block xor8(
			const block& x0, const block& x1, const block& x2, const block& x3,
			const block& x4, const block& x5, const block& x6, const block& x7) noexcept
		{
			const auto a = _mm_xor_si128(x0.mData, x1.mData);
			const auto b = _mm_xor_si128(x2.mData, x3.mData);
			const auto c = _mm_xor_si128(x4.mData, x5.mData);
			const auto d = _mm_xor_si128(x6.mData, x7.mData);
			return block(_mm_xor_si128(_mm_xor_si128(a, b), _mm_xor_si128(c, d)));
		}

		OC_FORCEINLINE void golayTranspose24x12(
			const block* __restrict input,
			block* __restrict output,
			u64 outputCount = 12) noexcept
		{
#define GOLAY_ROW(I, A, B, C, D, E, F, G) \
			if (outputCount > I) output[I] = xor8(input[23], input[(I) + A], \
				input[(I) + B], input[(I) + C], input[(I) + D], \
				input[(I) + E], input[(I) + F], input[(I) + G])
			GOLAY_ROW(0, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(1, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(2, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(3, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(4, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(5, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(6, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(7, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(8, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(9, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(10, 0, 1, 5, 6, 7, 9, 11);
			GOLAY_ROW(11, 0, 1, 5, 6, 7, 9, 11);
#undef GOLAY_ROW
		}

		OC_FORCEINLINE void suffixAccumulate(block* values) noexcept
		{
			for (u64 index = Code::outerLength - 1; index-- > 0;)
				values[index] ^= values[index + 1];
		}

		OC_FORCEINLINE void transposePermutation(
			const block* __restrict input,
			block* __restrict output,
			const u8* __restrict sourceToDestination) noexcept
		{
			for (u64 source = 0; source < Code::outerLength; ++source)
				output[source] = input[sourceToDestination[source]];
		}

		OC_FORCEINLINE void outerTransposeRow(
			const Implementation& state,
			u64 row,
			const block* __restrict baWord,
			block* __restrict output,
			u64 outputCount) noexcept
		{
			alignas(32) block stage0[Code::outerLength];
			alignas(32) block stage1[Code::outerLength];
			for (u64 coordinate = 0; coordinate < Code::outerLength; ++coordinate)
				stage0[coordinate] = baWord[coordinate];

			suffixAccumulate(stage0);
			const u64 permutationBase = state.mode == GolayBa3OuterMode::Reused
				? 0 : 2 * row * Code::outerLength;
			const u8* __restrict permutation1 = state.baPermutations.data() + permutationBase;
			const u8* __restrict permutation2 = permutation1 + Code::outerLength;
			transposePermutation(stage0, stage1, permutation2);
			suffixAccumulate(stage1);
			transposePermutation(stage1, stage0, permutation1);

			for (u64 constituent = 0; outputCount != 0; ++constituent)
			{
				const u64 count = std::min<u64>(12, outputCount);
				golayTranspose24x12(
					stage0 + 24 * constituent,
					output + 12 * constituent,
					count);
				outputCount -= count;
			}
		}

		void scalarOuterTransposeRow(
			const Implementation& state,
			u64 row,
			const block* innerWord,
			const u32* baToInner,
			block* output,
			u64 outputCount) noexcept
		{
			std::array<block, Code::outerLength> first;
			std::array<block, Code::outerLength> second;
			const u32* route = baToInner + row * Code::outerLength;
			for (u64 coordinate = 0; coordinate < Code::outerLength; ++coordinate)
				first[coordinate] = innerWord[route[coordinate]];
			for (u64 coordinate = Code::outerLength - 1; coordinate-- > 0;)
				first[coordinate] ^= first[coordinate + 1];

			const u64 permutationBase = state.mode == GolayBa3OuterMode::Reused
				? 0 : 2 * row * Code::outerLength;
			const u8* permutation1 = state.baPermutations.data() + permutationBase;
			const u8* permutation2 = permutation1 + Code::outerLength;
			for (u64 source = 0; source < Code::outerLength; ++source)
				second[source] = first[permutation2[source]];
			for (u64 coordinate = Code::outerLength - 1; coordinate-- > 0;)
				second[coordinate] ^= second[coordinate + 1];
			for (u64 source = 0; source < Code::outerLength; ++source)
				first[source] = second[permutation1[source]];

			for (u64 message = 0; message < outputCount; ++message)
			{
				const u64 constituent = message / 12;
				const u64 localMessage = message % 12;
				u32 mask = GolayGeneratorRows[localMessage];
				block value(0, 0);
				while (mask)
				{
					const unsigned coordinate = std::countr_zero(mask);
					value ^= first[24 * constituent + coordinate];
					mask &= mask - 1;
				}
				output[message] = value;
			}
		}

		void storePacked24(u8* destination, u32 value) noexcept
		{
			destination[0] = static_cast<u8>(value);
			destination[1] = static_cast<u8>(value >> 8);
			destination[2] = static_cast<u8>(value >> 16);
		}

		OC_FORCEINLINE u32 loadPacked24(const u8* source) noexcept
		{
			u32 value;
			std::memcpy(&value, source, sizeof(value));
			return value & 0x00ffffffU;
		}

		void initializeRoutes(Implementation& state, u64 routeSeed)
		{
			std::vector<u32> baToInner(Code::codeBlocks);
			std::array<u8, Code::outerLength> coordinateToRegion;
			std::vector<u16> rowToPosition(Code::outerRows);
			std::vector<u16> regionRoutes(Code::outerLength * Code::outerRows);

			for (u64 region = 0; region < Code::outerLength; ++region)
			{
				std::iota(rowToPosition.begin(), rowToPosition.end(), u16{ 0 });
				for (u64 last = Code::outerRows; last > 1; --last)
					std::swap(rowToPosition[last - 1], rowToPosition[uniformBelow(routeSeed, last)]);
				std::memcpy(
					regionRoutes.data() + region * Code::outerRows,
					rowToPosition.data(),
					Code::outerRows * sizeof(u16));
			}

			for (u64 row = 0; row < Code::outerRows; ++row)
			{
				samplePermutation(coordinateToRegion, routeSeed);
				for (u64 coordinate = 0; coordinate < Code::outerLength; ++coordinate)
				{
					const u64 region = coordinateToRegion[coordinate];
					const u64 position = regionRoutes[region * Code::outerRows + row];
					baToInner[row * Code::outerLength + coordinate] =
						static_cast<u32>(region * Code::outerRows + position);
				}
			}

			std::vector<u32> innerToBa(Code::codeBlocks);
			for (u64 ba = 0; ba < Code::codeBlocks; ++ba)
				innerToBa[baToInner[ba]] = static_cast<u32>(ba);
			state.innerToSlot24.resize(3 * Code::codeBlocks + sizeof(u32));
			state.bucketOffsets24.resize(3 * Code::codeBlocks + sizeof(u32));
			std::array<u32, RouteBuckets> counts;
			counts.fill(static_cast<u32>(TileBlocks));
			for (u64 inner = Code::codeBlocks; inner-- > 0;)
			{
				const u32 ba = innerToBa[inner];
				const u32 bucket = ba / static_cast<u32>(TileBlocks);
				const u32 local = ba - bucket * static_cast<u32>(TileBlocks);
				const u32 slot = bucket * static_cast<u32>(TileBlocks) + --counts[bucket];
				storePacked24(state.innerToSlot24.data() + 3 * inner, slot);
				storePacked24(state.bucketOffsets24.data() + 3 * slot, local);
			}
		}

		std::vector<u32> reconstructBaToInner(const Implementation& state)
		{
			std::vector<u32> baToInner(Code::codeBlocks);
			for (u64 inner = 0; inner < Code::codeBlocks; ++inner)
			{
				const u32 slot = loadPacked24(state.innerToSlot24.data() + 3 * inner);
				const u32 bucket = slot / static_cast<u32>(TileBlocks);
				const u32 local = loadPacked24(state.bucketOffsets24.data() + 3 * slot);
				baToInner[bucket * TileBlocks + local] = static_cast<u32>(inner);
			}
			return baToInner;
		}

		void initializeBaPermutations(
			Implementation& state,
			u64 baSeed,
			GolayBa3OuterMode mode)
		{
			const u64 draws = mode == GolayBa3OuterMode::Reused ? 1 : Code::outerRows;
			state.baPermutations.resize(2 * draws * Code::outerLength);
			std::array<u8, Code::outerLength> permutation;
			for (u64 draw = 0; draw < draws; ++draw)
				for (u64 layer = 0; layer < 2; ++layer)
				{
					samplePermutation(permutation, baSeed);
					std::memcpy(
						state.baPermutations.data() +
							(2 * draw + layer) * Code::outerLength,
						permutation.data(),
						Code::outerLength);
				}
		}

		void scalarInnerTranspose(
			const Implementation& state,
			const block* input,
			block* output) noexcept
		{
			std::array<block, 19> recurrence;
			recurrence.fill(block(0, 0));
			for (u64 epoch = Code::innerEpochs; epoch-- > 0;)
			{
				const block* node = input + 128 * epoch;
				for (u64 point = 0; point < 128; ++point)
				{
					block addend(0, 0);
					u32 column = detail::Rm2Sub19Columns[point];
					while (column)
					{
						const unsigned bit = std::countr_zero(column);
						addend ^= recurrence[bit];
						column &= column - 1;
					}
					output[128 * epoch + point] = node[point] ^ addend;
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
				if (epoch + 1 == Code::innerEpochs)
					recurrence = syndrome;
				else
				{
					std::array<block, 19> multiplied;
					const auto& masks = state.inner.schedules()[epoch].rowMasks;
					for (unsigned row = 0; row < 19; ++row)
					{
						block value(0, 0);
						u32 mask = masks[row];
						while (mask)
						{
							const unsigned bit = std::countr_zero(mask);
							value ^= recurrence[bit];
							mask &= mask - 1;
						}
						multiplied[row] = value ^ syndrome[row];
					}
					recurrence = multiplied;
				}
			}
		}
	}

	GolayBa3Rm2SubB240K20::Workspace::Workspace()
	{
		::new (static_cast<void*>(mStorage)) EncodeWorkspace;
	}

	GolayBa3Rm2SubB240K20::Workspace::~Workspace()
	{
		std::destroy_at(&encodeWorkspace(mStorage));
	}

	GolayBa3Rm2SubB240K20::GolayBa3Rm2SubB240K20()
	{
		::new (static_cast<void*>(mStorage)) Implementation;
	}

	GolayBa3Rm2SubB240K20::~GolayBa3Rm2SubB240K20()
	{
		std::destroy_at(&implementation(mStorage));
	}

	void GolayBa3Rm2SubB240K20::init(
		u64 routeSeed,
		u64 baSeed,
		u64 innerCoefficientSeed,
		GolayBa3OuterMode outerMode)
	{
		auto& state = implementation(mStorage);
		state.initialized = false;
		state.mode = outerMode;
		initializeRoutes(state, routeSeed);
		initializeBaPermutations(state, baSeed, outerMode);
		state.inner.init(innerEpochs, innerCoefficientSeed);
		state.initialized = true;
	}

	GolayBa3OuterMode GolayBa3Rm2SubB240K20::outerMode() const
	{
		const auto& state = implementation(mStorage);
		if (!state.initialized)
			throw std::logic_error("GolayBa3Rm2SubB240K20 is not initialized");
		return state.mode;
	}

	u64 GolayBa3Rm2SubB240K20::setupBytes() const noexcept
	{
		const auto& state = implementation(mStorage);
		return state.innerToSlot24.capacity() * sizeof(u8) +
			state.bucketOffsets24.capacity() * sizeof(u8) +
			state.baPermutations.capacity() * sizeof(u8) +
			state.inner.schedules().capacity() *
				sizeof(detail::Rm2Sub19FieldTransposeSchedule);
	}

	u64 GolayBa3Rm2SubB240K20::workspaceBytes() const noexcept
	{
		return (codeBlocks + TileBlocks) * sizeof(block);
	}

	void GolayBa3Rm2SubB240K20::dualEncodeTo(
		const block* __restrict input,
		u64 inputSize,
		block* __restrict output,
		u64 outputSize,
		Workspace& workspace) const
	{
		if (!implementation(mStorage).initialized)
			throw std::logic_error("GolayBa3Rm2SubB240K20 is not initialized");
		if (inputSize != codeBlocks || outputSize != messageBlocks)
			throw std::invalid_argument("Golay-BA-3/RM2Sub span size mismatch");
		dualEncodeUnchecked(input, output, workspace);
	}

	void GolayBa3Rm2SubB240K20::dualEncodeUnchecked(
		const block* __restrict input,
		block* __restrict output,
		Workspace& workspace) const noexcept
	{
		const auto& state = implementation(mStorage);
		auto& scratch = encodeWorkspace(workspace.mStorage);
		block* __restrict values = scratch.bucketValues.data();
		const u8* slotBytes = state.innerToSlot24.data() + 3 * codeBlocks;
		state.inner.emitReverse(input, codeBlocks,
			[&](u64, block value) {
				slotBytes -= 3;
				values[loadPacked24(slotBytes)] = value;
			});

		constexpr u64 fullRows = messageBlocks / outerDimension;
		constexpr u64 finalCount = messageBlocks % outerDimension;
		constexpr u64 activeRows = fullRows + (finalCount != 0);
		block* __restrict tile = scratch.tile.data();
		for (u64 bucket = 0; bucket < RouteBuckets; ++bucket)
		{
			const u64 bucketBase = bucket * TileBlocks;
			const u8* offsetBytes = state.bucketOffsets24.data() + 3 * bucketBase;
			for (u64 offset = 0; offset < TileBlocks; ++offset)
			{
				if (offset + 32 < TileBlocks)
				{
					const u32 future = loadPacked24(offsetBytes + 3 * 32);
					_mm_prefetch(reinterpret_cast<const char*>(tile + future), _MM_HINT_T0);
				}
				const u32 local = loadPacked24(offsetBytes);
				offsetBytes += 3;
				tile[local] = values[bucketBase + offset];
			}

			const u64 rowBase = bucket * TileRows;
			if (rowBase >= activeRows)
				break;
			const u64 rows = std::min<u64>(TileRows, activeRows - rowBase);
			for (u64 localRow = 0; localRow < rows; ++localRow)
			{
				const u64 row = rowBase + localRow;
				const u64 count = row < fullRows ? outerDimension : finalCount;
				outerTransposeRow(
					state, row, tile + localRow * outerLength,
					output + row * outerDimension, count);
			}
		}
	}

	void GolayBa3Rm2SubTestAccess::validateSetup(
		const GolayBa3Rm2SubB240K20& code)
	{
		const auto& state = implementation(code.mStorage);
		if (!state.initialized)
			throw std::logic_error("test setup is not initialized");

		const auto baToInner = reconstructBaToInner(state);
		std::vector<u8> seen(Code::codeBlocks, 0);
		for (u32 destination : baToInner)
		{
			if (destination >= Code::codeBlocks || seen[destination])
				throw std::runtime_error("factored route is not a permutation");
			seen[destination] = 1;
		}

		const u64 draws = state.mode == GolayBa3OuterMode::Reused ? 1 : Code::outerRows;
		std::array<u8, Code::outerLength> permutationSeen;
		for (u64 draw = 0; draw < draws; ++draw)
			for (u64 layer = 0; layer < 2; ++layer)
			{
				permutationSeen.fill(0);
				const u8* permutation = state.baPermutations.data() +
					(2 * draw + layer) * Code::outerLength;
				for (u64 coordinate = 0; coordinate < Code::outerLength; ++coordinate)
				{
					const u8 destination = permutation[coordinate];
					if (destination >= Code::outerLength || permutationSeen[destination])
						throw std::runtime_error("BA schedule is not a permutation");
					permutationSeen[destination] = 1;
				}
			}

		std::array<u64, 25> spectrum{};
		for (u32 message = 0; message < (1U << 12); ++message)
		{
			u32 word = 0;
			for (unsigned bit = 0; bit < 12; ++bit)
				if ((message >> bit) & 1U)
					word ^= GolayGeneratorRows[bit];
			++spectrum[std::popcount(word)];
		}
		if (spectrum[0] != 1 || spectrum[8] != 759 || spectrum[12] != 2576 ||
			spectrum[16] != 759 || spectrum[24] != 1 ||
			std::accumulate(spectrum.begin(), spectrum.end(), u64{ 0 }) != (1U << 12))
			throw std::runtime_error("Golay generator spectrum mismatch");
	}

	void GolayBa3Rm2SubTestAccess::dualEncodeReference(
		const GolayBa3Rm2SubB240K20& code,
		const block* input,
		block* output,
		GolayBa3Rm2SubB240K20::Workspace& workspace)
	{
		const auto& state = implementation(code.mStorage);
		auto& scratch = encodeWorkspace(workspace.mStorage).bucketValues;
		scalarInnerTranspose(state, input, scratch.data());
		const auto baToInner = reconstructBaToInner(state);
		constexpr u64 fullRows = Code::messageBlocks / Code::outerDimension;
		constexpr u64 finalCount = Code::messageBlocks % Code::outerDimension;
		for (u64 row = 0; row < fullRows; ++row)
			scalarOuterTransposeRow(
				state, row, scratch.data(), baToInner.data(),
				output + row * Code::outerDimension,
				Code::outerDimension);
		if constexpr (finalCount != 0)
			scalarOuterTransposeRow(
				state, fullRows, scratch.data(), baToInner.data(),
				output + fullRows * Code::outerDimension,
				finalCount);
	}
}
