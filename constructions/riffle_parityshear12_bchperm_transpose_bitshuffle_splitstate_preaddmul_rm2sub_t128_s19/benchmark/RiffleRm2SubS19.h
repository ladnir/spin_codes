#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

#include <immintrin.h>
#include <cryptoTools/Common/Defines.h>

namespace osuCrypto
{
	namespace detail
	{
		inline constexpr std::array<u32, 11> Rm2Sub19QuadraticMasks{
			0x11d97b, 0x18ef1a, 0x11ae86, 0x17ce2d, 0x063684, 0x001a59,
			0x14a88c, 0x13af70, 0x1b188e, 0x1eee71, 0x0ee30d,
		};

		inline constexpr std::array<unsigned, 21> Rm2Sub19MonomialPositions{
			0x03, 0x05, 0x09, 0x11, 0x21, 0x41, 0x06, 0x0a, 0x12, 0x22, 0x42, 0x0c, 0x14, 0x24, 0x44, 0x18, 0x28, 0x48, 0x30, 0x50, 0x60,
		};

		inline constexpr std::array<u32, 128> Rm2Sub19Columns{
			0x00001, 0x00003, 0x00005, 0x62907, 0x00009, 0x1070b, 0x2a10d, 0x58f0f,
			0x00011, 0x55c13, 0x15415, 0x22117, 0x3ef19, 0x7b41b, 0x01a1d, 0x2681f,
			0x00021, 0x56b23, 0x48325, 0x7c127, 0x13129, 0x55d2b, 0x7132d, 0x5562f,
			0x6cf31, 0x6f833, 0x31835, 0x50637, 0x41139, 0x5213b, 0x3673d, 0x47e3f,
			0x00041, 0x2a343, 0x6be45, 0x23447, 0x6d649, 0x5724b, 0x2c94d, 0x7444f,
			0x18d51, 0x67253, 0x66755, 0x7b157, 0x4b459, 0x24c5b, 0x1ff5d, 0x12e5f,
			0x65861, 0x19063, 0x46565, 0x58467, 0x1bf69, 0x7706b, 0x1236d, 0x1c56f,
			0x11a71, 0x38e73, 0x27375, 0x6ce77, 0x51279, 0x6817b, 0x4da7d, 0x1607f,
			0x00081, 0x28983, 0x29e85, 0x63e87, 0x60b89, 0x5858b, 0x6348d, 0x3938f,
			0x79891, 0x04d93, 0x45295, 0x5ae97, 0x27c99, 0x4ae9b, 0x3179d, 0x3ec9f,
			0x702a1, 0x0e0a3, 0x11fa5, 0x0d4a7, 0x038a9, 0x6ddab, 0x484ad, 0x448af,
			0x655b1, 0x4ebb3, 0x11cb5, 0x58bb7, 0x280b9, 0x139bb, 0x768bd, 0x2f8bf,
			0x3cfc1, 0x3e5c3, 0x7efc5, 0x1ecc7, 0x312c9, 0x23fcb, 0x593cd, 0x297cf,
			0x5dad1, 0x0acd3, 0x0aed5, 0x3f1d7, 0x6e8d9, 0x299db, 0x13ddd, 0x365df,
			0x295e1, 0x7d4e3, 0x236e5, 0x15ee7, 0x379e9, 0x73feb, 0x17bed, 0x314ef,
			0x24ff1, 0x252f3, 0x3b8f5, 0x58cf7, 0x04cf9, 0x156fb, 0x31afd, 0x429ff,
		};

		inline constexpr std::array<unsigned, 4> Rm2Sub19Group0{0, 16, 6, 2};
		inline constexpr std::array<unsigned, 4> Rm2Sub19Group1{14, 18, 1, 10};
		inline constexpr std::array<unsigned, 4> Rm2Sub19Group2{9, 12, 3, 5};
		inline constexpr std::array<unsigned, 4> Rm2Sub19Group3{17, 13, 15, 8};
		inline constexpr std::array<unsigned, 3> Rm2Sub19Group4{11, 7, 4};

		struct Rm2Sub19FieldTransposeSchedule
		{
			alignas(32) std::array<u32, 19> rowMasks{};
		};

		OC_FORCEINLINE __m128i rm2s19Xor(__m128i left, __m128i right) noexcept
		{
			return _mm_xor_si128(left, right);
		}

		inline u32 rm2s19FieldMultiplyScalar(u32 left, u32 right) noexcept
		{
			constexpr u32 modulusLow = 0x27; // x^19+x^5+x^2+x+1.
			u32 result = 0;
			for (unsigned bit = 0; bit < 19; ++bit)
			{
				result ^= left & (0U - (right & 1U));
				right >>= 1;
				const u32 carry = left >> 18;
				left = (left << 1) & 0x7ffffU;
				left ^= modulusLow & (0U - carry);
			}
			return result;
		}

		inline Rm2Sub19FieldTransposeSchedule rm2s19MakeFieldTransposeSchedule(
			u32 coefficient) noexcept
		{
			Rm2Sub19FieldTransposeSchedule schedule{};
			for (unsigned row = 0; row < 19; ++row)
				schedule.rowMasks[row] = rm2s19FieldMultiplyScalar(1U << row, coefficient);
			return schedule;
		}

		OC_FORCEINLINE void rm2s19BuildNibbleTable(
			const __m128i* input, __m128i* table) noexcept
		{
			table[0] = _mm_setzero_si128();
			table[1] = input[0];
			table[2] = input[1];
			table[3] = rm2s19Xor(input[0], input[1]);
			table[4] = input[2];
			table[5] = rm2s19Xor(input[2], table[1]);
			table[6] = rm2s19Xor(input[2], table[2]);
			table[7] = rm2s19Xor(input[2], table[3]);
			table[8] = input[3];
			table[9] = rm2s19Xor(input[3], table[1]);
			table[10] = rm2s19Xor(input[3], table[2]);
			table[11] = rm2s19Xor(input[3], table[3]);
			table[12] = rm2s19Xor(input[3], table[4]);
			table[13] = rm2s19Xor(input[3], table[5]);
			table[14] = rm2s19Xor(input[3], table[6]);
			table[15] = rm2s19Xor(input[3], table[7]);
		}

		OC_FORCEINLINE void rm2s19BuildThreeBitTable(
			const __m128i* input, __m128i* table) noexcept
		{
			table[0] = _mm_setzero_si128();
			table[1] = input[0];
			table[2] = input[1];
			table[3] = rm2s19Xor(input[0], input[1]);
			table[4] = input[2];
			table[5] = rm2s19Xor(input[2], table[1]);
			table[6] = rm2s19Xor(input[2], table[2]);
			table[7] = rm2s19Xor(input[2], table[3]);
		}

		OC_FORCEINLINE void rm2s19FieldMultiplyTranspose(
			__m128i* state, const Rm2Sub19FieldTransposeSchedule& schedule) noexcept
		{
			alignas(32) __m128i table[4][16];
			for (unsigned group = 0; group < 4; ++group)
				rm2s19BuildNibbleTable(state + 4 * group, table[group]);
			alignas(32) __m128i finalTable[8];
			rm2s19BuildThreeBitTable(state + 16, finalTable);
			alignas(32) __m128i output[19];
			for (unsigned row = 0; row < 19; ++row)
			{
				const u32 mask = schedule.rowMasks[row];
				auto value = rm2s19Xor(table[0][mask & 15], table[1][(mask >> 4) & 15]);
				value = rm2s19Xor(value, table[2][(mask >> 8) & 15]);
				value = rm2s19Xor(value, table[3][(mask >> 12) & 15]);
				output[row] = rm2s19Xor(value, finalTable[(mask >> 16) & 7]);
			}
			std::memcpy(state, output, sizeof(output));
		}

		template<unsigned Distance>
		OC_FORCEINLINE void rm2s19ZetaPair(__m128i* values, unsigned low) noexcept
		{
			auto lhs = _mm256_load_si256(reinterpret_cast<const __m256i*>(values + low));
			const auto rhs = _mm256_load_si256(
				reinterpret_cast<const __m256i*>(values + low + Distance));
			_mm256_store_si256(reinterpret_cast<__m256i*>(values + low), _mm256_xor_si256(lhs, rhs));
		}

#if defined(__AVX512F__)
		template<unsigned Distance>
		OC_FORCEINLINE void rm2s19ZetaQuad(__m128i* values, unsigned low) noexcept
		{
			auto lhs = _mm512_loadu_si512(values + low);
			const auto rhs = _mm512_loadu_si512(values + low + Distance);
			_mm512_storeu_si512(values + low, _mm512_xor_si512(lhs, rhs));
		}
#endif

		template<unsigned Distance>
		OC_FORCEINLINE void rm2s19ZetaStage(__m128i* values) noexcept
		{
			for (unsigned base = 0; base < 128; base += 2 * Distance)
			{
				if constexpr (Distance == 1)
					values[base] = rm2s19Xor(values[base], values[base + 1]);
				else
					for (unsigned offset = 0; offset < Distance; offset += 2)
						rm2s19ZetaPair<Distance>(values, base + offset);
			}
		}

		OC_FORCEINLINE void rm2s19ZetaTransposePruned(__m128i* values) noexcept
		{
#if defined(__AVX512F__)
			for (unsigned offset = 0; offset < 64; offset += 4) rm2s19ZetaQuad<64>(values, offset);
			for (unsigned base = 0; base < 128; base += 64)
				for (unsigned offset = 0; offset < 32; offset += 4) rm2s19ZetaQuad<32>(values, base + offset);
			for (unsigned base = 0; base < 128; base += 32)
				for (unsigned offset = 0; offset < 16; offset += 4) rm2s19ZetaQuad<16>(values, base + offset);
			for (unsigned base = 0; base <= 96; base += 16)
				for (unsigned offset = 0; offset < 8; offset += 4) rm2s19ZetaQuad<8>(values, base + offset);
			for (unsigned base = 0; base <= 48; base += 8) rm2s19ZetaQuad<4>(values, base);
			for (unsigned base = 64; base <= 80; base += 8) rm2s19ZetaQuad<4>(values, base);
			rm2s19ZetaQuad<4>(values, 96);
#else
			rm2s19ZetaStage<64>(values);
			rm2s19ZetaStage<32>(values);
			rm2s19ZetaStage<16>(values);
			for (unsigned base = 0; base <= 96; base += 16)
				for (unsigned offset = 0; offset < 8; offset += 2) rm2s19ZetaPair<8>(values, base + offset);
			for (unsigned base = 0; base <= 48; base += 8)
				for (unsigned offset = 0; offset < 4; offset += 2) rm2s19ZetaPair<4>(values, base + offset);
			for (unsigned base = 64; base <= 80; base += 8)
				for (unsigned offset = 0; offset < 4; offset += 2) rm2s19ZetaPair<4>(values, base + offset);
			rm2s19ZetaPair<4>(values, 96); rm2s19ZetaPair<4>(values, 98);
#endif
			for (unsigned base = 0; base <= 24; base += 4) rm2s19ZetaPair<2>(values, base);
			for (unsigned base = 32; base <= 40; base += 4) rm2s19ZetaPair<2>(values, base);
			rm2s19ZetaPair<2>(values, 48);
			for (unsigned base = 64; base <= 72; base += 4) rm2s19ZetaPair<2>(values, base);
			rm2s19ZetaPair<2>(values, 80); rm2s19ZetaPair<2>(values, 96);
#define RIFFLE_RM2S19_D1(Low) values[Low] = rm2s19Xor(values[Low], values[(Low) + 1])
			RIFFLE_RM2S19_D1(0);  RIFFLE_RM2S19_D1(2);  RIFFLE_RM2S19_D1(4);  RIFFLE_RM2S19_D1(6);
			RIFFLE_RM2S19_D1(8);  RIFFLE_RM2S19_D1(10); RIFFLE_RM2S19_D1(12); RIFFLE_RM2S19_D1(16);
			RIFFLE_RM2S19_D1(18); RIFFLE_RM2S19_D1(20); RIFFLE_RM2S19_D1(24); RIFFLE_RM2S19_D1(32);
			RIFFLE_RM2S19_D1(34); RIFFLE_RM2S19_D1(36); RIFFLE_RM2S19_D1(40); RIFFLE_RM2S19_D1(48);
			RIFFLE_RM2S19_D1(64); RIFFLE_RM2S19_D1(66); RIFFLE_RM2S19_D1(68); RIFFLE_RM2S19_D1(72);
			RIFFLE_RM2S19_D1(80); RIFFLE_RM2S19_D1(96);
#undef RIFFLE_RM2S19_D1
		}

		OC_FORCEINLINE void rm2s19FinishB(__m128i* values, __m128i* output) noexcept
		{
			rm2s19ZetaTransposePruned(values);
			output[0] = values[0];
			for (unsigned linear = 0; linear < 7; ++linear) output[1 + linear] = values[1U << linear];
			alignas(32) __m128i monomials[21];
			for (unsigned monomial = 0; monomial < 21; ++monomial)
				monomials[monomial] = values[Rm2Sub19MonomialPositions[monomial]];
			const auto q21 = rm2s19Xor(monomials[11], monomials[20]);
			const auto q22 = rm2s19Xor(monomials[15], q21);
			const auto q23 = rm2s19Xor(monomials[9], monomials[13]);
			const auto q24 = rm2s19Xor(monomials[10], q22);
			const auto q25 = rm2s19Xor(monomials[0], monomials[14]);
			const auto q26 = rm2s19Xor(monomials[2], monomials[7]);
			const auto q27 = rm2s19Xor(monomials[4], monomials[6]);
			const auto q28 = rm2s19Xor(monomials[17], q23);
			const auto q29 = rm2s19Xor(monomials[5], q24);
			const auto q30 = rm2s19Xor(monomials[3], monomials[19]);
			const auto q31 = rm2s19Xor(monomials[18], q25);
			const auto q32 = rm2s19Xor(monomials[1], monomials[16]);
			const auto q33 = rm2s19Xor(q26, q32);
			const auto q34 = rm2s19Xor(monomials[18], q26);
			const auto q35 = rm2s19Xor(monomials[3], q22);
			const auto q36 = rm2s19Xor(monomials[12], q27);
			const auto q37 = rm2s19Xor(q27, q29);
			const auto q38 = rm2s19Xor(q28, q37);
			const auto q39 = rm2s19Xor(monomials[8], q30);
			const auto q40 = rm2s19Xor(monomials[2], q31);
			const auto q41 = rm2s19Xor(monomials[3], monomials[9]);
			const auto q42 = rm2s19Xor(q23, q24);
			output[8] = rm2s19Xor(rm2s19Xor(rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[5], monomials[8]), q25), q32), q35), q36);
			output[9] = rm2s19Xor(rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[1], monomials[4]), monomials[14]), q39), q42);
			output[10] = rm2s19Xor(q33, q42);
			output[11] = rm2s19Xor(rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[16], monomials[17]), q29), q40), q41);
			output[12] = rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[10], monomials[12]), q28), q34);
			output[13] = rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[0], monomials[11]), q36), q41);
			output[14] = rm2s19Xor(rm2s19Xor(monomials[13], q34), q35);
			output[15] = rm2s19Xor(rm2s19Xor(monomials[8], monomials[16]), q38);
			output[16] = rm2s19Xor(rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[12], monomials[17]), q21), q30), q33);
			output[17] = rm2s19Xor(rm2s19Xor(monomials[19], q31), q38);
			output[18] = rm2s19Xor(rm2s19Xor(rm2s19Xor(monomials[15], q28), q39), q40);
		}

		template<std::size_t Point>
		OC_FORCEINLINE __m128i rm2s19GroupedAAddend(
			const __m128i table[4][16], const __m128i finalTable[8]) noexcept
		{
			constexpr u32 column = Rm2Sub19Columns[Point];
			constexpr unsigned i0 = ((column >> 0) & 1U) | (((column >> 16) & 1U) << 1) |
				(((column >> 6) & 1U) << 2) | (((column >> 2) & 1U) << 3);
			constexpr unsigned i1 = ((column >> 14) & 1U) | (((column >> 18) & 1U) << 1) |
				(((column >> 1) & 1U) << 2) | (((column >> 10) & 1U) << 3);
			constexpr unsigned i2 = ((column >> 9) & 1U) | (((column >> 12) & 1U) << 1) |
				(((column >> 3) & 1U) << 2) | (((column >> 5) & 1U) << 3);
			constexpr unsigned i3 = ((column >> 17) & 1U) | (((column >> 13) & 1U) << 1) |
				(((column >> 15) & 1U) << 2) | (((column >> 8) & 1U) << 3);
			constexpr unsigned i4 = ((column >> 11) & 1U) | (((column >> 7) & 1U) << 1) |
				(((column >> 4) & 1U) << 2);
			static_assert(i0 != 0);
			auto addend = table[0][i0];
			if constexpr (i1 != 0) addend = rm2s19Xor(addend, table[1][i1]);
			if constexpr (i2 != 0) addend = rm2s19Xor(addend, table[2][i2]);
			if constexpr (i3 != 0) addend = rm2s19Xor(addend, table[3][i3]);
			if constexpr (i4 != 0) addend = rm2s19Xor(addend, finalTable[i4]);
			return addend;
		}

		template<std::size_t ReversePoint, typename Emit>
		OC_FORCEINLINE void rm2s19EmitPoint(
			const block* input, __m128i* values, const __m128i table[4][16],
			const __m128i finalTable[8], u64 innerBase, Emit& emit) noexcept
		{
			constexpr std::size_t point = 127 - ReversePoint;
			const auto value = rm2s19Xor(input[point].mData, rm2s19GroupedAAddend<point>(table, finalTable));
			values[point] = value;
			emit(innerBase + point, block(value));
		}

		template<typename Emit, std::size_t... ReversePoints>
		OC_FORCEINLINE void rm2s19EmitPoints(
			const block* input, __m128i* values, const __m128i table[4][16],
			const __m128i finalTable[8], u64 innerBase, Emit& emit,
			std::index_sequence<ReversePoints...>) noexcept
		{
			(rm2s19EmitPoint<ReversePoints>(input, values, table, finalTable, innerBase, emit), ...);
		}

		template<typename Emit>
		OC_FORCEINLINE void rm2s19AddAGroupedAndEmit(
			const block* input, const __m128i* state, __m128i* values,
			u64 innerBase, Emit& emit) noexcept
		{
			alignas(32) __m128i groupedState[4][4]{
				{state[0], state[16], state[6], state[2]},
				{state[14], state[18], state[1], state[10]},
				{state[9], state[12], state[3], state[5]},
				{state[17], state[13], state[15], state[8]},
			};
			alignas(32) __m128i table[4][16];
			for (unsigned group = 0; group < 4; ++group) rm2s19BuildNibbleTable(groupedState[group], table[group]);
			alignas(32) __m128i finalState[3]{state[11], state[7], state[4]};
			alignas(32) __m128i finalTable[8];
			rm2s19BuildThreeBitTable(finalState, finalTable);
			rm2s19EmitPoints(input, values, table, finalTable, innerBase, emit, std::make_index_sequence<128>{});
		}

		template<typename Emit>
		OC_FORCEINLINE void rm2s19CopyAndEmit(
			const block* input, __m128i* values, u64 innerBase, Emit& emit) noexcept
		{
			for (u64 offset = 128; offset-- > 0;)
			{
				const auto value = input[offset].mData;
				values[offset] = value;
				emit(innerBase + offset, block(value));
			}
		}
	}

	class RiffleRm2SubS19Transpose
	{
	public:
		static constexpr u64 stepBlocks = 128;
		static constexpr u64 stateBlocks = 19;

		void init(u64 epochs, u64 coefficientSeed)
		{
			if (epochs == 0) throw std::invalid_argument("RM2Sub epoch count must be positive");
			mSchedules.resize(epochs);
			for (u64 epoch = 0; epoch < epochs; ++epoch)
			{
				u32 coefficient;
				do coefficient = static_cast<u32>(splitmix64(coefficientSeed)) & 0x7ffffU;
				while (coefficient == 0);
				mSchedules[epoch] = detail::rm2s19MakeFieldTransposeSchedule(coefficient);
			}
		}

		template<typename Emit>
		OC_FORCEINLINE void emitReverse(
			const block* __restrict input, u64 codeBlocks, Emit&& emit) const noexcept
		{
			const u64 epochs = codeBlocks / stepBlocks;
			alignas(32) __m128i state[19]{};
			alignas(32) __m128i syndrome[19];
			alignas(32) __m128i values[128];
			for (u64 epoch = epochs; epoch-- > 0;)
			{
				const block* node = input + epoch * stepBlocks;
				const u64 innerBase = epoch * stepBlocks;
				if (epoch + 1 == epochs) detail::rm2s19CopyAndEmit(node, values, innerBase, emit);
				else detail::rm2s19AddAGroupedAndEmit(node, state, values, innerBase, emit);
				if (epoch == 0) break;
				detail::rm2s19FinishB(values, syndrome);
				if (epoch + 1 == epochs) std::memcpy(state, syndrome, sizeof(state));
				else
				{
					detail::rm2s19FieldMultiplyTranspose(state, mSchedules[epoch]);
					for (unsigned bit = 0; bit < 19; ++bit) state[bit] = detail::rm2s19Xor(state[bit], syndrome[bit]);
				}
			}
		}

		const std::vector<detail::Rm2Sub19FieldTransposeSchedule>& schedules() const noexcept
		{
			return mSchedules;
		}

	private:
		std::vector<detail::Rm2Sub19FieldTransposeSchedule> mSchedules;

		static u64 splitmix64(u64& state) noexcept
		{
			u64 value = (state += 0x9e3779b97f4a7c15ULL);
			value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
			value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
			return value ^ (value >> 31);
		}
	};
}
