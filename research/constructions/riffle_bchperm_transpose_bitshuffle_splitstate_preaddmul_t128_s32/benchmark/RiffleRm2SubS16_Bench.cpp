#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <immintrin.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <pthread.h>
#include <sched.h>
#endif

namespace
{
    using Block = __m128i;

    constexpr std::size_t WordBlocks = std::size_t{1} << 21;
    constexpr std::size_t StepBlocks = 128;
    constexpr std::size_t StateBlocks = 16;
    constexpr std::size_t Epochs = WordBlocks / StepBlocks;
    constexpr std::uint16_t FieldModulusLow = 0x100b; // x^16+x^12+x^3+x+1.

    constexpr std::array<std::uint32_t, 8> QuadraticMasks{
        0x1e6718, 0x04f9b0, 0x0b7e1e, 0x0ee8a9,
        0x09389b, 0x0eb578, 0x116008, 0x041fce,
    };

    constexpr std::array<unsigned, 21> MonomialPositions{
        0x03, 0x05, 0x09, 0x11, 0x21, 0x41,
        0x06, 0x0a, 0x12, 0x22, 0x42,
        0x0c, 0x14, 0x24, 0x44,
        0x18, 0x28, 0x48,
        0x30, 0x50,
        0x60,
    };

    // A's coordinate columns, also the columns of B=A^T.  These are copied
    // verbatim from receipts/smaller_state/s16_rm2sub_selection.json.
    constexpr std::array<std::uint16_t, 128> Columns{
        0x0001, 0x0003, 0x0005, 0x1807, 0x0009, 0x940b, 0xa00d, 0x2c0f,
        0x0011, 0x8413, 0x9a15, 0x0617, 0x9e19, 0x8e1b, 0xa41d, 0xac1f,
        0x0021, 0xfd23, 0xa325, 0x4627, 0xb629, 0xdf2b, 0xb52d, 0xc42f,
        0x2a31, 0x5333, 0x1335, 0x7237, 0x0239, 0xef3b, 0x9b3d, 0x6e3f,
        0x0041, 0x3743, 0x8545, 0xaa47, 0x7f49, 0xdc4b, 0x5a4d, 0xe14f,
        0x5451, 0xe753, 0x4b55, 0xe057, 0xb559, 0x925b, 0x0a5d, 0x355f,
        0xab61, 0x6163, 0x8d65, 0x5f67, 0x6269, 0x3c6b, 0xe46d, 0xa26f,
        0xd571, 0x9b73, 0x6975, 0x3f77, 0x8279, 0x587b, 0x9e7d, 0x5c7f,
        0x0081, 0x2a83, 0xa585, 0x9787, 0x4f89, 0xf18b, 0x4a8d, 0xec8f,
        0x2d91, 0x8393, 0x1295, 0xa497, 0xfc99, 0xc69b, 0x639d, 0x419f,
        0x3da1, 0xeaa3, 0x3ba5, 0xf4a7, 0xc4a9, 0x87ab, 0x62ad, 0x39af,
        0x3ab1, 0x69b3, 0xa6b5, 0xedb7, 0x5db9, 0x9abb, 0x61bd, 0xbebf,
        0x41c1, 0x5cc3, 0x61c5, 0x64c7, 0x71c9, 0xf8cb, 0xf1cd, 0x60cf,
        0x38d1, 0xa1d3, 0x82d5, 0x03d7, 0x96d9, 0x9bdb, 0x8cdd, 0x99df,
        0xd7e1, 0x37e3, 0x54e5, 0xace7, 0x51e9, 0x25eb, 0x72ed, 0x1eef,
        0x84f1, 0xe0f3, 0x9df5, 0xe1f7, 0x9cf9, 0x6cfb, 0x25fd, 0xcdff,
    };

    struct FieldSchedule
    {
        alignas(32) std::array<std::uint16_t, StateBlocks> rowMasks{};
    };

    volatile std::uint64_t BenchmarkSink = 0;

    inline Block xorBlock(Block left, Block right) noexcept
    {
        return _mm_xor_si128(left, right);
    }

    inline bool equalBlock(Block left, Block right) noexcept
    {
        const Block difference = xorBlock(left, right);
        return _mm_testz_si128(difference, difference) != 0;
    }

    std::uint64_t splitmix64(std::uint64_t& state) noexcept
    {
        std::uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    std::uint16_t fieldMultiplyScalar(std::uint16_t left, std::uint16_t right) noexcept
    {
        std::uint16_t result = 0;
        for (unsigned bit = 0; bit < 16; ++bit)
        {
            result ^= left & static_cast<std::uint16_t>(
                0U - static_cast<unsigned>(right & 1));
            right >>= 1;
            const std::uint16_t carry = left >> 15;
            left = static_cast<std::uint16_t>(left << 1);
            left ^= FieldModulusLow & static_cast<std::uint16_t>(0U - carry);
        }
        return result;
    }

    FieldSchedule makeFieldSchedule(std::uint16_t coefficient) noexcept
    {
        FieldSchedule schedule{};
        for (unsigned input = 0; input < 16; ++input)
        {
            const std::uint16_t product = fieldMultiplyScalar(
                static_cast<std::uint16_t>(std::uint16_t{1} << input), coefficient);
            for (unsigned output = 0; output < 16; ++output)
                schedule.rowMasks[output] |= static_cast<std::uint16_t>(
                    ((product >> output) & 1U) << input);
        }
        return schedule;
    }

    inline void buildNibbleTable(const Block* input, Block* table) noexcept
    {
        const Block zero = _mm_setzero_si128();
        table[0] = zero;
        table[1] = input[0];
        table[2] = input[1];
        table[3] = xorBlock(input[0], input[1]);
        table[4] = input[2];
        table[5] = xorBlock(input[2], table[1]);
        table[6] = xorBlock(input[2], table[2]);
        table[7] = xorBlock(input[2], table[3]);
        table[8] = input[3];
        table[9] = xorBlock(input[3], table[1]);
        table[10] = xorBlock(input[3], table[2]);
        table[11] = xorBlock(input[3], table[3]);
        table[12] = xorBlock(input[3], table[4]);
        table[13] = xorBlock(input[3], table[5]);
        table[14] = xorBlock(input[3], table[6]);
        table[15] = xorBlock(input[3], table[7]);
    }

    inline void fieldMultiplyBitsliced(
        Block* state, const FieldSchedule& schedule) noexcept
    {
        alignas(32) Block table[4][16];
        for (unsigned group = 0; group < 4; ++group)
            buildNibbleTable(state + 4 * group, table[group]);

        alignas(32) Block output[16];
        for (unsigned row = 0; row < 16; ++row)
        {
            const std::uint16_t mask = schedule.rowMasks[row];
            Block value = xorBlock(table[0][mask & 15], table[1][(mask >> 4) & 15]);
            value = xorBlock(value, table[2][(mask >> 8) & 15]);
            output[row] = xorBlock(value, table[3][mask >> 12]);
        }
        std::memcpy(state, output, sizeof(output));
    }

    template<unsigned Distance>
    inline void zetaForwardStage(Block* values) noexcept
    {
        for (unsigned base = 0; base < 128; base += 2 * Distance)
        {
            if constexpr (Distance == 1)
            {
                values[base + 1] = xorBlock(values[base + 1], values[base]);
            }
            else
            {
                for (unsigned offset = 0; offset < Distance; offset += 2)
                {
                    const auto low = _mm256_load_si256(
                        reinterpret_cast<const __m256i*>(values + base + offset));
                    auto high = _mm256_load_si256(
                        reinterpret_cast<const __m256i*>(values + base + Distance + offset));
                    high = _mm256_xor_si256(high, low);
                    _mm256_store_si256(
                        reinterpret_cast<__m256i*>(values + base + Distance + offset), high);
                }
            }
        }
    }

    inline void zetaForward(Block* values) noexcept
    {
        zetaForwardStage<1>(values);
        zetaForwardStage<2>(values);
        zetaForwardStage<4>(values);
        zetaForwardStage<8>(values);
        zetaForwardStage<16>(values);
        zetaForwardStage<32>(values);
        zetaForwardStage<64>(values);
    }

    template<unsigned Distance>
    inline void zetaTransposeStage(Block* values) noexcept
    {
        for (unsigned base = 0; base < 128; base += 2 * Distance)
        {
            if constexpr (Distance == 1)
            {
                values[base] = xorBlock(values[base], values[base + 1]);
            }
            else
            {
                for (unsigned offset = 0; offset < Distance; offset += 2)
                {
                    auto low = _mm256_load_si256(
                        reinterpret_cast<const __m256i*>(values + base + offset));
                    const auto high = _mm256_load_si256(
                        reinterpret_cast<const __m256i*>(values + base + Distance + offset));
                    low = _mm256_xor_si256(low, high);
                    _mm256_store_si256(
                        reinterpret_cast<__m256i*>(values + base + offset), low);
                }
            }
        }
    }

    inline void zetaTranspose(Block* values) noexcept
    {
        zetaTransposeStage<64>(values);
        zetaTransposeStage<32>(values);
        zetaTransposeStage<16>(values);
        zetaTransposeStage<8>(values);
        zetaTransposeStage<4>(values);
        zetaTransposeStage<2>(values);
        zetaTransposeStage<1>(values);
    }

    template<unsigned Distance>
    inline void zetaTransposePrunedPair(Block* values, unsigned lowPosition) noexcept
    {
        auto low = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(values + lowPosition));
        const auto high = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(values + lowPosition + Distance));
        low = _mm256_xor_si256(low, high);
        _mm256_store_si256(reinterpret_cast<__m256i*>(values + lowPosition), low);
    }

#if defined(__AVX512F__)
    template<unsigned Distance>
    inline void zetaTransposePrunedQuad(Block* values, unsigned lowPosition) noexcept
    {
        auto low = _mm512_loadu_si512(values + lowPosition);
        const auto high = _mm512_loadu_si512(values + lowPosition + Distance);
        low = _mm512_xor_si512(low, high);
        _mm512_storeu_si512(values + lowPosition, low);
    }
#endif

    inline void zetaTransposePruned(Block* values) noexcept
    {
        // This is the exact backward slice needed for the constant, seven
        // linear, and 21 quadratic RM correlations.  Counts by stage are
        // 64,64,64,56,44,32,22 block XORs, versus 64 at every full stage.
#if defined(__AVX512F__)
        for (unsigned offset = 0; offset < 64; offset += 4)
            zetaTransposePrunedQuad<64>(values, offset);
        for (unsigned base = 0; base < 128; base += 64)
            for (unsigned offset = 0; offset < 32; offset += 4)
                zetaTransposePrunedQuad<32>(values, base + offset);
        for (unsigned base = 0; base < 128; base += 32)
            for (unsigned offset = 0; offset < 16; offset += 4)
                zetaTransposePrunedQuad<16>(values, base + offset);
        for (unsigned base = 0; base <= 96; base += 16)
            for (unsigned offset = 0; offset < 8; offset += 4)
                zetaTransposePrunedQuad<8>(values, base + offset);
        for (unsigned base = 0; base <= 48; base += 8)
            zetaTransposePrunedQuad<4>(values, base);
        for (unsigned base = 64; base <= 80; base += 8)
            zetaTransposePrunedQuad<4>(values, base);
        zetaTransposePrunedQuad<4>(values, 96);
#else
        zetaTransposeStage<64>(values);
        zetaTransposeStage<32>(values);
        zetaTransposeStage<16>(values);
        for (unsigned base = 0; base <= 96; base += 16)
            for (unsigned offset = 0; offset < 8; offset += 2)
                zetaTransposePrunedPair<8>(values, base + offset);
        for (unsigned base = 0; base <= 48; base += 8)
            for (unsigned offset = 0; offset < 4; offset += 2)
                zetaTransposePrunedPair<4>(values, base + offset);
        for (unsigned base = 64; base <= 80; base += 8)
            for (unsigned offset = 0; offset < 4; offset += 2)
                zetaTransposePrunedPair<4>(values, base + offset);
        zetaTransposePrunedPair<4>(values, 96);
        zetaTransposePrunedPair<4>(values, 98);
#endif

        for (unsigned base = 0; base <= 24; base += 4)
            zetaTransposePrunedPair<2>(values, base);
        for (unsigned base = 32; base <= 40; base += 4)
            zetaTransposePrunedPair<2>(values, base);
        zetaTransposePrunedPair<2>(values, 48);
        for (unsigned base = 64; base <= 72; base += 4)
            zetaTransposePrunedPair<2>(values, base);
        zetaTransposePrunedPair<2>(values, 80);
        zetaTransposePrunedPair<2>(values, 96);

#define RIFFLE_PRUNED_D1(Low) \
        values[Low] = xorBlock(values[Low], values[(Low) + 1])
        RIFFLE_PRUNED_D1(0);  RIFFLE_PRUNED_D1(2);
        RIFFLE_PRUNED_D1(4);  RIFFLE_PRUNED_D1(6);
        RIFFLE_PRUNED_D1(8);  RIFFLE_PRUNED_D1(10);
        RIFFLE_PRUNED_D1(12); RIFFLE_PRUNED_D1(16);
        RIFFLE_PRUNED_D1(18); RIFFLE_PRUNED_D1(20);
        RIFFLE_PRUNED_D1(24); RIFFLE_PRUNED_D1(32);
        RIFFLE_PRUNED_D1(34); RIFFLE_PRUNED_D1(36);
        RIFFLE_PRUNED_D1(40); RIFFLE_PRUNED_D1(48);
        RIFFLE_PRUNED_D1(64); RIFFLE_PRUNED_D1(66);
        RIFFLE_PRUNED_D1(68); RIFFLE_PRUNED_D1(72);
        RIFFLE_PRUNED_D1(80); RIFFLE_PRUNED_D1(96);
#undef RIFFLE_PRUNED_D1
    }

    template<std::uint32_t Mask>
    inline Block xorSelectedFixed(const Block* values) noexcept
    {
        std::uint32_t mask = Mask;
        constexpr unsigned first = std::countr_zero(Mask);
        Block result = values[first];
        mask &= mask - 1;
        while (mask)
        {
            const unsigned index = std::countr_zero(mask);
            result = xorBlock(result, values[index]);
            mask &= mask - 1;
        }
        return result;
    }

    inline Block xorSelected(const Block* values, std::uint32_t mask) noexcept
    {
        const unsigned first = std::countr_zero(mask);
        Block result = values[first];
        mask &= mask - 1;
        while (mask)
        {
            const unsigned index = std::countr_zero(mask);
            result = xorBlock(result, values[index]);
            mask &= mask - 1;
        }
        return result;
    }

    inline void addAWithRm(Block* word, const Block* state) noexcept
    {
        alignas(32) Block values[128];
        std::memset(values, 0, sizeof(values));
        values[0] = state[0];
        for (unsigned linear = 0; linear < 7; ++linear)
            values[1U << linear] = state[1 + linear];

        alignas(32) Block quadraticTable[2][16];
        buildNibbleTable(state + 8, quadraticTable[0]);
        buildNibbleTable(state + 12, quadraticTable[1]);
        for (unsigned monomial = 0; monomial < 21; ++monomial)
        {
            unsigned mask = 0;
            for (unsigned quadratic = 0; quadratic < 8; ++quadratic)
                mask |= ((QuadraticMasks[quadratic] >> monomial) & 1U) << quadratic;
            values[MonomialPositions[monomial]] = xorBlock(
                quadraticTable[0][mask & 15], quadraticTable[1][mask >> 4]);
        }

        zetaForward(values);
        for (unsigned point = 0; point < 128; point += 2)
        {
            auto output = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(word + point));
            const auto addend = _mm256_load_si256(
                reinterpret_cast<const __m256i*>(values + point));
            output = _mm256_xor_si256(output, addend);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(word + point), output);
        }
    }

    inline void addAWithTables(Block* word, const Block* state) noexcept
    {
        alignas(32) Block table[4][16];
        for (unsigned group = 0; group < 4; ++group)
            buildNibbleTable(state + 4 * group, table[group]);
        for (unsigned point = 0; point < 128; ++point)
        {
            const std::uint16_t column = Columns[point];
            Block value = xorBlock(table[0][column & 15], table[1][(column >> 4) & 15]);
            value = xorBlock(value, table[2][(column >> 8) & 15]);
            value = xorBlock(value, table[3][column >> 12]);
            word[point] = xorBlock(word[point], value);
        }
    }

    template<unsigned Bits>
    inline void buildFullTable(const Block* input, Block* table) noexcept
    {
        static_assert(Bits >= 1 && Bits <= 6);
        table[0] = _mm_setzero_si128();
        for (unsigned mask = 1; mask < (1U << Bits); ++mask)
        {
            const unsigned bit = std::countr_zero(mask);
            const unsigned remainder = mask & (mask - 1);
            table[mask] = remainder == 0 ? input[bit] : xorBlock(table[remainder], input[bit]);
        }
    }

    inline void addAWithTables556(Block* word, const Block* state) noexcept
    {
        alignas(32) Block low[32];
        alignas(32) Block middle[32];
        alignas(32) Block high[64];
        buildFullTable<5>(state, low);
        buildFullTable<5>(state + 5, middle);
        buildFullTable<6>(state + 10, high);
        for (unsigned point = 0; point < 128; ++point)
        {
            const std::uint16_t column = Columns[point];
            Block value = xorBlock(low[column & 31], middle[(column >> 5) & 31]);
            value = xorBlock(value, high[column >> 10]);
            word[point] = xorBlock(word[point], value);
        }
    }

    inline void computeBWithRm(const Block* word, Block* output) noexcept
    {
        alignas(32) Block values[128];
        for (unsigned point = 0; point < 128; point += 2)
        {
            const auto input = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(word + point));
            _mm256_store_si256(reinterpret_cast<__m256i*>(values + point), input);
        }
        zetaTranspose(values);
        output[0] = values[0];
        for (unsigned linear = 0; linear < 7; ++linear)
            output[1 + linear] = values[1U << linear];

        alignas(32) Block monomials[21];
        for (unsigned monomial = 0; monomial < 21; ++monomial)
            monomials[monomial] = values[MonomialPositions[monomial]];
        output[8] = xorSelectedFixed<QuadraticMasks[0]>(monomials);
        output[9] = xorSelectedFixed<QuadraticMasks[1]>(monomials);
        output[10] = xorSelectedFixed<QuadraticMasks[2]>(monomials);
        output[11] = xorSelectedFixed<QuadraticMasks[3]>(monomials);
        output[12] = xorSelectedFixed<QuadraticMasks[4]>(monomials);
        output[13] = xorSelectedFixed<QuadraticMasks[5]>(monomials);
        output[14] = xorSelectedFixed<QuadraticMasks[6]>(monomials);
        output[15] = xorSelectedFixed<QuadraticMasks[7]>(monomials);
    }

    inline void finishBWithPrunedRm(Block* values, Block* output) noexcept
    {
        zetaTransposePruned(values);
        output[0] = values[0];
        for (unsigned linear = 0; linear < 7; ++linear)
            output[1 + linear] = values[1U << linear];

        alignas(32) Block monomials[21];
        for (unsigned monomial = 0; monomial < 21; ++monomial)
            monomials[monomial] = values[MonomialPositions[monomial]];

        // Exact 42-XOR circuit for the eight selected quadratic masks.
        // Direct evaluation costs 75 XORs.  Signal numbers match the
        // 21 inputs followed by the 15 shared gates below.
        const Block a21 = xorBlock(monomials[3], monomials[13]);
        const Block a22 = xorBlock(monomials[19], a21);
        const Block a23 = xorBlock(monomials[4], monomials[12]);
        const Block a24 = xorBlock(monomials[7], monomials[11]);
        const Block a25 = xorBlock(monomials[8], monomials[18]);
        const Block a26 = xorBlock(monomials[17], a22);
        const Block a27 = xorBlock(monomials[9], monomials[10]);
        const Block a28 = xorBlock(monomials[5], monomials[15]);
        const Block a29 = xorBlock(monomials[14], a26);
        const Block a30 = xorBlock(a23, a28);
        const Block a31 = xorBlock(monomials[0], a24);
        const Block a32 = xorBlock(monomials[6], a25);
        const Block a33 = xorBlock(monomials[2], a27);
        const Block a34 = xorBlock(monomials[1], monomials[16]);
        const Block a35 = xorBlock(a23, a34);

        output[8] = xorBlock(xorBlock(xorBlock(xorBlock(
            monomials[4], monomials[20]), a25), a27), a29);
        output[9] = xorBlock(xorBlock(xorBlock(xorBlock(
            monomials[13], monomials[14]), a24), a25), a30);
        output[10] = xorBlock(xorBlock(xorBlock(
            monomials[11], a29), a33), a35);
        output[11] = xorBlock(xorBlock(xorBlock(
            monomials[18], a28), a29), a31);
        output[12] = xorBlock(xorBlock(a22, a31), a35);
        output[13] = xorBlock(xorBlock(xorBlock(
            monomials[10], a26), a30), a32);
        output[14] = xorBlock(xorBlock(xorBlock(
            monomials[14], monomials[16]), monomials[20]), a21);
        output[15] = xorBlock(xorBlock(xorBlock(xorBlock(xorBlock(
            monomials[1], monomials[3]), monomials[12]), a24), a32), a33);
    }

    inline void computeBWithPrunedRm(const Block* word, Block* output) noexcept
    {
        alignas(32) Block values[128];
        for (unsigned point = 0; point < 128; point += 2)
        {
            const auto input = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(word + point));
            _mm256_store_si256(reinterpret_cast<__m256i*>(values + point), input);
        }
        finishBWithPrunedRm(values, output);
    }

    inline void addAAndComputeBWithTables(
        Block* word, const Block* state, Block* output) noexcept
    {
        alignas(32) Block table[4][16];
        alignas(32) Block values[128];
        for (unsigned group = 0; group < 4; ++group)
            buildNibbleTable(state + 4 * group, table[group]);
        for (unsigned point = 0; point < 128; ++point)
        {
            const std::uint16_t column = Columns[point];
            Block addend = xorBlock(
                table[0][column & 15], table[1][(column >> 4) & 15]);
            addend = xorBlock(addend, table[2][(column >> 8) & 15]);
            addend = xorBlock(addend, table[3][column >> 12]);
            const Block encoded = xorBlock(word[point], addend);
            word[point] = encoded;
            values[point] = encoded;
        }
        finishBWithPrunedRm(values, output);
    }

    template<std::size_t Point>
    inline Block groupedAAddend(const Block table[4][16]) noexcept
    {
        constexpr std::uint16_t column = Columns[Point];
        constexpr unsigned index0 =
            ((column >> 0) & 1U) | (((column >> 2) & 1U) << 1) |
            (((column >> 3) & 1U) << 2) | (((column >> 4) & 1U) << 3);
        constexpr unsigned index1 =
            ((column >> 1) & 1U) | (((column >> 10) & 1U) << 1) |
            (((column >> 13) & 1U) << 2) | (((column >> 15) & 1U) << 3);
        constexpr unsigned index2 =
            ((column >> 5) & 1U) | (((column >> 9) & 1U) << 1) |
            (((column >> 11) & 1U) << 2) | (((column >> 12) & 1U) << 3);
        constexpr unsigned index3 =
            ((column >> 6) & 1U) | (((column >> 7) & 1U) << 1) |
            (((column >> 8) & 1U) << 2) | (((column >> 14) & 1U) << 3);
        static_assert(index0 != 0); // The constant RM coordinate is in group zero.

        Block addend = table[0][index0];
        if constexpr (index1 != 0)
            addend = xorBlock(addend, table[1][index1]);
        if constexpr (index2 != 0)
            addend = xorBlock(addend, table[2][index2]);
        if constexpr (index3 != 0)
            addend = xorBlock(addend, table[3][index3]);
        return addend;
    }

    template<std::size_t Point>
    inline void addAGroupedPoint(
        Block* word, Block* values, const Block table[4][16]) noexcept
    {
        const Block encoded = xorBlock(word[Point], groupedAAddend<Point>(table));
        word[Point] = encoded;
        values[Point] = encoded;
    }

    template<std::size_t... Points>
    inline void addAGroupedPoints(
        Block* word,
        Block* values,
        const Block table[4][16],
        std::index_sequence<Points...>) noexcept
    {
        (addAGroupedPoint<Points>(word, values, table), ...);
    }

    template<std::size_t... Points>
    inline void addAGroupedPointsOnly(
        Block* word,
        const Block table[4][16],
        std::index_sequence<Points...>) noexcept
    {
        ((word[Points] = xorBlock(word[Points], groupedAAddend<Points>(table))), ...);
    }

    inline void addAWithGroupedTables(Block* word, const Block* state) noexcept
    {
        alignas(32) Block groupedState[4][4]{
            {state[0], state[2], state[3], state[4]},
            {state[1], state[10], state[13], state[15]},
            {state[5], state[9], state[11], state[12]},
            {state[6], state[7], state[8], state[14]},
        };
        alignas(32) Block table[4][16];
        for (unsigned group = 0; group < 4; ++group)
            buildNibbleTable(groupedState[group], table[group]);
        addAGroupedPointsOnly(word, table, std::make_index_sequence<128>{});
    }

    inline void addAAndComputeBWithGroupedTables(
        Block* word, const Block* state, Block* output) noexcept
    {
        alignas(32) Block groupedState[4][4]{
            {state[0], state[2], state[3], state[4]},
            {state[1], state[10], state[13], state[15]},
            {state[5], state[9], state[11], state[12]},
            {state[6], state[7], state[8], state[14]},
        };
        alignas(32) Block table[4][16];
        alignas(32) Block values[128];
        for (unsigned group = 0; group < 4; ++group)
            buildNibbleTable(groupedState[group], table[group]);
        addAGroupedPoints(
            word, values, table, std::make_index_sequence<128>{});
        finishBWithPrunedRm(values, output);
    }

#if defined(_MSC_VER)
#define RIFFLE_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define RIFFLE_NOINLINE __attribute__((noinline))
#else
#define RIFFLE_NOINLINE
#endif
    RIFFLE_NOINLINE void addAAndComputeBWithPaar(
        Block* word, const Block* state, Block* output) noexcept
    {
        alignas(32) Block values[128];
#include "RiffleRm2SubS16PaarA.inc"
        finishBWithPrunedRm(values, output);
    }
#undef RIFFLE_NOINLINE

    inline void addAAndComputeBWithRm(
        Block* word, const Block* state, Block* output) noexcept
    {
        alignas(32) Block values[128];
        std::memset(values, 0, sizeof(values));
        values[0] = state[0];
        for (unsigned linear = 0; linear < 7; ++linear)
            values[1U << linear] = state[1 + linear];

        alignas(32) Block quadraticTable[2][16];
        buildNibbleTable(state + 8, quadraticTable[0]);
        buildNibbleTable(state + 12, quadraticTable[1]);
        for (unsigned monomial = 0; monomial < 21; ++monomial)
        {
            unsigned mask = 0;
            for (unsigned quadratic = 0; quadratic < 8; ++quadratic)
                mask |= ((QuadraticMasks[quadratic] >> monomial) & 1U) << quadratic;
            values[MonomialPositions[monomial]] = xorBlock(
                quadraticTable[0][mask & 15], quadraticTable[1][mask >> 4]);
        }
        zetaForward(values);

        for (unsigned point = 0; point < 128; point += 2)
        {
            auto encoded = _mm256_xor_si256(
                _mm256_loadu_si256(reinterpret_cast<const __m256i*>(word + point)),
                _mm256_load_si256(reinterpret_cast<const __m256i*>(values + point)));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(word + point), encoded);
            _mm256_store_si256(reinterpret_cast<__m256i*>(values + point), encoded);
        }
        finishBWithPrunedRm(values, output);
    }

    inline void addADense(Block* word, const Block* state) noexcept
    {
        for (unsigned point = 0; point < 128; ++point)
        {
            Block value = _mm_setzero_si128();
            std::uint16_t column = Columns[point];
            while (column)
            {
                const unsigned index = std::countr_zero(column);
                value = xorBlock(value, state[index]);
                column &= column - 1;
            }
            word[point] = xorBlock(word[point], value);
        }
    }

    inline void computeBDense(const Block* word, Block* output) noexcept
    {
        std::fill_n(output, 16, _mm_setzero_si128());
        for (unsigned point = 0; point < 128; ++point)
        {
            std::uint16_t column = Columns[point];
            while (column)
            {
                const unsigned index = std::countr_zero(column);
                output[index] = xorBlock(output[index], word[point]);
                column &= column - 1;
            }
        }
    }

    inline void fieldMultiplyDense(Block* state, const FieldSchedule& schedule) noexcept
    {
        alignas(32) Block output[16];
        for (unsigned row = 0; row < 16; ++row)
            output[row] = xorSelected(state, schedule.rowMasks[row]);
        std::memcpy(state, output, sizeof(output));
    }

    enum class AImplementation {
        rm, table, table556, fusedPruned, paarPruned, rmFusedPruned,
        groupedFusedPruned
    };

    void applyInnerFusedPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept;
    void applyInnerPaarPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept;
    void applyInnerRmFusedPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept;
    void applyInnerGroupedFusedPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept;

    template<AImplementation AImpl>
    void applyInner(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept
    {
        if constexpr (AImpl == AImplementation::fusedPruned)
        {
            applyInnerFusedPruned(word, epochs, schedule);
            return;
        }
        if constexpr (AImpl == AImplementation::paarPruned)
        {
            applyInnerPaarPruned(word, epochs, schedule);
            return;
        }
        if constexpr (AImpl == AImplementation::rmFusedPruned)
        {
            applyInnerRmFusedPruned(word, epochs, schedule);
            return;
        }
        if constexpr (AImpl == AImplementation::groupedFusedPruned)
        {
            applyInnerGroupedFusedPruned(word, epochs, schedule);
            return;
        }
        alignas(32) Block state[16]{};
        alignas(32) Block syndrome[16];
        for (std::size_t epoch = 0; epoch < epochs; ++epoch)
        {
            Block* node = word + epoch * StepBlocks;
            if (epoch != 0)
            {
                if constexpr (AImpl == AImplementation::rm)
                    addAWithRm(node, state);
                else if constexpr (AImpl == AImplementation::table556)
                    addAWithTables556(node, state);
                else
                    addAWithTables(node, state);
            }
            if (epoch + 1 == epochs)
                break;
            computeBWithRm(node, syndrome);
            fieldMultiplyBitsliced(state, schedule[epoch]);
            for (unsigned bit = 0; bit < 16; ++bit)
                state[bit] = xorBlock(state[bit], syndrome[bit]);
        }
        BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
    }

    void applyInnerDense(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept
    {
        alignas(32) Block state[16]{};
        alignas(32) Block syndrome[16];
        for (std::size_t epoch = 0; epoch < epochs; ++epoch)
        {
            Block* node = word + epoch * StepBlocks;
            if (epoch != 0)
                addADense(node, state);
            if (epoch + 1 == epochs)
                break;
            computeBDense(node, syndrome);
            fieldMultiplyDense(state, schedule[epoch]);
            for (unsigned bit = 0; bit < 16; ++bit)
                state[bit] = xorBlock(state[bit], syndrome[bit]);
        }
        BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
    }

    void applyInnerFusedPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept
    {
        alignas(32) Block state[16]{};
        alignas(32) Block syndrome[16];
        if (epochs == 1)
            return;

        computeBWithPrunedRm(word, syndrome);
        fieldMultiplyBitsliced(state, schedule[0]);
        for (unsigned bit = 0; bit < 16; ++bit)
            state[bit] = xorBlock(state[bit], syndrome[bit]);

        for (std::size_t epoch = 1; epoch + 1 < epochs; ++epoch)
        {
            addAAndComputeBWithTables(
                word + epoch * StepBlocks, state, syndrome);
            fieldMultiplyBitsliced(state, schedule[epoch]);
            for (unsigned bit = 0; bit < 16; ++bit)
                state[bit] = xorBlock(state[bit], syndrome[bit]);
        }
        addAWithGroupedTables(word + (epochs - 1) * StepBlocks, state);
        BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
    }

    void applyInnerPaarPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept
    {
        alignas(32) Block state[16]{};
        alignas(32) Block syndrome[16];
        if (epochs == 1)
            return;

        computeBWithPrunedRm(word, syndrome);
        fieldMultiplyBitsliced(state, schedule[0]);
        for (unsigned bit = 0; bit < 16; ++bit)
            state[bit] = xorBlock(state[bit], syndrome[bit]);

        for (std::size_t epoch = 1; epoch + 1 < epochs; ++epoch)
        {
            addAAndComputeBWithPaar(
                word + epoch * StepBlocks, state, syndrome);
            fieldMultiplyBitsliced(state, schedule[epoch]);
            for (unsigned bit = 0; bit < 16; ++bit)
                state[bit] = xorBlock(state[bit], syndrome[bit]);
        }
        addAWithTables(word + (epochs - 1) * StepBlocks, state);
        BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
    }

    void applyInnerRmFusedPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept
    {
        alignas(32) Block state[16]{};
        alignas(32) Block syndrome[16];
        if (epochs == 1)
            return;

        computeBWithPrunedRm(word, syndrome);
        fieldMultiplyBitsliced(state, schedule[0]);
        for (unsigned bit = 0; bit < 16; ++bit)
            state[bit] = xorBlock(state[bit], syndrome[bit]);

        for (std::size_t epoch = 1; epoch + 1 < epochs; ++epoch)
        {
            addAAndComputeBWithRm(
                word + epoch * StepBlocks, state, syndrome);
            fieldMultiplyBitsliced(state, schedule[epoch]);
            for (unsigned bit = 0; bit < 16; ++bit)
                state[bit] = xorBlock(state[bit], syndrome[bit]);
        }
        addAWithRm(word + (epochs - 1) * StepBlocks, state);
        BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
    }

    void applyInnerGroupedFusedPruned(
        Block* word, std::size_t epochs, const FieldSchedule* schedule) noexcept
    {
        alignas(32) Block state[16]{};
        alignas(32) Block syndrome[16];
        if (epochs == 1)
            return;

        computeBWithPrunedRm(word, syndrome);
        std::memcpy(state, syndrome, sizeof(state));

        for (std::size_t epoch = 1; epoch + 1 < epochs; ++epoch)
        {
            addAAndComputeBWithGroupedTables(
                word + epoch * StepBlocks, state, syndrome);
            fieldMultiplyBitsliced(state, schedule[epoch]);
            for (unsigned bit = 0; bit < 16; ++bit)
                state[bit] = xorBlock(state[bit], syndrome[bit]);
        }
        addAWithGroupedTables(word + (epochs - 1) * StepBlocks, state);
        BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
    }

    void pinThread(unsigned cpu)
    {
#ifdef _WIN32
        if (!SetThreadAffinityMask(GetCurrentThread(), DWORD_PTR{1} << cpu))
            throw std::runtime_error("SetThreadAffinityMask failed");
#else
        cpu_set_t set;
        CPU_ZERO(&set);
        CPU_SET(cpu, &set);
        if (pthread_setaffinity_np(pthread_self(), sizeof(set), &set) != 0)
            throw std::runtime_error("pthread_setaffinity_np failed");
#endif
    }

    void verify(const std::vector<FieldSchedule>& schedule)
    {
        constexpr std::size_t TestEpochs = 5;
        std::vector<Block> source(TestEpochs * StepBlocks);
        std::uint64_t randomState = 0x682130ff903a17b5ULL;
        for (Block& value : source)
        {
            const std::uint64_t low = splitmix64(randomState);
            const std::uint64_t high = splitmix64(randomState);
            value = _mm_set_epi64x(static_cast<long long>(high), static_cast<long long>(low));
        }
        auto dense = source;
        auto rm = source;
        auto table = source;
        auto table556 = source;
        auto fusedPruned = source;
        auto paarPruned = source;
        auto rmFusedPruned = source;
        auto groupedFusedPruned = source;
        applyInnerDense(dense.data(), TestEpochs, schedule.data());
        applyInner<AImplementation::rm>(rm.data(), TestEpochs, schedule.data());
        applyInner<AImplementation::table>(table.data(), TestEpochs, schedule.data());
        applyInner<AImplementation::table556>(
            table556.data(), TestEpochs, schedule.data());
        applyInnerFusedPruned(fusedPruned.data(), TestEpochs, schedule.data());
        applyInnerPaarPruned(paarPruned.data(), TestEpochs, schedule.data());
        applyInnerRmFusedPruned(
            rmFusedPruned.data(), TestEpochs, schedule.data());
        applyInnerGroupedFusedPruned(
            groupedFusedPruned.data(), TestEpochs, schedule.data());
        for (std::size_t index = 0; index < source.size(); ++index)
            if (!equalBlock(dense[index], rm[index]) ||
                !equalBlock(dense[index], table[index]) ||
                !equalBlock(dense[index], table556[index]) ||
                !equalBlock(dense[index], fusedPruned[index]) ||
                !equalBlock(dense[index], paarPruned[index]) ||
                !equalBlock(dense[index], rmFusedPruned[index]) ||
                !equalBlock(dense[index], groupedFusedPruned[index]))
                throw std::runtime_error("optimized inner disagrees with dense reference");
    }

    std::uint64_t checksum(const Block* word, std::size_t blocks) noexcept
    {
        std::uint64_t result = 0;
        for (std::size_t index = 0; index < blocks; index += 4096)
        {
            result ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(word[index]));
            result ^= static_cast<std::uint64_t>(
                _mm_extract_epi64(word[std::min(index + 127, blocks - 1)], 1));
        }
        return result;
    }

    template<AImplementation AImpl>
    void benchmark(unsigned trials, unsigned cpu)
    {
        pinThread(cpu);
        std::uint64_t randomState = 0x10389ca7b615f24dULL;
        std::vector<Block> source(WordBlocks);
        std::vector<Block> working(WordBlocks);
        for (Block& value : source)
        {
            const std::uint64_t low = splitmix64(randomState);
            const std::uint64_t high = splitmix64(randomState);
            value = _mm_set_epi64x(static_cast<long long>(high), static_cast<long long>(low));
        }

        std::vector<FieldSchedule> schedule(Epochs - 1);
        for (FieldSchedule& entry : schedule)
        {
            std::uint16_t coefficient = static_cast<std::uint16_t>(splitmix64(randomState));
            if (coefficient == 0)
                coefficient = 1;
            entry = makeFieldSchedule(coefficient);
        }
        verify(schedule);

        working = source;
        applyInner<AImpl>(working.data(), Epochs, schedule.data());

        std::vector<double> samples;
        samples.reserve(trials);
        for (unsigned trial = 0; trial < trials; ++trial)
        {
            working = source; // Deliberately outside the timed interval.
            const auto start = std::chrono::steady_clock::now();
            applyInner<AImpl>(working.data(), Epochs, schedule.data());
            const auto stop = std::chrono::steady_clock::now();
            samples.push_back(std::chrono::duration<double, std::milli>(stop - start).count());
            BenchmarkSink ^= checksum(working.data(), working.size());
        }
        std::sort(samples.begin(), samples.end());
        const double median = samples[samples.size() / 2];
        const double minimum = samples.front();
        const double gib = static_cast<double>(WordBlocks * sizeof(Block)) /
            static_cast<double>(std::uint64_t{1} << 30);

        std::cout << std::fixed << std::setprecision(6)
                  << "{\n"
                  << "  \"schema\": \"riffle-rm2sub-s16-inner-benchmark-v1\",\n"
                  << "  \"construction\": \"SplitState-PreAddMul-RM2Sub t=128 s=16\",\n"
                  << "  \"a_implementation\": \""
                  << (AImpl == AImplementation::groupedFusedPruned ?
                      "optimal_grouped_4bit_tables_fused_with_pruned_b" :
                      AImpl == AImplementation::rmFusedPruned ?
                      "rm_zeta_avx2_inplace_fused_with_pruned_b" :
                      AImpl == AImplementation::paarPruned ?
                      "paar_299_xors_fused_with_pruned_b" :
                      AImpl == AImplementation::fusedPruned ?
                      "nibble_tables_fused_with_pruned_b" :
                      AImpl == AImplementation::rm ? "rm_zeta_avx2" :
                      AImpl == AImplementation::table556 ? "tables_5_5_6" :
                      "nibble_tables")
                  << "\",\n"
                  << "  \"b_implementation\": \""
                  << ((AImpl == AImplementation::fusedPruned ||
                       AImpl == AImplementation::paarPruned ||
                       AImpl == AImplementation::rmFusedPruned ||
                       AImpl == AImplementation::groupedFusedPruned) ?
                      "pruned_rm_transpose_zeta_avx2" : "rm_transpose_zeta_avx2")
                  << "\",\n"
                  << "  \"field_implementation\": \"bitsliced_nibble_tables\",\n"
                  << "  \"word_blocks\": " << WordBlocks << ",\n"
                  << "  \"step_blocks\": " << StepBlocks << ",\n"
                  << "  \"state_blocks\": " << StateBlocks << ",\n"
                  << "  \"epochs\": " << Epochs << ",\n"
                  << "  \"trials\": " << trials << ",\n"
                  << "  \"median_ms\": " << median << ",\n"
                  << "  \"minimum_ms\": " << minimum << ",\n"
                  << "  \"median_input_gib_per_second\": " << gib / (median / 1000.0) << ",\n"
                  << "  \"checksum\": \"0x" << std::hex << BenchmarkSink << std::dec << "\",\n"
                  << "  \"correctness\": \"PASS: dense, RM, 4x4-table, and 5+5+6-table implementations agree\"\n"
                  << "}\n";
    }

    enum class Phase { aTable, bRm, field };

    template<Phase Selected>
    void benchmarkPhase(unsigned trials, unsigned cpu)
    {
        pinThread(cpu);
        std::uint64_t randomState = 0x10389ca7b615f24dULL;
        std::vector<Block> source(WordBlocks);
        std::vector<Block> working(WordBlocks);
        for (Block& value : source)
        {
            const std::uint64_t low = splitmix64(randomState);
            const std::uint64_t high = splitmix64(randomState);
            value = _mm_set_epi64x(static_cast<long long>(high), static_cast<long long>(low));
        }
        std::vector<FieldSchedule> schedule(Epochs - 1);
        for (FieldSchedule& entry : schedule)
        {
            std::uint16_t coefficient = static_cast<std::uint16_t>(splitmix64(randomState));
            if (coefficient == 0)
                coefficient = 1;
            entry = makeFieldSchedule(coefficient);
        }
        verify(schedule);

        std::vector<double> samples;
        samples.reserve(trials);
        for (unsigned trial = 0; trial < trials; ++trial)
        {
            working = source;
            alignas(32) Block state[16];
            alignas(32) Block syndrome[16];
            std::memcpy(state, source.data(), sizeof(state));
            const auto start = std::chrono::steady_clock::now();
            for (std::size_t epoch = 0; epoch + 1 < Epochs; ++epoch)
            {
                Block* node = working.data() + epoch * StepBlocks;
                if constexpr (Selected == Phase::aTable)
                {
                    addAWithTables(node, state);
                    state[epoch & 15] = xorBlock(state[epoch & 15], node[epoch & 127]);
                }
                else if constexpr (Selected == Phase::bRm)
                {
                    computeBWithRm(node, syndrome);
                    for (unsigned bit = 0; bit < 16; ++bit)
                        state[bit] = xorBlock(state[bit], syndrome[bit]);
                }
                else
                {
                    fieldMultiplyBitsliced(state, schedule[epoch]);
                }
            }
            const auto stop = std::chrono::steady_clock::now();
            samples.push_back(std::chrono::duration<double, std::milli>(stop - start).count());
            BenchmarkSink ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(state[0]));
            if constexpr (Selected == Phase::aTable)
                BenchmarkSink ^= checksum(working.data(), working.size());
        }
        std::sort(samples.begin(), samples.end());
        std::cout << std::fixed << std::setprecision(6)
                  << "{\n"
                  << "  \"schema\": \"riffle-rm2sub-s16-phase-benchmark-v1\",\n"
                  << "  \"phase\": \""
                  << (Selected == Phase::aTable ? "A_nibble_tables" :
                      Selected == Phase::bRm ? "B_rm_transpose_zeta_avx2" :
                      "GF16_bitsliced_nibble_tables")
                  << "\",\n"
                  << "  \"calls\": " << Epochs - 1 << ",\n"
                  << "  \"trials\": " << trials << ",\n"
                  << "  \"median_ms\": " << samples[samples.size() / 2] << ",\n"
                  << "  \"minimum_ms\": " << samples.front() << ",\n"
                  << "  \"checksum\": \"0x" << std::hex << BenchmarkSink << std::dec << "\"\n"
                  << "}\n";
    }
}

int main(int argc, char** argv)
{
    try
    {
        std::string implementation = "grouped";
        unsigned trials = 15;
        unsigned cpu = 15;
        for (int index = 1; index < argc; ++index)
        {
            const std::string option = argv[index];
            auto value = [&]() -> std::string {
                if (++index == argc)
                    throw std::invalid_argument("missing option value");
                return argv[index];
            };
            if (option == "--implementation")
                implementation = value();
            else if (option == "--trials")
                trials = static_cast<unsigned>(std::stoul(value()));
            else if (option == "--cpu")
                cpu = static_cast<unsigned>(std::stoul(value()));
            else
                throw std::invalid_argument("unknown option: " + option);
        }
        if (trials == 0)
            throw std::invalid_argument("trials must be positive");
        if (implementation == "rm")
            benchmark<AImplementation::rm>(trials, cpu);
        else if (implementation == "table")
            benchmark<AImplementation::table>(trials, cpu);
        else if (implementation == "table556")
            benchmark<AImplementation::table556>(trials, cpu);
        else if (implementation == "fused")
            benchmark<AImplementation::fusedPruned>(trials, cpu);
        else if (implementation == "paar")
            benchmark<AImplementation::paarPruned>(trials, cpu);
        else if (implementation == "rm-fused")
            benchmark<AImplementation::rmFusedPruned>(trials, cpu);
        else if (implementation == "grouped")
            benchmark<AImplementation::groupedFusedPruned>(trials, cpu);
        else if (implementation == "a-table")
            benchmarkPhase<Phase::aTable>(trials, cpu);
        else if (implementation == "b-rm")
            benchmarkPhase<Phase::bRm>(trials, cpu);
        else if (implementation == "field")
            benchmarkPhase<Phase::field>(trials, cpu);
        else
            throw std::invalid_argument(
                "implementation must be rm, table, table556, fused, paar, rm-fused, grouped, a-table, b-rm, or field");
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
