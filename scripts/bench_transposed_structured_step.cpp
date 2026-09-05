#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <intrin.h>
#include <immintrin.h>
#include <windows.h>

namespace {

struct alignas(16) Block {
    __m128i value;
};

inline Block zero_block() noexcept { return {_mm_setzero_si128()}; }

inline Block xor_block(Block left, Block right) noexcept {
    return {_mm_xor_si128(left.value, right.value)};
}

inline bool equal_block(Block left, Block right) noexcept {
    const __m128i difference = _mm_xor_si128(left.value, right.value);
    return _mm_testz_si128(difference, difference) != 0;
}

struct Product128 {
    std::uint64_t low;
    std::uint64_t high;
};

inline Product128 carryless_multiply(std::uint64_t left, std::uint64_t right) noexcept {
    const __m128i product = _mm_clmulepi64_si128(
        _mm_cvtsi64_si128(static_cast<long long>(left)),
        _mm_cvtsi64_si128(static_cast<long long>(right)), 0x00);
    return {
        static_cast<std::uint64_t>(_mm_cvtsi128_si64(product)),
        static_cast<std::uint64_t>(_mm_extract_epi64(product, 1)),
    };
}

inline Product128 xor_product(Product128 left, Product128 right) noexcept {
    return {left.low ^ right.low, left.high ^ right.high};
}

inline Product128 shifted_word(std::uint64_t value, unsigned shift) noexcept {
    if (shift == 0) {
        return {value, 0};
    }
    return {value << shift, value >> (64 - shift)};
}

template <unsigned D>
constexpr std::uint64_t width_mask() noexcept {
    static_assert(D > 0 && D <= 64);
    if constexpr (D == 64) {
        return ~std::uint64_t{0};
    } else {
        return (std::uint64_t{1} << D) - 1;
    }
}

template <unsigned D>
struct FieldModulus;

template <>
struct FieldModulus<31> {
    static constexpr std::array<unsigned, 4> exponents{0, 1, 2, 3};
    static constexpr std::uint64_t low = 0x0f;
};
template <>
struct FieldModulus<32> {
    static constexpr std::array<unsigned, 4> exponents{0, 1, 3, 16};
    static constexpr std::uint64_t low = 0x1000b;
};
template <>
struct FieldModulus<46> {
    static constexpr std::array<unsigned, 4> exponents{0, 1, 2, 11};
    static constexpr std::uint64_t low = 0x807;
};
template <>
struct FieldModulus<47> {
    static constexpr std::array<unsigned, 4> exponents{0, 1, 3, 7};
    static constexpr std::uint64_t low = 0x8b;
};
template <>
struct FieldModulus<48> {
    static constexpr std::array<unsigned, 4> exponents{0, 1, 2, 17};
    static constexpr std::uint64_t low = 0x20007;
};
template <>
struct FieldModulus<64> {
    static constexpr std::array<unsigned, 4> exponents{0, 1, 2, 11};
    static constexpr std::uint64_t low = 0x807;
};

template <unsigned D>
inline Product128 multiply_by_modulus_low(std::uint64_t value) noexcept {
    constexpr auto exponents = FieldModulus<D>::exponents;
    Product128 result = shifted_word(value, exponents[0]);
    result = xor_product(result, shifted_word(value, exponents[1]));
    result = xor_product(result, shifted_word(value, exponents[2]));
    result = xor_product(result, shifted_word(value, exponents[3]));
    return result;
}

template <unsigned D>
inline std::uint64_t high_above_degree(Product128 value) noexcept {
    if constexpr (D == 64) {
        return value.high;
    } else {
        return (value.low >> D) | (value.high << (64 - D));
    }
}

template <unsigned D>
inline std::uint64_t field_multiply_scalar(
    std::uint64_t value, std::uint64_t coefficient) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    Product128 product = carryless_multiply(value & mask, coefficient & mask);

    // For p(x)=x^D+r(x), replace the high part H*x^D by H*r(x).
    const std::uint64_t first_high = high_above_degree<D>(product);
    Product128 reduced{product.low & mask, 0};
    reduced = xor_product(
        reduced, carryless_multiply(first_high, FieldModulus<D>::low));

    // Every selected r has degree at most 17.  One more fold therefore
    // completes reduction for all widths instantiated by this benchmark.
    const std::uint64_t second_high = high_above_degree<D>(reduced);
    Product128 final_value{reduced.low & mask, 0};
    final_value = xor_product(
        final_value, carryless_multiply(second_high, FieldModulus<D>::low));
    return final_value.low & mask;
}

template <unsigned D>
inline __m256i vector_lane_high(__m256i value) noexcept {
    if constexpr (D == 64) {
        return _mm256_srli_si256(value, 8);
    } else {
        const __m256i from_low = _mm256_srli_epi64(value, D);
        const __m256i from_high = _mm256_slli_epi64(
            _mm256_srli_si256(value, 8), 64 - D);
        return _mm256_or_si256(from_low, from_high);
    }
}

template <unsigned D>
inline __m256i field_multiply_vector(
    __m256i value, std::uint64_t coefficient) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    const __m256i lane_mask = _mm256_set_epi64x(
        0, static_cast<long long>(mask), 0, static_cast<long long>(mask));
    const __m256i coefficient_vector = _mm256_set_epi64x(
        0, static_cast<long long>(coefficient & mask),
        0, static_cast<long long>(coefficient & mask));
    const __m256i modulus_vector = _mm256_set_epi64x(
        0, static_cast<long long>(FieldModulus<D>::low),
        0, static_cast<long long>(FieldModulus<D>::low));
    __m256i product = _mm256_clmulepi64_epi128(
        _mm256_and_si256(value, lane_mask), coefficient_vector, 0x00);
    const __m256i first_high =
        _mm256_and_si256(vector_lane_high<D>(product), lane_mask);
    __m256i reduced = _mm256_xor_si256(
        _mm256_and_si256(product, lane_mask),
        _mm256_clmulepi64_epi128(first_high, modulus_vector, 0x00));
    const __m256i second_high =
        _mm256_and_si256(vector_lane_high<D>(reduced), lane_mask);
    reduced = _mm256_xor_si256(
        _mm256_and_si256(reduced, lane_mask),
        _mm256_clmulepi64_epi128(second_high, modulus_vector, 0x00));
    return _mm256_and_si256(reduced, lane_mask);
}

template <unsigned D>
inline std::uint64_t field_multiply_reference(
    std::uint64_t value, std::uint64_t coefficient) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    value &= mask;
    coefficient &= mask;
    std::uint64_t result = 0;
    for (unsigned bit = 0; bit < D; ++bit) {
        const std::uint64_t selected = std::uint64_t{0} - (coefficient & 1);
        result ^= value & selected;
        coefficient >>= 1;
        const bool carry = (value >> (D - 1)) != 0;
        value = (value << 1) & mask;
        value ^= FieldModulus<D>::low & (std::uint64_t{0} - carry);
    }
    return result;
}

template <unsigned D>
struct ToeplitzCoefficient {
    std::uint64_t reversed_low;
    std::uint64_t reversed_high;
};

template <unsigned D>
inline std::uint64_t toeplitz_transpose_scalar(
    std::uint64_t value, ToeplitzCoefficient<D> coefficient) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    const Product128 lower = carryless_multiply(value & mask, coefficient.reversed_low);
    std::uint64_t middle = lower.high;
    if constexpr (2 * D - 1 > 64) {
        const Product128 upper =
            carryless_multiply(value & mask, coefficient.reversed_high);
        middle ^= upper.low;
    }
    constexpr unsigned start = D - 1;
    std::uint64_t shifted;
    if constexpr (start == 0) {
        shifted = lower.low;
    } else {
        shifted = (lower.low >> start) | (middle << (64 - start));
    }
    return shifted & mask;
}

template <unsigned D>
inline __m256i toeplitz_transpose_vector(
    __m256i value, ToeplitzCoefficient<D> coefficient) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    const __m256i lane_mask = _mm256_set_epi64x(
        0, static_cast<long long>(mask), 0, static_cast<long long>(mask));
    const __m256i lower_coefficient = _mm256_set_epi64x(
        0, static_cast<long long>(coefficient.reversed_low),
        0, static_cast<long long>(coefficient.reversed_low));
    const __m256i input = _mm256_and_si256(value, lane_mask);
    const __m256i lower =
        _mm256_clmulepi64_epi128(input, lower_coefficient, 0x00);
    __m256i middle = _mm256_srli_si256(lower, 8);
    if constexpr (2 * D - 1 > 64) {
        const __m256i upper_coefficient = _mm256_set_epi64x(
            0, static_cast<long long>(coefficient.reversed_high),
            0, static_cast<long long>(coefficient.reversed_high));
        middle = _mm256_xor_si256(
            middle,
            _mm256_clmulepi64_epi128(input, upper_coefficient, 0x00));
    }
    constexpr int start = D - 1;
    __m256i shifted;
    if constexpr (start == 0) {
        shifted = lower;
    } else {
        shifted = _mm256_or_si256(
            _mm256_srli_epi64(lower, start),
            _mm256_slli_epi64(middle, 64 - start));
    }
    return _mm256_and_si256(shifted, lane_mask);
}

template <unsigned D>
inline std::uint64_t toeplitz_transpose_reference(
    std::uint64_t value, ToeplitzCoefficient<D> coefficient) noexcept {
    std::uint64_t output = 0;
    for (unsigned column = 0; column < D; ++column) {
        unsigned parity = 0;
        for (unsigned row = 0; row < D; ++row) {
            const unsigned kernel_index = D - 1 + column - row;
            const std::uint64_t word = kernel_index < 64
                ? coefficient.reversed_low
                : coefficient.reversed_high;
            parity ^= static_cast<unsigned>((value >> row) & 1) &
                static_cast<unsigned>((word >> (kernel_index & 63)) & 1);
        }
        output |= static_cast<std::uint64_t>(parity & 1) << column;
    }
    return output;
}

// In-place transpose of a 64 by 64 binary matrix.  The words are rows.
inline void transpose64(std::uint64_t* rows) noexcept {
    std::uint64_t mask = 0x00000000ffffffffULL;
    for (unsigned shift = 32; shift != 0; shift >>= 1) {
        for (unsigned base = 0; base < 64; base = (base + shift + 1) & ~shift) {
            const std::uint64_t difference =
                (rows[base] ^ (rows[base + shift] >> shift)) & mask;
            rows[base] ^= difference;
            rows[base + shift] ^= difference << shift;
        }
        mask ^= mask << (shift >> 1);
    }
}

inline std::uint64_t reverse_bits64(std::uint64_t value) noexcept {
    value = ((value & 0x5555555555555555ULL) << 1) |
            ((value >> 1) & 0x5555555555555555ULL);
    value = ((value & 0x3333333333333333ULL) << 2) |
            ((value >> 2) & 0x3333333333333333ULL);
    value = ((value & 0x0f0f0f0f0f0f0f0fULL) << 4) |
            ((value >> 4) & 0x0f0f0f0f0f0f0f0fULL);
    return _byteswap_uint64(value);
}

template <unsigned D>
constexpr unsigned four_russians_groups() noexcept {
    return (D + 3) / 4;
}

template <unsigned D>
using FourRussiansTable =
    std::array<Block, four_russians_groups<D>() * 16>;

template <unsigned D>
inline void build_four_russians_table(
    const Block* values, FourRussiansTable<D>& table) noexcept {
    static_assert(D > 0 && D <= 64);
    constexpr unsigned group_bits = 4;
    constexpr unsigned group_size = 1u << group_bits;
    constexpr unsigned groups = four_russians_groups<D>();

    // Four-Russians application of a binary D by D matrix to D 128-bit
    // elements.  Tables retain all inputs, so values may be overwritten.
    for (unsigned group = 0; group < groups; ++group) {
        Block* const entries = table.data() + group * group_size;
        entries[0] = zero_block();
        for (unsigned bit = 0; bit < group_bits; ++bit) {
            const unsigned input_index = group * group_bits + bit;
            entries[1u << bit] = input_index < D ? values[input_index] : zero_block();
        }
        entries[3] = xor_block(entries[1], entries[2]);
        entries[5] = xor_block(entries[1], entries[4]);
        entries[6] = xor_block(entries[2], entries[4]);
        entries[7] = xor_block(entries[3], entries[4]);
        entries[9] = xor_block(entries[1], entries[8]);
        entries[10] = xor_block(entries[2], entries[8]);
        entries[11] = xor_block(entries[3], entries[8]);
        entries[12] = xor_block(entries[4], entries[8]);
        entries[13] = xor_block(entries[5], entries[8]);
        entries[14] = xor_block(entries[6], entries[8]);
        entries[15] = xor_block(entries[7], entries[8]);
    }
}

template <unsigned D>
inline void apply_four_russians_table(
    Block* values, const std::uint64_t* rows,
    const FourRussiansTable<D>& table) noexcept {
    constexpr unsigned group_bits = 4;
    constexpr unsigned group_size = 1u << group_bits;
    constexpr unsigned groups = four_russians_groups<D>();

    unsigned row = 0;
    if constexpr (D >= 40) {
        for (; row + 4 <= D; row += 4) {
            std::uint64_t selections0 = rows[row];
            std::uint64_t selections1 = rows[row + 1];
            std::uint64_t selections2 = rows[row + 2];
            std::uint64_t selections3 = rows[row + 3];
            Block output0 = table[selections0 & (group_size - 1)];
            Block output1 = table[selections1 & (group_size - 1)];
            Block output2 = table[selections2 & (group_size - 1)];
            Block output3 = table[selections3 & (group_size - 1)];
            selections0 >>= group_bits;
            selections1 >>= group_bits;
            selections2 >>= group_bits;
            selections3 >>= group_bits;
            for (unsigned group = 1; group < groups; ++group) {
                const Block* const entries = table.data() + group * group_size;
                output0 = xor_block(
                    output0, entries[selections0 & (group_size - 1)]);
                output1 = xor_block(
                    output1, entries[selections1 & (group_size - 1)]);
                output2 = xor_block(
                    output2, entries[selections2 & (group_size - 1)]);
                output3 = xor_block(
                    output3, entries[selections3 & (group_size - 1)]);
                selections0 >>= group_bits;
                selections1 >>= group_bits;
                selections2 >>= group_bits;
                selections3 >>= group_bits;
            }
            values[row] = output0;
            values[row + 1] = output1;
            values[row + 2] = output2;
            values[row + 3] = output3;
        }
    }
    for (; row < D; ++row) {
        std::uint64_t selections = rows[row];
        Block output = table[selections & (group_size - 1)];
        selections >>= group_bits;
        for (unsigned group = 1; group < groups; ++group) {
            output = xor_block(
                output, table[group * group_size + (selections & (group_size - 1))]);
            selections >>= group_bits;
        }
        values[row] = output;
    }
}

template <unsigned D>
inline void apply_direct_rows(Block* values, const std::uint64_t* rows) noexcept {
    alignas(64) FourRussiansTable<D> table;
    build_four_russians_table<D>(values, table);
    apply_four_russians_table<D>(values, rows, table);
}

template <unsigned D>
inline void apply_direct_field(Block* values, std::uint64_t coefficient) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    alignas(64) std::array<std::uint64_t, 64> matrix{};
    std::uint64_t column = coefficient & mask;
    for (unsigned index = 0; index < D; ++index) {
        // transpose64 reverses both coordinate axes.
        matrix[63 - index] = column;
        const std::uint64_t carry = column >> (D - 1);
        column = ((column << 1) & mask) ^
            (FieldModulus<D>::low & (std::uint64_t{0} - carry));
    }
    transpose64(matrix.data());
    std::array<std::uint64_t, D> rows;
    for (unsigned row = 0; row < D; ++row) {
        rows[row] = matrix[63 - row] & mask;
    }
    apply_direct_rows<D>(values, rows.data());
}

template <unsigned D>
inline void build_toeplitz_rows(
    ToeplitzCoefficient<D> coefficient, std::uint64_t* rows) noexcept {
    constexpr std::uint64_t mask = width_mask<D>();
    std::uint64_t row_mask = reverse_bits64(coefficient.reversed_low) >> (64 - D);
    std::uint64_t tail_bits;
    if constexpr (D == 64) {
        tail_bits = coefficient.reversed_high;
    } else {
        tail_bits = (coefficient.reversed_low >> D) |
            (coefficient.reversed_high << (64 - D));
    }
    for (unsigned row = 0; row < D; ++row) {
        rows[row] = row_mask;
        if (row + 1 < D) {
            row_mask = ((row_mask << 1) & mask) | (tail_bits & 1);
            tail_bits >>= 1;
        }
    }
}

template <unsigned D>
inline void apply_direct_toeplitz(
    Block* values, ToeplitzCoefficient<D> coefficient) noexcept {
    std::array<std::uint64_t, D> rows;
    build_toeplitz_rows<D>(coefficient, rows.data());
    apply_direct_rows<D>(values, rows.data());
}

inline void xor_block_range(
    Block* destination, const Block* source, unsigned count) noexcept {
    unsigned index = 0;
    for (; index + 2 <= count; index += 2) {
        const __m256i left = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination + index));
        const __m256i right = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + index));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination + index),
            _mm256_xor_si256(left, right));
    }
    if (index < count) {
        destination[index] = xor_block(destination[index], source[index]);
    }
}

template <unsigned D>
inline void xor_toeplitz_diagonal(
    Block* output, const Block* input, unsigned kernel_index) noexcept {
    const int shift = static_cast<int>(kernel_index) - static_cast<int>(D - 1);
    const unsigned displacement = static_cast<unsigned>(shift < 0 ? -shift : shift);
    const unsigned count = D - displacement;
    const unsigned input_start = shift < 0 ? displacement : 0;
    const unsigned output_start = shift > 0 ? displacement : 0;
    xor_block_range(output + output_start, input + input_start, count);
}

template <unsigned D>
inline void apply_diagonal_toeplitz(
    Block* values, ToeplitzCoefficient<D> coefficient) noexcept {
    std::array<Block, D> output;
    output.fill(zero_block());

    std::uint64_t bits = coefficient.reversed_low;
    while (bits != 0) {
        const unsigned kernel_index = std::countr_zero(bits);
        xor_toeplitz_diagonal<D>(output.data(), values, kernel_index);
        bits &= bits - 1;
    }
    bits = coefficient.reversed_high;
    while (bits != 0) {
        const unsigned kernel_index = 64 + std::countr_zero(bits);
        if (kernel_index >= 2 * D - 1) {
            break;
        }
        xor_toeplitz_diagonal<D>(output.data(), values, kernel_index);
        bits &= bits - 1;
    }
    std::copy(output.begin(), output.end(), values);
}

template <unsigned N>
constexpr std::uint64_t low_bits_mask() noexcept {
    static_assert(N > 0 && N <= 64);
    if constexpr (N == 64) {
        return ~std::uint64_t{0};
    } else {
        return (std::uint64_t{1} << N) - 1;
    }
}

template <unsigned N>
inline void karatsuba_block_by_binary(
    Block* output, const Block* input, std::uint64_t coefficient) noexcept {
    static_assert(std::has_single_bit(N) && N <= 64);
    constexpr unsigned product_size = 2 * N - 1;
    std::fill_n(output, product_size, zero_block());

    if constexpr (N <= 32) {
        std::uint64_t bits = coefficient & low_bits_mask<N>();
        while (bits != 0) {
            const unsigned shift = std::countr_zero(bits);
            xor_block_range(output + shift, input, N);
            bits &= bits - 1;
        }
    } else {
        constexpr unsigned half = N / 2;
        constexpr unsigned half_product_size = 2 * half - 1;
        std::array<Block, half> input_sum;
        for (unsigned index = 0; index < half; ++index) {
            input_sum[index] = xor_block(input[index], input[index + half]);
        }

        const std::uint64_t low = coefficient & low_bits_mask<half>();
        const std::uint64_t high = coefficient >> half;
        alignas(32) std::array<Block, half_product_size> product_low;
        alignas(32) std::array<Block, half_product_size> product_sum;
        alignas(32) std::array<Block, half_product_size> product_high;
        karatsuba_block_by_binary<half>(product_low.data(), input, low);
        karatsuba_block_by_binary<half>(
            product_sum.data(), input_sum.data(), low ^ high);
        karatsuba_block_by_binary<half>(
            product_high.data(), input + half, high);

        for (unsigned index = 0; index < half_product_size; ++index) {
            output[index] = xor_block(output[index], product_low[index]);
            const Block middle = xor_block(
                xor_block(product_sum[index], product_low[index]),
                product_high[index]);
            output[half + index] = xor_block(output[half + index], middle);
            output[2 * half + index] =
                xor_block(output[2 * half + index], product_high[index]);
        }
    }
}

template <unsigned D>
constexpr unsigned karatsuba_width() noexcept {
    if constexpr (D <= 8) {
        return 8;
    } else if constexpr (D <= 16) {
        return 16;
    } else if constexpr (D <= 32) {
        return 32;
    } else {
        return 64;
    }
}

template <unsigned D>
inline void apply_karatsuba_toeplitz(
    Block* values, ToeplitzCoefficient<D> coefficient) noexcept {
    constexpr unsigned N = karatsuba_width<D>();
    constexpr unsigned product_size = 2 * N - 1;
    std::array<Block, N> padded_input;
    padded_input.fill(zero_block());
    std::copy_n(values, D, padded_input.begin());

    std::uint64_t coefficient_low;
    std::uint64_t coefficient_high;
    if constexpr (N == 64) {
        coefficient_low = coefficient.reversed_low;
        coefficient_high = coefficient.reversed_high;
    } else {
        coefficient_low = coefficient.reversed_low & low_bits_mask<N>();
        coefficient_high =
            (coefficient.reversed_low >> N) & low_bits_mask<N>();
    }

    alignas(32) std::array<Block, product_size> product_low;
    alignas(32) std::array<Block, product_size> product_high;
    karatsuba_block_by_binary<N>(
        product_low.data(), padded_input.data(), coefficient_low);
    karatsuba_block_by_binary<N>(
        product_high.data(), padded_input.data(), coefficient_high);

    for (unsigned row = 0; row < D; ++row) {
        const unsigned product_index = D - 1 + row;
        Block result = product_low[product_index];
        if (product_index >= N) {
            result = xor_block(result, product_high[product_index - N]);
        }
        values[row] = result;
    }
}

template <std::size_t BlockShift, std::size_t RowShift, std::size_t J = 0>
inline std::enable_if_t<J == (std::size_t{1} << BlockShift)>
avx_transpose_block_iter1(__m256i*) noexcept {}

template <std::size_t BlockShift, std::size_t RowShift, std::size_t J = 0>
inline std::enable_if_t<
    (J < (std::size_t{1} << BlockShift)) && (BlockShift > 0) &&
    (BlockShift < 6) && (RowShift >= 1)>
avx_transpose_block_iter1(__m256i* values) noexcept {
    avx_transpose_block_iter1<BlockShift, RowShift, J + (std::size_t{1} << RowShift)>(values);
    std::uint64_t mask = ~std::uint64_t{0} << 32;
    for (int bit = 4; bit >= static_cast<int>(BlockShift); --bit) {
        mask ^= mask >> (std::size_t{1} << bit);
    }
    __m256i& left = values[J / 2];
    __m256i& right = values[J / 2 + (std::size_t{1} << (BlockShift - 1))];
    if constexpr (BlockShift == 1) {
        __m256i low = _mm256_permute2x128_si256(left, right, 0x20);
        __m256i high = _mm256_permute2x128_si256(left, right, 0x31);
        __m256i difference = _mm256_xor_si256(low, _mm256_slli_epi16(high, 1));
        difference = _mm256_and_si256(difference, _mm256_set1_epi16(static_cast<short>(0xaaaa)));
        low = _mm256_xor_si256(low, difference);
        high = _mm256_xor_si256(high, _mm256_srli_epi16(difference, 1));
        left = _mm256_permute2x128_si256(low, high, 0x20);
        right = _mm256_permute2x128_si256(low, high, 0x31);
    }
    constexpr int shift = static_cast<int>(std::size_t{1} << BlockShift);
    __m256i difference = _mm256_xor_si256(left, _mm256_slli_epi64(right, shift));
    difference = _mm256_and_si256(
        difference, _mm256_set1_epi64x(static_cast<long long>(mask)));
    left = _mm256_xor_si256(left, difference);
    right = _mm256_xor_si256(right, _mm256_srli_epi64(difference, shift));
}

template <std::size_t BlockShift, std::size_t RowShift, std::size_t J = 0>
inline std::enable_if_t<(J < (std::size_t{1} << BlockShift)) && (BlockShift == 6)>
avx_transpose_block_iter1(__m256i* values) noexcept {
    avx_transpose_block_iter1<BlockShift, RowShift, J + (std::size_t{1} << RowShift)>(values);
    __m256i& left = values[J / 2];
    __m256i& right = values[J / 2 + (std::size_t{1} << (BlockShift - 1))];
    const __m256i output_left = _mm256_unpacklo_epi64(left, right);
    const __m256i output_right = _mm256_unpackhi_epi64(left, right);
    left = output_left;
    right = output_right;
}

template <std::size_t BlockShift, std::size_t RowShift, std::size_t Rows>
inline std::enable_if_t<Rows == 0>
avx_transpose_block_iter2(__m256i*) noexcept {}

template <std::size_t BlockShift, std::size_t RowShift, std::size_t Rows>
inline std::enable_if_t<(Rows > 0)>
avx_transpose_block_iter2(__m256i* values) noexcept {
    constexpr std::size_t matrix_size = std::size_t{1} << (BlockShift + 1);
    constexpr std::size_t offset = Rows - matrix_size;
    avx_transpose_block_iter2<BlockShift, RowShift, offset>(values);
    avx_transpose_block_iter1<BlockShift, RowShift>(values + offset / 2);
}

template <std::size_t BlockShift, std::size_t MatrixShift,
          std::size_t RowShift, std::size_t MatrixRowsShift>
inline std::enable_if_t<BlockShift == MatrixShift>
avx_transpose_block(__m256i*) noexcept {}

template <std::size_t BlockShift, std::size_t MatrixShift,
          std::size_t RowShift, std::size_t MatrixRowsShift>
inline std::enable_if_t<(BlockShift < MatrixShift)>
avx_transpose_block(__m256i* values) noexcept {
    avx_transpose_block_iter2<
        BlockShift, RowShift, (std::size_t{1} << (MatrixRowsShift + MatrixShift))>(values);
    avx_transpose_block<BlockShift + 1, MatrixShift, RowShift, MatrixRowsShift>(values);
}

constexpr std::size_t avx_block_shift = 4;
constexpr std::size_t avx_block_size = std::size_t{1} << avx_block_shift;

template <std::size_t Iteration = 7>
inline std::enable_if_t<(Iteration <= avx_block_shift + 1)>
avx_transpose(__m256i* values) noexcept {
    for (std::size_t offset = 0; offset < 64; offset += avx_block_size) {
        avx_transpose_block<1, Iteration, 1, avx_block_shift + 1 - Iteration>(
            values + offset);
    }
}

template <std::size_t Iteration = 7>
inline std::enable_if_t<(Iteration > avx_block_shift + 1)>
avx_transpose(__m256i* values) noexcept {
    avx_transpose<Iteration - avx_block_shift>(values);
    constexpr std::size_t block_shift = Iteration - avx_block_shift;
    constexpr std::size_t mask =
        (std::size_t{1} << (Iteration - 1)) - (std::size_t{1} << (block_shift - 1));
    if constexpr (Iteration == 7) {
        for (std::size_t index = 0; index < (std::size_t{1} << (block_shift - 1)); ++index) {
            avx_transpose_block<block_shift, Iteration, block_shift, 0>(values + index);
        }
    } else {
        for (std::size_t index = 0; index < 64; index = (index + mask + 1) & ~mask) {
            avx_transpose_block<block_shift, Iteration, block_shift, 0>(values + index);
        }
    }
}

inline void transpose128_avx(Block* values) noexcept {
    avx_transpose(reinterpret_cast<__m256i*>(values));
}

template <unsigned D>
inline void transform_field_square(Block* square, std::uint64_t coefficient) noexcept {
    transpose128_avx(square);
    auto* lanes = reinterpret_cast<__m256i*>(square);
    for (unsigned index = 0; index < 64; ++index) {
        lanes[index] = field_multiply_vector<D>(lanes[index], coefficient);
    }
    transpose128_avx(square);
}

template <unsigned D>
inline void transform_toeplitz_square(
    Block* square, ToeplitzCoefficient<D> coefficient) noexcept {
    transpose128_avx(square);
    auto* lanes = reinterpret_cast<__m256i*>(square);
    for (unsigned index = 0; index < 64; ++index) {
        lanes[index] = toeplitz_transpose_vector<D>(lanes[index], coefficient);
    }
    transpose128_avx(square);
}

template <unsigned D, typename ScalarMap>
inline void transform_bitsliced_square(Block* square, ScalarMap&& map) noexcept {
    static_assert(D <= 64);
    transpose128_avx(square);
    for (unsigned lane_index = 0; lane_index < 128; ++lane_index) {
        Block& lane = square[lane_index];
        const std::uint64_t input =
            static_cast<std::uint64_t>(_mm_cvtsi128_si64(lane.value));
        lane.value = _mm_cvtsi64_si128(static_cast<long long>(map(input)));
    }
    transpose128_avx(square);
}

template <unsigned D, typename ScalarMap>
inline void apply_bitsliced(Block* values, ScalarMap&& map) noexcept {
    alignas(64) std::array<Block, 128> square{};
    for (unsigned row = 0; row < D; ++row) {
        square[row] = values[row];
    }
    transform_bitsliced_square<D>(square.data(), std::forward<ScalarMap>(map));
    for (unsigned row = 0; row < D; ++row) {
        values[row] = square[row];
    }
}

template <unsigned D, typename ScalarMap>
std::array<Block, D> dense_reference(
    const std::array<Block, D>& input, ScalarMap&& scalar_map) {
    std::array<Block, D> output;
    output.fill(zero_block());
    for (unsigned column = 0; column < D; ++column) {
        const std::uint64_t image = scalar_map(std::uint64_t{1} << column);
        for (unsigned row = 0; row < D; ++row) {
            if ((image >> row) & 1) {
                output[row] = xor_block(output[row], input[column]);
            }
        }
    }
    return output;
}

std::uint64_t splitmix64(std::uint64_t& state) noexcept {
    std::uint64_t value = (state += 0x9e3779b97f4a7c15ULL);
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

Block random_block(std::uint64_t& state) noexcept {
    const std::uint64_t low = splitmix64(state);
    const std::uint64_t high = splitmix64(state);
    return {_mm_set_epi64x(static_cast<long long>(high), static_cast<long long>(low))};
}

template <unsigned D>
ToeplitzCoefficient<D> random_toeplitz(std::uint64_t& state) noexcept {
    constexpr unsigned bits = 2 * D - 1;
    ToeplitzCoefficient<D> result{splitmix64(state), splitmix64(state)};
    if constexpr (bits < 128) {
        constexpr unsigned high_bits = bits > 64 ? bits - 64 : 0;
        if constexpr (high_bits == 0) {
            result.reversed_high = 0;
            result.reversed_low &= width_mask<bits>();
        } else if constexpr (high_bits < 64) {
            result.reversed_high &= (std::uint64_t{1} << high_bits) - 1;
        }
    }
    return result;
}

template <unsigned D>
void test_width() {
    std::uint64_t random_state = 0x746573742d776964ULL ^ D;
    constexpr std::uint64_t mask = width_mask<D>();
    for (unsigned trial = 0; trial < 80; ++trial) {
        const std::uint64_t left = splitmix64(random_state) & mask;
        const std::uint64_t right = splitmix64(random_state) & mask;
        if (field_multiply_scalar<D>(left, right) !=
            field_multiply_reference<D>(left, right)) {
            throw std::runtime_error("optimized field reduction disagrees with reference");
        }

        std::array<Block, D> input;
        for (Block& value : input) {
            value = random_block(random_state);
        }
        const std::uint64_t coefficient = splitmix64(random_state) & mask;
        auto expected_field = dense_reference<D>(input, [coefficient](std::uint64_t x) {
            return field_multiply_reference<D>(x, coefficient);
        });
        auto actual_field = input;
        apply_bitsliced<D>(actual_field.data(), [coefficient](std::uint64_t x) {
            return field_multiply_scalar<D>(x, coefficient);
        });
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_field[index], actual_field[index])) {
                throw std::runtime_error("bitsliced field map disagrees with dense reference");
            }
        }
        auto direct_field = input;
        apply_direct_field<D>(direct_field.data(), coefficient);
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_field[index], direct_field[index])) {
                throw std::runtime_error("direct field map disagrees with dense reference");
            }
        }
        alignas(64) std::array<Block, 128> vector_field_square{};
        std::copy(input.begin(), input.end(), vector_field_square.begin());
        transform_field_square<D>(vector_field_square.data(), coefficient);
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_field[index], vector_field_square[index])) {
                throw std::runtime_error("vector field map disagrees with dense reference");
            }
        }

        const auto toeplitz = random_toeplitz<D>(random_state);
        const std::uint64_t toeplitz_input = splitmix64(random_state) & mask;
        if (toeplitz_transpose_scalar<D>(toeplitz_input, toeplitz) !=
            toeplitz_transpose_reference<D>(toeplitz_input, toeplitz)) {
            throw std::runtime_error("PCLMUL Toeplitz map disagrees with direct convolution");
        }
        auto expected_toeplitz = dense_reference<D>(input, [toeplitz](std::uint64_t x) {
            return toeplitz_transpose_scalar<D>(x, toeplitz);
        });
        auto actual_toeplitz = input;
        apply_bitsliced<D>(actual_toeplitz.data(), [toeplitz](std::uint64_t x) {
            return toeplitz_transpose_scalar<D>(x, toeplitz);
        });
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_toeplitz[index], actual_toeplitz[index])) {
                throw std::runtime_error("bitsliced Toeplitz map disagrees with dense reference");
            }
        }
        auto direct_toeplitz = input;
        apply_direct_toeplitz<D>(direct_toeplitz.data(), toeplitz);
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_toeplitz[index], direct_toeplitz[index])) {
                throw std::runtime_error("direct Toeplitz map disagrees with dense reference");
            }
        }
        auto diagonal_toeplitz = input;
        apply_diagonal_toeplitz<D>(diagonal_toeplitz.data(), toeplitz);
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_toeplitz[index], diagonal_toeplitz[index])) {
                throw std::runtime_error("diagonal Toeplitz map disagrees with dense reference");
            }
        }
        auto karatsuba_toeplitz = input;
        apply_karatsuba_toeplitz<D>(karatsuba_toeplitz.data(), toeplitz);
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_toeplitz[index], karatsuba_toeplitz[index])) {
                throw std::runtime_error("Karatsuba Toeplitz map disagrees with dense reference");
            }
        }
        alignas(64) std::array<Block, 128> vector_toeplitz_square{};
        std::copy(input.begin(), input.end(), vector_toeplitz_square.begin());
        transform_toeplitz_square<D>(vector_toeplitz_square.data(), toeplitz);
        for (unsigned index = 0; index < D; ++index) {
            if (!equal_block(expected_toeplitz[index], vector_toeplitz_square[index])) {
                throw std::runtime_error("vector Toeplitz map disagrees with dense reference");
            }
        }
    }
}

enum class Kernel { field, toeplitz };
enum class Implementation {
    four_russians,
    preexpanded_rows,
    diagonal,
    karatsuba,
    movement_floor,
    bitsliced,
};

template <unsigned D>
struct ToeplitzRows {
    std::array<std::uint64_t, D> values;
};

template <unsigned T, unsigned S>
void apply_field_chain_bitsliced(
    Block* word, std::size_t word_blocks, const std::uint64_t* coefficients) noexcept {
    constexpr unsigned D = T + S;
    constexpr std::uint64_t mask = width_mask<D>();
    std::array<Block, S> state;
    state.fill(zero_block());
    alignas(64) std::array<Block, 128> square{};
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, square.begin());
        std::copy_n(state.begin(), S, square.begin() + T);
        const std::uint64_t coefficient = coefficients[step] & mask;
        transform_field_square<D>(square.data(), coefficient);
        std::copy_n(square.begin(), T, position);
        std::copy_n(square.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void apply_toeplitz_chain_bitsliced(
    Block* word, std::size_t word_blocks,
    const ToeplitzCoefficient<T + S>* coefficients) noexcept {
    constexpr unsigned D = T + S;
    std::array<Block, S> state;
    state.fill(zero_block());
    alignas(64) std::array<Block, 128> square{};
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, square.begin());
        std::copy_n(state.begin(), S, square.begin() + T);
        const ToeplitzCoefficient<D> coefficient = coefficients[step];
        transform_toeplitz_square<D>(square.data(), coefficient);
        std::copy_n(square.begin(), T, position);
        std::copy_n(square.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void apply_field_chain_direct(
    Block* word, std::size_t word_blocks, const std::uint64_t* coefficients) noexcept {
    constexpr unsigned D = T + S;
    std::array<Block, S> state;
    state.fill(zero_block());
    std::array<Block, D> node;
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, node.begin());
        std::copy_n(state.begin(), S, node.begin() + T);
        apply_direct_field<D>(node.data(), coefficients[step]);
        std::copy_n(node.begin(), T, position);
        std::copy_n(node.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void apply_toeplitz_chain_direct(
    Block* word, std::size_t word_blocks,
    const ToeplitzCoefficient<T + S>* coefficients) noexcept {
    constexpr unsigned D = T + S;
    std::array<Block, S> state;
    state.fill(zero_block());
    std::array<Block, D> node;
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, node.begin());
        std::copy_n(state.begin(), S, node.begin() + T);
        apply_direct_toeplitz<D>(node.data(), coefficients[step]);
        std::copy_n(node.begin(), T, position);
        std::copy_n(node.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void apply_toeplitz_chain_preexpanded(
    Block* word, std::size_t word_blocks,
    const ToeplitzRows<T + S>* rows) noexcept {
    constexpr unsigned D = T + S;
    std::array<Block, S> state;
    state.fill(zero_block());
    std::array<Block, D> node;
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, node.begin());
        std::copy_n(state.begin(), S, node.begin() + T);
        apply_direct_rows<D>(node.data(), rows[step].values.data());
        std::copy_n(node.begin(), T, position);
        std::copy_n(node.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void apply_movement_floor_chain(Block* word, std::size_t word_blocks) noexcept {
    static_assert(S <= T);
    std::array<Block, S> state;
    state.fill(zero_block());
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        for (unsigned index = 0; index < S; ++index) {
            const Block input = position[index];
            position[index] = xor_block(input, state[index]);
            state[index] = input;
        }
        for (unsigned index = S; index < T; ++index) {
            position[index] = xor_block(position[index], state[index % S]);
        }
    }
}

template <unsigned T, unsigned S>
void apply_toeplitz_chain_diagonal(
    Block* word, std::size_t word_blocks,
    const ToeplitzCoefficient<T + S>* coefficients) noexcept {
    constexpr unsigned D = T + S;
    std::array<Block, S> state;
    state.fill(zero_block());
    std::array<Block, D> node;
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, node.begin());
        std::copy_n(state.begin(), S, node.begin() + T);
        apply_diagonal_toeplitz<D>(node.data(), coefficients[step]);
        std::copy_n(node.begin(), T, position);
        std::copy_n(node.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void apply_toeplitz_chain_karatsuba(
    Block* word, std::size_t word_blocks,
    const ToeplitzCoefficient<T + S>* coefficients) noexcept {
    constexpr unsigned D = T + S;
    std::array<Block, S> state;
    state.fill(zero_block());
    std::array<Block, D> node;
    const std::size_t steps = word_blocks / T;
    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        std::copy_n(position, T, node.begin());
        std::copy_n(state.begin(), S, node.begin() + T);
        apply_karatsuba_toeplitz<D>(node.data(), coefficients[step]);
        std::copy_n(node.begin(), T, position);
        std::copy_n(node.begin() + T, S, state.begin());
    }
}

template <unsigned T, unsigned S>
void test_reverse_chains() {
    constexpr unsigned D = T + S;
    constexpr std::size_t steps = 5;
    constexpr std::size_t word_blocks = T * steps;
    std::uint64_t random_state = 0x636861696e2d7465ULL ^ D;
    std::array<Block, word_blocks> source;
    for (Block& value : source) {
        value = random_block(random_state);
    }

    std::array<std::uint64_t, steps> field_coefficients;
    std::array<ToeplitzCoefficient<D>, steps> toeplitz_coefficients;
    for (std::size_t step = 0; step < steps; ++step) {
        field_coefficients[step] = splitmix64(random_state) & width_mask<D>();
        toeplitz_coefficients[step] = random_toeplitz<D>(random_state);
    }

    auto expected_field = source;
    std::array<Block, S> field_state;
    field_state.fill(zero_block());
    for (std::size_t step = steps; step-- > 0;) {
        std::array<Block, D> node;
        std::copy_n(expected_field.begin() + step * T, T, node.begin());
        std::copy_n(field_state.begin(), S, node.begin() + T);
        const std::uint64_t coefficient = field_coefficients[step];
        const auto transformed = dense_reference<D>(node, [coefficient](std::uint64_t x) {
            return field_multiply_reference<D>(x, coefficient);
        });
        std::copy_n(transformed.begin(), T, expected_field.begin() + step * T);
        std::copy_n(transformed.begin() + T, S, field_state.begin());
    }
    auto actual_field = source;
    apply_field_chain_bitsliced<T, S>(
        actual_field.data(), word_blocks, field_coefficients.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_field[index], actual_field[index])) {
            throw std::runtime_error("optimized reverse field chain disagrees with reference");
        }
    }
    auto direct_field = source;
    apply_field_chain_direct<T, S>(
        direct_field.data(), word_blocks, field_coefficients.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_field[index], direct_field[index])) {
            throw std::runtime_error("direct reverse field chain disagrees with reference");
        }
    }

    auto expected_toeplitz = source;
    std::array<Block, S> toeplitz_state;
    toeplitz_state.fill(zero_block());
    for (std::size_t step = steps; step-- > 0;) {
        std::array<Block, D> node;
        std::copy_n(expected_toeplitz.begin() + step * T, T, node.begin());
        std::copy_n(toeplitz_state.begin(), S, node.begin() + T);
        const auto coefficient = toeplitz_coefficients[step];
        const auto transformed = dense_reference<D>(node, [coefficient](std::uint64_t x) {
            return toeplitz_transpose_reference<D>(x, coefficient);
        });
        std::copy_n(transformed.begin(), T, expected_toeplitz.begin() + step * T);
        std::copy_n(transformed.begin() + T, S, toeplitz_state.begin());
    }
    auto actual_toeplitz = source;
    apply_toeplitz_chain_bitsliced<T, S>(
        actual_toeplitz.data(), word_blocks, toeplitz_coefficients.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_toeplitz[index], actual_toeplitz[index])) {
            throw std::runtime_error("optimized reverse Toeplitz chain disagrees with reference");
        }
    }
    auto direct_toeplitz = source;
    apply_toeplitz_chain_direct<T, S>(
        direct_toeplitz.data(), word_blocks, toeplitz_coefficients.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_toeplitz[index], direct_toeplitz[index])) {
            throw std::runtime_error("direct reverse Toeplitz chain disagrees with reference");
        }
    }
    std::array<ToeplitzRows<D>, steps> expanded_rows;
    for (std::size_t step = 0; step < steps; ++step) {
        build_toeplitz_rows<D>(
            toeplitz_coefficients[step], expanded_rows[step].values.data());
    }
    auto preexpanded_toeplitz = source;
    apply_toeplitz_chain_preexpanded<T, S>(
        preexpanded_toeplitz.data(), word_blocks, expanded_rows.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_toeplitz[index], preexpanded_toeplitz[index])) {
            throw std::runtime_error(
                "preexpanded-row Toeplitz chain disagrees with reference");
        }
    }
    auto diagonal_toeplitz = source;
    apply_toeplitz_chain_diagonal<T, S>(
        diagonal_toeplitz.data(), word_blocks, toeplitz_coefficients.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_toeplitz[index], diagonal_toeplitz[index])) {
            throw std::runtime_error(
                "diagonal reverse Toeplitz chain disagrees with reference");
        }
    }
    auto karatsuba_toeplitz = source;
    apply_toeplitz_chain_karatsuba<T, S>(
        karatsuba_toeplitz.data(), word_blocks, toeplitz_coefficients.data());
    for (std::size_t index = 0; index < word_blocks; ++index) {
        if (!equal_block(expected_toeplitz[index], karatsuba_toeplitz[index])) {
            throw std::runtime_error(
                "Karatsuba reverse Toeplitz chain disagrees with reference");
        }
    }
}

std::uint64_t checksum(const std::vector<Block>& values) noexcept {
    std::uint64_t result = 0;
    for (std::size_t index = 0; index < values.size(); index += 4093) {
        result ^= static_cast<std::uint64_t>(_mm_cvtsi128_si64(values[index].value));
        result = std::rotl(result, 17);
        result ^= static_cast<std::uint64_t>(_mm_extract_epi64(values[index].value, 1));
    }
    return result;
}

inline std::uint64_t cycle_clock_begin() noexcept {
    _mm_lfence();
    const std::uint64_t value = __rdtsc();
    _mm_lfence();
    return value;
}

inline std::uint64_t cycle_clock_end() noexcept {
    unsigned auxiliary = 0;
    const std::uint64_t value = __rdtscp(&auxiliary);
    _mm_lfence();
    return value;
}

inline std::uint64_t cycle_clock_overhead() noexcept {
    std::uint64_t result = ~std::uint64_t{0};
    for (unsigned trial = 0; trial < 1000; ++trial) {
        const std::uint64_t begin = cycle_clock_begin();
        const std::uint64_t end = cycle_clock_end();
        result = (std::min)(result, end - begin);
    }
    return result;
}

struct PhaseCycles {
    std::uint64_t assemble = 0;
    std::uint64_t rows = 0;
    std::uint64_t table = 0;
    std::uint64_t apply = 0;
    std::uint64_t commit = 0;
    std::uint64_t samples = 0;
};

inline std::uint64_t adjusted_cycles(
    std::uint64_t begin, std::uint64_t end, std::uint64_t overhead) noexcept {
    const std::uint64_t elapsed = end - begin;
    return elapsed > overhead ? elapsed - overhead : 0;
}

template <unsigned T, unsigned S>
PhaseCycles profile_four_russians_chain(
    Block* word, std::size_t word_blocks,
    const ToeplitzCoefficient<T + S>* coefficients,
    std::uint64_t overhead) noexcept {
    constexpr unsigned D = T + S;
    constexpr std::size_t sample_mask = 255;
    std::array<Block, S> state;
    state.fill(zero_block());
    std::array<Block, D> node;
    std::array<std::uint64_t, D> rows;
    alignas(64) FourRussiansTable<D> table;
    PhaseCycles result;
    const std::size_t steps = word_blocks / T;

    for (std::size_t step = steps; step-- > 0;) {
        Block* const position = word + step * T;
        if ((step & sample_mask) == 0) {
            ++result.samples;
            std::uint64_t begin = cycle_clock_begin();
            std::copy_n(position, T, node.begin());
            std::copy_n(state.begin(), S, node.begin() + T);
            std::uint64_t end = cycle_clock_end();
            result.assemble += adjusted_cycles(begin, end, overhead);

            begin = cycle_clock_begin();
            build_toeplitz_rows<D>(coefficients[step], rows.data());
            end = cycle_clock_end();
            result.rows += adjusted_cycles(begin, end, overhead);

            begin = cycle_clock_begin();
            build_four_russians_table<D>(node.data(), table);
            end = cycle_clock_end();
            result.table += adjusted_cycles(begin, end, overhead);

            begin = cycle_clock_begin();
            apply_four_russians_table<D>(node.data(), rows.data(), table);
            end = cycle_clock_end();
            result.apply += adjusted_cycles(begin, end, overhead);

            begin = cycle_clock_begin();
            std::copy_n(node.begin(), T, position);
            std::copy_n(node.begin() + T, S, state.begin());
            end = cycle_clock_end();
            result.commit += adjusted_cycles(begin, end, overhead);
        } else {
            std::copy_n(position, T, node.begin());
            std::copy_n(state.begin(), S, node.begin() + T);
            build_toeplitz_rows<D>(coefficients[step], rows.data());
            build_four_russians_table<D>(node.data(), table);
            apply_four_russians_table<D>(node.data(), rows.data(), table);
            std::copy_n(node.begin(), T, position);
            std::copy_n(node.begin() + T, S, state.begin());
        }
    }
    return result;
}

template <unsigned T, unsigned S>
void profile_kernel_phases(std::size_t word_blocks, std::uint64_t seed) {
    constexpr unsigned D = T + S;
    if (word_blocks % T) {
        throw std::invalid_argument("word length must be divisible by t");
    }
    const std::size_t steps = word_blocks / T;
    std::vector<Block> source(word_blocks);
    std::vector<Block> working(word_blocks);
    std::vector<ToeplitzCoefficient<D>> coefficients(steps);
    std::uint64_t random_state = seed;
    for (Block& value : source) {
        value = random_block(random_state);
    }
    for (auto& coefficient : coefficients) {
        coefficient = random_toeplitz<D>(random_state);
    }

    std::memcpy(working.data(), source.data(), source.size() * sizeof(Block));
    apply_toeplitz_chain_direct<T, S>(
        working.data(), word_blocks, coefficients.data());
    std::memcpy(working.data(), source.data(), source.size() * sizeof(Block));
    const std::uint64_t overhead = cycle_clock_overhead();
    const PhaseCycles cycles = profile_four_russians_chain<T, S>(
        working.data(), word_blocks, coefficients.data(), overhead);
    const double denominator = static_cast<double>(cycles.samples);
    const double movement = static_cast<double>(cycles.assemble + cycles.commit) /
        denominator;
    const double total = movement + static_cast<double>(
        cycles.rows + cycles.table + cycles.apply) / denominator;

    std::cout << std::fixed << std::setprecision(3)
              << "{\n"
              << "  \"schema\": \"riffle-four-russians-phase-profile-v1\",\n"
              << "  \"step_bits\": " << T << ",\n"
              << "  \"state_bits\": " << S << ",\n"
              << "  \"map_width_bits\": " << D << ",\n"
              << "  \"samples\": " << cycles.samples << ",\n"
              << "  \"sample_period_steps\": 256,\n"
              << "  \"timestamp_overhead_ticks\": " << overhead << ",\n"
              << "  \"movement_ticks_per_step\": " << movement << ",\n"
              << "  \"row_generation_ticks_per_step\": "
              << static_cast<double>(cycles.rows) / denominator << ",\n"
              << "  \"table_construction_ticks_per_step\": "
              << static_cast<double>(cycles.table) / denominator << ",\n"
              << "  \"table_application_ticks_per_step\": "
              << static_cast<double>(cycles.apply) / denominator << ",\n"
              << "  \"accounted_ticks_per_step\": " << total << ",\n"
              << "  \"checksum\": \"0x" << std::hex << checksum(working)
              << std::dec << "\"\n"
              << "}\n";
}

template <unsigned T, unsigned S>
void benchmark_kernel(
    Kernel kernel, Implementation implementation, std::size_t word_blocks,
    unsigned trials, std::uint64_t seed) {
    constexpr unsigned D = T + S;
    if (word_blocks % T) {
        throw std::invalid_argument("word length must be divisible by t");
    }
    if (kernel == Kernel::field &&
        implementation != Implementation::four_russians &&
        implementation != Implementation::movement_floor &&
        implementation != Implementation::bitsliced) {
        throw std::invalid_argument(
            "diagonal and Karatsuba implementations apply only to Toeplitz");
    }
    const std::size_t steps = word_blocks / T;
    std::vector<Block> source(word_blocks);
    std::vector<Block> working(word_blocks);
    std::uint64_t random_state = seed;
    for (Block& value : source) {
        value = random_block(random_state);
    }
    std::vector<std::uint64_t> field_coefficients;
    std::vector<ToeplitzCoefficient<D>> toeplitz_coefficients;
    std::vector<ToeplitzRows<D>> preexpanded_rows;
    if (kernel == Kernel::field) {
        field_coefficients.resize(steps);
        for (std::uint64_t& coefficient : field_coefficients) {
            coefficient = splitmix64(random_state) & width_mask<D>();
        }
    } else {
        toeplitz_coefficients.resize(steps);
        for (auto& coefficient : toeplitz_coefficients) {
            coefficient = random_toeplitz<D>(random_state);
        }
        if (implementation == Implementation::preexpanded_rows) {
            preexpanded_rows.resize(steps);
            for (std::size_t step = 0; step < steps; ++step) {
                build_toeplitz_rows<D>(
                    toeplitz_coefficients[step],
                    preexpanded_rows[step].values.data());
            }
        }
    }

    auto invoke = [&] {
        if (implementation == Implementation::movement_floor) {
            apply_movement_floor_chain<T, S>(working.data(), word_blocks);
            return;
        }
        if (kernel == Kernel::field) {
            if (implementation == Implementation::four_russians) {
                apply_field_chain_direct<T, S>(
                    working.data(), word_blocks, field_coefficients.data());
            } else {
                apply_field_chain_bitsliced<T, S>(
                    working.data(), word_blocks, field_coefficients.data());
            }
        } else {
            if (implementation == Implementation::four_russians) {
                apply_toeplitz_chain_direct<T, S>(
                    working.data(), word_blocks, toeplitz_coefficients.data());
            } else if (implementation == Implementation::preexpanded_rows) {
                apply_toeplitz_chain_preexpanded<T, S>(
                    working.data(), word_blocks, preexpanded_rows.data());
            } else if (implementation == Implementation::diagonal) {
                apply_toeplitz_chain_diagonal<T, S>(
                    working.data(), word_blocks, toeplitz_coefficients.data());
            } else if (implementation == Implementation::karatsuba) {
                apply_toeplitz_chain_karatsuba<T, S>(
                    working.data(), word_blocks, toeplitz_coefficients.data());
            } else {
                apply_toeplitz_chain_bitsliced<T, S>(
                    working.data(), word_blocks, toeplitz_coefficients.data());
            }
        }
    };

    std::memcpy(working.data(), source.data(), source.size() * sizeof(Block));
    invoke();

    std::vector<double> milliseconds;
    milliseconds.reserve(trials);
    std::uint64_t final_checksum = 0;
    for (unsigned trial = 0; trial < trials; ++trial) {
        std::memcpy(working.data(), source.data(), source.size() * sizeof(Block));
        const auto begin = std::chrono::steady_clock::now();
        invoke();
        const auto end = std::chrono::steady_clock::now();
        milliseconds.push_back(
            std::chrono::duration<double, std::milli>(end - begin).count());
        final_checksum ^= std::rotl(checksum(working), trial & 63);
    }
    std::sort(milliseconds.begin(), milliseconds.end());
    const double median = milliseconds[milliseconds.size() / 2];
    const double minimum = milliseconds.front();
    const double gib = static_cast<double>(word_blocks * sizeof(Block)) /
        static_cast<double>(std::uint64_t{1} << 30);

    std::cout << std::fixed << std::setprecision(6)
              << "{\n"
              << "  \"schema\": \"riffle-transposed-structured-step-benchmark-v1\",\n"
              << "  \"kernel\": \"" << (kernel == Kernel::field ? "field" : "toeplitz") << "\",\n"
              << "  \"implementation\": \""
              << (implementation == Implementation::four_russians
                    ? "block_four_russians_compact_coefficients"
                    : implementation == Implementation::preexpanded_rows
                        ? "block_four_russians_preexpanded_rows"
                        : implementation == Implementation::diagonal
                        ? "block_avx2_nonzero_diagonals"
                        : implementation == Implementation::karatsuba
                            ? "block_compile_time_karatsuba"
                            : implementation == Implementation::movement_floor
                                ? "block_recurrent_movement_floor"
                            : "avx2_transpose128_vpclmulqdq_two_lanes")
              << "\",\n"
              << "  \"coefficient_bits_per_step\": "
              << (kernel == Kernel::field ? D : 2 * D - 1) << ",\n"
              << "  \"step_bits\": " << T << ",\n"
              << "  \"state_bits\": " << S << ",\n"
              << "  \"map_width_bits\": " << D << ",\n"
              << "  \"word_blocks\": " << word_blocks << ",\n"
              << "  \"inner_steps\": " << steps << ",\n"
              << "  \"trials\": " << trials << ",\n"
              << "  \"median_ms\": " << median << ",\n"
              << "  \"minimum_ms\": " << minimum << ",\n"
              << "  \"median_input_gib_per_second\": " << gib / (median / 1000.0) << ",\n"
              << "  \"checksum\": \"0x" << std::hex << final_checksum << std::dec << "\",\n"
              << "  \"scope\": \"Reverse-order in-place inner chain on 128-bit elements; setup and input reset are outside timing.\"\n"
              << "}\n";
}

struct Arguments {
    bool self_test = false;
    bool profile_phases = false;
    Kernel kernel = Kernel::field;
    Implementation implementation = Implementation::four_russians;
    unsigned step_bits = 16;
    unsigned state_bits = 16;
    std::size_t word_blocks = std::size_t{1} << 21;
    unsigned trials = 11;
    std::uint64_t seed = 0x726966666c652d31ULL;
};

Arguments parse_arguments(int argc, char** argv) {
    Arguments result;
    for (int index = 1; index < argc; ++index) {
        const std::string option = argv[index];
        auto value = [&]() -> std::string {
            if (++index >= argc) {
                throw std::invalid_argument("missing value after " + option);
            }
            return argv[index];
        };
        if (option == "--self-test") {
            result.self_test = true;
        } else if (option == "--profile-phases") {
            result.profile_phases = true;
        } else if (option == "--kernel") {
            const std::string name = value();
            if (name == "field") {
                result.kernel = Kernel::field;
            } else if (name == "toeplitz") {
                result.kernel = Kernel::toeplitz;
            } else {
                throw std::invalid_argument("kernel must be field or toeplitz");
            }
        } else if (option == "--implementation") {
            const std::string name = value();
            if (name == "direct" || name == "four-russians") {
                result.implementation = Implementation::four_russians;
            } else if (name == "preexpanded-rows") {
                result.implementation = Implementation::preexpanded_rows;
            } else if (name == "diagonal") {
                result.implementation = Implementation::diagonal;
            } else if (name == "karatsuba") {
                result.implementation = Implementation::karatsuba;
            } else if (name == "movement-floor") {
                result.implementation = Implementation::movement_floor;
            } else if (name == "bitsliced") {
                result.implementation = Implementation::bitsliced;
            } else {
                throw std::invalid_argument(
                    "implementation must be four-russians, preexpanded-rows, diagonal, karatsuba, movement-floor, or bitsliced");
            }
        } else if (option == "--t") {
            result.step_bits = static_cast<unsigned>(std::stoul(value()));
        } else if (option == "--s") {
            result.state_bits = static_cast<unsigned>(std::stoul(value()));
        } else if (option == "--word-blocks") {
            result.word_blocks = static_cast<std::size_t>(std::stoull(value()));
        } else if (option == "--trials") {
            result.trials = static_cast<unsigned>(std::stoul(value()));
        } else if (option == "--seed") {
            result.seed = std::stoull(value(), nullptr, 0);
        } else {
            throw std::invalid_argument("unknown option " + option);
        }
    }
    if (result.trials == 0) {
        throw std::invalid_argument("trials must be positive");
    }
    return result;
}

template <unsigned T, unsigned S>
void dispatch_if(const Arguments& arguments) {
    if (arguments.profile_phases) {
        profile_kernel_phases<T, S>(arguments.word_blocks, arguments.seed);
    } else {
        benchmark_kernel<T, S>(
            arguments.kernel, arguments.implementation, arguments.word_blocks,
            arguments.trials, arguments.seed);
    }
}

void dispatch(const Arguments& arguments) {
    if (arguments.step_bits == 16 && arguments.state_bits == 15) {
        dispatch_if<16, 15>(arguments);
    } else if (arguments.step_bits == 16 && arguments.state_bits == 16) {
        dispatch_if<16, 16>(arguments);
    } else if (arguments.step_bits == 32 && arguments.state_bits == 14) {
        dispatch_if<32, 14>(arguments);
    } else if (arguments.step_bits == 32 && arguments.state_bits == 15) {
        dispatch_if<32, 15>(arguments);
    } else if (arguments.step_bits == 32 && arguments.state_bits == 16) {
        dispatch_if<32, 16>(arguments);
    } else if (arguments.step_bits == 32 && arguments.state_bits == 32) {
        dispatch_if<32, 32>(arguments);
    } else {
        throw std::invalid_argument(
            "supported (t,s): (16,15), (16,16), (32,14), (32,15), (32,16), (32,32)");
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Arguments arguments = parse_arguments(argc, argv);
        int cpu_features[4]{};
        __cpuidex(cpu_features, 7, 0);
        if ((arguments.self_test ||
             arguments.implementation == Implementation::bitsliced) &&
            (cpu_features[2] & (1 << 10)) == 0) {
            throw std::runtime_error("the optimized kernel requires VPCLMULQDQ");
        }
        SetThreadAffinityMask(GetCurrentThread(), 1);
        SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST);
        if (arguments.self_test) {
            test_width<31>();
            test_width<32>();
            test_width<46>();
            test_width<47>();
            test_width<48>();
            test_width<64>();
            test_reverse_chains<16, 15>();
            test_reverse_chains<16, 16>();
            test_reverse_chains<32, 14>();
            test_reverse_chains<32, 15>();
            test_reverse_chains<32, 16>();
            test_reverse_chains<32, 32>();
            std::cout << "self_test,pass\n";
            return 0;
        }
        dispatch(arguments);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error," << error.what() << '\n';
        return 1;
    }
}
