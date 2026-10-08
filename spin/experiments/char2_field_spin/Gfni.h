#pragma once
#include "SmallField.h"
#include <immintrin.h>
#include <cstring>

#ifndef SPIN_CHAR2_CLMUL_REDUCTION
#define SPIN_CHAR2_CLMUL_REDUCTION 0
#endif
#ifndef SPIN_CHAR2_SCALAR_LOOKUP
#define SPIN_CHAR2_SCALAR_LOOKUP 0
#endif

namespace spin::research::char2 {

inline __m128i smallmul128(__m128i input, Coeff128 coefficient) noexcept {
    const auto low = _mm_gf2p8affine_epi64_epi8(input,
        _mm_set1_epi64x(static_cast<long long>(coefficient.low)), 0);
    const auto high = _mm_gf2p8affine_epi64_epi8(input,
        _mm_set1_epi64x(static_cast<long long>(coefficient.high)), 0);
    const auto overflow = _mm_srli_si128(high, 15);
#if SPIN_CHAR2_CLMUL_REDUCTION
    // overflow has at most seven bits, so this unreduced product fits 14 bits.
    const auto fold = _mm_clmulepi64_si128(overflow, _mm_set1_epi64x(0x87), 0x00);
    return _mm_xor_si128(_mm_xor_si128(low, _mm_slli_si128(high, 1)), fold);
#else
    const auto reduceLow = _mm_gf2p8affine_epi64_epi8(overflow,
        _mm_set1_epi64x(static_cast<long long>(affineMatrix(0x87, false))), 0);
    const auto reduceHigh = _mm_gf2p8affine_epi64_epi8(overflow,
        _mm_set1_epi64x(static_cast<long long>(affineMatrix(0x87, true))), 0);
    return _mm_xor_si128(_mm_xor_si128(low, _mm_slli_si128(high, 1)),
        _mm_xor_si128(reduceLow, _mm_slli_si128(reduceHigh, 1)));
#endif
}

inline __m128i smallmul128(__m128i input, std::uint8_t coefficient) noexcept {
    return smallmul128(input, coefficient128(coefficient));
}

inline __m512i smallmul4(__m512i input, __m512i lowMatrix,
                        __m512i highMatrix) noexcept {
    const auto low = _mm512_gf2p8affine_epi64_epi8(input, lowMatrix, 0);
    const auto high = _mm512_gf2p8affine_epi64_epi8(input, highMatrix, 0);
    const auto overflow = _mm512_bsrli_epi128(high, 15);
#if SPIN_CHAR2_CLMUL_REDUCTION
    const auto fold = _mm512_clmulepi64_epi128(overflow, _mm512_set1_epi64(0x87), 0x00);
    return _mm512_xor_si512(_mm512_xor_si512(low, _mm512_bslli_epi128(high, 1)), fold);
#else
    const auto reduceLow = _mm512_gf2p8affine_epi64_epi8(overflow,
        _mm512_set1_epi64(static_cast<long long>(affineMatrix(0x87, false))), 0);
    const auto reduceHigh = _mm512_gf2p8affine_epi64_epi8(overflow,
        _mm512_set1_epi64(static_cast<long long>(affineMatrix(0x87, true))), 0);
    return _mm512_xor_si512(_mm512_xor_si512(low, _mm512_bslli_epi128(high, 1)),
        _mm512_xor_si512(reduceLow, _mm512_bslli_epi128(reduceHigh, 1)));
#endif
}

inline __m512i smallmul4(__m512i input, const Coeff4& coefficients) noexcept {
    return smallmul4(input, _mm512_load_si512(coefficients.low.data()),
                    _mm512_load_si512(coefficients.high.data()));
}

inline __m512i smallmul4(__m512i input, Coeff128 coefficient) noexcept {
    return smallmul4(input,
        _mm512_set1_epi64(static_cast<long long>(coefficient.low)),
        _mm512_set1_epi64(static_cast<long long>(coefficient.high)));
}

inline __m512i smallmul4(__m512i input, std::uint8_t coefficient) noexcept {
    return smallmul4(input, coefficient128(coefficient));
}

inline __m512i smallmul4(__m512i input, const std::uint8_t* coefficients) noexcept {
    std::uint32_t packed;
    std::memcpy(&packed, coefficients, sizeof(packed));
#if SPIN_CHAR2_SCALAR_LOOKUP
    const auto c0 = coefficient128(std::uint8_t(packed));
    const auto c1 = coefficient128(std::uint8_t(packed >> 8));
    const auto c2 = coefficient128(std::uint8_t(packed >> 16));
    const auto c3 = coefficient128(std::uint8_t(packed >> 24));
    return smallmul4(input,
        _mm512_setr_epi64(static_cast<long long>(c0.low), static_cast<long long>(c0.low),
                         static_cast<long long>(c1.low), static_cast<long long>(c1.low),
                         static_cast<long long>(c2.low), static_cast<long long>(c2.low),
                         static_cast<long long>(c3.low), static_cast<long long>(c3.low)),
        _mm512_setr_epi64(static_cast<long long>(c0.high), static_cast<long long>(c0.high),
                         static_cast<long long>(c1.high), static_cast<long long>(c1.high),
                         static_cast<long long>(c2.high), static_cast<long long>(c2.high),
                         static_cast<long long>(c3.high), static_cast<long long>(c3.high)));
#else
    auto indices = _mm256_cvtepu8_epi32(_mm_cvtsi32_si128(static_cast<int>(packed)));
    indices = _mm256_permutevar8x32_epi32(indices, _mm256_setr_epi32(0,0,1,1,2,2,3,3));
    return smallmul4(input, _mm512_i32gather_epi64(indices, affineTable.low.data(), 8),
                    _mm512_i32gather_epi64(indices, affineTable.high.data(), 8));
#endif
}

// Eight independently weighted GF(2^64) records, using the polynomial x^64+
// x^4+x^3+x+1. The high product byte crosses all byte boundaries in a record.
inline __m512i smallmul8x64(__m512i input, __m512i lowMatrix,
                           __m512i highMatrix) noexcept {
    const auto low = _mm512_gf2p8affine_epi64_epi8(input, lowMatrix, 0);
    const auto high = _mm512_gf2p8affine_epi64_epi8(input, highMatrix, 0);
    const auto overflow = _mm512_srli_epi64(high, 56);
    const auto reduceLow = _mm512_gf2p8affine_epi64_epi8(overflow,
        _mm512_set1_epi64(static_cast<long long>(affineMatrix(0x1b, false))), 0);
    const auto reduceHigh = _mm512_gf2p8affine_epi64_epi8(overflow,
        _mm512_set1_epi64(static_cast<long long>(affineMatrix(0x1b, true))), 0);
    return _mm512_xor_si512(_mm512_xor_si512(low, _mm512_slli_epi64(high, 8)),
        _mm512_xor_si512(reduceLow, _mm512_slli_epi64(reduceHigh, 8)));
}

// Eight GF(2^32) records in low dwords of separate 64-bit lanes. This deliberately
// avoids pretending a GFNI matrix can vary between the two dwords of one lane.
inline __m512i smallmul8x32Expanded(__m512i input, __m512i lowMatrix,
                                   __m512i highMatrix) noexcept {
    const auto low = _mm512_gf2p8affine_epi64_epi8(input, lowMatrix, 0);
    const auto high = _mm512_gf2p8affine_epi64_epi8(input, highMatrix, 0);
    const auto overflow = _mm512_srli_epi64(high, 24);
    const auto reduceLow = _mm512_gf2p8affine_epi64_epi8(overflow,
        _mm512_set1_epi64(static_cast<long long>(affineMatrix(0x8d, false))), 0);
    const auto reduceHigh = _mm512_gf2p8affine_epi64_epi8(overflow,
        _mm512_set1_epi64(static_cast<long long>(affineMatrix(0x8d, true))), 0);
    const auto result = _mm512_xor_si512(_mm512_xor_si512(low, _mm512_slli_epi64(high, 8)),
        _mm512_xor_si512(reduceLow, _mm512_slli_epi64(reduceHigh, 8)));
    return _mm512_and_si512(result, _mm512_set1_epi64(0xffffffffULL));
}

} // namespace spin::research::char2
