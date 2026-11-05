#pragma once
#include "Field128.h"
#include <immintrin.h>

namespace spin::research::char2::field128_lift {
#if defined(_MSC_VER)
#define FIELD128_LIFT_INLINE __forceinline
#else
#define FIELD128_LIFT_INLINE inline __attribute__((always_inline))
#endif

// Four independent native polynomial elements in the four 128-bit lanes.
// Requires AVX512F/BW and VPCLMULQDQ. No allocation, coefficient expansion,
// payload transpose, or change of basis is hidden in this arithmetic kernel.
FIELD128_LIFT_INLINE __m512i reduce4(__m512i low,__m512i high) noexcept {
    const auto polynomial=_mm512_set1_epi64(0x87);
    // If H=h0+X^64 h1 and h1*r=a0+X^64 a1, fold as
    // L+(h0+a1)*r+X^64*a0, r=X^7+X^2+X+1.
    const auto a=_mm512_clmulepi64_epi128(high,polynomial,0x01);
    const auto q=_mm512_clmulepi64_epi128(
        _mm512_xor_si512(high,_mm512_bsrli_epi128(a,8)),polynomial,0x00);
    return _mm512_ternarylogic_epi64(low,q,_mm512_bslli_epi128(a,8),0x96);
}

FIELD128_LIFT_INLINE __m512i multiply4(__m512i value,__m512i coefficient) noexcept {
    const auto p0=_mm512_clmulepi64_epi128(value,coefficient,0x00);
    const auto p1=_mm512_clmulepi64_epi128(value,coefficient,0x11);
    const auto a=_mm512_xor_si512(value,_mm512_shuffle_epi32(value,static_cast<_MM_PERM_ENUM>(0x4e)));
    const auto b=_mm512_xor_si512(coefficient,_mm512_shuffle_epi32(coefficient,static_cast<_MM_PERM_ENUM>(0x4e)));
    const auto p2=_mm512_clmulepi64_epi128(a,b,0x00);
    const auto middle=_mm512_ternarylogic_epi64(p0,p1,p2,0x96);
    return reduce4(_mm512_xor_si512(p0,_mm512_bslli_epi128(middle,8)),
                   _mm512_xor_si512(p1,_mm512_bsrli_epi128(middle,8)));
}

FIELD128_LIFT_INLINE __m512i multiply4(__m512i value,const Scalar* coefficient) noexcept {
    return multiply4(value,_mm512_loadu_si512(coefficient));
}

// Exact four-record loads/stores. Input and output may be identical; each
// lane has its own full 128-bit coefficient. No masked tail or overread.
FIELD128_LIFT_INLINE void multiply4Records(const Scalar* input,const Scalar* coefficients,
                                         Scalar* output) noexcept {
    const auto value=_mm512_loadu_si512(input);
    const auto coefficient=_mm512_loadu_si512(coefficients);
    _mm512_storeu_si512(output,multiply4(value,coefficient));
}

#undef FIELD128_LIFT_INLINE
} // namespace spin::research::char2::field128_lift
