#pragma once
#include "SparseByteField128Simd.h"

// Experimental arithmetic control only. The measured CLMUL primitive and
// LiteralEncoder.h remain unchanged. This header uses the same field and
// raw sixteen-byte coefficients; no setup expansion or basis conversion.
namespace spin::research::char2::field128_lift::sparse_byte {
#if defined(_MSC_VER)
#define SPARSE128_GFNI_INLINE __forceinline
#else
#define SPARSE128_GFNI_INLINE inline __attribute__((always_inline))
#endif

template<unsigned I>
SPARSE128_GFNI_INLINE __m512i gfniCoefficientProduct(__m512i a,__m512i b) noexcept {
    static_assert(I<16);
    // Each independent128-bit record broadcasts its own coefficient b[I].
    return _mm512_gf2p8mul_epi8(a,_mm512_shuffle_epi8(b,_mm512_set1_epi8(static_cast<char>(I))));
}
template<unsigned I>
SPARSE128_GFNI_INLINE void gfniConvolutionPair(__m512i a,__m512i b,__m512i& low,__m512i& high) noexcept {
    static_assert(I>=2&&I<=14&&I%2==0);
    const auto p=gfniCoefficientProduct<I>(a,b),q=gfniCoefficientProduct<I+1>(a,b);
    low=_mm512_ternarylogic_epi64(low,_mm512_bslli_epi128(p,I),_mm512_bslli_epi128(q,I+1),0x96);
    high=_mm512_ternarylogic_epi64(high,_mm512_bsrli_epi128(p,16-I),_mm512_bsrli_epi128(q,15-I),0x96);
}

// Sixteen coefficientwise GF256 products form the unreduced polynomial
// convolution. Its low/high halves contain degrees0..15 and16..30.
// Static source schedule per four records:16 variable GFNI,16 broadcasts,
// 30 byte shifts, and15 XOR/ternary accumulations. The shared sparse-field
// reduction adds2 constant GFNI and its shifts/XORs. There are no CLMULs.
SPARSE128_GFNI_INLINE __m512i multiply4Gfni(__m512i a,__m512i b) noexcept {
    const auto p0=gfniCoefficientProduct<0>(a,b),p1=gfniCoefficientProduct<1>(a,b);
    auto low=_mm512_xor_si512(p0,_mm512_bslli_epi128(p1,1));
    auto high=_mm512_bsrli_epi128(p1,15);
    gfniConvolutionPair<2>(a,b,low,high);
    gfniConvolutionPair<4>(a,b,low,high);
    gfniConvolutionPair<6>(a,b,low,high);
    gfniConvolutionPair<8>(a,b,low,high);
    gfniConvolutionPair<10>(a,b,low,high);
    gfniConvolutionPair<12>(a,b,low,high);
    gfniConvolutionPair<14>(a,b,low,high);
    return reduceExtension(low,high);
}
SPARSE128_GFNI_INLINE __m512i multiply4Gfni(__m512i value,const Element* coefficients) noexcept {
    return multiply4Gfni(value,_mm512_loadu_si512(coefficients));
}
SPARSE128_GFNI_INLINE void multiply4GfniRecords(const Element* input,const Element* coefficients,Element* output) noexcept {
    const auto a=_mm512_loadu_si512(input),b=_mm512_loadu_si512(coefficients);
    _mm512_storeu_si512(output,multiply4Gfni(a,b));
}

#undef SPARSE128_GFNI_INLINE
} // namespace spin::research::char2::field128_lift::sparse_byte

// Optional, deliberately scoped compiler-force-include adapter. Defining
// FIELD128_LIFT_GFNI_ENCODER and force-including THIS header compiles the
// otherwise unmodified literal_encoder.cpp with the alternate primitive.
// SparseByteField128Simd.h has already been parsed, so its measured functions
// keep their original names. The identifier replacement affects only the
// subsequent LiteralEncoder.h body, then is removed. That header's existing
// dependencies contain no other multiply4 identifier. Do not combine binary
// objects built with different adapters: encoder function names are shared.
// Example from the char2_field_spin directory:
// g++ -O3 -std=c++20 -march=znver4 -DFIELD128_LIFT_GFNI_ENCODER=1 \
//   -include field128_lift/SparseByteField128Gfni.h \
//   field128_lift/literal_encoder.cpp -o literal_encoder_gfni16
#if defined(FIELD128_LIFT_GFNI_ENCODER) && FIELD128_LIFT_GFNI_ENCODER
#ifdef multiply4
#error "GFNI encoder adapter requires the multiply4 identifier to be unmodified"
#endif
#define multiply4 multiply4Gfni
#include "LiteralEncoder.h"
#undef multiply4
#endif
