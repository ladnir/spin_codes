#pragma once
#include "SparseByteField128.h"
#include <immintrin.h>

namespace spin::research::char2::field128_lift::sparse_byte {
#if defined(_MSC_VER)
#define SPARSE128_INLINE __forceinline
#else
#define SPARSE128_INLINE inline __attribute__((always_inline))
#endif

struct RawProduct {__m512i low,high;};

// Four independent, unreduced 128x128 carryless products. The caller's
// 16-bit coefficient spacing prevents adjacent AES products from overlapping.
SPARSE128_INLINE RawProduct carryless128(__m512i a,__m512i b) noexcept {
    const auto p=_mm512_clmulepi64_epi128(a,b,0x00);
    const auto q=_mm512_clmulepi64_epi128(a,b,0x11);
    const auto aa=_mm512_xor_si512(a,_mm512_shuffle_epi32(a,static_cast<_MM_PERM_ENUM>(0x4e)));
    const auto bb=_mm512_xor_si512(b,_mm512_shuffle_epi32(b,static_cast<_MM_PERM_ENUM>(0x4e)));
    const auto r=_mm512_clmulepi64_epi128(aa,bb,0x00);
    const auto middle=_mm512_ternarylogic_epi64(p,q,r,0x96);
    return {_mm512_xor_si512(p,_mm512_bslli_epi128(middle,8)),
            _mm512_xor_si512(q,_mm512_bsrli_epi128(middle,8))};
}

// Each 16-bit slot is a polynomial of degree at most14 over GF2.
// X^8=0x1b in AES GF256, so the low byte after this operation is its
// field remainder. The shuffle discards each unreduced high byte.
SPARSE128_INLINE __m512i reduceAndPack(__m512i low,__m512i high) noexcept {
    const auto r=_mm512_set1_epi8(0x1b);
    low=_mm512_xor_si512(low,_mm512_gf2p8mul_epi8(_mm512_srli_epi16(low,8),r));
    high=_mm512_xor_si512(high,_mm512_gf2p8mul_epi8(_mm512_srli_epi16(high,8),r));
    const auto lowBytes=_mm512_setr_epi64(0x0e0c0a0806040200LL,static_cast<long long>(0x8080808080808080ULL),
        0x0e0c0a0806040200LL,static_cast<long long>(0x8080808080808080ULL),
        0x0e0c0a0806040200LL,static_cast<long long>(0x8080808080808080ULL),
        0x0e0c0a0806040200LL,static_cast<long long>(0x8080808080808080ULL));
    const auto highBytes=_mm512_setr_epi64(static_cast<long long>(0x8080808080808080ULL),0x0e0c0a0806040200LL,
        static_cast<long long>(0x8080808080808080ULL),0x0e0c0a0806040200LL,
        static_cast<long long>(0x8080808080808080ULL),0x0e0c0a0806040200LL,
        static_cast<long long>(0x8080808080808080ULL),0x0e0c0a0806040200LL);
    return _mm512_xor_si512(_mm512_shuffle_epi8(low,lowBytes),_mm512_shuffle_epi8(high,highBytes));
}

// Reduce L+Y^16 H using Y^16=Y^5+Y^3+0x32. The first fold's overflow
// J has degree at most4 (at most3 for a product), so its second fold fits.
SPARSE128_INLINE __m512i reduceExtension(__m512i low,__m512i high) noexcept {
    const auto c=_mm512_set1_epi8(0x32);
    const auto overflow=_mm512_xor_si512(_mm512_bsrli_epi128(high,13),_mm512_bsrli_epi128(high,11));
    const auto first=_mm512_ternarylogic_epi64(_mm512_gf2p8mul_epi8(high,c),
        _mm512_bslli_epi128(high,3),_mm512_bslli_epi128(high,5),0x96);
    const auto second=_mm512_ternarylogic_epi64(_mm512_gf2p8mul_epi8(overflow,c),
        _mm512_bslli_epi128(overflow,3),_mm512_bslli_epi128(overflow,5),0x96);
    return _mm512_ternarylogic_epi64(low,first,second,0x96);
}

// Genuine GF128 multiplication, with four independent 16-byte keys.
// Raw keys are used directly: no per-key matrix or expanded leaf table.
// Static arithmetic count: nine VPCLMUL and six GFNI multiplications per
// four records, plus the explicit shifts, shuffles, and Boolean operations.
// Requires AVX512F/BW, GFNI, and VPCLMULQDQ. No basis conversion/allocation.
SPARSE128_INLINE __m512i multiply4(__m512i value,__m512i coefficient) noexcept {
    const auto zero=_mm512_setzero_si512();
    const auto a0=_mm512_unpacklo_epi8(value,zero),a1=_mm512_unpackhi_epi8(value,zero);
    const auto b0=_mm512_unpacklo_epi8(coefficient,zero),b1=_mm512_unpackhi_epi8(coefficient,zero);
    const auto p=carryless128(a0,b0),q=carryless128(a1,b1);
    const auto r=carryless128(_mm512_xor_si512(a0,a1),_mm512_xor_si512(b0,b1));
    const auto middle0=_mm512_ternarylogic_epi64(p.low,q.low,r.low,0x96);
    const auto middle1=_mm512_ternarylogic_epi64(p.high,q.high,r.high,0x96);
    const auto low=reduceAndPack(p.low,_mm512_xor_si512(p.high,middle0));
    const auto high=reduceAndPack(_mm512_xor_si512(q.low,middle1),q.high);
    return reduceExtension(low,high);
}
SPARSE128_INLINE __m512i multiply4(__m512i value,const Element* coefficients) noexcept {
    return multiply4(value,_mm512_loadu_si512(coefficients));
}
SPARSE128_INLINE void multiply4Records(const Element* input,const Element* coefficients,Element* output) noexcept {
    const auto value=_mm512_loadu_si512(input),c=_mm512_loadu_si512(coefficients);
    _mm512_storeu_si512(output,multiply4(value,c));
}

#undef SPARSE128_INLINE
} // namespace spin::research::char2::field128_lift::sparse_byte
