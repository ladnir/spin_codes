#pragma once
#include "SparseByteField128MulY.h"
#include <immintrin.h>

namespace spin::research::char2::field128_lift::sparse_byte {
#if defined(_MSC_VER)
#define SPARSE128_MULY_INLINE __forceinline
#else
#define SPARSE128_MULY_INLINE inline __attribute__((always_inline))
#endif

// Four independent 16-byte records, with no lane crossing. Static work:
// two byte shifts, one byte shuffle, one GFNI, and one ternary XOR.
// Requires AVX512F/BW and GFNI; VPCLMULQDQ is not required.
SPARSE128_MULY_INLINE __m512i multiplyY4(__m512i value) noexcept {
    const auto duplicateTop=_mm512_setr_epi64(
        static_cast<long long>(0x80800f800f808080ULL),static_cast<long long>(0x8080808080808080ULL),
        static_cast<long long>(0x80800f800f808080ULL),static_cast<long long>(0x8080808080808080ULL),
        static_cast<long long>(0x80800f800f808080ULL),static_cast<long long>(0x8080808080808080ULL),
        static_cast<long long>(0x80800f800f808080ULL),static_cast<long long>(0x8080808080808080ULL));
    const auto shifted=_mm512_bslli_epi128(value,1);
    const auto top=_mm512_bsrli_epi128(value,15);
    const auto constant=_mm512_gf2p8mul_epi8(top,_mm512_set1_epi8(0x32));
    const auto taps=_mm512_shuffle_epi8(value,duplicateTop); // a_15 at byte3 and byte5.
    return _mm512_ternarylogic_epi64(shifted,constant,taps,0x96);
}

// Unaligned byte addresses are permitted. Loading all64 bytes before the
// store also permits exact aliasing and arbitrary partial overlap.
SPARSE128_MULY_INLINE void multiplyY4Records(const void* input,void* output) noexcept {
    const auto value=_mm512_loadu_si512(input);
    _mm512_storeu_si512(output,multiplyY4(value));
}

#undef SPARSE128_MULY_INLINE
} // namespace spin::research::char2::field128_lift::sparse_byte
