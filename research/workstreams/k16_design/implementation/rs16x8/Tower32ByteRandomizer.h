// Choose forward symbol maps as ADJOINTS of field multiplications. Their
// transposes are ordinary multiplications, implemented with byte GFNI MUL.
// This family is exactly transitive on nonzero binary symbols: for x != 0,
// (M_c^T-M_d^T)x = M_(c+d)^T x != 0 when c != d. Independent uniform scalars
// therefore preserve the proof's fixed-input symbol law.
#pragma once
#include "Tower32Randomizer.h"

namespace spin::research::rs::tower32byte {
inline void coefficients16(std::uint16_t scalar, std::uint8_t* output) {
    const auto c0 = std::uint8_t(scalar), c1 = std::uint8_t(scalar >> 8);
    output[0] = c0;
    output[1] = tower32::multiply8(tower32::Mu,c1);
    output[2] = c0 ^ c1;
}

inline std::array<std::uint8_t,9> coefficients(std::uint32_t scalar) {
    if(!scalar) throw std::invalid_argument("tower32byte requires a nonzero scalar");
    const auto c0 = std::uint16_t(scalar), c1 = std::uint16_t(scalar >> 16);
    std::array<std::uint8_t,9> result{};
    coefficients16(c0,result.data());
    coefficients16(tower32::multiply16(tower32::Nu,c1),result.data()+3);
    coefficients16(c0 ^ c1,result.data()+6);
    return result;
}

// Rows of ordinary field multiplication, which is the transposed encoder's
// sampled symbol map for this construction. Not tower32::adjointRows.
inline std::array<std::uint32_t,32> multiplyRows(std::uint32_t scalar) {
    if(!scalar) throw std::invalid_argument("tower32byte requires a nonzero scalar");
    std::array<std::uint32_t,32> result{};
    for(unsigned column = 0; column < 32; ++column) {
        const auto image = tower32::multiply32(scalar, std::uint32_t{1} << column);
        for(unsigned row = 0; row < 32; ++row)
            result[row] |= ((image >> row) & 1U) << column;
    }
    return result;
}

#if defined(__AVX512F__) || defined(SPIN_PACKET_ENABLE_VBMI)
struct Pair { __m512i lo,hi; };
static SPIN_FORCEINLINE Pair quadraticMultiply(__m512i x0, __m512i x1, const std::uint8_t* c) {
    const auto m0 = _mm512_set1_epi8(char(c[0]));
    const auto m1 = _mm512_set1_epi8(char(c[1]));
    const auto m2 = _mm512_set1_epi8(char(c[2]));
    const auto p0 = _mm512_gf2p8mul_epi8(x0,m0);
    const auto p1 = _mm512_gf2p8mul_epi8(x1,m1);
    const auto p2 = _mm512_gf2p8mul_epi8(_mm512_xor_si512(x0,x1),m2);
    return {_mm512_xor_si512(p0,p1),_mm512_xor_si512(p0,p2)};
}

// Four byteplanes, one 64-bit payload half; 9 byte multiplications + 15 XORs.
// Tables contain 9 bytes per symbol, not 9 eight-by-eight affine matrices.
static SPIN_FORCEINLINE void applyMultiply(__m512i& y0, __m512i& y1,
    __m512i& y2, __m512i& y3, const std::uint8_t* coeff) {
    const auto p = quadraticMultiply(y0,y1,coeff);
    const auto q = quadraticMultiply(y2,y3,coeff+3);
    const auto r = quadraticMultiply(_mm512_xor_si512(y0,y2),_mm512_xor_si512(y1,y3),coeff+6);
    y0 = _mm512_xor_si512(p.lo,q.lo);
    y1 = _mm512_xor_si512(p.hi,q.hi);
    y2 = _mm512_xor_si512(p.lo,r.lo);
    y3 = _mm512_xor_si512(p.hi,r.hi);
}
#endif
} // namespace spin::research::rs::tower32byte
