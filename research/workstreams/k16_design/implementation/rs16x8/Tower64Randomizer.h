// Third quadratic tower over Tower32Randomizer's exact GF(2^32).
// Adjoin w with w^2+w+0x20000000=0. The reduction constant has absolute
// trace one. Byte basis (1,u,v,uv,w,uw,vw,uvw), least significant first.
// Uniform nonzero scalar multiplication has the same fixed-input law as
// uniform GL64, not the same multi-message correlations or setup seeds.
#pragma once
#include "Tower32Randomizer.h"

namespace spin::research::rs::tower64 {
inline constexpr std::uint32_t Theta = 0x20000000;

inline std::uint64_t multiply64(std::uint64_t a, std::uint64_t b) {
    const auto a0 = std::uint32_t(a), a1 = std::uint32_t(a >> 32);
    const auto b0 = std::uint32_t(b), b1 = std::uint32_t(b >> 32);
    const auto lo = tower32::multiply32(a0,b0) ^ tower32::multiply32(Theta,tower32::multiply32(a1,b1));
    const auto hi = tower32::multiply32(a0,b1) ^ tower32::multiply32(a1,b0) ^ tower32::multiply32(a1,b1);
    return std::uint64_t(lo) | (std::uint64_t(hi) << 32);
}

// Unlike the public randomizer factory, the recursive coefficient circuit
// accepts zero intermediate constants. No branch enters the encoding path.
inline void coefficients32(std::uint32_t scalar, std::uint64_t* output) {
    const auto c0 = std::uint16_t(scalar), c1 = std::uint16_t(scalar >> 16);
    tower32::coefficients16(c0,output);
    tower32::coefficients16(tower32::multiply16(tower32::Nu,c1),output+3);
    tower32::coefficients16(c0 ^ c1,output+6);
}

inline std::array<std::uint64_t,27> adjointCoefficients(std::uint64_t scalar) {
    if(!scalar) throw std::invalid_argument("tower64 requires a nonzero scalar");
    const auto c0 = std::uint32_t(scalar), c1 = std::uint32_t(scalar >> 32);
    std::array<std::uint64_t,27> result{};
    coefficients32(c0,result.data());
    coefficients32(tower32::multiply32(Theta,c1),result.data()+9);
    coefficients32(c0 ^ c1,result.data()+18);
    return result;
}

inline std::array<std::uint64_t,64> adjointRows(std::uint64_t scalar) {
    if(!scalar) throw std::invalid_argument("tower64 requires a nonzero scalar");
    std::array<std::uint64_t,64> result{};
    for(unsigned bit = 0; bit < 64; ++bit) result[bit] = multiply64(scalar, std::uint64_t{1} << bit);
    return result;
}

// Words must supply independent uniform 64-bit words in the ideal setup
// model. A 32-bit source requires two words and must not use this helper.
template<class Words>
inline std::uint64_t sampleNonzero(Words& words) {
    static_assert(sizeof(decltype(words())) >= sizeof(std::uint64_t), "64-bit setup words required");
    std::uint64_t scalar;
    do { scalar = std::uint64_t(words()); } while(!scalar);
    return scalar;
}

#if defined(__AVX512F__) || defined(SPIN_PACKET_ENABLE_VBMI)
// Eight byteplanes, one 64-bit payload half: 27 GFNI and 57 XOR operations.
// A 128-bit payload uses two independent calls to this same circuit.
static SPIN_FORCEINLINE void applyTranspose(__m512i& y0, __m512i& y1,
    __m512i& y2, __m512i& y3, __m512i& y4, __m512i& y5,
    __m512i& y6, __m512i& y7, const std::uint64_t* coefficients) {
    auto p0 = _mm512_xor_si512(y0,y4), p1 = _mm512_xor_si512(y1,y5);
    auto p2 = _mm512_xor_si512(y2,y6), p3 = _mm512_xor_si512(y3,y7);
    auto q0 = y0, q1 = y1, q2 = y2, q3 = y3;
    auto r0 = y4, r1 = y5, r2 = y6, r3 = y7;
    tower32::applyTranspose(p0,p1,p2,p3,coefficients);
    tower32::applyTranspose(q0,q1,q2,q3,coefficients+9);
    tower32::applyTranspose(r0,r1,r2,r3,coefficients+18);
    y0 = _mm512_xor_si512(p0,r0); y1 = _mm512_xor_si512(p1,r1);
    y2 = _mm512_xor_si512(p2,r2); y3 = _mm512_xor_si512(p3,r3);
    y4 = _mm512_xor_si512(q0,r0); y5 = _mm512_xor_si512(q1,r1);
    y6 = _mm512_xor_si512(q2,r2); y7 = _mm512_xor_si512(q3,r3);
}
#endif
} // namespace spin::research::rs::tower64
