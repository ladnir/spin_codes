// Exact GF(2^32) nonzero-scalar randomizer. Research candidate, not GL32.
// Basis: bytes (1,u,v,uv), least significant first; GF256 modulus0x11b,
// u^2+u+0x20=0 and v^2+v+0x2000=0. Both constants have absolute trace1.
// Independently uniform nonzero scalars preserve the fixed-symbol law used
// by the outer first-moment proof. Setup and seeds differ from dense GL32.
#pragma once
#include "../../../../../spin/src/kernels/Block.h"
#include <array>
#include <cstdint>
#include <stdexcept>

namespace spin::research::rs::tower32 {
inline constexpr std::uint8_t Mu = 0x20;
inline constexpr std::uint16_t Nu = 0x2000;

inline std::uint8_t multiply8(std::uint8_t a, std::uint8_t b) {
    std::uint8_t result = 0;
    for(unsigned bit = 0; bit < 8; ++bit) {
        if(b & 1U) result ^= a;
        a = std::uint8_t((unsigned(a) << 1) ^ ((a & 128U) ? 0x11bU : 0U));
        b >>= 1;
    }
    return result;
}

inline std::uint16_t multiply16(std::uint16_t a, std::uint16_t b) {
    const auto a0 = std::uint8_t(a), a1 = std::uint8_t(a >> 8);
    const auto b0 = std::uint8_t(b), b1 = std::uint8_t(b >> 8);
    const auto lo = multiply8(a0,b0) ^ multiply8(Mu,multiply8(a1,b1));
    const auto hi = multiply8(a0,b1) ^ multiply8(a1,b0) ^ multiply8(a1,b1);
    return std::uint16_t(lo | (unsigned(hi) << 8));
}

inline std::uint32_t multiply32(std::uint32_t a, std::uint32_t b) {
    const auto a0 = std::uint16_t(a), a1 = std::uint16_t(a >> 16);
    const auto b0 = std::uint16_t(b), b1 = std::uint16_t(b >> 16);
    const auto lo = multiply16(a0,b0) ^ multiply16(Nu,multiply16(a1,b1));
    const auto hi = multiply16(a0,b1) ^ multiply16(a1,b0) ^ multiply16(a1,b1);
    return std::uint32_t(lo) | (std::uint32_t(hi) << 16);
}

inline std::uint64_t adjointMatrix8(std::uint8_t scalar) {
    std::uint64_t matrix = 0;
    for(unsigned bit = 0; bit < 8; ++bit)
        matrix |= std::uint64_t(multiply8(scalar,std::uint8_t(1U << bit))) << (8 * (7 - bit));
    return matrix;
}

inline void coefficients16(std::uint16_t scalar, std::uint64_t* output) {
    const auto c0 = std::uint8_t(scalar), c1 = std::uint8_t(scalar >> 8);
    output[0] = adjointMatrix8(c0);
    output[1] = adjointMatrix8(multiply8(Mu,c1));
    output[2] = adjointMatrix8(c0 ^ c1);
}

inline std::array<std::uint64_t,9> adjointCoefficients(std::uint32_t scalar) {
    if(!scalar) throw std::invalid_argument("tower32 requires a nonzero scalar");
    const auto c0 = std::uint16_t(scalar), c1 = std::uint16_t(scalar >> 16);
    std::array<std::uint64_t,9> result{};
    coefficients16(c0,result.data());
    coefficients16(multiply16(Nu,c1),result.data()+3);
    coefficients16(c0 ^ c1,result.data()+6);
    return result;
}

// Rows of the binary transpose are columns of forward multiplication.
inline std::array<std::uint32_t,32> adjointRows(std::uint32_t scalar) {
    if(!scalar) throw std::invalid_argument("tower32 requires a nonzero scalar");
    std::array<std::uint32_t,32> result{};
    for(unsigned bit = 0; bit < 32; ++bit) result[bit] = multiply32(scalar, std::uint32_t{1} << bit);
    return result;
}

// Setup only. With ideal independent uniform input words, rejecting zero
// produces an exactly uniform nonzero scalar; there is no modulo bias.
template<class Words>
inline std::uint32_t sampleNonzero(Words& words) {
    std::uint32_t scalar;
    do { scalar = std::uint32_t(words()); } while(!scalar);
    return scalar;
}

#if defined(__AVX512F__) || defined(SPIN_PACKET_ENABLE_VBMI)
struct Pair { __m512i lo,hi; };
static SPIN_FORCEINLINE Pair quadraticTranspose(__m512i y0, __m512i y1, const std::uint64_t* c) {
    const auto m0 = _mm512_set1_epi64(c[0]);
    const auto m1 = _mm512_set1_epi64(c[1]);
    const auto m2 = _mm512_set1_epi64(c[2]);
    const auto p0 = _mm512_gf2p8affine_epi64_epi8(_mm512_xor_si512(y0,y1),m0,0);
    const auto p1 = _mm512_gf2p8affine_epi64_epi8(y0,m1,0);
    const auto p2 = _mm512_gf2p8affine_epi64_epi8(y1,m2,0);
    return {_mm512_xor_si512(p0,p2),_mm512_xor_si512(p1,p2)};
}

// Four byteplanes, one64-bit payload half. Call once per payload half.
// Fully explicit9-GFNI/15-XOR circuit; no allocations or runtime dispatch.
static SPIN_FORCEINLINE void applyTranspose(__m512i& y0, __m512i& y1,
    __m512i& y2, __m512i& y3, const std::uint64_t* coefficients) {
    const auto p = quadraticTranspose(_mm512_xor_si512(y0,y2),_mm512_xor_si512(y1,y3),coefficients);
    const auto q = quadraticTranspose(y0,y1,coefficients+3);
    const auto r = quadraticTranspose(y2,y3,coefficients+6);
    y0 = _mm512_xor_si512(p.lo,r.lo);
    y1 = _mm512_xor_si512(p.hi,r.hi);
    y2 = _mm512_xor_si512(q.lo,r.lo);
    y3 = _mm512_xor_si512(q.hi,r.hi);
}
#endif
} // namespace spin::research::rs::tower32
