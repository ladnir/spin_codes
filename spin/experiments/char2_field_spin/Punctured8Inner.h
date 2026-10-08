#pragma once
#include "Cubic8Inner.h"

namespace spin::research::char2 {

// Exact normalized full-Cauchy maps; A here includes eta=2. All matrices
// are symmetric XOR circulants. F=242*A*C, so no separate F circuit is used.
inline constexpr std::uint8_t p8A[8]={172,59,114,102,161,248,77,150};
inline constexpr std::uint8_t p8C[8]={52,172,158,17,209,32,250,29};
inline constexpr std::uint8_t p8F[8]={184,226,103,44,151,1,216,173};
template<bool Emission,unsigned Index>
constexpr std::uint8_t p8ZetaCoefficient() noexcept {
    std::uint8_t value=0;
    for(unsigned i=0;i!=8;++i)if((i&Index)==Index)value^=Emission?p8A[i]:p8C[i];
    return value;
}
inline void p8Zeta(__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3,
                   __m512i& z4,__m512i& z5,__m512i& z6,__m512i& z7) noexcept {
    c8Zeta(z0,z1,z2,z3);c8Zeta(z4,z5,z6,z7);
    z0=_mm512_xor_si512(z0,z4);z1=_mm512_xor_si512(z1,z5);
    z2=_mm512_xor_si512(z2,z6);z3=_mm512_xor_si512(z3,z7);
}
// Four packets at once, one coordinate per register. Nineteen GFNI sites;
// descending subset convolution keeps eight work registers and seven constants.
template<bool Emission>
inline void p8Circulant(__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3,
                        __m512i& z4,__m512i& z5,__m512i& z6,__m512i& z7) noexcept {
    static_assert(p8ZetaCoefficient<Emission,0>()==1);
    p8Zeta(z0,z1,z2,z3,z4,z5,z6,z7);
    const auto g1=_mm512_set1_epi8(p8ZetaCoefficient<Emission,1>());
    const auto g2=_mm512_set1_epi8(p8ZetaCoefficient<Emission,2>());
    const auto g3=_mm512_set1_epi8(p8ZetaCoefficient<Emission,3>());
    const auto g4=_mm512_set1_epi8(p8ZetaCoefficient<Emission,4>());
    const auto g5=_mm512_set1_epi8(p8ZetaCoefficient<Emission,5>());
    const auto g6=_mm512_set1_epi8(p8ZetaCoefficient<Emission,6>());
    const auto g7=_mm512_set1_epi8(p8ZetaCoefficient<Emission,7>());
    z7=_mm512_ternarylogic_epi32(z7,_mm512_gf2p8mul_epi8(z6,g1),_mm512_gf2p8mul_epi8(z5,g2),0x96);
    z7=_mm512_ternarylogic_epi32(z7,_mm512_gf2p8mul_epi8(z4,g3),_mm512_gf2p8mul_epi8(z3,g4),0x96);
    z7=_mm512_ternarylogic_epi32(z7,_mm512_gf2p8mul_epi8(z2,g5),_mm512_gf2p8mul_epi8(z1,g6),0x96);
    z7=_mm512_xor_si512(z7,_mm512_gf2p8mul_epi8(z0,g7));
    z6=_mm512_ternarylogic_epi32(z6,_mm512_gf2p8mul_epi8(z4,g2),_mm512_gf2p8mul_epi8(z2,g4),0x96);
    z6=_mm512_xor_si512(z6,_mm512_gf2p8mul_epi8(z0,g6));
    z5=_mm512_ternarylogic_epi32(z5,_mm512_gf2p8mul_epi8(z4,g1),_mm512_gf2p8mul_epi8(z1,g4),0x96);
    z5=_mm512_xor_si512(z5,_mm512_gf2p8mul_epi8(z0,g5));
    z4=_mm512_xor_si512(z4,_mm512_gf2p8mul_epi8(z0,g4));
    z3=_mm512_ternarylogic_epi32(z3,_mm512_gf2p8mul_epi8(z2,g1),_mm512_gf2p8mul_epi8(z1,g2),0x96);
    z3=_mm512_xor_si512(z3,_mm512_gf2p8mul_epi8(z0,g3));
    z2=_mm512_xor_si512(z2,_mm512_gf2p8mul_epi8(z0,g2));
    z1=_mm512_xor_si512(z1,_mm512_gf2p8mul_epi8(z0,g1));
    p8Zeta(z0,z1,z2,z3,z4,z5,z6,z7);
}
inline void p8History(__m512i& moment0,__m512i& moment1,__m512i& h0,__m512i& h1,
                       const std::uint8_t* diagonal) noexcept {
    const auto old0=h0,old1=h1;
    h0=_mm512_gf2p8mul_epi8(_mm512_xor_si512(old0,moment0),sfBroadcastFour(diagonal));
    h1=_mm512_gf2p8mul_epi8(_mm512_xor_si512(old1,moment1),sfBroadcastFour(diagonal+4));
    const auto scale=_mm512_set1_epi8(static_cast<char>(242));
    moment0=_mm512_xor_si512(old0,_mm512_gf2p8mul_epi8(moment0,scale));
    moment1=_mm512_xor_si512(old1,_mm512_gf2p8mul_epi8(moment1,scale));
}
// Input is coordinate-major within four consecutive chronological packets:
// [j0:p0,p1,p2,p3; ...; j6:p0,p1,p2,p3]. Coordinate7 is punctured.
// Output is packet-major (32 records), before the input diagonal and route.
#ifndef SPIN_P8_FORCE_INLINE
#define SPIN_P8_FORCE_INLINE 1
#endif
#if SPIN_P8_FORCE_INLINE && defined(_MSC_VER)
#define P8_PACKET_INLINE __forceinline
#elif SPIN_P8_FORCE_INLINE
#define P8_PACKET_INLINE inline __attribute__((always_inline))
#else
#define P8_PACKET_INLINE inline
#endif
P8_PACKET_INLINE void p8InnerTransposeFour(const Block* input,Block* output,
                                 __m512i& h0,__m512i& h1,const std::uint8_t* diagonal) noexcept {
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_loadu_si512(input+12);
    auto z4=_mm512_loadu_si512(input+16),z5=_mm512_loadu_si512(input+20);
    auto z6=_mm512_loadu_si512(input+24),z7=_mm512_setzero_si512();
    p8Circulant<true>(z0,z1,z2,z3,z4,z5,z6,z7);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    p8History(z3,z7,h0,h1,diagonal+24);
    p8History(z2,z6,h0,h1,diagonal+16);
    p8History(z1,z5,h0,h1,diagonal+8);
    p8History(z0,z4,h0,h1,diagonal);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    p8Circulant<false>(z0,z1,z2,z3,z4,z5,z6,z7);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z4);
    _mm512_storeu_si512(output+8,z1);_mm512_storeu_si512(output+12,z5);
    _mm512_storeu_si512(output+16,z2);_mm512_storeu_si512(output+20,z6);
    _mm512_storeu_si512(output+24,z3);_mm512_storeu_si512(output+28,z7);
}
#undef P8_PACKET_INLINE
} // namespace spin::research::char2
