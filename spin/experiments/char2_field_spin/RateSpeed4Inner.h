#pragma once
#include "RateSpeed8Inner.h"

namespace spin::research::char2 {

// Normalized full8 Cauchy packet, rows0..7 and columns128..135.
// A below includes eta2. F=107*A*C; both A and C have unit row sum.
// Four input/state symbols, first3 emitted symbols retained, initial state
// zero and terminal state discarded by the surrounding code.
inline constexpr std::uint8_t rs4A[4]={30,235,104,156};
inline constexpr std::uint8_t rs4C[4]={27,28,18,20};
inline constexpr std::uint8_t rs4F[4]={86,90,226,133};
inline constexpr std::uint8_t rs4Scale=107;

template<bool Emission,unsigned Index>
constexpr std::uint8_t rs4ZetaCoefficient()noexcept{
    static_assert(Index<4);
    std::uint8_t value=0;
    for(unsigned i=0;i!=4;++i)if((i&Index)==Index)value^=Emission?rs4A[i]:rs4C[i];
    return value;
}
constexpr bool rs4CheckMapConstants()noexcept{
    if(rs4ZetaCoefficient<true,0>()!=1||rs4ZetaCoefficient<false,0>()!=1)return false;
    for(unsigned row=0;row!=4;++row)for(unsigned column=0;column!=4;++column){
        std::uint8_t value=0;
        for(unsigned k=0;k!=4;++k)value^=sfMultiply(rs4A[row^k],rs4C[k^column]);
        if(sfMultiply(rs4Scale,value)!=rs4F[row^column])return false;
    }
    return true;
}
static_assert(rs4CheckMapConstants());

#if defined(_MSC_VER)
#define RS4_INLINE __forceinline
#else
#define RS4_INLINE inline __attribute__((always_inline))
#endif

// Four independent four-coordinate XOR circulants, one packet per record
// lane. The normalized subset circuit uses exactly five GFNI sites.
template<bool Emission,bool Affine>
RS4_INLINE void rs4Circulant(__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3)noexcept{
    c8Zeta(z0,z1,z2,z3);
    z3=_mm512_ternarylogic_epi32(z3,rs8ConstantMultiply<Affine,rs4ZetaCoefficient<Emission,1>()>(z2),
        rs8ConstantMultiply<Affine,rs4ZetaCoefficient<Emission,2>()>(z1),0x96);
    z3=_mm512_xor_si512(z3,rs8ConstantMultiply<Affine,rs4ZetaCoefficient<Emission,3>()>(z0));
    z2=_mm512_xor_si512(z2,rs8ConstantMultiply<Affine,rs4ZetaCoefficient<Emission,2>()>(z0));
    z1=_mm512_xor_si512(z1,rs8ConstantMultiply<Affine,rs4ZetaCoefficient<Emission,1>()>(z0));
    c8Zeta(z0,z1,z2,z3);
}
template<bool AffineScale>
RS4_INLINE void rs4HistoryDiagonal(__m512i& moment,__m512i& history,__m512i diagonal)noexcept{
    const auto old=history;
    history=_mm512_gf2p8mul_epi8(_mm512_xor_si512(old,moment),diagonal);
    moment=_mm512_xor_si512(old,rs8ConstantMultiply<AffineScale,rs4Scale>(moment));
}
template<bool AffineScale>
RS4_INLINE void rs4History(__m512i& moment,__m512i& history,const std::uint8_t* diagonal)noexcept{
    rs4HistoryDiagonal<AffineScale>(moment,history,sfBroadcastFour(diagonal));
}

// input[4*j+phase], j=0,1,2, holds twelve retained records from four
// consecutive chronological packets. The missing coordinate3 is zero.
// output[4*j+phase], j=0..3, contains sixteen folded raw-input adjoints,
// before input diagonals and inverse routing. history is the four-coordinate
// terminal dual state on entry and preceding dual state on return.
// Eighteen GFNI sites here plus four external Dx sites per16 raw records.
// PackedDiagonals reuses one compact16-byte load for all four state updates;
// it changes neither the coefficient layout nor the chronological recurrence.
template<bool AffineMaps=false,bool AffineScale=false,bool PackedDiagonals=false>
RS4_INLINE void rs4InnerTransposeFour(const Block* input,Block* output,__m512i& history,
                                      const std::uint8_t* diagonal)noexcept{
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_setzero_si512();
    rs4Circulant<true,AffineMaps>(z0,z1,z2,z3);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    if constexpr(PackedDiagonals){
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal)));
        rs4HistoryDiagonal<AffineScale>(z3,history,sf16ExpandRow<3>(packed));
        rs4HistoryDiagonal<AffineScale>(z2,history,sf16ExpandRow<2>(packed));
        rs4HistoryDiagonal<AffineScale>(z1,history,sf16ExpandRow<1>(packed));
        rs4HistoryDiagonal<AffineScale>(z0,history,sf16ExpandRow<0>(packed));
    }else{
        rs4History<AffineScale>(z3,history,diagonal+12);
        rs4History<AffineScale>(z2,history,diagonal+8);
        rs4History<AffineScale>(z1,history,diagonal+4);
        rs4History<AffineScale>(z0,history,diagonal);
    }
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    rs4Circulant<false,AffineMaps>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
    _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
}

// Two consecutive four-packet tiles: input[0..11] is the earlier tile,
// input[12..23] the later tile. Output keeps the same two folded16-record
// tiles, not one coordinate-major eight-packet tile. Both emission maps are
// independent of history and are prepared before its eight reverse updates.
// The later injection map is independent of the earlier history updates.
// Arithmetic is unchanged: thirty-six GFNI sites plus eight external Dx
// sites per32 raw records. This variant changes scheduling, not the map.
template<bool AffineMaps=false,bool AffineScale=false,bool PackedDiagonals=true>
RS4_INLINE void rs4InnerTransposeEight(const Block* input,Block* output,__m512i& history,
                                       const std::uint8_t* diagonal)noexcept{
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_setzero_si512();
    auto w0=_mm512_loadu_si512(input+12),w1=_mm512_loadu_si512(input+16);
    auto w2=_mm512_loadu_si512(input+20),w3=_mm512_setzero_si512();
    rs4Circulant<true,AffineMaps>(z0,z1,z2,z3);
    rs4Circulant<true,AffineMaps>(w0,w1,w2,w3);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(w0,w1,w2,w3,w0,w1,w2,w3);
    if constexpr(PackedDiagonals){
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal+16)));
        rs4HistoryDiagonal<AffineScale>(w3,history,sf16ExpandRow<3>(packed));
        rs4HistoryDiagonal<AffineScale>(w2,history,sf16ExpandRow<2>(packed));
        rs4HistoryDiagonal<AffineScale>(w1,history,sf16ExpandRow<1>(packed));
        rs4HistoryDiagonal<AffineScale>(w0,history,sf16ExpandRow<0>(packed));
    }else{
        rs4History<AffineScale>(w3,history,diagonal+28);
        rs4History<AffineScale>(w2,history,diagonal+24);
        rs4History<AffineScale>(w1,history,diagonal+20);
        rs4History<AffineScale>(w0,history,diagonal+16);
    }
    c8TransposeRecords(w0,w1,w2,w3,w0,w1,w2,w3);
    rs4Circulant<false,AffineMaps>(w0,w1,w2,w3);
    if constexpr(PackedDiagonals){
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal)));
        rs4HistoryDiagonal<AffineScale>(z3,history,sf16ExpandRow<3>(packed));
        rs4HistoryDiagonal<AffineScale>(z2,history,sf16ExpandRow<2>(packed));
        rs4HistoryDiagonal<AffineScale>(z1,history,sf16ExpandRow<1>(packed));
        rs4HistoryDiagonal<AffineScale>(z0,history,sf16ExpandRow<0>(packed));
    }else{
        rs4History<AffineScale>(z3,history,diagonal+12);
        rs4History<AffineScale>(z2,history,diagonal+8);
        rs4History<AffineScale>(z1,history,diagonal+4);
        rs4History<AffineScale>(z0,history,diagonal);
    }
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    rs4Circulant<false,AffineMaps>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
    _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
    _mm512_storeu_si512(output+16,w0);_mm512_storeu_si512(output+20,w1);
    _mm512_storeu_si512(output+24,w2);_mm512_storeu_si512(output+28,w3);
}

#undef RS4_INLINE
} // namespace spin::research::char2
