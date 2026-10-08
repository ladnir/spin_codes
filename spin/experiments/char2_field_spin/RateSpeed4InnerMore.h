#pragma once
#include "RateSpeed4InnerNext.h"

namespace spin::research::char2 {

#if defined(_MSC_VER)
#define RS4_MORE_INLINE __forceinline
#else
#define RS4_MORE_INLINE inline __attribute__((always_inline))
#endif

// The fourth retained emission input is zero. Write b=A1+A3, c=A2+A3,
// d=A3 and S=x0+x1+x2. With T=b*x2+c*x1+d*S, the four outputs are
// [x0+b*S+c*S+T, x1+b*S+T, x2+c*S+T, T]. This is the normalized subset
// circuit with its input-zero and cancelling-XOR identities made explicit.
// Five GFNI and six Boolean sites, versus five and ten in the generic
// source circuit after eliminating operations whose operand is zero.
template<bool Affine=true>
RS4_MORE_INLINE void rs4CompactEmission(__m512i& x0,__m512i& x1,
    __m512i& x2,__m512i& x3) noexcept {
    constexpr auto b=rs4ZetaCoefficient<true,1>();
    constexpr auto c=rs4ZetaCoefficient<true,2>();
    constexpr auto d=rs4ZetaCoefficient<true,3>();
    const auto sum=_mm512_ternarylogic_epi32(x0,x1,x2,0x96);
    const auto bs=rs8ConstantMultiply<Affine,b>(sum);
    const auto cs=rs8ConstantMultiply<Affine,c>(sum);
    const auto tail=_mm512_ternarylogic_epi32(rs8ConstantMultiply<Affine,b>(x2),
        rs8ConstantMultiply<Affine,c>(x1),rs8ConstantMultiply<Affine,d>(sum),0x96);
    x0=_mm512_xor_si512(_mm512_ternarylogic_epi32(x0,bs,cs,0x96),tail);
    x1=_mm512_ternarylogic_epi32(x1,bs,tail,0x96);
    x2=_mm512_ternarylogic_epi32(x2,cs,tail,0x96);
    x3=tail;
}

// The corresponding general-input identity has
// T=b*(x2+x3)+c*(x1+x3)+d*(x0+x1+x2+x3), and last output x3+T.
// Five GFNI and ten Boolean sites, versus five and twelve in the generic
// subset source circuit. Both maps rely on the already-checked row sum1.
template<bool Affine=true>
RS4_MORE_INLINE void rs4CompactInjection(__m512i& x0,__m512i& x1,
    __m512i& x2,__m512i& x3) noexcept {
    constexpr auto b=rs4ZetaCoefficient<false,1>();
    constexpr auto c=rs4ZetaCoefficient<false,2>();
    constexpr auto d=rs4ZetaCoefficient<false,3>();
    const auto sum=_mm512_xor_si512(_mm512_ternarylogic_epi32(x0,x1,x2,0x96),x3);
    const auto bs=rs8ConstantMultiply<Affine,b>(sum);
    const auto cs=rs8ConstantMultiply<Affine,c>(sum);
    const auto tail=_mm512_ternarylogic_epi32(
        rs8ConstantMultiply<Affine,b>(_mm512_xor_si512(x2,x3)),
        rs8ConstantMultiply<Affine,c>(_mm512_xor_si512(x1,x3)),
        rs8ConstantMultiply<Affine,d>(sum),0x96);
    x0=_mm512_xor_si512(_mm512_ternarylogic_epi32(x0,bs,cs,0x96),tail);
    x1=_mm512_ternarylogic_epi32(x1,bs,tail,0x96);
    x2=_mm512_ternarylogic_epi32(x2,cs,tail,0x96);
    x3=_mm512_xor_si512(x3,tail);
}

template<bool Compact>
RS4_MORE_INLINE void rs4MorePrepareFour(const Block* input,
    __m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3) noexcept {
    if constexpr(Compact){
        z0=_mm512_loadu_si512(input);z1=_mm512_loadu_si512(input+4);
        z2=_mm512_loadu_si512(input+8);
        rs4CompactEmission<true>(z0,z1,z2,z3);
        c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    }else rs4PrepareEmissionFour<true>(input,z0,z1,z2,z3);
}

template<bool Compact>
RS4_MORE_INLINE void rs4MoreFinishFour(__m512i z0,__m512i z1,__m512i z2,__m512i z3,
    Block* output,__m512i& history,const std::uint8_t* diagonal) noexcept {
    if constexpr(Compact){
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal)));
        rs4HistoryDiagonal<true>(z3,history,sf16ExpandRow<3>(packed));
        rs4HistoryDiagonal<true>(z2,history,sf16ExpandRow<2>(packed));
        rs4HistoryDiagonal<true>(z1,history,sf16ExpandRow<1>(packed));
        rs4HistoryDiagonal<true>(z0,history,sf16ExpandRow<0>(packed));
        c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
        rs4CompactInjection<true>(z0,z1,z2,z3);
        _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
        _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
    }else rs4FinishStoreFour<true,true>(z0,z1,z2,z3,output,history,diagonal);
}

// Same interface and exact schedule as the selected16-packet checkpoint.
// Only the two constant maps change their Boolean circuits. No setup and
// no extra coefficient or state storage;72 GFNI sites remain unchanged.
RS4_MORE_INLINE void rs4InnerTransposeSixteenCompact(const Block* input,Block* output,
    __m512i& history,const std::uint8_t* diagonal) noexcept {
    __m512i z0,z1,z2,z3,w0,w1,w2,w3;
    rs4MorePrepareFour<true>(input+36,w0,w1,w2,w3);
    rs4MorePrepareFour<true>(input+24,z0,z1,z2,z3);
    rs4MoreFinishFour<true>(w0,w1,w2,w3,output+48,history,diagonal+48);
    rs4MorePrepareFour<true>(input+12,w0,w1,w2,w3);
    rs4MoreFinishFour<true>(z0,z1,z2,z3,output+32,history,diagonal+32);
    rs4MorePrepareFour<true>(input,z0,z1,z2,z3);
    rs4MoreFinishFour<true>(w0,w1,w2,w3,output+16,history,diagonal+16);
    rs4MoreFinishFour<true>(z0,z1,z2,z3,output,history,diagonal);
}

// Whole-region entry point for the compact-map candidate. packetCount must
// be a positive multiple of16. Input has3*packetCount retained records;
// output and diagonal have4*packetCount records/bytes, in ordinary folded
// four-packet tiles. History retains its original representation.
inline void rs4InnerTransposeCompact(const Block* input,Block* output,
    __m512i& history,const std::uint8_t* diagonal,std::size_t packetCount) noexcept {
    for(auto packets=packetCount;packets;){
        packets-=16;
        rs4InnerTransposeSixteenCompact(input+3*packets,output+4*packets,history,diagonal+4*packets);
    }
}

// Continuous two-tile lookahead across16-packet body boundaries. Same
// arithmetic and live-vector budget as the selected checkpoint. Each loop
// body is explicitly unrolled for16 packets and leaves the next two tiles
// prepared, avoiding the checkpoint's drain/refill at each body boundary.
// packetCount must be a positive multiple of16; there is no padded access.
// Compact=false isolates scheduling; true combines the compact maps.
template<bool Compact=false>
inline void rs4InnerTransposeRolling(const Block* input,Block* output,
    __m512i& history,const std::uint8_t* diagonal,std::size_t packetCount) noexcept {
    auto records=4*packetCount;
    __m512i z0,z1,z2,z3,w0,w1,w2,w3;
    rs4MorePrepareFour<Compact>(input+3*((records-16)/4),w0,w1,w2,w3);
    rs4MorePrepareFour<Compact>(input+3*((records-32)/4),z0,z1,z2,z3);
    for(;records>64;records-=64){
        rs4MoreFinishFour<Compact>(w0,w1,w2,w3,output+records-16,history,diagonal+records-16);
        rs4MorePrepareFour<Compact>(input+3*((records-48)/4),w0,w1,w2,w3);
        rs4MoreFinishFour<Compact>(z0,z1,z2,z3,output+records-32,history,diagonal+records-32);
        rs4MorePrepareFour<Compact>(input+3*((records-64)/4),z0,z1,z2,z3);
        rs4MoreFinishFour<Compact>(w0,w1,w2,w3,output+records-48,history,diagonal+records-48);
        rs4MorePrepareFour<Compact>(input+3*((records-80)/4),w0,w1,w2,w3);
        rs4MoreFinishFour<Compact>(z0,z1,z2,z3,output+records-64,history,diagonal+records-64);
        rs4MorePrepareFour<Compact>(input+3*((records-96)/4),z0,z1,z2,z3);
    }
    rs4MoreFinishFour<Compact>(w0,w1,w2,w3,output+48,history,diagonal+48);
    rs4MorePrepareFour<Compact>(input+12,w0,w1,w2,w3);
    rs4MoreFinishFour<Compact>(z0,z1,z2,z3,output+32,history,diagonal+32);
    rs4MorePrepareFour<Compact>(input,z0,z1,z2,z3);
    rs4MoreFinishFour<Compact>(w0,w1,w2,w3,output+16,history,diagonal+16);
    rs4MoreFinishFour<Compact>(z0,z1,z2,z3,output,history,diagonal);
}

#undef RS4_MORE_INLINE
} // namespace spin::research::char2
