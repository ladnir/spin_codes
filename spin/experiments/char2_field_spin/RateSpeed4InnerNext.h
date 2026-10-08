#pragma once
#include "RateSpeed4Inner.h"

namespace spin::research::char2 {

#if defined(_MSC_VER)
#define RS4_NEXT_INLINE __forceinline
#else
#define RS4_NEXT_INLINE inline __attribute__((always_inline))
#endif

// All routines in this header preserve the RateSpeed4 sampled linear map.
// The gauge variant changes the representation of the carried dual state;
// the pair-lookahead variant keeps the original state representation.
constexpr std::uint8_t rs4NextInverse(std::uint8_t x) noexcept {
    std::uint8_t result=1;
    for(unsigned power=254;power;power>>=1){
        if(power&1)result=sfMultiply(result,x);
        x=sfMultiply(x,x);
    }
    return result;
}
inline constexpr auto rs4NextInverses=[]() constexpr {
    std::array<std::uint8_t,256> table{};
    for(unsigned i=1;i!=256;++i)table[i]=rs4NextInverse(static_cast<std::uint8_t>(i));
    return table;
}();

// Setup only. records must be a multiple of16. Allocate2*records coefficient
// bytes; for each four-packet tile its32 bytes are:
//   [inverse P_(t+1), coordinate-major16][P_(t+1), coordinate-major16].
// Here P_t[j]=d_t[j]*P_(t+1)[j], with P_end[j]=1. P must span the entire
// chronological chain, including region boundaries where history survives.
// The kernel carries H_t[j]=h_t[j]/P_t[j], starting with the usual zero.
// The returned P_0 converts the final carried state back to ordinary h_0.
inline std::array<std::uint8_t,4> rs4PrepareGauge(
    const std::uint8_t* diagonal,std::size_t records,std::uint8_t* coefficients){
    if(records%16)throw std::invalid_argument("RateSpeed4 gauge tile geometry");
    std::array<std::uint8_t,4> suffix{1,1,1,1},inverseSuffix{1,1,1,1};
    for(std::size_t step=records/4;step--;){
        const auto tile=step/4,phase=step%4;
        for(unsigned coordinate=0;coordinate!=4;++coordinate){
            const auto d=diagonal[4*step+coordinate];
            if(!d)throw std::invalid_argument("RateSpeed4 zero state diagonal");
            coefficients[32*tile+4*coordinate+phase]=inverseSuffix[coordinate];
            coefficients[32*tile+16+4*coordinate+phase]=suffix[coordinate];
            suffix[coordinate]=sfMultiply(d,suffix[coordinate]);
            inverseSuffix[coordinate]=sfMultiply(rs4NextInverses[d],inverseSuffix[coordinate]);
        }
    }
    return suffix;
}

// Each128-bit record lane is one chronological packet. Return
// [x0+x1+x2+x3, x1+x2+x3, x2+x3, x3] in two shuffle/XOR stages.
RS4_NEXT_INLINE __m512i rs4SuffixFour(__m512i x) noexcept {
    x=_mm512_xor_si512(x,_mm512_maskz_shuffle_i32x4(0x00ff,x,x,0xee));
    return _mm512_xor_si512(x,_mm512_maskz_shuffle_i32x4(0x0fff,x,x,0x39));
}
RS4_NEXT_INLINE __m512i rs4PackFirstRecords(
    __m512i x0,__m512i x1,__m512i x2,__m512i x3) noexcept {
    const auto low=_mm512_shuffle_i32x4(x0,x1,0x00);
    const auto high=_mm512_shuffle_i32x4(x2,x3,0x00);
    return _mm512_shuffle_i32x4(low,high,0x88);
}

// m=A*y, u=m/P_(t+1). The inclusive suffix plus incoming H is H_t.
// The required pre-injection moment is
//   h_(t+1)+107*m = P_(t+1)*H_t + (107+1)*m.
// This form lets the affine constant multiplies execute independently of
// the normalization/suffix path. The carried-history dependence path is
// one XOR per tile: suffix packing is independent of the incoming history.
template<bool AffineScale>
RS4_NEXT_INLINE void rs4GaugeFour(__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3,
    __m512i& history,const std::uint8_t* coefficients) noexcept {
    const auto inverse=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(coefficients)));
    const auto direct=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(coefficients+16)));
    const auto u0=rs4SuffixFour(_mm512_gf2p8mul_epi8(z0,sf16ExpandRow<0>(inverse)));
    const auto u1=rs4SuffixFour(_mm512_gf2p8mul_epi8(z1,sf16ExpandRow<1>(inverse)));
    const auto u2=rs4SuffixFour(_mm512_gf2p8mul_epi8(z2,sf16ExpandRow<2>(inverse)));
    const auto u3=rs4SuffixFour(_mm512_gf2p8mul_epi8(z3,sf16ExpandRow<3>(inverse)));
    const auto h0=_mm512_xor_si512(u0,_mm512_shuffle_i32x4(history,history,0x00));
    const auto h1=_mm512_xor_si512(u1,_mm512_shuffle_i32x4(history,history,0x55));
    const auto h2=_mm512_xor_si512(u2,_mm512_shuffle_i32x4(history,history,0xaa));
    const auto h3=_mm512_xor_si512(u3,_mm512_shuffle_i32x4(history,history,0xff));
    history=_mm512_xor_si512(history,rs4PackFirstRecords(u0,u1,u2,u3));
    z0=_mm512_xor_si512(_mm512_gf2p8mul_epi8(h0,sf16ExpandRow<0>(direct)),rs8ConstantMultiply<AffineScale,rs4Scale^1>(z0));
    z1=_mm512_xor_si512(_mm512_gf2p8mul_epi8(h1,sf16ExpandRow<1>(direct)),rs8ConstantMultiply<AffineScale,rs4Scale^1>(z1));
    z2=_mm512_xor_si512(_mm512_gf2p8mul_epi8(h2,sf16ExpandRow<2>(direct)),rs8ConstantMultiply<AffineScale,rs4Scale^1>(z2));
    z3=_mm512_xor_si512(_mm512_gf2p8mul_epi8(h3,sf16ExpandRow<3>(direct)),rs8ConstantMultiply<AffineScale,rs4Scale^1>(z3));
}

// input24/output32 have exactly the ordinary folded8-packet layout.
// coefficients points to64 bytes produced by rs4PrepareGauge, for these
// two tiles. history is normalized H, not ordinary h. Forty-four GFNI sites
// per32 raw records, versus36 baseline; no serial GFNI history chain and
// no coordinate/packet record transposes. Fifteen lane shuffles per tile
// replace16 transpose shuffles, with eight coefficient-expansion shuffles
// versus the baseline four. The change trades throughput for latency.
template<bool AffineMaps=true,bool AffineScale=true>
RS4_NEXT_INLINE void rs4InnerTransposeEightGauge(const Block* input,Block* output,
    __m512i& history,const std::uint8_t* coefficients) noexcept {
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_setzero_si512();
    auto w0=_mm512_loadu_si512(input+12),w1=_mm512_loadu_si512(input+16);
    auto w2=_mm512_loadu_si512(input+20),w3=_mm512_setzero_si512();
    rs4Circulant<true,AffineMaps>(z0,z1,z2,z3);
    rs4Circulant<true,AffineMaps>(w0,w1,w2,w3);
    rs4GaugeFour<AffineScale>(w0,w1,w2,w3,history,coefficients+32);
    rs4Circulant<false,AffineMaps>(w0,w1,w2,w3);
    rs4GaugeFour<AffineScale>(z0,z1,z2,z3,history,coefficients);
    rs4Circulant<false,AffineMaps>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
    _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
    _mm512_storeu_si512(output+16,w0);_mm512_storeu_si512(output+20,w1);
    _mm512_storeu_si512(output+24,w2);_mm512_storeu_si512(output+28,w3);
}

// Two reverse steps, with dLow*dHigh already formed. The next pair can
// start after one variable GFNI rather than two. Low's independent product
// is deliberately separated from the history-dependent product.
template<bool AffineScale>
RS4_NEXT_INLINE void rs4HistoryPair(__m512i& low,__m512i& high,__m512i& history,
    __m512i dLow,__m512i dHigh,__m512i dPair) noexcept {
    const auto old=history;
    const auto highInput=_mm512_xor_si512(old,high);
    const auto lowInjection=_mm512_gf2p8mul_epi8(low,dLow);
    const auto middle=_mm512_gf2p8mul_epi8(highInput,dHigh);
    history=_mm512_xor_si512(_mm512_gf2p8mul_epi8(highInput,dPair),lowInjection);
    high=_mm512_xor_si512(old,rs8ConstantMultiply<AffineScale,rs4Scale>(high));
    low=_mm512_xor_si512(middle,rs8ConstantMultiply<AffineScale,rs4Scale>(low));
}
RS4_NEXT_INLINE __m512i rs4AdjacentDiagonalProducts(__m512i packed) noexcept {
    const auto select=_mm512_setr_epi32(0x07060504,0x03020100,0x0f0e0d0c,0x0b0a0908,
        0x07060504,0x03020100,0x0f0e0d0c,0x0b0a0908,
        0x07060504,0x03020100,0x0f0e0d0c,0x0b0a0908,
        0x07060504,0x03020100,0x0f0e0d0c,0x0b0a0908);
    return _mm512_gf2p8mul_epi8(packed,_mm512_shuffle_epi8(packed,select));
}

// Optional pair setup. Allocate48 bytes per32 raw records. Each chunk has
// the original chronological32 diagonal bytes followed by16 products:
// d0*d1, d2*d3, d4*d5, d6*d7, each containing four state coordinates.
inline void rs4PreparePairs(const std::uint8_t* diagonal,std::size_t records,
    std::uint8_t* coefficients){
    if(records%32)throw std::invalid_argument("RateSpeed4 pair tile geometry");
    for(std::size_t base=0;base!=records;base+=32){
        auto* target=coefficients+(base/32)*48;
        std::memcpy(target,diagonal+base,32);
        for(unsigned pair=0;pair!=4;++pair)for(unsigned j=0;j!=4;++j)
            target[32+4*pair+j]=sfMultiply(diagonal[base+8*pair+j],diagonal[base+8*pair+4+j]);
    }
}

// Precomputed=false accepts the original32 diagonal bytes and ordinary
// history; it is a direct replacement for rs4InnerTransposeEight.
// Precomputed=true accepts48 coefficient bytes from rs4PreparePairs.
// The arithmetic counts are42 GFNI sites (raw) and40 (precomputed), versus
//36 baseline. The serial history path has four variable GFNI, versus eight.
template<bool AffineMaps=true,bool AffineScale=true,bool Precomputed=false>
RS4_NEXT_INLINE void rs4InnerTransposeEightPairs(const Block* input,Block* output,
    __m512i& history,const std::uint8_t* diagonal) noexcept {
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_setzero_si512();
    auto w0=_mm512_loadu_si512(input+12),w1=_mm512_loadu_si512(input+16);
    auto w2=_mm512_loadu_si512(input+20),w3=_mm512_setzero_si512();
    rs4Circulant<true,AffineMaps>(z0,z1,z2,z3);
    rs4Circulant<true,AffineMaps>(w0,w1,w2,w3);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(w0,w1,w2,w3,w0,w1,w2,w3);
    const auto early=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal)));
    const auto late=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal+16)));
    __m512i earlyPairs,latePairs;
    if constexpr(Precomputed){
        earlyPairs=latePairs=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal+32)));
    }else{
        earlyPairs=rs4AdjacentDiagonalProducts(early);
        latePairs=rs4AdjacentDiagonalProducts(late);
    }
    rs4HistoryPair<AffineScale>(w2,w3,history,sf16ExpandRow<2>(late),sf16ExpandRow<3>(late),
        sf16ExpandRow<Precomputed?3:2>(latePairs));
    rs4HistoryPair<AffineScale>(w0,w1,history,sf16ExpandRow<0>(late),sf16ExpandRow<1>(late),
        sf16ExpandRow<Precomputed?2:0>(latePairs));
    c8TransposeRecords(w0,w1,w2,w3,w0,w1,w2,w3);
    rs4Circulant<false,AffineMaps>(w0,w1,w2,w3);
    rs4HistoryPair<AffineScale>(z2,z3,history,sf16ExpandRow<2>(early),sf16ExpandRow<3>(early),
        sf16ExpandRow<Precomputed?1:2>(earlyPairs));
    rs4HistoryPair<AffineScale>(z0,z1,history,sf16ExpandRow<0>(early),sf16ExpandRow<1>(early),
        sf16ExpandRow<0>(earlyPairs));
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    rs4Circulant<false,AffineMaps>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
    _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
    _mm512_storeu_si512(output+16,w0);_mm512_storeu_si512(output+20,w1);
    _mm512_storeu_si512(output+24,w2);_mm512_storeu_si512(output+28,w3);
}

template<bool AffineMaps>
RS4_NEXT_INLINE void rs4PrepareEmissionFour(const Block* input,
    __m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3) noexcept {
    z0=_mm512_loadu_si512(input);z1=_mm512_loadu_si512(input+4);
    z2=_mm512_loadu_si512(input+8);z3=_mm512_setzero_si512();
    rs4Circulant<true,AffineMaps>(z0,z1,z2,z3);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
}

template<bool AffineMaps,bool AffineScale>
RS4_NEXT_INLINE void rs4FinishStoreFour(__m512i z0,__m512i z1,__m512i z2,__m512i z3,
    Block* output,__m512i& history,const std::uint8_t* diagonal) noexcept {
    const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(diagonal)));
    rs4HistoryDiagonal<AffineScale>(z3,history,sf16ExpandRow<3>(packed));
    rs4HistoryDiagonal<AffineScale>(z2,history,sf16ExpandRow<2>(packed));
    rs4HistoryDiagonal<AffineScale>(z1,history,sf16ExpandRow<1>(packed));
    rs4HistoryDiagonal<AffineScale>(z0,history,sf16ExpandRow<0>(packed));
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    rs4Circulant<false,AffineMaps>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
    _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
}

// Four folded four-packet tiles:48 retained input records,64 raw output
// records,64 ordinary state-diagonal bytes, and ordinary carried history.
// Identical72 GFNI sites to two baseline eight-packet calls. Only two tiles
// are live at once: completed late tiles are stored before their registers
// prepare another early tile. This extends emission lookahead across the
// usual eight-packet boundary without holding all16 packets in registers.
// Most useful when output points directly into the regional temporary and
// Dx is applied during routing. No allocation or separate setup is needed.
template<bool AffineMaps=true,bool AffineScale=true>
RS4_NEXT_INLINE void rs4InnerTransposeSixteenPipeline(const Block* input,Block* output,
    __m512i& history,const std::uint8_t* diagonal) noexcept {
    __m512i z0,z1,z2,z3,w0,w1,w2,w3;
    rs4PrepareEmissionFour<AffineMaps>(input+36,w0,w1,w2,w3);
    rs4PrepareEmissionFour<AffineMaps>(input+24,z0,z1,z2,z3);
    rs4FinishStoreFour<AffineMaps,AffineScale>(w0,w1,w2,w3,output+48,history,diagonal+48);
    rs4PrepareEmissionFour<AffineMaps>(input+12,w0,w1,w2,w3);
    rs4FinishStoreFour<AffineMaps,AffineScale>(z0,z1,z2,z3,output+32,history,diagonal+32);
    rs4PrepareEmissionFour<AffineMaps>(input,z0,z1,z2,z3);
    rs4FinishStoreFour<AffineMaps,AffineScale>(w0,w1,w2,w3,output+16,history,diagonal+16);
    rs4FinishStoreFour<AffineMaps,AffineScale>(z0,z1,z2,z3,output,history,diagonal);
}

#undef RS4_NEXT_INLINE
} // namespace spin::research::char2
