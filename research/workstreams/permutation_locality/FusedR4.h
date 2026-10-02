#pragma once
// Research-only implementation of the ORIGINAL four-transvection map.
// No new randomness or change to the distribution of the encoder.
#include "Inner.h"

namespace spin::research::fused_r4 {
using namespace spin::detail::kernel;

struct ScalarRow {
    unsigned dot[4];
    std::uint8_t select[20];
};
struct WideRow {
    unsigned dot[4];
    std::uint8_t index01[40],index23[40];
};
static_assert(sizeof(ScalarRow)==36 && sizeof(WideRow)==96);

// Setup only. If M=T4 T3 T2 T1 and Ti=I+ui vi^T, then
// M=I+sum_i ai vi^T, where ai=T4 ... T(i+1) ui.
// Reverse/transpose therefore reads ai and writes vi. The forward view below
// is used only by independent adjoint checks, outside the measured encoder.
inline ScalarRow makeScalar(const unsigned* masks,bool transpose=true) {
    unsigned a[4],v[4];
    for(unsigned i=0;i<4;++i) {
        a[i]=masks[2*i];v[i]=masks[2*i+1];
        for(unsigned r=i+1;r<4;++r)
            if(std::popcount(a[i]&masks[2*r+1])&1)a[i]^=masks[2*r];
    }
    ScalarRow result{};
    for(unsigned i=0;i<4;++i)result.dot[i]=transpose?a[i]:v[i];
    for(unsigned j=0;j<19;++j)for(unsigned i=0;i<4;++i)
        result.select[j]|=std::uint8_t((((transpose?v[i]:a[i])>>j)&1)<<i);
    return result;
}
inline WideRow makeWide(const ScalarRow& scalar) {
    WideRow result{};
    for(unsigned i=0;i<4;++i)result.dot[i]=scalar.dot[i];
    for(unsigned j=0;j<20;++j)for(unsigned half=0;half<2;++half) {
        result.index01[2*j+half]=std::uint8_t(2*(scalar.select[j]&3)+half);
        result.index23[2*j+half]=std::uint8_t(2*(scalar.select[j]>>2)+half);
    }
    return result;
}

// Four independent reductions: state is immutable until all four finish.
// Retain the baseline's sparse ctz traversal; no generic callback or allocation.
SPIN_FORCEINLINE __m128i dot(const __m128i* state,unsigned mask) {
    auto value=_mm_setzero_si128();
    while(mask) {
        const auto j=std::countr_zero(mask);mask&=mask-1;
        value=_mm_xor_si128(value,state[j]);
    }
    return value;
}
template<bool Add,unsigned J=0> SPIN_FORCEINLINE void scalarUpdate(
    __m128i* state,const __m128i* table,const std::uint8_t* select,const __m128i* syndrome) {
    if constexpr(Add)state[J]=_mm_ternarylogic_epi64(state[J],table[select[J]],syndrome[J],0x96);
    else state[J]=_mm_xor_si128(state[J],table[select[J]]);
    if constexpr(J+1<19)scalarUpdate<Add,J+1>(state,table,select,syndrome);
}
template<bool Add=false> SPIN_FORCEINLINE void step(
    __m128i* state,const ScalarRow& row,const __m128i* syndrome=nullptr) {
    alignas(32) __m128i table[16];
    table[0]=_mm_setzero_si128();
    table[1]=dot(state,row.dot[0]);
    table[2]=dot(state,row.dot[1]);
    table[4]=dot(state,row.dot[2]);
    table[8]=dot(state,row.dot[3]);
    table[3]=_mm_xor_si128(table[1],table[2]);
    table[5]=_mm_xor_si128(table[1],table[4]);
    table[6]=_mm_xor_si128(table[2],table[4]);
    table[7]=_mm_xor_si128(table[3],table[4]);
    table[9]=_mm_xor_si128(table[1],table[8]);
    table[10]=_mm_xor_si128(table[2],table[8]);
    table[11]=_mm_xor_si128(table[3],table[8]);
    table[12]=_mm_xor_si128(table[4],table[8]);
    table[13]=_mm_xor_si128(table[5],table[8]);
    table[14]=_mm_xor_si128(table[6],table[8]);
    table[15]=_mm_xor_si128(table[7],table[8]);
    scalarUpdate<Add>(state,table,row.select,syndrome);
}
SPIN_FORCEINLINE __m512i pairTable(__m128i d0,__m128i d1) {
    auto table=_mm512_castsi128_si512(_mm_setzero_si128());
    table=_mm512_inserti32x4(table,d0,1);
    table=_mm512_inserti32x4(table,d1,2);
    return _mm512_inserti32x4(table,_mm_xor_si128(d0,d1),3);
}
template<bool Add,unsigned G=0> SPIN_FORCEINLINE void wideUpdate(
    __m128i* state,const WideRow& row,__m512i table01,__m512i table23,const __m128i* syndrome) {
    const auto i01=_mm512_cvtepu8_epi64(_mm_loadl_epi64(
        reinterpret_cast<const __m128i*>(row.index01+8*G)));
    const auto i23=_mm512_cvtepu8_epi64(_mm_loadl_epi64(
        reinterpret_cast<const __m128i*>(row.index23+8*G)));
    const auto a=_mm512_permutexvar_epi64(i01,table01);
    const auto b=_mm512_permutexvar_epi64(i23,table23);
    if constexpr(G==4) {
        // Only three 128-bit state words remain; never touch a twentieth word.
        auto x=_mm512_maskz_loadu_epi64(0x3f,state+16);
        if constexpr(Add)x=_mm512_xor_si512(x,_mm512_maskz_loadu_epi64(0x3f,syndrome+16));
        _mm512_mask_storeu_epi64(state+16,0x3f,_mm512_ternarylogic_epi64(x,a,b,0x96));
    } else {
        auto x=_mm512_loadu_si512(state+4*G);
        if constexpr(Add)x=_mm512_xor_si512(x,_mm512_loadu_si512(syndrome+4*G));
        _mm512_storeu_si512(state+4*G,_mm512_ternarylogic_epi64(x,a,b,0x96));
        wideUpdate<Add,G+1>(state,row,table01,table23,syndrome);
    }
}
template<bool Add=false> SPIN_FORCEINLINE void step(
    __m128i* state,const WideRow& row,const __m128i* syndrome=nullptr) {
    const auto d0=dot(state,row.dot[0]),d1=dot(state,row.dot[1]);
    const auto d2=dot(state,row.dot[2]),d3=dot(state,row.dot[3]);
    wideUpdate<Add>(state,row,pairTable(d0,d1),pairTable(d2,d3),syndrome);
}
}
