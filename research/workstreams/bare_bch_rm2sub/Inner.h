#pragma once
#include "Spin.h"
#include "generated/SelectedMaps.h"
#include <array>
#include <bit>
#include <cstring>
#include <immintrin.h>
#include <utility>

namespace bare_spin {
#ifndef SPIN_GROUPED_A
#define SPIN_GROUPED_A 1
#endif
template<class Map> inline u32 fieldMultiply(u32 a,u32 b) noexcept {
    u32 r=0;
    for(unsigned i=0;i<Map::S;++i) {
        r^=a & (0U-(b&1)); b>>=1;
        const u32 carry=a>>(Map::S-1);
        a=((a<<1)&((1U<<Map::S)-1)) ^ (Map::modulusLow & (0U-carry));
    }
    return r;
}

// Fixed small tables, no per-epoch allocations or callback type erasure.
template<unsigned Count> OC_FORCEINLINE void nibble(const __m128i* in,__m128i* t) {
    t[0]=_mm_setzero_si128(); t[1]=in[0];
    if constexpr(Count>=2) { t[2]=in[1]; t[3]=_mm_xor_si128(t[1],t[2]); }
    if constexpr(Count>=3) {
        t[4]=in[2]; t[5]=_mm_xor_si128(t[4],t[1]);
        t[6]=_mm_xor_si128(t[4],t[2]); t[7]=_mm_xor_si128(t[4],t[3]);
    }
    if constexpr(Count==4) {
        t[8]=in[3];
        t[9]=_mm_xor_si128(t[8],t[1]); t[10]=_mm_xor_si128(t[8],t[2]);
        t[11]=_mm_xor_si128(t[8],t[3]); t[12]=_mm_xor_si128(t[8],t[4]);
        t[13]=_mm_xor_si128(t[8],t[5]); t[14]=_mm_xor_si128(t[8],t[6]);
        t[15]=_mm_xor_si128(t[8],t[7]);
    }
}
template<unsigned S,unsigned G=0> OC_FORCEINLINE void tables(const __m128i* state,__m128i t[][16]) {
    constexpr unsigned count=(S-4*G<4?S-4*G:4);
    nibble<count>(state+4*G,t[G]);
    if constexpr(4*(G+1)<S) tables<S,G+1>(state,t);
}
template<u32 Column,unsigned G=0> OC_FORCEINLINE __m128i fixedSum(const __m128i t[][16]) {
    constexpr unsigned index=Column&15;
    if constexpr((Column>>4)==0) return t[G][index];
    else if constexpr(index==0) return fixedSum<(Column>>4),G+1>(t);
    else return _mm_xor_si128(t[G][index],fixedSum<(Column>>4),G+1>(t));
}
template<unsigned S,unsigned G=0> OC_FORCEINLINE __m128i variableSum(const __m128i t[][16],u32 mask) {
    constexpr unsigned bits=(S-4*G<4?S-4*G:4);
    auto x=t[G][mask&((1U<<bits)-1)];
    if constexpr(4*(G+1)<S) x=_mm_xor_si128(x,variableSum<S,G+1>(t,mask>>4));
    return x;
}

// High-to-low zeta stages. Only branches needed by degree <=2 monomials remain.
// All indices and branch tests are compile-time constants, retaining SIMD pairs.
template<unsigned T,unsigned D,unsigned Base=0> OC_FORCEINLINE void zetaStage(__m128i* z) {
    if constexpr(std::popcount(Base/(2*D))<=2) {
        if constexpr(D==1) z[Base]=_mm_xor_si128(z[Base],z[Base+1]);
        else for(unsigned j=0;j<D;j+=2) {
            const auto a=_mm256_load_si256(reinterpret_cast<const __m256i*>(z+Base+j));
            const auto b=_mm256_load_si256(reinterpret_cast<const __m256i*>(z+Base+j+D));
            _mm256_store_si256(reinterpret_cast<__m256i*>(z+Base+j),_mm256_xor_si256(a,b));
        }
    }
    if constexpr(Base+2*D<T) zetaStage<T,D,Base+2*D>(z);
}
template<unsigned T,unsigned D=T/2> OC_FORCEINLINE void zeta(__m128i* z) {
    zetaStage<T,D>(z);
    if constexpr(D>1) zeta<T,D/2>(z);
}
template<class Map,std::size_t R,class Emit> OC_FORCEINLINE void emitPoint(
    const block* in,__m128i* values,const __m128i table[][16],std::size_t base,Emit& emit) {
    constexpr unsigned p=Map::T-1-R;
    constexpr auto column=(SPIN_GROUPED_A && Map::S==19)?Map::groupedColumns[p]:Map::columns[p];
    const auto v=_mm_xor_si128(in[p].mData,fixedSum<column>(table));
    values[p]=v; emit(base+p,block(v));
}
template<class Map,class Emit,std::size_t... R> OC_FORCEINLINE void emitPoints(
    const block* in,__m128i* values,const __m128i table[][16],std::size_t base,Emit& emit,std::index_sequence<R...>) {
    (emitPoint<Map,R>(in,values,table,base,emit),...);
}
template<class Map,class Emit> OC_FORCEINLINE void innerReverse(
    const block* in,std::size_t n,const u32* fieldRows,Emit&& emit) {
    alignas(32) __m128i state[Map::S]{};
    alignas(32) __m128i syndrome[Map::S],grouped[Map::S],values[Map::T],table[(Map::S+3)/4][16];
    const auto epochs=n/Map::T;
    for(std::size_t epoch=epochs;epoch-->0;) {
        const auto base=epoch*Map::T;
        if(epoch+1==epochs) {
            for(unsigned p=Map::T;p-->0;) { values[p]=in[base+p].mData; emit(base+p,in[base+p]); }
        } else {
            if constexpr(SPIN_GROUPED_A && Map::S==19) {
                for(unsigned j=0;j<Map::S;++j) grouped[j]=state[Map::groupOrder[j]];
                tables<Map::S>(grouped,table);
            } else tables<Map::S>(state,table);
            emitPoints<Map>(in+base,values,table,base,emit,std::make_index_sequence<Map::T>{});
        }
        if(epoch==0) break;
        zeta<Map::T>(values); Map::finish(values,syndrome);
        if(epoch+1==epochs) std::memcpy(state,syndrome,sizeof(state));
        else {
            tables<Map::S>(state,table);
            for(unsigned j=0;j<Map::S;++j)
                state[j]=_mm_xor_si128(variableSum<Map::S>(table,fieldRows[epoch*Map::S+j]),syndrome[j]);
        }
    }
}
}
