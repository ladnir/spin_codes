#pragma once
#include "Inner.h"

namespace bare_spin {
template<unsigned S,std::size_t... J> OC_FORCEINLINE void maskedRankOne(
    __m128i* state,const __m128i* syndrome,__m128i dot,u32 v,std::index_sequence<J...>) {
    ((state[J]=_mm_xor_si128(_mm_xor_si128(state[J],syndrome[J]),
         _mm_and_si128(dot,_mm_set1_epi32(-int((v>>J)&1))))),...);
}

// The existing emission table already contains all four-bit state sums.
// Setup permutes u into that table's grouped basis; v stays in the state basis.
// No allocation, per-block dispatch, or reconstruction of a dense state matrix.
template<class Map,bool Quarter,bool Masked,class Emit> OC_FORCEINLINE void mixerReverse(
    const block* in,std::size_t n,const u32* masks,Emit&& emit) {
    if constexpr(!Quarter) {
        innerReverse<Map>(in,n,masks,emit);
    } else {
        static_assert(Map::T==128 && Map::S==19);
        alignas(32) __m128i state[Map::S]{};
        alignas(32) __m128i syndrome[Map::S],grouped[Map::S],values[Map::T],table[(Map::S+3)/4][16];
        const auto epochs=n/Map::T;
        for(std::size_t epoch=epochs;epoch-->0;) {
            const auto base=epoch*Map::T;
            if(epoch+1==epochs) {
                for(unsigned p=Map::T;p-->0;) { values[p]=in[base+p].mData; emit(base+p,in[base+p]); }
            } else {
                if constexpr(SPIN_GROUPED_A) {
                    for(unsigned j=0;j<Map::S;++j) grouped[j]=state[Map::groupOrder[j]];
                    tables<Map::S>(grouped,table);
                } else tables<Map::S>(state,table);
                emitPoints<Map>(in+base,values,table,base,emit,std::make_index_sequence<Map::T>{});
            }
            if(epoch==0) break;
            zeta<Map::T>(values);Map::finish(values,syndrome);
            if(epoch+1==epochs) std::memcpy(state,syndrome,sizeof(state));
            else {
                const auto dot=variableSum<Map::S>(table,masks[2*epoch]);
                auto v=masks[2*epoch+1];
                if constexpr(Masked) {
                    maskedRankOne<Map::S>(state,syndrome,dot,v,std::make_index_sequence<Map::S>{});
                } else {
                    while(v) {
                        const auto j=std::countr_zero(v);v&=v-1;
                        state[j]=_mm_xor_si128(state[j],dot);
                    }
                    for(unsigned j=0;j<Map::S;++j) state[j]=_mm_xor_si128(state[j],syndrome[j]);
                }
            }
        }
    }
}
}
