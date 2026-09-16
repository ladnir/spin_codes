#pragma once
#include "../asymmetric/AsymmetricInner.h"
namespace bare_spin {
template<class Map,bool First,std::size_t R,class Emit> OC_FORCEINLINE void macroPoint(
    const block* in,__m128i* folded,const __m128i* state,std::size_t base,Emit& emit) {
    constexpr unsigned p=Map::T-1-R;
    const auto input=in[p].mData;
    if constexpr(p>=Map::T-128) folded[p%128]=input;
    else folded[p%128]=_mm_xor_si128(folded[p%128],input);
    if constexpr(First) emit(base+p,block(input));
    else emit(base+p,block(_mm_xor_si128(input,sparseFixedSum<Map::columns[p]>(state))));
}
template<class Map,bool First,class Emit,std::size_t... R> OC_FORCEINLINE void macroPoints(
    const block* in,__m128i* folded,const __m128i* state,std::size_t base,Emit& emit,std::index_sequence<R...>) {
    (macroPoint<Map,First,R>(in,folded,state,base,emit),...);
}
template<class Map,bool Quarter,class Emit> OC_FORCEINLINE void macroReverse(
    const block* in,std::size_t n,const u32* masks,Emit&& emit) {
    if constexpr(!Quarter) innerReverse<Map>(in,n,masks,emit);
    else {
        alignas(32) __m128i state[19]{},syndrome[19],folded[128];
        const auto epochs=n/Map::T;
        for(std::size_t e=epochs;e-->0;) {
            if(e+1==epochs) macroPoints<Map,true>(in+e*Map::T,folded,state,e*Map::T,emit,std::make_index_sequence<Map::T>{});
            else macroPoints<Map,false>(in+e*Map::T,folded,state,e*Map::T,emit,std::make_index_sequence<Map::T>{});
            if(e==0) break;
            zeta<128>(folded);Map::finish(folded,syndrome);
            for(unsigned r=Map::Rounds;r-->0;) {
                auto u=masks[2*(e*Map::Rounds+r)],v=masks[2*(e*Map::Rounds+r)+1];
                auto dot=_mm_setzero_si128();
                while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
                while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
            }
            for(unsigned j=0;j<19;++j) state[j]=_mm_xor_si128(state[j],syndrome[j]);
        }
    }
}
}
