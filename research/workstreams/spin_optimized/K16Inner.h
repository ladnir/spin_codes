#pragma once
// Include after Inner.h. Keep the measured reverse kernel's shared emission,
// pruned zeta transform, and sparse transvection update unchanged.
#include "Map64S12.h"
namespace bare_spin {
template<class Emit> OC_FORCEINLINE void inner64Reverse(
    const block* in,std::size_t n,const u32* masks,Emit&& emit) {
    using Map=Map64S12;
    alignas(32) __m128i state[Map::S]{},syndrome[Map::S],raw[Map::T];
    const auto epochs=n/Map::T;
    for(std::size_t epoch=epochs;epoch-->0;) {
        const auto base=epoch*Map::T;
        if(epoch+1==epochs) {
            for(unsigned p=Map::T;p-->0;) {raw[p]=in[base+p].mData;emit(base+p,in[base+p]);}
        } else Map::emitShared(in+base,raw,state,base,emit);
        if(epoch==0) break;
        zeta<Map::T>(raw);Map::finish(raw,syndrome);
        if(epoch+1==epochs) std::memcpy(state,syndrome,sizeof(state));
        else {
            auto u=masks[2*epoch];auto dot=_mm_setzero_si128();
            while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
            auto v=masks[2*epoch+1];
            while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
            for(unsigned j=0;j<Map::S;++j) state[j]=_mm_xor_si128(state[j],syndrome[j]);
        }
    }
}
}
