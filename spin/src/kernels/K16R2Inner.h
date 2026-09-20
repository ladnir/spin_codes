#pragma once
#include "ImtRounds.h"
namespace spin::detail::kernel {
template<class Emit> SPIN_FORCEINLINE void inner64R2Reverse(
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
            imtStep(state,masks[4*epoch+2],masks[4*epoch+3]);
            imtStep(state,masks[4*epoch],masks[4*epoch+1]);
            for(unsigned j=0;j<Map::S;++j) state[j]=_mm_xor_si128(state[j],syndrome[j]);
        }
    }
}
}
