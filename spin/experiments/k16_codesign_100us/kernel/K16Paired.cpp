#include "K16Paired.h"
#include "T64PairedMap.h"
#include <bit>
#include <cstring>

namespace spin::research::k16codesign {
namespace {
using namespace detail::packet::fast;
struct CachedRoutePaired {
    Block* scratch;const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+bases[packet]),value);
    }
};
template<bool Forward>
static void scalar(const Block* input,Block* output,const rs::Plan& plan) {
    alignas(16) __m128i state[16]{},next[16],feedback[16];
    for(std::size_t step=0;step<plan.n/64;++step) {
        const auto epoch=Forward?step:plan.n/64-1-step;
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<64;++p) {
            const auto i=64*epoch+p,routed=plan.route[i/4]+(i&3);
            const auto raw=input[Forward?routed:i].mData;
            auto value=raw;
            for(unsigned mask=paired::columns[p];mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[j]);feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            output[Forward?i:routed]=Block(value);
        }
        if(step+1==plan.n/64)break;
        for(unsigned j=0;j<16;++j) {
            auto value=feedback[j];
            if constexpr(Forward) {
                for(unsigned c=0;c<16;++c)
                    if((plan.reverseMatrices[epoch][c]>>j)&1U)value=_mm_xor_si128(value,state[c]);
            } else {
                for(unsigned mask=plan.reverseMatrices[epoch][j];mask;mask&=mask-1)
                    value=_mm_xor_si128(value,state[std::countr_zero(mask)]);
            }
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
void preparePaired(const rs::Plan& plan,PairedTables& tables) {
    tables.updates.resize(4*(plan.n/64));
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch)
        for(unsigned out=0;out<2;++out)
            for(unsigned in=0;in<2;++in) {
                std::uint64_t value=0;
                for(unsigned j=0;j<8;++j)
                    value|=std::uint64_t((plan.reverseMatrices[epoch][8*out+j]>>(8*in))&255U)<<(8*j);
                tables.updates[4*epoch+2*out+in]=value;
            }
}
void reverseRoutePaired(const Block* input,Block* scratch,const rs::Plan& plan,const PairedTables& tables) {
    CachedRoutePaired emit{scratch,plan.route.data()};
    alignas(64) __m128i moments[64],syndrome[16];
    alignas(64) __m512i packets[16];
    PackedState before;WordPackets16 state;
    for(std::size_t epoch=plan.n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;
        const bool first=epoch+1==plan.n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
        } else {
            if(epoch){packWordGroup<0>(state,before);packWordGroup<1>(state,before);}
            paired::stream(state,raw,16*epoch,emit,moments);
        }
        if(!epoch)break;
        if(first)packetMoments(packets,moments);
        paired::finishMono(moments,syndrome);
        if(first)wordPackets(syndrome,state);
        else applyPackedToPackets(before,tables.updates.data()+4*epoch,syndrome,state);
    }
}
void reverseRoutePairedScalar(const Block* input,Block* scratch,const rs::Plan& plan) {scalar<false>(input,scratch,plan);}
void forwardInnerPairedScalar(const Block* input,Block* output,const rs::Plan& plan) {scalar<true>(input,output,plan);}
}


