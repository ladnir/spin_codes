#include "K18Certified22.h"
#include "K18Map.h"
#include <bit>
#include <cstring>

namespace spin::research::k18cert22 {
namespace {
template<bool Forward>
void scalar(const Block* input,Block* output,const rs::Plan& plan,const Tables& tables) {
    // Deliberately literal certified coordinates, independent of the fast
    // monomial emission circuit, feedback transform, and GFNI update tables.
    alignas(16) __m128i state[22]{},next[22],feedback[22];
    const auto epochs=plan.n/64;
    for(std::size_t step=0; step<epochs; ++step) {
        const auto epoch=Forward?step:epochs-1-step;
        for(auto& value:feedback) value=_mm_setzero_si128();
        for(unsigned p=0; p<64; ++p) {
            const auto index=64*epoch+p;
            const auto routed=plan.route[index/4]+(index&3);
            const auto raw=input[Forward?routed:index].mData;
            auto value=raw;
            for(auto mask=map::oldColumns[p]; mask; mask&=mask-1) {
                const auto j=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            output[Forward?index:routed]=Block(value);
        }
        if(step+1==epochs) break;
        const auto& matrix=tables.reverseMatrices[epoch];
        for(unsigned j=0; j<22; ++j) {
            auto value=feedback[j];
            if constexpr(Forward) {
                for(unsigned c=0; c<22; ++c)
                    if((matrix[c]>>j)&1U) value=_mm_xor_si128(value,state[c]);
            } else {
                for(auto mask=matrix[j]; mask; mask&=mask-1)
                    value=_mm_xor_si128(value,state[std::countr_zero(mask)]);
            }
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
void reverseRouteScalar(const Block* input,Block* scratch,const rs::Plan& plan,const Tables& tables) {
    scalar<false>(input,scratch,plan,tables);
}
void forwardInnerScalar(const Block* scratch,Block* output,const rs::Plan& plan,const Tables& tables) {
    scalar<true>(scratch,output,plan,tables);
}
}
