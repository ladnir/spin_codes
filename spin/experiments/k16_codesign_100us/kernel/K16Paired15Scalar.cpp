#include "K16Paired15.h"
#include "T64Paired15Map.h"
#include <bit>
#include <cstring>
namespace spin::research::k16codesign {
namespace {
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
            for(unsigned mask=paired15::columns[p];mask;mask&=mask-1) {
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
void reverseRoutePaired15Scalar(const Block* in,Block* scratch,const rs::Plan& plan){scalar<false>(in,scratch,plan);}
void forwardInnerPaired15Scalar(const Block* in,Block* out,const rs::Plan& plan){scalar<true>(in,out,plan);}
}
