#include "K16T128.h"
#include "T128Map.h"
#include <bit>
#include <cstring>
#include <stdexcept>

namespace spin::research::k16codesign {
namespace {
using namespace detail::packet::large;
struct CachedRoute128 {
    Block* scratch; const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+bases[packet]),value);
    }
};
}

void prepareT128(const rs::Plan& plan,T128Tables& tables) {
    if(plan.n%128 || plan.reverseMatrices.size()<plan.n/128)
        throw std::invalid_argument("t128 requires whole 128-output steps and GL16 rows");
    tables.updates.resize(plan.n/128);
    for(std::size_t epoch=0;epoch<plan.n/128;++epoch)
        for(unsigned out=0;out<2;++out)
            for(unsigned in=0;in<2;++in) {
                std::uint64_t value=0;
                for(unsigned j=0;j<8;++j)
                    value|=std::uint64_t((plan.reverseMatrices[epoch][8*out+j]>>(8*in))&255U)<<(8*(7-j));
                tables.updates[epoch].matrix[2*out+in]=value;
            }
}

void reverseRouteT128(const Block* input,Block* scratch,const rs::Plan& plan,const T128Tables& tables) {
    CachedRoute128 emit{scratch,plan.route.data()};
    alignas(64) __m128i words[16],moments[128],syndrome[16];
    PackedState state{},feedback;
    for(std::size_t epoch=plan.n/128;epoch-->0;) {
        wideUnpackVbmi(state,words);
        t128::step(words,input+128*epoch,32*epoch,emit,moments);
        if(!epoch)break;
        t128::finish128(moments,syndrome);
        widePackVbmi(syndrome,feedback);
        denseStep16(state,tables.updates[epoch],&feedback);
    }
}

void reverseRouteT128Scalar(const Block* input,Block* scratch,const rs::Plan& plan) {
    alignas(16) __m128i state[16]{},next[16],feedback[16];
    for(std::size_t epoch=plan.n/128;epoch-->0;) {
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<128;++p) {
            const auto i=128*epoch+p;
            const auto raw=input[i].mData;
            auto value=raw;
            for(unsigned mask=t128::columns[p];mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            scratch[plan.route[i/4]+(i&3)]=Block(value);
        }
        if(!epoch)break;
        for(unsigned j=0;j<16;++j) {
            auto value=feedback[j];
            for(unsigned mask=plan.reverseMatrices[epoch][j];mask;mask&=mask-1)
                value=_mm_xor_si128(value,state[std::countr_zero(mask)]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}

void forwardInnerT128Scalar(const Block* routed,Block* encoded,const rs::Plan& plan) {
    alignas(16) __m128i state[16]{},next[16],feedback[16];
    for(std::size_t epoch=0;epoch<plan.n/128;++epoch) {
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<128;++p) {
            const auto i=128*epoch+p;
            const auto raw=routed[plan.route[i/4]+(i&3)].mData;
            auto value=raw;
            for(unsigned mask=t128::columns[p];mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            encoded[i]=Block(value);
        }
        if(epoch+1==plan.n/128)break;
        for(unsigned j=0;j<16;++j) {
            auto value=feedback[j];
            for(unsigned c=0;c<16;++c)
                if((plan.reverseMatrices[epoch][c]>>j)&1)
                    value=_mm_xor_si128(value,state[c]);
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
