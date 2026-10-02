#include "K16T128Border.h"
#include "T128BorderMap.h"
#include <bit>
#include <cstring>

namespace spin::research::k16codesign {
namespace {
using namespace detail::packet::large;
struct ExtraState { __m512i lo,hi; };
template<unsigned Extra>
static SPIN_FORCEINLINE ExtraState packExtra(const __m128i* words) {
    const auto zero=_mm_setzero_si128();
    const auto a=join(words[0],words[1],words[2],Extra==4?words[3]:zero);
    const auto b=_mm512_setzero_si512();
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    return {_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0),
            _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0)};
}
template<unsigned Extra>
static SPIN_FORCEINLINE void unpackExtra(const ExtraState& extra,__m128i* words) {
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto lo=_mm512_gf2p8affine_epi64_epi8(basis,extra.lo,0);
    const auto hi=_mm512_gf2p8affine_epi64_epi8(basis,extra.hi,0);
    const auto joint=_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi);
    words[0]=_mm512_castsi512_si128(joint);
    words[1]=_mm512_extracti32x4_epi32(joint,1);
    words[2]=_mm512_extracti32x4_epi32(joint,2);
    if constexpr(Extra==4)words[3]=_mm512_extracti32x4_epi32(joint,3);
}
static SPIN_FORCEINLINE void update(PackedState& state,ExtraState& extra,
    const PackedState& feedback,const ExtraState& extraFeedback,const std::uint64_t* row) {
    const auto r00=_mm512_set1_epi64(row[0]),r01=_mm512_set1_epi64(row[1]),r02=_mm512_set1_epi64(row[2]);
    const auto n0=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[0],r00,0),_mm512_gf2p8affine_epi64_epi8(state.v[2],r01,0),
        _mm512_gf2p8affine_epi64_epi8(extra.lo,r02,0),0x96),feedback.v[0]);
    const auto n1=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[1],r00,0),_mm512_gf2p8affine_epi64_epi8(state.v[3],r01,0),
        _mm512_gf2p8affine_epi64_epi8(extra.hi,r02,0),0x96),feedback.v[1]);
    const auto r10=_mm512_set1_epi64(row[3]),r11=_mm512_set1_epi64(row[4]),r12=_mm512_set1_epi64(row[5]);
    const auto n2=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[0],r10,0),_mm512_gf2p8affine_epi64_epi8(state.v[2],r11,0),
        _mm512_gf2p8affine_epi64_epi8(extra.lo,r12,0),0x96),feedback.v[2]);
    const auto n3=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[1],r10,0),_mm512_gf2p8affine_epi64_epi8(state.v[3],r11,0),
        _mm512_gf2p8affine_epi64_epi8(extra.hi,r12,0),0x96),feedback.v[3]);
    const auto r20=_mm512_set1_epi64(row[6]),r21=_mm512_set1_epi64(row[7]),r22=_mm512_set1_epi64(row[8]);
    const auto n4=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[0],r20,0),_mm512_gf2p8affine_epi64_epi8(state.v[2],r21,0),
        _mm512_gf2p8affine_epi64_epi8(extra.lo,r22,0),0x96),extraFeedback.lo);
    const auto n5=_mm512_xor_si512(_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(state.v[1],r20,0),_mm512_gf2p8affine_epi64_epi8(state.v[3],r21,0),
        _mm512_gf2p8affine_epi64_epi8(extra.hi,r22,0),0x96),extraFeedback.hi);
    state.v[0]=n0;state.v[1]=n1;state.v[2]=n2;state.v[3]=n3;extra.lo=n4;extra.hi=n5;
}
struct CachedRouteBorder {
    Block* scratch;const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+bases[packet]),value);
    }
};
template<unsigned S>
static SPIN_NOINLINE void reverse(const Block* input,Block* scratch,const rs::Plan& plan,const T128BorderTables& tables) {
    CachedRouteBorder emit{scratch,plan.route.data()};
    alignas(64) __m128i words[20],moments[128],syndrome[20];
    PackedState state{},feedback;ExtraState extra{};
    for(std::size_t epoch=plan.n/128;epoch-->0;) {
        wideUnpackVbmi(state,words);unpackExtra<S-16>(extra,words+16);
        t128border::Map<S>::step(words,input+128*epoch,32*epoch,emit,moments);
        if(!epoch)break;
        t128border::Map<S>::finish128(moments,syndrome);
        widePackVbmi(syndrome,feedback);
        const auto extraFeedback=packExtra<S-16>(syndrome+16);
        update(state,extra,feedback,extraFeedback,tables.packedUpdates[epoch].data());
    }
}
template<unsigned S,bool Forward>
static void scalar(const Block* input,Block* output,const rs::Plan& plan,const T128BorderTables& tables) {
    alignas(16) __m128i state[S]{},next[S],feedback[S];
    for(std::size_t step=0;step<plan.n/128;++step) {
        const auto epoch=Forward?step:plan.n/128-1-step;
        for(auto& value:feedback)value=_mm_setzero_si128();
        for(unsigned p=0;p<128;++p) {
            const auto i=128*epoch+p;
            const auto routed=plan.route[i/4]+(i&3);
            const auto raw=input[Forward?routed:i].mData;
            auto value=raw;
            for(unsigned mask=t128border::Map<S>::columns[p];mask;mask&=mask-1) {
                const auto j=std::countr_zero(mask);
                value=_mm_xor_si128(value,state[j]);
                feedback[j]=_mm_xor_si128(feedback[j],raw);
            }
            output[Forward?i:routed]=Block(value);
        }
        if(step+1==plan.n/128)break;
        for(unsigned j=0;j<S;++j) {
            auto value=feedback[j];
            if constexpr(Forward) {
                for(unsigned c=0;c<S;++c)
                    if((tables.reverseMatrices[epoch][c]>>j)&1U)value=_mm_xor_si128(value,state[c]);
            } else {
                for(unsigned mask=tables.reverseMatrices[epoch][j];mask;mask&=mask-1)
                    value=_mm_xor_si128(value,state[std::countr_zero(mask)]);
            }
            next[j]=value;
        }
        std::memcpy(state,next,sizeof(state));
    }
}
}
void reverseRouteT128Border(const Block* in,Block* scratch,const rs::Plan& plan,const T128BorderTables& tables) {
    if(tables.stateBits==19)reverse<19>(in,scratch,plan,tables);
    else reverse<20>(in,scratch,plan,tables);
}
void reverseRouteT128BorderScalar(const Block* in,Block* scratch,const rs::Plan& plan,const T128BorderTables& tables) {
    if(tables.stateBits==19)scalar<19,false>(in,scratch,plan,tables);
    else scalar<20,false>(in,scratch,plan,tables);
}
void forwardInnerT128BorderScalar(const Block* routed,Block* out,const rs::Plan& plan,const T128BorderTables& tables) {
    if(tables.stateBits==19)scalar<19,true>(routed,out,plan,tables);
    else scalar<20,true>(routed,out,plan,tables);
}
}
