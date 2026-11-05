#pragma once
#include "Gamma30DeferredEncoder.h"

// Exact completion/adjoint alternative on top of the deferred outer. Plan,
// FFT butterflies, field products, routes, stores and prefetches are frozen.
namespace spin::research::char2::field128_lift::gamma30_completion {
namespace baseline=spin::research::char2::field128_lift::gamma30;
namespace deferred=spin::research::char2::field128_lift::gamma30_deferred;
namespace literal=baseline::literal;
namespace field=baseline::field;
namespace sr=baseline::sr;
using Plan=baseline::Plan;
using Element=baseline::Element;
using Block=baseline::Block;
using State=baseline::State;
using Packet=baseline::Packet;
using Refresh=baseline::Refresh;
using Kind=baseline::Kind;
using deferred::outer;
using baseline::loadPacket;
using baseline::storePacket;

#if defined(_MSC_VER)
#define GAMMA30_COMPLETE_INLINE __forceinline
#define GAMMA30_COMPLETE_NOINLINE __declspec(noinline)
#else
#define GAMMA30_COMPLETE_INLINE inline __attribute__((always_inline))
#define GAMMA30_COMPLETE_NOINLINE inline __attribute__((noinline))
#endif

namespace detail {
using baseline::detail::clearPadding;

// For degree <=29 on the additive 32-node set, sum(q_j)=sum(j*q_j)=0.
// With S=sum_{j<30}q_j and T=sum_{j<30}j*q_j, the missing values are
// q30=T+31*S and q31=T+30*S. Fold 30 into the packed constants:
// B=sum_{j<30}(j XOR30)*q_j, so q31=B and q30=B+S. Translating every node
// by 128 does not change these completion coefficients. Ignored storage
// coordinates 30/31 have zero weights and are cleared before the forward map.
constexpr unsigned weight(unsigned j) noexcept {return j<30?(j^30U):0U;}
static_assert([]{
    for(unsigned j=0;j!=30;++j)
        if(baseline::detail::completion[0][j]!=(j^31U)||
           baseline::detail::completion[1][j]!=(j^30U))return false;
    return true;
}(),"completion identity must match the frozen Lagrange coefficients");

template<unsigned Vector>
GAMMA30_COMPLETE_INLINE __m512i weighted(__m512i x) noexcept {
    return sr::byteProduct4<true,weight(4*Vector),weight(4*Vector+1),
        weight(4*Vector+2),weight(4*Vector+3)>(x);
}

GAMMA30_COMPLETE_INLINE void complete(__m512i* x) noexcept {
    clearPadding(x);
    auto sum=_mm512_setzero_si512(),b=sum;
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        sum=_mm512_xor_si512(sum,x[i]);
        b=_mm512_xor_si512(b,weighted<i>(x[i]));
    });
    sum=hdbw::detail::repeatRecordSum(sum);
    b=hdbw::detail::repeatRecordSum(b);
    x[7]=_mm512_mask_mov_epi64(x[7],0x30,_mm512_xor_si512(b,sum));
    x[7]=_mm512_mask_mov_epi64(x[7],0xc0,b);
}

// Exact transpose of E: F^30 -> F^32, not another forward interpolation.
// For a=z30,b=z31, (E^T z)_j=z_j+a+(j XOR30)*(a+b), j<30.
// All 32 input coordinates matter; only the two output padding slots vanish.
GAMMA30_COMPLETE_INLINE void completeTranspose(__m512i* x) noexcept {
    const auto a=sr::repeat<2>(x[7]);
    const auto sum=_mm512_xor_si512(a,sr::repeat<3>(x[7]));
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        x[i]=_mm512_ternarylogic_epi64(x[i],a,weighted<i>(sum),0x96);
    });
    clearPadding(x);
}
} // namespace detail

// The following fixed-width wrappers preserve the frozen schedules and
// replace only the calls to complete/completeTranspose. Explicit qualified
// calls avoid ADL selecting a baseline function for the shared state types.
GAMMA30_COMPLETE_INLINE void emitAdd(const State& state,Packet& output) noexcept {
    State coefficients=state;
    detail::complete(coefficients.values);
    baseline::detail::fft32<128,Kind::Inverse>(coefficients.values);
    baseline::detail::emitCoset<0>(coefficients,output);
    baseline::detail::emitCoset<32>(coefficients,output);
}

GAMMA30_COMPLETE_INLINE State feedback(Packet& input) noexcept {
    baseline::detail::fft32<0,Kind::ForwardTranspose>(input.values);
    baseline::detail::fft32<32,Kind::ForwardTranspose>(input.values+8);
    State result;
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        result.values[i]=_mm512_xor_si512(input.values[i],input.values[i+8]);
    });
    baseline::detail::fft32<128,Kind::InverseTranspose>(result.values);
    detail::completeTranspose(result.values);
    return result;
}

GAMMA30_COMPLETE_INLINE void refreshMap(State& state) noexcept {
    detail::complete(state.values);
    const auto tail=state.values[7];
    baseline::detail::fft32<0,Kind::Inverse>(state.values);
    baseline::detail::fft32<32,Kind::Forward>(state.values);
    c8Unroll(std::make_index_sequence<7>{},[&](auto index){
        constexpr unsigned i=7-index;
        state.values[i]=_mm512_shuffle_i32x4(state.values[i-1],state.values[i],0x4e);
    });
    state.values[0]=_mm512_shuffle_i32x4(tail,state.values[0],0x4e);
    detail::clearPadding(state.values);
}

GAMMA30_COMPLETE_INLINE void finish(State& state,State syndrome,const Refresh& diagonal) noexcept {
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        syndrome.values[i]=field::multiply4(
            _mm512_xor_si512(state.values[i],field::multiplyY4(syndrome.values[i])),
            diagonal.before.data()+4*i);
    });
    detail::clearPadding(syndrome.values);
    gamma30_completion::refreshMap(syndrome);
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        state.values[i]=field::multiply4(syndrome.values[i],diagonal.after.data()+4*i);
    });
    detail::clearPadding(state.values);
}

template<bool Stream>
GAMMA30_COMPLETE_INLINE void innerPacket(const Block* input,const std::uint16_t* route,Block* output,
                                        State& state,const Refresh& diagonal) noexcept {
    Packet y=baseline::loadPacket(input,route);
    gamma30_completion::emitAdd(state,y);
    baseline::storePacket<Stream>(y,output);
    gamma30_completion::finish(state,gamma30_completion::feedback(y),diagonal);
}

template<bool Stream>
GAMMA30_COMPLETE_NOINLINE void inner(const Plan& plan,const Block* planes,Block* output){
    const auto& p=plan.outerPlan;
    if constexpr(Stream)if(reinterpret_cast<std::uintptr_t>(output)&63U){gamma30_completion::inner<false>(plan,planes,output);return;}
    State state{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            gamma30_completion::innerPacket<Stream>(source,p.packetRoute.data()+(begin+j)/2,output+begin+j,
                                                   state,plan.refresh[(begin+j)/64]);
        }
    }
    if constexpr(Stream)_mm_sfence();
}

template<bool Stream=true>
inline void forward(const Plan& plan,const Block* input,Block* output,Block* scratch){
    deferred::outer<true>(plan.outerPlan,input,scratch);
    gamma30_completion::inner<Stream>(plan,scratch,output);
}

#undef GAMMA30_COMPLETE_INLINE
#undef GAMMA30_COMPLETE_NOINLINE
} // namespace spin::research::char2::field128_lift::gamma30_completion
