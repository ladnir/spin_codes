#pragma once
#include "Gamma30CompletionEncoder.h"

// Exact route-load fusion only. All maps, diagonals, butterfly layouts,
// output stores, outer arithmetic and setup are the frozen completion path.
namespace spin::research::char2::field128_lift::gamma30_route {
namespace previous=spin::research::char2::field128_lift::gamma30_completion;
namespace baseline=previous::baseline;
namespace literal=previous::literal;
using Plan=previous::Plan;
using Element=previous::Element;
using Block=previous::Block;
using State=previous::State;
using Packet=previous::Packet;
using Refresh=previous::Refresh;
using previous::outer;
using previous::storePacket;
using previous::emitAdd;
using previous::feedback;
using previous::refreshMap;
using previous::finish;

#if defined(_MSC_VER)
#define GAMMA30_ROUTE_INLINE __forceinline
#define GAMMA30_ROUTE_NOINLINE __declspec(noinline)
#else
#define GAMMA30_ROUTE_INLINE inline __attribute__((always_inline))
#define GAMMA30_ROUTE_NOINLINE inline __attribute__((noinline))
#endif

namespace detail {
// Identical two pair loads, insert and independent pair flips as the frozen
// innerInput4, but compile-time route indices directly build ab0/cd0/ab1/cd1.
template<unsigned Lo,unsigned Hi>
GAMMA30_ROUTE_INLINE __m512i pairPair(const Block* input,const std::uint16_t* route) noexcept {
    const unsigned r0=route[Lo],r1=route[Hi];
    const auto low=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+(r0&0x7fffU)));
    const auto high=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+(r1&0x7fffU)));
    const auto value=_mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
    const auto swaps=__mmask8(((0U-(r0>>15))&0x0fU)|((0U-(r1>>15))&0xf0U));
    return _mm512_mask_shuffle_i64x2(value,swaps,value,value,0xb1);
}
} // namespace detail

GAMMA30_ROUTE_INLINE Packet loadPacket(const Block* input,const std::uint16_t* route) noexcept {
    Packet result;
    c8Unroll(std::make_index_sequence<4>{},[&](auto group){
        constexpr unsigned j=8*group;
        const auto ab0=detail::pairPair<j,j+2>(input,route);
        const auto cd0=detail::pairPair<j+4,j+6>(input,route);
        const auto ab1=detail::pairPair<j+1,j+3>(input,route);
        const auto cd1=detail::pairPair<j+5,j+7>(input,route);
        result.values[group]=_mm512_shuffle_i32x4(ab0,cd0,0x88);
        result.values[4+group]=_mm512_shuffle_i32x4(ab0,cd0,0xdd);
        result.values[8+group]=_mm512_shuffle_i32x4(ab1,cd1,0x88);
        result.values[12+group]=_mm512_shuffle_i32x4(ab1,cd1,0xdd);
    });
    return result;
}

template<bool Stream>
GAMMA30_ROUTE_INLINE void innerPacket(const Block* input,const std::uint16_t* route,Block* output,
                                     State& state,const Refresh& diagonal) noexcept {
    Packet y=gamma30_route::loadPacket(input,route);
    previous::emitAdd(state,y);
    previous::storePacket<Stream>(y,output);
    previous::finish(state,previous::feedback(y),diagonal);
}

template<bool Stream>
GAMMA30_ROUTE_NOINLINE void inner(const Plan& plan,const Block* planes,Block* output){
    const auto& p=plan.outerPlan;
    if constexpr(Stream)if(reinterpret_cast<std::uintptr_t>(output)&63U){gamma30_route::inner<false>(plan,planes,output);return;}
    State state{};
    for(unsigned region=0;region!=256;++region){
        const auto begin=region*p.rows;const auto* source=planes+region*p.stride;
        for(std::size_t j=0;j!=p.rows;j+=64){
            if(region+1<256){
                const auto* next=reinterpret_cast<const char*>(planes+(region+1)*p.stride+j);
                for(unsigned h=0;h!=16;++h)_mm_prefetch(next+64*h,_MM_HINT_T1);
            }
            gamma30_route::innerPacket<Stream>(source,p.packetRoute.data()+(begin+j)/2,output+begin+j,
                                               state,plan.refresh[(begin+j)/64]);
        }
    }
    if constexpr(Stream)_mm_sfence();
}

template<bool Stream=true>
inline void forward(const Plan& plan,const Block* input,Block* output,Block* scratch){
    previous::outer<true>(plan.outerPlan,input,scratch);
    gamma30_route::inner<Stream>(plan,scratch,output);
}

#undef GAMMA30_ROUTE_INLINE
#undef GAMMA30_ROUTE_NOINLINE
} // namespace spin::research::char2::field128_lift::gamma30_route
