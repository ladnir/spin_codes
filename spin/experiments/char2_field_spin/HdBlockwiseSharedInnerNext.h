#pragma once
#include "HdBlockwiseSharedInner.h"

namespace spin::research::char2::hdbwsharednext {
// Same exact beta64 map as HdBlockwiseSharedInner.h.
// This control uses63 rather than59 GFNI calls, trading four products for
// whole-vector feedback and emission. No changed map or certificate claim.
using State=hdbwshared::State;
using Mixer=hdbwshared::Mixer;
using Scalar=hdbwshared::Scalar;
inline constexpr unsigned packetWidth=64,stateSize=4,beta=64;
inline Mixer prepareMixer(const Scalar* u,const Scalar* v){return hdbwshared::prepareMixer(u,v);}
inline void loadState(State& s,const Block* p) noexcept {hdbwshared::loadState(s,p);}
inline void storeState(Block* p,const State& s) noexcept {hdbwshared::storeState(p,s);}
#if defined(_MSC_VER)
#define HDBSHN_INLINE __forceinline
#else
#define HDBSHN_INLINE inline __attribute__((always_inline))
#endif
namespace detail {
using hdbwshared::detail::Moments;
using hdbwshared::detail::permute;
using hdbwshared::detail::select;
using hdbwshared::detail::product;

// Nine fixed GFNI calls, with no 128-bit extract or insert. Each reduction
// stays in the existing four-record layout, then broadcasts its record sum.
HDBSHN_INLINE __m512i feedback(const Moments& m) noexcept {
    const auto p0=product<0,224,193,100>(m.v0);
    const auto p1=product<54,180,171,0>(m.v1);
    const auto p2=product<0,28,120,235>(m.v2);
    const auto p3=product<2,48,224,237>(m.v3);
    const auto p4=product<24,64,155,65>(m.v4);
    const auto p5=product<1,6,8,0>(m.v5);
    const auto z3=hdbw::detail::repeatRecordSum(_mm512_xor_si512(
        _mm512_xor_si512(_mm512_xor_si512(p0,p1),_mm512_xor_si512(p2,p3)),
        _mm512_xor_si512(p4,p5)));
    const auto a1=product<0,8,16,32>(m.v0);
    const auto a2=product<0,28,120,235>(m.v0);
    const auto tail=product<2,4,6,0>(select<0,4,4,4>(m.v3,m.v4));
    const auto z1=_mm512_xor_si512(hdbw::detail::repeatRecordSum(a1),
        _mm512_xor_si512(permute<0,0,0,0>(m.v2),
        _mm512_xor_si512(permute<0,0,0,0>(tail),permute<1,1,1,1>(tail))));
    const auto z2=_mm512_xor_si512(hdbw::detail::repeatRecordSum(a2),
        _mm512_xor_si512(permute<0,0,0,0>(m.v3),permute<2,2,2,2>(tail)));
    auto h=permute<0,0,0,0>(m.v0);
    h=_mm512_mask_xor_epi64(h,0x33,h,z1);
    h=_mm512_mask_xor_epi64(h,0x0f,h,z2);
    return _mm512_mask_xor_epi64(h,0x03,h,z3);
}

// The low two coordinate bits occupy the four ZMM record lanes. The high
// four bits select16 whole-vector evaluations of a degree-two Boolean
// polynomial. Ten GFNI calls build11 vector coefficients; the pruned zeta
// evaluation uses26 ZMM XORs, not120 individual-record XORs and repacking.
HDBSHN_INLINE void emission(__m512i q,Block* output) noexcept {
    const auto q0=permute<0,0,0,0>(q);
    const auto q1=permute<1,1,1,1>(q);
    const auto q2=permute<2,2,2,2>(q);
    const auto q3=permute<3,3,3,3>(q);
    const auto v1=product<64,65,66,67>(q1);
    const auto v2=product<161,161,160,160>(q2);
    const auto v3=product<213,116,206,110>(q3);
    const auto a0=_mm512_xor_si512(_mm512_xor_si512(q0,v1),_mm512_xor_si512(v2,v3));
    const auto b1=product<4,8,16,32>(q1);
    const auto b2=product<6,28,120,235>(q2);
    const auto b3=product<49,222,22,75>(q3);
    const auto base=_mm512_xor_si512(_mm512_xor_si512(b1,b2),b3);
    const auto low1=product<6,28,120,235>(q3);
    const auto low2=product<8,48,224,237>(q3);
    const auto pairs0=product<64,155,65,54>(q3);
    const auto pairs1=product<180,171,0,0>(q3);
    auto a1=permute<0,0,0,0>(base);
    a1=_mm512_mask_xor_epi64(a1,0xcc,a1,permute<0,0,0,0>(low1));
    a1=_mm512_mask_xor_epi64(a1,0xf0,a1,permute<0,0,0,0>(low2));
    auto a2=permute<1,1,1,1>(base);
    a2=_mm512_mask_xor_epi64(a2,0xcc,a2,permute<1,1,1,1>(low1));
    a2=_mm512_mask_xor_epi64(a2,0xf0,a2,permute<1,1,1,1>(low2));
    auto a4=permute<2,2,2,2>(base);
    a4=_mm512_mask_xor_epi64(a4,0xcc,a4,permute<2,2,2,2>(low1));
    a4=_mm512_mask_xor_epi64(a4,0xf0,a4,permute<2,2,2,2>(low2));
    auto a8=permute<3,3,3,3>(base);
    a8=_mm512_mask_xor_epi64(a8,0xcc,a8,permute<3,3,3,3>(low1));
    a8=_mm512_mask_xor_epi64(a8,0xf0,a8,permute<3,3,3,3>(low2));
    const auto a3=permute<0,0,0,0>(pairs0);
    const auto a5=permute<1,1,1,1>(pairs0);
    const auto a9=permute<2,2,2,2>(pairs0);
    const auto a6=permute<3,3,3,3>(pairs0);
    const auto a10=permute<0,0,0,0>(pairs1);
    const auto a12=permute<1,1,1,1>(pairs1);
    const auto z0_1=_mm512_xor_si512(a1,a0);
    const auto z0_3=_mm512_xor_si512(a3,a2);
    const auto z0_5=_mm512_xor_si512(a5,a4);
    const auto z0_9=_mm512_xor_si512(a9,a8);
    const auto z1_2=_mm512_xor_si512(a2,a0);
    const auto z1_3=_mm512_xor_si512(z0_3,z0_1);
    const auto z1_6=_mm512_xor_si512(a6,a4);
    const auto z1_7=_mm512_xor_si512(a6,z0_5);
    const auto z1_10=_mm512_xor_si512(a10,a8);
    const auto z1_11=_mm512_xor_si512(a10,z0_9);
    const auto z2_4=_mm512_xor_si512(a4,a0);
    const auto z2_5=_mm512_xor_si512(z0_5,z0_1);
    const auto z2_6=_mm512_xor_si512(z1_6,z1_2);
    const auto z2_7=_mm512_xor_si512(z1_7,z1_3);
    const auto z2_12=_mm512_xor_si512(a12,a8);
    const auto z2_13=_mm512_xor_si512(a12,z0_9);
    const auto z2_14=_mm512_xor_si512(a12,z1_10);
    const auto z2_15=_mm512_xor_si512(a12,z1_11);
    const auto z3_8=_mm512_xor_si512(a8,a0);
    const auto z3_9=_mm512_xor_si512(z0_9,z0_1);
    const auto z3_10=_mm512_xor_si512(z1_10,z1_2);
    const auto z3_11=_mm512_xor_si512(z1_11,z1_3);
    const auto z3_12=_mm512_xor_si512(z2_12,z2_4);
    const auto z3_13=_mm512_xor_si512(z2_13,z2_5);
    const auto z3_14=_mm512_xor_si512(z2_14,z2_6);
    const auto z3_15=_mm512_xor_si512(z2_15,z2_7);
    _mm512_store_si512(output+0,a0);
    _mm512_store_si512(output+4,z0_1);
    _mm512_store_si512(output+8,z1_2);
    _mm512_store_si512(output+12,z1_3);
    _mm512_store_si512(output+16,z2_4);
    _mm512_store_si512(output+20,z2_5);
    _mm512_store_si512(output+24,z2_6);
    _mm512_store_si512(output+28,z2_7);
    _mm512_store_si512(output+32,z3_8);
    _mm512_store_si512(output+36,z3_9);
    _mm512_store_si512(output+40,z3_10);
    _mm512_store_si512(output+44,z3_11);
    _mm512_store_si512(output+48,z3_12);
    _mm512_store_si512(output+52,z3_13);
    _mm512_store_si512(output+56,z3_14);
    _mm512_store_si512(output+60,z3_15);
}

template<bool Routed,class Index>
HDBSHN_INLINE void packet(const Block* input,const Index* route,Block* output,
                          State& state,const Mixer& mixer) noexcept {
    const auto q=state.values[0];
    hdbwsmall::detail::mix(state,mixer);
    alignas(64) Block emitted[64];
    emission(q,emitted);
    Moments m{};
#define HDBSHN_TILE(I) hdbwshared::detail::cubeTile<I,Routed>(input,route,output,emitted,m)
    HDBSHN_TILE(0);HDBSHN_TILE(1);HDBSHN_TILE(2);HDBSHN_TILE(3);
    HDBSHN_TILE(4);HDBSHN_TILE(5);HDBSHN_TILE(6);HDBSHN_TILE(7);
#undef HDBSHN_TILE
    state.values[0]=_mm512_xor_si512(state.values[0],hdbwsharednext::detail::feedback(m));
}
} // namespace detail

template<bool Fft=(HDBSH_FFT_CONTROL!=0)>
HDBSHN_INLINE void packetForwardRaw(const Block* input,Block* output,
                                    State& state,const Mixer& mixer) noexcept {
    if constexpr(Fft)hdbwshared::packetForwardRaw<true>(input,output,state,mixer);
    else detail::packet<false,std::uint16_t>(input,nullptr,output,state,mixer);
}
template<bool Fft=(HDBSH_FFT_CONTROL!=0),class Index>
HDBSHN_INLINE void packetForwardRouted(const Block* input,const Index* route,Block* output,
                                       State& state,const Mixer& mixer) noexcept {
    if constexpr(Fft)hdbwshared::packetForwardRouted<true>(input,route,output,state,mixer);
    else detail::packet<true>(input,route,output,state,mixer);
}
#undef HDBSHN_INLINE
} // namespace spin::research::char2::hdbwsharednext
