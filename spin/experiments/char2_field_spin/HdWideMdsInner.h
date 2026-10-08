#pragma once
#include "HdSumrankInner.h"

namespace spin::research::char2::hdwide {
using hdsumrank::Coefficient;
using hdsumrank::Scalar;
using Kind=hdbw::detail::Kind;
#if defined(_MSC_VER)
#define HDW_INLINE __forceinline
#else
#define HDW_INLINE inline __attribute__((always_inline))
#endif

template<unsigned S> struct alignas(64) State {
    static_assert(S==20||S==24||S==32);
    __m512i v[S/4]{};
};
template<unsigned S> struct Setup {
    std::array<Coefficient,S> before{},after{},v{};
};
struct alignas(64) Padded {__m512i v[8]{};};

template<unsigned Offset,Kind K>
HDW_INLINE void fft32(Padded& x){
    if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){hdsumrank::affinePair<hdbw::detail::normalized(4,Offset),K>(x.v[i],x.v[i+4]);});
        hdsumrank::naturalFft16<Offset,K,true>(x.v[0],x.v[1],x.v[2],x.v[3]);
        hdsumrank::naturalFft16<Offset^16,K,true>(x.v[4],x.v[5],x.v[6],x.v[7]);
    }else{
        hdsumrank::naturalFft16<Offset,K,true>(x.v[0],x.v[1],x.v[2],x.v[3]);
        hdsumrank::naturalFft16<Offset^16,K,true>(x.v[4],x.v[5],x.v[6],x.v[7]);
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){hdsumrank::affinePair<hdbw::detail::normalized(4,Offset),K>(x.v[i],x.v[i+4]);});
    }
}

template<unsigned Offset,unsigned S,bool Routed,class Index>
HDW_INLINE void tile(const Block* input,const Index* route,Block* output,
                      const State<S>& state,State<S>& feedback){
    Padded emit;
    c8Unroll(std::make_index_sequence<S/4>{},[&](auto i){emit.v[i]=state.v[i];});
    fft32<Offset,Kind::Forward>(emit);
    Padded inputWords;
    c8Unroll(std::make_index_sequence<8>{},[&](auto i){
        inputWords.v[i]=hdbw::detail::input4<Routed>(input,route,Offset+4*i);
        _mm512_storeu_si512(output+Offset+4*i,_mm512_xor_si512(inputWords.v[i],emit.v[i]));
    });
    fft32<Offset,Kind::ForwardTranspose>(inputWords);
    c8Unroll(std::make_index_sequence<S/4>{},[&](auto i){feedback.v[i]=_mm512_xor_si512(feedback.v[i],inputWords.v[i]);});
}

template<unsigned S>
HDW_INLINE void finish(State<S>& state,const State<S>& feedback,const Setup<S>& setup){
    __m512i dot=_mm512_setzero_si512();Padded refreshed;
    c8Unroll(std::make_index_sequence<S/4>{},[&](auto i){
        dot=_mm512_xor_si512(dot,wide16::multiply512(state.v[i],setup.v.data()+4*i));
        refreshed.v[i]=wide16::multiply512(feedback.v[i],setup.before.data()+4*i);
    });
    dot=hdbw::detail::repeatRecordSum(dot);
    // H is the leading S-by-S submatrix of the Cauchy interpolation
    // from 0..31 to 32..63. Pad inputs, evaluate, then retain S outputs.
    fft32<0,Kind::Inverse>(refreshed);fft32<32,Kind::Forward>(refreshed);
    c8Unroll(std::make_index_sequence<S/4>{},[&](auto i){
        state.v[i]=_mm512_xor_si512(state.v[i],wide16::multiply512(_mm512_xor_si512(refreshed.v[i],dot),setup.after.data()+4*i));
    });
}

template<unsigned S,bool Routed,class Index>
HDW_INLINE void packet(const Block* input,const Index* route,Block* output,
                        State<S>& state,const Setup<S>& setup){
    State<S> feedback;
    tile<0,S,Routed>(input,route,output,state,feedback);
    tile<32,S,Routed>(input,route,output,state,feedback);
    tile<64,S,Routed>(input,route,output,state,feedback);
    tile<96,S,Routed>(input,route,output,state,feedback);
    finish(state,feedback,setup);
}
#undef HDW_INLINE
} // namespace spin::research::char2::hdwide
