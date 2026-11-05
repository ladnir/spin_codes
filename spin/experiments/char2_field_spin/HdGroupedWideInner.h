#pragma once
#include "HdGroupedInner.h"

#ifndef HDGROUP_WIDE_INNER
#define HDGROUP_WIDE_INNER 0
#endif
static_assert(HDGROUP_WIDE_INNER==0||HDGROUP_WIDE_INNER==1);

namespace spin::research::char2::hdgroup {
#if defined(_MSC_VER)
#define HDGW_INLINE __forceinline
#else
#define HDGW_INLINE inline __attribute__((always_inline))
#endif

// Global A=RS[64,16], C=A^T. The state is always striped: coordinate i
// occupies vector i%4, record lane i/4, including at the refresh interface.
// Input uses innerInput4's record/packet route representation; the four lanes
// of vector j are logical records 16*lane+j. All schedules below are unrolled.
template<bool Stream,class MixerData>
HDGW_INLINE void innerPacket16(const Block* input,const std::uint16_t* route,Block* output,
                              hdsumrank::State& state,const MixerData&,
                              const hdsumrank::ByteDiagonalRefresh& refresh){
    namespace sr=hdsumrank;
    // One vector per coefficient/evaluation position. The full packet is
    // local to this non-coroutine helper; no extra input or output arrays.
    __m512i x[16];
    c8Unroll(std::make_index_sequence<16>{},[&](auto i){
        x[i]=sr::repeat<i/4>(state.values[i%4]);
    });
    sr::stripedStage<3,false>(x);
    sr::stripedStage<2,false>(x);
    sr::stripedStage<1,false>(x);
    sr::stripedStage<0,false>(x);

    c8Unroll(std::make_index_sequence<4>{},[&](auto g){
        constexpr unsigned j=4*g;
        x[j]=_mm512_xor_si512(x[j],innerInput4<j>(input,route));
        x[j+1]=_mm512_xor_si512(x[j+1],innerInput4<j+1>(input,route));
        x[j+2]=_mm512_xor_si512(x[j+2],innerInput4<j+2>(input,route));
        x[j+3]=_mm512_xor_si512(x[j+3],innerInput4<j+3>(input,route));
        auto a=x[j],b=x[j+1],c=x[j+2],d=x[j+3];
        if constexpr(!HDGROUP_NATIVE_OUTPUT)sr::transposeRecords4(a,b,c,d);
        constexpr unsigned out0=HDGROUP_NATIVE_OUTPUT?4*j:j;
        constexpr unsigned out1=HDGROUP_NATIVE_OUTPUT?4*j+4:16+j;
        constexpr unsigned out2=HDGROUP_NATIVE_OUTPUT?4*j+8:32+j;
        constexpr unsigned out3=HDGROUP_NATIVE_OUTPUT?4*j+12:48+j;
        if constexpr(Stream){
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out0),a);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out1),b);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out2),c);
            _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out3),d);
        }else{
            _mm512_storeu_si512(output+out0,a);_mm512_storeu_si512(output+out1,b);
            _mm512_storeu_si512(output+out2,c);_mm512_storeu_si512(output+out3,d);
        }
    });

    // A^T A=0: products of its degree-at-most-15 basis polynomials have
    // degree at most 30, whose sum over the 64-point additive space is zero.
    // Thus C Y=C X, and the output vectors can be reused as feedback input.
    sr::stripedStage<0,true>(x);
    sr::stripedStage<1,true>(x);
    sr::stripedStage<2,true>(x);
    sr::stripedStage<3,true>(x);
    sr::State feedback;
    c8Unroll(std::make_index_sequence<4>{},[&](auto r){
        // Each x[i] holds the four point-coset contributions to coefficient
        // i. Reduce the lanes directly into the persistent striped layout.
        auto a=x[r],b=x[r+4],c=x[r+8],d=x[r+12];
        sr::transposeRecords4(a,b,c,d);
        feedback.values[r]=_mm512_xor_si512(_mm512_xor_si512(a,b),_mm512_xor_si512(c,d));
    });
    finishStriped(state,feedback,refresh);
}

#undef HDGW_INLINE
} // namespace spin::research::char2::hdgroup
