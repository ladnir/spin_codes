#pragma once
#include "HdSumrankInner.h"

#ifndef HDGROUP_INNER
#define HDGROUP_INNER 0
#endif
#ifndef HDGROUP_FUSE_REFRESH
#define HDGROUP_FUSE_REFRESH 0
#endif
#ifndef HDGROUP_XOR_REFRESH
#define HDGROUP_XOR_REFRESH 0
#endif
#ifndef HDGROUP_PACKET_ROUTE
#define HDGROUP_PACKET_ROUTE 0
#endif
#ifndef HDGROUP_ROUTE_WIDTH
#define HDGROUP_ROUTE_WIDTH 4
#endif
#ifndef HDGROUP_SHIFTED_PAIRS
#define HDGROUP_SHIFTED_PAIRS 0
#endif
#ifndef HDGROUP_PAIR_SWAPS
#define HDGROUP_PAIR_SWAPS 0
#endif
#ifndef HDGROUP_DECODED_ROUTE
#define HDGROUP_DECODED_ROUTE 0
#endif
static_assert(!HDGROUP_DECODED_ROUTE || (HDGROUP_PAIR_SWAPS && HDGROUP_ROUTE_WIDTH==2));
#ifndef HDGROUP_NATIVE_OUTPUT
#define HDGROUP_NATIVE_OUTPUT 0
#endif
static_assert(HDGROUP_INNER>=0&&HDGROUP_INNER<=5);
static_assert(HDGROUP_XOR_REFRESH==0||HDGROUP_XOR_REFRESH==1);
static_assert(!HDGROUP_XOR_REFRESH||HDGROUP_INNER>=4,
              "The experimental XOR refresh requires persistent striped state");
static_assert(HDGROUP_PACKET_ROUTE==0||HDGROUP_PACKET_ROUTE==1);
static_assert(HDGROUP_ROUTE_WIDTH==2||HDGROUP_ROUTE_WIDTH==4);
// 0: even pairs; 1: independent random parity; 2: parity=region&1.
static_assert(HDGROUP_SHIFTED_PAIRS>=0&&HDGROUP_SHIFTED_PAIRS<=2);
static_assert(!HDGROUP_SHIFTED_PAIRS||(HDGROUP_PACKET_ROUTE&&HDGROUP_ROUTE_WIDTH==2),
              "Shifted pairs require width-two packet routing");
static_assert(HDGROUP_PAIR_SWAPS==0||HDGROUP_PAIR_SWAPS==1);
static_assert(!HDGROUP_PAIR_SWAPS||(HDGROUP_PACKET_ROUTE&&HDGROUP_ROUTE_WIDTH==2),
              "Pair swaps require width-two packet routing");
static_assert(HDGROUP_NATIVE_OUTPUT==0||HDGROUP_NATIVE_OUTPUT==1);
static_assert(!(HDGROUP_PACKET_ROUTE||HDGROUP_NATIVE_OUTPUT)||HDGROUP_INNER>=2,
              "Packet routes and native output require a striped-input inner variant");

namespace spin::research::char2::hdgroup {
#if defined(_MSC_VER)
#define HDGI_INLINE __forceinline
#else
#define HDGI_INLINE inline __attribute__((always_inline))
#endif

// Physical position of logical record i in one 64-record inner step. This
// fixed output permutation is independent of the sampled regional route.
constexpr unsigned nativeOutputIndex(unsigned i){return 4*(i%16)+i/16;}

template<unsigned J>
HDGI_INLINE __m512i innerInput4(const Block* input,const std::uint16_t* route){
    static_assert(J<16);
    if constexpr(HDGROUP_PACKET_ROUTE){
        if constexpr(HDGROUP_ROUTE_WIDTH==4){
            // One entry names four contiguous outer-group records, with
            // packet lanes assigned to the four independent RS[16,4] blocks.
            return _mm512_loadu_si512(input+route[J]);
        }else if constexpr(HDGROUP_PAIR_SWAPS){
            // Bit15 independently reverses each pair's two 128-bit records.
            // The address still names a contiguous pair (including padding
            // for an odd pair wrapping from the last record to the first).
            const unsigned r0=route[2*J],r1=route[2*J+1];
#if HDGROUP_DECODED_ROUTE
            // Setup-derived packet layout: 32 clean bases, then 16 mask bytes.
            // No fresh randomness and no per-encode route expansion.
            const auto low=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+r0));
            const auto high=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+r1));
            const auto swaps=__mmask8(reinterpret_cast<const std::uint8_t*>(route+32)[J]);
            const auto value=_mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
#else
            const auto low=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+(r0&0x7fffU)));
            const auto high=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+(r1&0x7fffU)));
            const auto value=_mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
            const auto swaps=__mmask8(((0U-(r0>>15))&0x0fU)|((0U-(r1>>15))&0xf0U));
#endif
            // 0xb1 selects 128-bit lanes [1,0,3,2]; each four-bit mask
            // selects one complete 256-bit pair, never half of a record.
            return _mm512_mask_shuffle_i64x2(value,swaps,value,value,0xb1);
        }else{
            // Two independently positioned source pairs fill one ZMM. The
            // low/high record of each pair keeps its packet-lane identity.
            // Shifted pairs require input[rows]=input[0], so the final odd
            // base wraps through padding without a branch in this kernel.
            const auto low=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+route[2*J]));
            const auto high=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input+route[2*J+1]));
            return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
        }
    }else return hdsumrank::stripedInput4<true>(input,route,4*J);
}

// One short RS coset at a time, striped across the four independent RS[16,4]
// blocks. No 16-vector packet array survives across cosets. A/C and the state
// recurrence do not depend on the input route or physical output layout.
template<unsigned Offset,bool Stream>
HDGI_INLINE void innerCoset(const Block* input,const std::uint16_t* route,Block* output,
                           const hdsumrank::State& q,hdsumrank::State& feedback){
    static_assert(Offset<16&&Offset%4==0);
    namespace sr=hdsumrank;
    using Kind=hdbw::detail::Kind;
    constexpr auto c2=hdbw::detail::normalized(1,Offset);
    constexpr auto c0=hdbw::detail::normalized(0,Offset);
    constexpr auto c1=hdbw::detail::normalized(0,Offset+2);
    auto a=q.values[0],b=q.values[1],c=q.values[2],d=q.values[3];
    sr::affinePair<c2,Kind::Forward>(a,c);sr::affinePair<c2,Kind::Forward>(b,d);
    sr::affinePair<c0,Kind::Forward>(a,b);sr::affinePair<c1,Kind::Forward>(c,d);
    __m512i v0,v1,v2,v3;
    if constexpr(HDGROUP_INNER>=2){
        v0=innerInput4<Offset>(input,route);
        v1=innerInput4<Offset+1>(input,route);
        v2=innerInput4<Offset+2>(input,route);
        v3=innerInput4<Offset+3>(input,route);
        a=_mm512_xor_si512(a,v0);b=_mm512_xor_si512(b,v1);
        c=_mm512_xor_si512(c,v2);d=_mm512_xor_si512(d,v3);
        // C=A^T and CA=0, so feedback may consume the output instead of X.
        if constexpr(HDGROUP_INNER==3||HDGROUP_INNER==5){v0=a;v1=b;v2=c;v3=d;}
        if constexpr(!HDGROUP_NATIVE_OUTPUT)sr::transposeRecords4(a,b,c,d);
    }else{
        sr::transposeRecords4(a,b,c,d);
        v0=sr::stripedInput4<true>(input,route,Offset);
        v1=sr::stripedInput4<true>(input,route,16+Offset);
        v2=sr::stripedInput4<true>(input,route,32+Offset);
        v3=sr::stripedInput4<true>(input,route,48+Offset);
        a=_mm512_xor_si512(a,v0);b=_mm512_xor_si512(b,v1);
        c=_mm512_xor_si512(c,v2);d=_mm512_xor_si512(d,v3);
    }
    constexpr unsigned out0=HDGROUP_NATIVE_OUTPUT?nativeOutputIndex(Offset):Offset;
    constexpr unsigned out1=HDGROUP_NATIVE_OUTPUT?nativeOutputIndex(Offset+1):16+Offset;
    constexpr unsigned out2=HDGROUP_NATIVE_OUTPUT?nativeOutputIndex(Offset+2):32+Offset;
    constexpr unsigned out3=HDGROUP_NATIVE_OUTPUT?nativeOutputIndex(Offset+3):48+Offset;
    if constexpr(Stream){
        _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out0),a);
        _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out1),b);
        _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out2),c);
        _mm512_stream_si512(reinterpret_cast<__m512i*>(output+out3),d);
    }else{
        _mm512_storeu_si512(output+out0,a);_mm512_storeu_si512(output+out1,b);
        _mm512_storeu_si512(output+out2,c);_mm512_storeu_si512(output+out3,d);
    }
    if constexpr(HDGROUP_INNER<2)sr::transposeRecords4(v0,v1,v2,v3);
    sr::affinePair<c0,Kind::ForwardTranspose>(v0,v1);sr::affinePair<c1,Kind::ForwardTranspose>(v2,v3);
    sr::affinePair<c2,Kind::ForwardTranspose>(v0,v2);sr::affinePair<c2,Kind::ForwardTranspose>(v1,v3);
    if constexpr(Offset==0){feedback.values[0]=v0;feedback.values[1]=v1;feedback.values[2]=v2;feedback.values[3]=v3;}
    else{
        feedback.values[0]=_mm512_xor_si512(feedback.values[0],v0);
        feedback.values[1]=_mm512_xor_si512(feedback.values[1],v1);
        feedback.values[2]=_mm512_xor_si512(feedback.values[2],v2);
        feedback.values[3]=_mm512_xor_si512(feedback.values[3],v3);
    }
}

template<unsigned Level,unsigned Offset,hdbw::detail::Kind K>
HDGI_INLINE void stripedPair(__m512i& a,__m512i& b){
    static_assert(Level<2&&(K==hdbw::detail::Kind::Forward||K==hdbw::detail::Kind::Inverse));
    constexpr auto c0=hdbw::detail::normalized(Level,Offset);
    constexpr auto c1=hdbw::detail::normalized(Level,Offset+4);
    constexpr auto c2=hdbw::detail::normalized(Level,Offset+8);
    constexpr auto c3=hdbw::detail::normalized(Level,Offset+12);
    const auto product=[](__m512i x){return hdsumrank::byteProduct4<true,c0,c1,c2,c3>(x);};
    if constexpr(K==hdbw::detail::Kind::Inverse)b=_mm512_xor_si512(a,b);
    if constexpr(c0==1&&c1==1&&c2==1&&c3==1)a=_mm512_xor_si512(a,b);
    else if constexpr(c0||c1||c2||c3)a=_mm512_xor_si512(a,product(b));
    if constexpr(K==hdbw::detail::Kind::Forward)b=_mm512_xor_si512(a,b);
}

template<unsigned Half,unsigned Offset,hdbw::detail::Kind K>
HDGI_INLINE __m512i stripedHigh(__m512i x){
    static_assert((Half==1||Half==2)&&(K==hdbw::detail::Kind::Forward||K==hdbw::detail::Kind::Inverse));
    // State coordinate i is in vector i%4, record lane i/4. These are
    // the size-eight/sixteen stages, not the natural layout's low stages.
    constexpr auto c0=hdbw::detail::normalized(Half==1?2:3,Offset);
    constexpr auto c1=hdbw::detail::normalized(Half==1?2:3,Offset^(Half==1?8:0));
    constexpr int lowShuffle=Half==1?0xa0:0x44,highShuffle=Half==1?0xf5:0xee;
    constexpr __mmask8 lowMask=Half==1?0x33:0x0f,highMask=static_cast<__mmask8>(~lowMask);
    const auto a=_mm512_shuffle_i32x4(x,x,lowShuffle);
    if constexpr(c0==0&&c1==0)return _mm512_mask_xor_epi64(x,highMask,x,a);
    else{
        const auto b=_mm512_shuffle_i32x4(x,x,highShuffle);
        if constexpr(K==hdbw::detail::Kind::Forward){
            const auto p=[&](){
                if constexpr(Half==1)return hdsumrank::byteProduct4<true,c0,c0^1,c1,c1^1>(b);
                else return hdsumrank::byteProduct4<true,c0,c0,c0^1,c0^1>(b);
            }();
            return _mm512_xor_si512(a,p);
        }else{
            const auto sum=_mm512_xor_si512(a,b);
            const auto p=[&](){
                if constexpr(Half==1)return hdsumrank::byteProduct4<true,c0,c0,c1,c1>(sum);
                else return hdsumrank::byteProduct4<true,c0>(sum);
            }();
            return _mm512_mask_xor_epi64(sum,lowMask,a,p);
        }
    }
}

template<unsigned Offset,hdbw::detail::Kind K,bool SkipHigh=false>
HDGI_INLINE void stripedFft16(hdsumrank::State& x){
    auto& a=x.values[0];auto& b=x.values[1];auto& c=x.values[2];auto& d=x.values[3];
    if constexpr(K==hdbw::detail::Kind::Inverse){
        stripedPair<0,Offset,K>(a,b);stripedPair<0,Offset+2,K>(c,d);
        stripedPair<1,Offset,K>(a,c);stripedPair<1,Offset,K>(b,d);
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            x.values[i]=stripedHigh<1,Offset,K>(x.values[i]);
            if constexpr(!SkipHigh)x.values[i]=stripedHigh<2,Offset,K>(x.values[i]);
        });
    }else{
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            if constexpr(!SkipHigh)x.values[i]=stripedHigh<2,Offset,K>(x.values[i]);
            x.values[i]=stripedHigh<1,Offset,K>(x.values[i]);
        });
        stripedPair<1,Offset,K>(a,c);stripedPair<1,Offset,K>(b,d);
        stripedPair<0,Offset,K>(a,b);stripedPair<0,Offset+2,K>(c,d);
    }
}

HDGI_INLINE void finishStriped(hdsumrank::State& state,hdsumrank::State feedback,
                              const hdsumrank::ByteDiagonalRefresh& t){
    namespace sr=hdsumrank;
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){
        feedback.values[i]=sr::byteRefreshProduct<i>(_mm512_xor_si512(state.values[i],feedback.values[i]),t.before);
    });
    if constexpr(HDGROUP_XOR_REFRESH){
        // A distinct setup ensemble: D2 (I+11^T) D1, not the certified
        // interpolation refresh. All sixteen 128-bit records contribute
        // once to the common sum; their physical striped order is irrelevant.
        auto sum=_mm512_xor_si512(_mm512_xor_si512(feedback.values[0],feedback.values[1]),
                                 _mm512_xor_si512(feedback.values[2],feedback.values[3]));
        sum=_mm512_xor_si512(sum,_mm512_shuffle_i64x2(sum,sum,0x4e));
        sum=_mm512_xor_si512(sum,_mm512_shuffle_i64x2(sum,sum,0xb1));
        feedback.values[0]=_mm512_xor_si512(feedback.values[0],sum);
        feedback.values[1]=_mm512_xor_si512(feedback.values[1],sum);
        feedback.values[2]=_mm512_xor_si512(feedback.values[2],sum);
        feedback.values[3]=_mm512_xor_si512(feedback.values[3],sum);
    }else if constexpr(HDGROUP_FUSE_REFRESH){
        stripedFft16<0,hdbw::detail::Kind::Inverse,true>(feedback);
        // The adjacent size-sixteen stages compose to (a,b) ->
        // (a+c(a+b), b+c(a+b)). One record shuffle replaces three.
        constexpr auto c=hdbw::detail::normalized(3,16);
        c8Unroll(std::make_index_sequence<4>{},[&](auto i){
            const auto x=feedback.values[i];
            const auto sum=_mm512_xor_si512(x,_mm512_shuffle_i32x4(x,x,0x4e));
            feedback.values[i]=_mm512_xor_si512(x,sr::byteProduct4<true,c>(sum));
        });
        stripedFft16<16,hdbw::detail::Kind::Forward,true>(feedback);
    }else{
        stripedFft16<0,hdbw::detail::Kind::Inverse>(feedback);
        stripedFft16<16,hdbw::detail::Kind::Forward>(feedback);
    }
    c8Unroll(std::make_index_sequence<4>{},[&](auto i){state.values[i]=sr::byteRefreshProduct<i>(feedback.values[i],t.after);});
}

template<bool Stream,class MixerData,class RefreshData>
HDGI_INLINE void innerPacket4(const Block* input,const std::uint16_t* route,Block* output,
                             hdsumrank::State& state,const MixerData& mixer,const RefreshData& refresh){
    namespace sr=hdsumrank;
    sr::State q=state,feedback;
    if constexpr(HDGROUP_INNER<4)sr::transposeRecords4(q.values[0],q.values[1],q.values[2],q.values[3]);
    innerCoset<0,Stream>(input,route,output,q,feedback);
    innerCoset<4,Stream>(input,route,output,q,feedback);
    innerCoset<8,Stream>(input,route,output,q,feedback);
    innerCoset<12,Stream>(input,route,output,q,feedback);
    if constexpr(HDGROUP_INNER>=4)finishStriped(state,feedback,refresh);
    else{
        sr::transposeRecords4(feedback.values[0],feedback.values[1],feedback.values[2],feedback.values[3]);
        sr::mdsFinish<true,true,2>(state,feedback,mixer,refresh);
    }
}
#undef HDGI_INLINE
}
