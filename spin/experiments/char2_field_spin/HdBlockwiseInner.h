#pragma once
#include "Cauchy256InnerFft.h"
#include "WideDiagonal16.h"
#include <array>
#include <cstdint>
#include <stdexcept>
#include <type_traits>

#ifndef HDBW_HW_GATHER
#define HDBW_HW_GATHER 0
#endif

namespace spin::research::char2::hdbw {

// Isolated t64/s16 forward COST CONTROL, not a distance certificate.
// All records contain eight independent GF65536 tower-field coordinates.
// U={0,...,63}, V={128,...,143} in AES GF256. A is interpolation from
// the sixteen V values to U, C=A^T. State stays in the V evaluation basis:
//   y=x+A*q; q_next=q+u*(v^T*q)+C*x.
// A has shape64x16, C16x64, and C*A=0. The fresh mixer is invertible iff
// 1+v^T*u!=0. Setup, including the required distribution, belongs to caller.
// No input diagonal, pre-prefix, puncture, or terminal-state flush is added.
inline constexpr unsigned packetWidth=64;
inline constexpr unsigned stateSize=16;
using Scalar=std::uint16_t;
using Coefficient=wide16::Coefficient;

struct alignas(64) State16 { __m512i values[4]{}; };
struct Mixer16 { std::array<Coefficient,16> u{},v{}; };
static_assert(sizeof(State16)==16*sizeof(Block));
static_assert(sizeof(Mixer16)==128);

#if defined(_MSC_VER)
#define HDBW_INLINE __forceinline
#else
#define HDBW_INLINE inline __attribute__((always_inline))
#endif

inline constexpr std::uint8_t gf8Inverse(std::uint8_t x) noexcept {
    return c256_inner_fft_detail::power(x,254);
}
inline constexpr std::uint8_t aCoefficient(unsigned output,unsigned state) noexcept {
    std::uint8_t numerator=1,denominator=1;
    for(unsigned j=0;j!=16;++j)if(j!=state){
        numerator=sfMultiply(numerator,static_cast<std::uint8_t>(output^(128+j)));
        denominator=sfMultiply(denominator,static_cast<std::uint8_t>(state^j));
    }
    return sfMultiply(numerator,gf8Inverse(denominator));
}
inline constexpr std::uint8_t cCoefficient(unsigned state,unsigned input) noexcept {
    return aCoefficient(input,state);
}

// Setup/scalar-reference interface. Zero entries of u or v are allowed;
// u as a vector must be nonzero. A uniform sampler must not independently
// reject the zero coordinates: it samples u uniformly from F^16\{0}, and
// v uniformly conditioned on dot(v,u)!=1.
inline bool validMixer16(const Scalar* u,const Scalar* v) noexcept {
    Scalar dot=0,nonzero=0;
    for(unsigned i=0;i!=16;++i){
        nonzero|=u[i];dot^=medium32::quadraticMultiply(u[i],v[i]);
    }
    return nonzero!=0&&dot!=1;
}
inline Mixer16 prepareMixer(const Scalar* u,const Scalar* v){
    if(!validMixer16(u,v))throw std::invalid_argument("blockwise mixer requires nonzero u and dot(v,u)!=1");
    Mixer16 result;
    wide16::prepare(u,result.u.data(),16);wide16::prepare(v,result.v.data(),16);
    return result;
}
HDBW_INLINE void loadState(State16& state,const Block* input16) noexcept {
    state.values[0]=_mm512_loadu_si512(input16);
    state.values[1]=_mm512_loadu_si512(input16+4);
    state.values[2]=_mm512_loadu_si512(input16+8);
    state.values[3]=_mm512_loadu_si512(input16+12);
}
HDBW_INLINE void storeState(Block* output16,const State16& state) noexcept {
    _mm512_storeu_si512(output16,state.values[0]);
    _mm512_storeu_si512(output16+4,state.values[1]);
    _mm512_storeu_si512(output16+8,state.values[2]);
    _mm512_storeu_si512(output16+12,state.values[3]);
}

namespace detail {
using Kind=c256_inner_fft_detail::Kind;

constexpr std::uint8_t normalized(unsigned level,unsigned offset) noexcept {
    std::uint8_t value=1,scale=1;
    for(unsigned j=0;j!=(1U<<level);++j){
        value=sfMultiply(value,static_cast<std::uint8_t>(offset^j));
        scale=sfMultiply(scale,static_cast<std::uint8_t>((1U<<level)^j));
    }
    return sfMultiply(value,gf8Inverse(scale));
}

template<unsigned Coefficient,Kind K>
HDBW_INLINE void pair(__m512i& a,__m512i& b) noexcept {
    if constexpr(K==Kind::Forward){
        if constexpr(Coefficient==1)a=_mm512_xor_si512(a,b);
        else if constexpr(Coefficient!=0)a=_mm512_xor_si512(a,sfConstant4<Coefficient>(b));
        b=_mm512_xor_si512(a,b);
    }else if constexpr(K==Kind::Inverse){
        b=_mm512_xor_si512(a,b);
        if constexpr(Coefficient==1)a=_mm512_xor_si512(a,b);
        else if constexpr(Coefficient!=0)a=_mm512_xor_si512(a,sfConstant4<Coefficient>(b));
    }else if constexpr(K==Kind::ForwardTranspose){
        a=_mm512_xor_si512(a,b);
        if constexpr(Coefficient==1)b=_mm512_xor_si512(b,a);
        else if constexpr(Coefficient!=0)b=_mm512_xor_si512(b,sfConstant4<Coefficient>(a));
    }else{
        if constexpr(Coefficient==1)b=_mm512_xor_si512(b,a);
        else if constexpr(Coefficient!=0)b=_mm512_xor_si512(b,sfConstant4<Coefficient>(a));
        a=_mm512_xor_si512(a,b);
    }
}

// Packed four-coordinate leaves, adapted from Cauchy256InnerFft to the
// natural binary basis0,1,2,4,8. The original eight-node helper uses a
// different ordered subspace and must not be reused with these offsets.
template<unsigned Half,unsigned Offset,Kind K>
HDBW_INLINE __m512i smallStage(__m512i x) noexcept {
    static_assert(Half==1||Half==2);
    constexpr auto c0=normalized(Half==1?0:1,Offset);
    constexpr auto c1=normalized(Half==1?0:1,Offset^(Half==1?2:0));
    constexpr int lowShuffle=Half==1?0xa0:0x44;
    constexpr int highShuffle=Half==1?0xf5:0xee;
    constexpr __mmask8 lowMask=Half==1?0x33:0x0f;
    constexpr __mmask8 highMask=static_cast<__mmask8>(~lowMask);
    if constexpr(c0==0&&c1==0){
        if constexpr(K==Kind::Forward||K==Kind::Inverse){
            const auto a=_mm512_shuffle_i32x4(x,x,lowShuffle);
            return _mm512_mask_xor_epi64(x,highMask,x,a);
        }else{
            const auto b=_mm512_shuffle_i32x4(x,x,highShuffle);
            return _mm512_mask_xor_epi64(x,lowMask,x,b);
        }
    }else{
        const auto a=_mm512_shuffle_i32x4(x,x,lowShuffle);
        const auto b=_mm512_shuffle_i32x4(x,x,highShuffle);
        if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
            constexpr unsigned flip=K==Kind::Forward?0:1;
            const auto coefficients=[](){
                if constexpr(Half==1)return c256_inner_fft_detail::coefficients<c0^flip,c0^flip^1,c1^flip,c1^flip^1>();
                else return c256_inner_fft_detail::coefficients<c0^flip,c0^flip,c0^flip^1,c0^flip^1>();
            }();
            if constexpr(K==Kind::Forward)return _mm512_xor_si512(a,_mm512_gf2p8mul_epi8(b,coefficients));
            else return _mm512_xor_si512(b,_mm512_gf2p8mul_epi8(a,coefficients));
        }else{
            const auto coefficients=[](){
                if constexpr(Half==1)return c256_inner_fft_detail::coefficients<c0,c0,c1,c1>();
                else return c256_inner_fft_detail::coefficients<c0,c0,c0,c0>();
            }();
            const auto sum=_mm512_xor_si512(a,b);
            const auto product=_mm512_gf2p8mul_epi8(sum,coefficients);
            if constexpr(K==Kind::Inverse)return _mm512_mask_xor_epi64(sum,lowMask,a,product);
            else return _mm512_mask_xor_epi64(sum,highMask,b,product);
        }
    }
}

template<unsigned Offset,Kind K>
HDBW_INLINE void transform16(__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3) noexcept {
    if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
        pair<normalized(3,Offset),K>(z0,z2);pair<normalized(3,Offset),K>(z1,z3);
        pair<normalized(2,Offset),K>(z0,z1);pair<normalized(2,Offset^8),K>(z2,z3);
        z0=smallStage<2,Offset,K>(z0);z1=smallStage<2,Offset^4,K>(z1);
        z2=smallStage<2,Offset^8,K>(z2);z3=smallStage<2,Offset^12,K>(z3);
        z0=smallStage<1,Offset,K>(z0);z1=smallStage<1,Offset^4,K>(z1);
        z2=smallStage<1,Offset^8,K>(z2);z3=smallStage<1,Offset^12,K>(z3);
    }else{
        z0=smallStage<1,Offset,K>(z0);z1=smallStage<1,Offset^4,K>(z1);
        z2=smallStage<1,Offset^8,K>(z2);z3=smallStage<1,Offset^12,K>(z3);
        z0=smallStage<2,Offset,K>(z0);z1=smallStage<2,Offset^4,K>(z1);
        z2=smallStage<2,Offset^8,K>(z2);z3=smallStage<2,Offset^12,K>(z3);
        pair<normalized(2,Offset),K>(z0,z1);pair<normalized(2,Offset^8),K>(z2,z3);
        pair<normalized(3,Offset),K>(z0,z2);pair<normalized(3,Offset),K>(z1,z3);
    }
}

HDBW_INLINE __m512i repeatRecordSum(__m512i x) noexcept {
    x=_mm512_xor_si512(x,_mm512_shuffle_i32x4(x,x,0x4e));
    return _mm512_xor_si512(x,_mm512_shuffle_i32x4(x,x,0xb1));
}
HDBW_INLINE void mix(State16& state,const Mixer16& mixer) noexcept {
    const auto a=wide16::multiply512(state.values[0],mixer.v.data());
    const auto b=wide16::multiply512(state.values[1],mixer.v.data()+4);
    const auto c=wide16::multiply512(state.values[2],mixer.v.data()+8);
    const auto d=wide16::multiply512(state.values[3],mixer.v.data()+12);
    const auto dot=repeatRecordSum(_mm512_xor_si512(_mm512_xor_si512(a,b),_mm512_xor_si512(c,d)));
    state.values[0]=_mm512_xor_si512(state.values[0],wide16::multiply512(dot,mixer.u.data()));
    state.values[1]=_mm512_xor_si512(state.values[1],wide16::multiply512(dot,mixer.u.data()+4));
    state.values[2]=_mm512_xor_si512(state.values[2],wide16::multiply512(dot,mixer.u.data()+8));
    state.values[3]=_mm512_xor_si512(state.values[3],wide16::multiply512(dot,mixer.u.data()+12));
}

template<bool Routed,class Index>
HDBW_INLINE __m512i input4(const Block* input,const Index* route,unsigned offset) noexcept {
    if constexpr(!Routed)return _mm512_loadu_si512(input+offset);
    else{
        static_assert(std::is_same_v<Index,std::uint16_t>||std::is_same_v<Index,std::uint32_t>);
        if constexpr(HDBW_HW_GATHER){
            if constexpr(std::is_same_v<Index,std::uint16_t>){
                const auto packed=_mm_loadl_epi64(reinterpret_cast<const __m128i*>(route+offset));
                const auto duplicate=_mm_unpacklo_epi16(packed,packed);
                const auto word=_mm256_slli_epi32(_mm256_cvtepu16_epi32(duplicate),1);
                const auto slots=_mm256_or_si256(word,_mm256_setr_epi32(0,1,0,1,0,1,0,1));
                return _mm512_i32gather_epi64(slots,input,8);
            }else{
                // Widen before doubling so the full uint32_t record-index
                // range remains positive and does not overflow signed slots.
                const auto packed=_mm_loadu_si128(reinterpret_cast<const __m128i*>(route+offset));
                const auto duplicate=_mm256_permutevar8x32_epi32(_mm256_broadcastsi128_si256(packed),
                    _mm256_setr_epi32(0,0,1,1,2,2,3,3));
                const auto word=_mm512_slli_epi64(_mm512_cvtepu32_epi64(duplicate),1);
                const auto slots=_mm512_or_si512(word,_mm512_setr_epi64(0,1,0,1,0,1,0,1));
                return _mm512_i64gather_epi64(slots,input,8);
            }
        }
        auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset])));
        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset+1])),1);
        x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset+2])),2);
        return _mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+route[offset+3])),3);
    }
}

template<unsigned Tile,bool Routed,class Index>
HDBW_INLINE void tile(const Block* input,const Index* route,Block* output,
    __m512i q0,__m512i q1,__m512i q2,__m512i q3,
    __m512i& c0,__m512i& c1,__m512i& c2,__m512i& c3) noexcept {
    static_assert(Tile<4);
    constexpr unsigned begin=16*Tile,offset=128+begin;
    auto x0=input4<Routed>(input,route,begin),x1=input4<Routed>(input,route,begin+4);
    auto x2=input4<Routed>(input,route,begin+8),x3=input4<Routed>(input,route,begin+12);
    transform16<offset,Kind::Forward>(q0,q1,q2,q3);
    _mm512_storeu_si512(output+begin,_mm512_xor_si512(x0,q0));
    _mm512_storeu_si512(output+begin+4,_mm512_xor_si512(x1,q1));
    _mm512_storeu_si512(output+begin+8,_mm512_xor_si512(x2,q2));
    _mm512_storeu_si512(output+begin+12,_mm512_xor_si512(x3,q3));
    transform16<offset,Kind::ForwardTranspose>(x0,x1,x2,x3);
    c0=_mm512_xor_si512(c0,x0);c1=_mm512_xor_si512(c1,x1);
    c2=_mm512_xor_si512(c2,x2);c3=_mm512_xor_si512(c3,x3);
}

template<bool Routed,class Index>
HDBW_INLINE void packet(const Block* input,const Index* route,Block* output,
                        State16& state,const Mixer16& mixer) noexcept {
    auto q0=state.values[0],q1=state.values[1],q2=state.values[2],q3=state.values[3];
    mix(state,mixer); // Old state remains in q0..q3 for emission, not M*q.
    transform16<0,Kind::Inverse>(q0,q1,q2,q3);
    auto c0=_mm512_setzero_si512(),c1=c0,c2=c0,c3=c0;
    tile<0,Routed>(input,route,output,q0,q1,q2,q3,c0,c1,c2,c3);
    tile<1,Routed>(input,route,output,q0,q1,q2,q3,c0,c1,c2,c3);
    tile<2,Routed>(input,route,output,q0,q1,q2,q3,c0,c1,c2,c3);
    tile<3,Routed>(input,route,output,q0,q1,q2,q3,c0,c1,c2,c3);
    transform16<0,Kind::InverseTranspose>(c0,c1,c2,c3);
    state.values[0]=_mm512_xor_si512(state.values[0],c0);
    state.values[1]=_mm512_xor_si512(state.values[1],c1);
    state.values[2]=_mm512_xor_si512(state.values[2],c2);
    state.values[3]=_mm512_xor_si512(state.values[3],c3);
}
} // namespace detail

// Standalone map hooks for basis/oracle checks. Input/output may overlap
// only when exact alias is specifically stated. Neither may overlap state.
HDBW_INLINE void emit64(const Block* inputState16,Block* output64) noexcept {
    auto q0=_mm512_loadu_si512(inputState16),q1=_mm512_loadu_si512(inputState16+4);
    auto q2=_mm512_loadu_si512(inputState16+8),q3=_mm512_loadu_si512(inputState16+12);
    detail::transform16<0,detail::Kind::Inverse>(q0,q1,q2,q3);
#define HDBW_EMIT(OFFSET) { auto a=q0,b=q1,c=q2,d=q3; \
    detail::transform16<128+OFFSET,detail::Kind::Forward>(a,b,c,d); \
    _mm512_storeu_si512(output64+OFFSET,a);_mm512_storeu_si512(output64+OFFSET+4,b); \
    _mm512_storeu_si512(output64+OFFSET+8,c);_mm512_storeu_si512(output64+OFFSET+12,d); }
    HDBW_EMIT(0);HDBW_EMIT(16);HDBW_EMIT(32);HDBW_EMIT(48);
#undef HDBW_EMIT
}
HDBW_INLINE void inject16(const Block* input64,Block* outputState16) noexcept {
    auto c0=_mm512_setzero_si512(),c1=c0,c2=c0,c3=c0;
#define HDBW_INJECT(OFFSET) { auto a=_mm512_loadu_si512(input64+OFFSET),b=_mm512_loadu_si512(input64+OFFSET+4); \
    auto c=_mm512_loadu_si512(input64+OFFSET+8),d=_mm512_loadu_si512(input64+OFFSET+12); \
    detail::transform16<128+OFFSET,detail::Kind::ForwardTranspose>(a,b,c,d); \
    c0=_mm512_xor_si512(c0,a);c1=_mm512_xor_si512(c1,b); \
    c2=_mm512_xor_si512(c2,c);c3=_mm512_xor_si512(c3,d); }
    HDBW_INJECT(0);HDBW_INJECT(16);HDBW_INJECT(32);HDBW_INJECT(48);
#undef HDBW_INJECT
    detail::transform16<0,detail::Kind::InverseTranspose>(c0,c1,c2,c3);
    _mm512_storeu_si512(outputState16,c0);_mm512_storeu_si512(outputState16+4,c1);
    _mm512_storeu_si512(outputState16+8,c2);_mm512_storeu_si512(outputState16+12,c3);
}

// Exactly64 records, no alignment requirement on payload or route. Raw
// exact input/output alias is supported; other overlap is unsupported.
// State and mixer must be disjoint from all payload buffers.
HDBW_INLINE void packetForwardRaw(const Block* input64,Block* output64,
                                  State16& state,const Mixer16& mixer) noexcept {
    detail::packet<false,std::uint16_t>(input64,nullptr,output64,state,mixer);
}
HDBW_INLINE void packetForward(const Block* input64,Block* output64,
                               State16& state,const Mixer16& mixer) noexcept {
    packetForwardRaw(input64,output64,state,mixer);
}
// Every route entry indexes a record relative to region. Routed output must
// be disjoint from the entire readable region (a permutation can read later).
template<class Index>
HDBW_INLINE void packetForwardRouted(const Block* region,const Index* route64,
                                     Block* output64,State16& state,const Mixer16& mixer) noexcept {
    detail::packet<true>(region,route64,output64,state,mixer);
}

#undef HDBW_INLINE
} // namespace spin::research::char2::hdbw
