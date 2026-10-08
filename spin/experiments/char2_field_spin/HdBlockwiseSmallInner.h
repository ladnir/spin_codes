#pragma once
#include "HdBlockwiseInner.h"

namespace spin::research::char2::hdbwsmall {

// Isolated small-state COST CONTROLS, not distance certificates.
// U={0,...,T-1}, V={128,...,128+S-1}; A interpolates V values to U,
// C=A^T, and the live state stays in the original V evaluation basis.
//   y=x+A*q; q_next=q+u*(v^T*q)+C*x.
// Fixed maps act over embedded AES GF256. The rank-one mixer acts over
// genuine GF65536, independently on the eight field lanes in each record.
// The frozen s16 header supplies only its exact packed butterfly/load leaves.
// T defaults to64 for compatibility; supported T are16,32,64 with T>=2*S.
// This last condition ensures C*A=0: the Lagrange products have degree
// at most2*S-2, below T-1, and their sum over the additive space U vanishes.
inline constexpr unsigned packetWidth=64;
using Scalar=hdbw::Scalar;
using Coefficient=hdbw::Coefficient;

template<unsigned S> struct alignas(64) State {
    static_assert(S==4||S==8,"small blockwise state must have 4 or 8 records");
    __m512i values[S/4]{};
};
template<unsigned S> struct Mixer {
    static_assert(S==4||S==8,"small blockwise state must have 4 or 8 records");
    std::array<Coefficient,S> u{},v{};
};
static_assert(sizeof(State<4>)==4*sizeof(Block));
static_assert(sizeof(State<8>)==8*sizeof(Block));
static_assert(sizeof(Mixer<4>)==32&&sizeof(Mixer<8>)==64);

#if defined(_MSC_VER)
#define HDBWS_INLINE __forceinline
#else
#define HDBWS_INLINE inline __attribute__((always_inline))
#endif

template<unsigned S,unsigned T=64>
inline constexpr std::uint8_t aCoefficient(unsigned output,unsigned state) noexcept {
    static_assert(S==4||S==8);
    static_assert((T==16||T==32||T==64)&&T>=2*S);
    std::uint8_t numerator=1,denominator=1;
    for(unsigned j=0;j!=S;++j)if(j!=state){
        numerator=sfMultiply(numerator,static_cast<std::uint8_t>(output^(128+j)));
        denominator=sfMultiply(denominator,static_cast<std::uint8_t>(state^j));
    }
    return sfMultiply(numerator,hdbw::gf8Inverse(denominator));
}
template<unsigned S,unsigned T=64>
inline constexpr std::uint8_t cCoefficient(unsigned state,unsigned input) noexcept {
    return aCoefficient<S,T>(input,state);
}

// Setup only. Sample u uniformly in F^S\{0} and v uniformly conditioned
// on dot(v,u)!=1; individual zero coordinates must remain permitted.
template<unsigned S>
inline bool validMixer(const Scalar* u,const Scalar* v) noexcept {
    static_assert(S==4||S==8);
    Scalar dot=0,nonzero=0;
    for(unsigned i=0;i!=S;++i){
        nonzero|=u[i];dot^=medium32::quadraticMultiply(u[i],v[i]);
    }
    return nonzero!=0&&dot!=1;
}
template<unsigned S>
inline Mixer<S> prepareMixer(const Scalar* u,const Scalar* v){
    if(!validMixer<S>(u,v))throw std::invalid_argument("blockwise mixer requires nonzero u and dot(v,u)!=1");
    Mixer<S> result;
    wide16::prepare(u,result.u.data(),S);wide16::prepare(v,result.v.data(),S);
    return result;
}
template<unsigned S>
HDBWS_INLINE void loadState(State<S>& state,const Block* input) noexcept {
    state.values[0]=_mm512_loadu_si512(input);
    if constexpr(S==8)state.values[1]=_mm512_loadu_si512(input+4);
}
template<unsigned S>
HDBWS_INLINE void storeState(Block* output,const State<S>& state) noexcept {
    _mm512_storeu_si512(output,state.values[0]);
    if constexpr(S==8)_mm512_storeu_si512(output+4,state.values[1]);
}

namespace detail {
using Kind=hdbw::detail::Kind;
enum class Operation { Emit,Inject,Packet };

// One or two ZMMs carry the entire S-point transform. The S=4 second
// argument is a compile-time dead value; it never loads a fifth record.
template<unsigned S,unsigned Offset,Kind K>
HDBWS_INLINE void transform(__m512i& z0,__m512i& z1) noexcept {
    static_assert(S==4||S==8);
    if constexpr(K==Kind::Forward||K==Kind::InverseTranspose){
        if constexpr(S==8)hdbw::detail::pair<hdbw::detail::normalized(2,Offset),K>(z0,z1);
        z0=hdbw::detail::smallStage<2,Offset,K>(z0);
        if constexpr(S==8)z1=hdbw::detail::smallStage<2,Offset^4,K>(z1);
        z0=hdbw::detail::smallStage<1,Offset,K>(z0);
        if constexpr(S==8)z1=hdbw::detail::smallStage<1,Offset^4,K>(z1);
    }else{
        z0=hdbw::detail::smallStage<1,Offset,K>(z0);
        if constexpr(S==8)z1=hdbw::detail::smallStage<1,Offset^4,K>(z1);
        z0=hdbw::detail::smallStage<2,Offset,K>(z0);
        if constexpr(S==8)z1=hdbw::detail::smallStage<2,Offset^4,K>(z1);
        if constexpr(S==8)hdbw::detail::pair<hdbw::detail::normalized(2,Offset),K>(z0,z1);
    }
}

template<unsigned S>
HDBWS_INLINE void mix(State<S>& state,const Mixer<S>& mixer) noexcept {
    auto dot=wide16::multiply512(state.values[0],mixer.v.data());
    if constexpr(S==8)dot=_mm512_xor_si512(dot,wide16::multiply512(state.values[1],mixer.v.data()+4));
    dot=hdbw::detail::repeatRecordSum(dot);
    state.values[0]=_mm512_xor_si512(state.values[0],wide16::multiply512(dot,mixer.u.data()));
    if constexpr(S==8)state.values[1]=_mm512_xor_si512(state.values[1],wide16::multiply512(dot,mixer.u.data()+4));
}

template<unsigned S,unsigned T,unsigned Tile,Operation Op,bool Routed,class Index>
HDBWS_INLINE void tile(const Block* input,const Index* route,Block* output,
    __m512i q0,__m512i q1,__m512i& c0,__m512i& c1) noexcept {
    static_assert(S==4||S==8);
    static_assert((T==16||T==32||T==64)&&T>=2*S);
    static_assert(Tile<T/S);
    constexpr unsigned begin=S*Tile,offset=128+begin;
    auto x0=_mm512_setzero_si512(),x1=x0;
    if constexpr(Op!=Operation::Emit){
        x0=hdbw::detail::input4<Routed>(input,route,begin);
        if constexpr(S==8)x1=hdbw::detail::input4<Routed>(input,route,begin+4);
    }
    if constexpr(Op!=Operation::Inject){
        transform<S,offset,Kind::Forward>(q0,q1);
        if constexpr(Op==Operation::Packet){
            _mm512_storeu_si512(output+begin,_mm512_xor_si512(x0,q0));
            if constexpr(S==8)_mm512_storeu_si512(output+begin+4,_mm512_xor_si512(x1,q1));
        }else{
            _mm512_storeu_si512(output+begin,q0);
            if constexpr(S==8)_mm512_storeu_si512(output+begin+4,q1);
        }
    }
    if constexpr(Op!=Operation::Emit){
        // Injection uses retained raw inputs, including with exact alias.
        transform<S,offset,Kind::ForwardTranspose>(x0,x1);
        c0=_mm512_xor_si512(c0,x0);
        if constexpr(S==8)c1=_mm512_xor_si512(c1,x1);
    }
}

template<unsigned S,unsigned T,Operation Op,bool Routed,class Index>
HDBWS_INLINE void tiles(const Block* input,const Index* route,Block* output,
    __m512i q0,__m512i q1,__m512i& c0,__m512i& c1) noexcept {
    // The cosets and their GFNI constants are explicitly compile-time.
#define HDBWS_TILE(I) tile<S,T,I,Op,Routed>(input,route,output,q0,q1,c0,c1)
    HDBWS_TILE(0);HDBWS_TILE(1);
    if constexpr(T/S>=4){HDBWS_TILE(2);HDBWS_TILE(3);}
    if constexpr(T/S>=8){
        HDBWS_TILE(4);HDBWS_TILE(5);HDBWS_TILE(6);HDBWS_TILE(7);
    }
    if constexpr(T/S>=16){
        HDBWS_TILE(8);HDBWS_TILE(9);HDBWS_TILE(10);HDBWS_TILE(11);
        HDBWS_TILE(12);HDBWS_TILE(13);HDBWS_TILE(14);HDBWS_TILE(15);
    }
#undef HDBWS_TILE
}

template<unsigned S,unsigned T,bool Routed,class Index>
HDBWS_INLINE void packet(const Block* input,const Index* route,Block* output,
                         State<S>& state,const Mixer<S>& mixer) noexcept {
    auto q0=state.values[0],q1=_mm512_setzero_si512();
    if constexpr(S==8)q1=state.values[1];
    mix(state,mixer);
    transform<S,0,Kind::Inverse>(q0,q1);
    auto c0=_mm512_setzero_si512(),c1=c0;
    tiles<S,T,Operation::Packet,Routed>(input,route,output,q0,q1,c0,c1);
    transform<S,0,Kind::InverseTranspose>(c0,c1);
    state.values[0]=_mm512_xor_si512(state.values[0],c0);
    if constexpr(S==8)state.values[1]=_mm512_xor_si512(state.values[1],c1);
}
} // namespace detail

// Basis/scalar-oracle hooks. emit reads exactly S records and writes T;
// inject reads T records and writes S. Both support exact input/output alias.
template<unsigned S,unsigned T=64>
HDBWS_INLINE void emit(const Block* inputState,Block* output) noexcept {
    static_assert(S==4||S==8);
    auto q0=_mm512_loadu_si512(inputState),q1=_mm512_setzero_si512();
    if constexpr(S==8)q1=_mm512_loadu_si512(inputState+4);
    detail::transform<S,0,detail::Kind::Inverse>(q0,q1);
    auto unused0=_mm512_setzero_si512(),unused1=unused0;
    detail::tiles<S,T,detail::Operation::Emit,false,std::uint16_t>(
        nullptr,nullptr,output,q0,q1,unused0,unused1);
}
template<unsigned S>
HDBWS_INLINE void emit64(const Block* inputState,Block* output64) noexcept {
    emit<S,64>(inputState,output64);
}
template<unsigned S,unsigned T=64>
HDBWS_INLINE void inject(const Block* input,Block* outputState) noexcept {
    static_assert(S==4||S==8);
    const auto zero=_mm512_setzero_si512();
    auto c0=zero,c1=zero;
    detail::tiles<S,T,detail::Operation::Inject,false,std::uint16_t>(
        input,nullptr,nullptr,zero,zero,c0,c1);
    detail::transform<S,0,detail::Kind::InverseTranspose>(c0,c1);
    _mm512_storeu_si512(outputState,c0);
    if constexpr(S==8)_mm512_storeu_si512(outputState+4,c1);
}

// Exactly T records. Payload and route need not be aligned. State/mixer
// must be disjoint from payload; raw exact input/output alias is supported.
template<unsigned S,unsigned T=64>
HDBWS_INLINE void packetForwardRaw(const Block* input,Block* output,
                                   State<S>& state,const Mixer<S>& mixer) noexcept {
    detail::packet<S,T,false,std::uint16_t>(input,nullptr,output,state,mixer);
}
template<unsigned S,unsigned T=64>
HDBWS_INLINE void packetForward(const Block* input,Block* output,
                                State<S>& state,const Mixer<S>& mixer) noexcept {
    packetForwardRaw<S,T>(input,output,state,mixer);
}
// Route entries index records relative to region. Routed output must be
// disjoint from the entire readable region. HDBW_HW_GATHER selects the same
// optional input4 implementation as the frozen s16 header; default is0.
template<unsigned S,unsigned T=64,class Index>
HDBWS_INLINE void packetForwardRouted(const Block* region,const Index* route,
                                      Block* output,State<S>& state,const Mixer<S>& mixer) noexcept {
    detail::packet<S,T,true>(region,route,output,state,mixer);
}

#undef HDBWS_INLINE
} // namespace spin::research::char2::hdbwsmall
