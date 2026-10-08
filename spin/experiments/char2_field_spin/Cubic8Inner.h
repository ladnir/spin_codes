#pragma once
#include "Subfield16.h"

namespace spin::research::char2 {

// Certified local geometry: actual AES GF4, eight halves of the old GF16
// phases. A=(u+v)^3+1, C=(u+v)^3+3, C*A=0, eta=2. C is NOT A^T.
inline constexpr std::uint8_t c8H4[4]={0,1,188,189};
inline constexpr std::uint8_t c8Representatives[4]={0,92,176,236};
inline constexpr std::uint8_t c8Offsets[4]={0,2,6,4};
constexpr std::uint8_t c8Label(unsigned phase,unsigned row) noexcept {
    return c8Offsets[phase/2]^c8Representatives[2*(phase%2)+row/4]^c8H4[row%4];
}
constexpr std::uint8_t c8Coefficient(unsigned phase,unsigned row,unsigned column,
                                    bool injection) noexcept {
    const auto x=std::uint8_t(c8Label(phase,row)^16^c8H4[column]);
    return sfMultiply(sfMultiply(x,x),x)^(injection?3:1);
}

template<std::uint8_t A,std::uint8_t B,std::uint8_t C,std::uint8_t D>
inline __m512i c8FourConstants() noexcept {
    constexpr auto repeat=[](std::uint8_t x){return static_cast<long long>(0x0101010101010101ULL*x);};
    return _mm512_setr_epi64(repeat(A),repeat(A),repeat(B),repeat(B),
                            repeat(C),repeat(C),repeat(D),repeat(D));
}

// Eight VSHUFI32X4 operations; each lane is one complete 128-bit record.
inline void c8TransposeRecords(__m512i a,__m512i b,__m512i c,__m512i d,
                                __m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3) noexcept {
    const auto ab0=_mm512_shuffle_i32x4(a,b,0x44),ab1=_mm512_shuffle_i32x4(a,b,0xee);
    const auto cd0=_mm512_shuffle_i32x4(c,d,0x44),cd1=_mm512_shuffle_i32x4(c,d,0xee);
    x0=_mm512_shuffle_i32x4(ab0,cd0,0x88);x1=_mm512_shuffle_i32x4(ab0,cd0,0xdd);
    x2=_mm512_shuffle_i32x4(ab1,cd1,0x88);x3=_mm512_shuffle_i32x4(ab1,cd1,0xdd);
}
inline void c8Zeta(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3) noexcept {
    const auto lo=_mm512_xor_si512(x0,x1),hi=_mm512_xor_si512(x2,x3);
    x0=_mm512_xor_si512(lo,hi);x1=_mm512_xor_si512(x1,x3);x2=hi;
}

template<unsigned Phase,bool Transpose,unsigned Lane,unsigned Zeta>
constexpr std::uint8_t c8PackedCoefficient() noexcept {
    static_assert(Phase<8&&Lane<4&&Zeta>=1&&Zeta<=3);
    constexpr bool injection=Lane<2?!Transpose:Transpose;
    constexpr unsigned first=4*(Lane%2);
    constexpr auto g3=c8Coefficient(Phase,first+3,0,injection);
    if constexpr(Zeta==3)return g3;
    else return c8Coefficient(Phase,first+Zeta,0,injection)^g3;
}
template<unsigned Phase,bool Transpose,unsigned Zeta>
inline __m512i c8PackedConstants() noexcept {
    return c8FourConstants<c8PackedCoefficient<Phase,Transpose,0,Zeta>(),
        c8PackedCoefficient<Phase,Transpose,1,Zeta>(),
        c8PackedCoefficient<Phase,Transpose,2,Zeta>(),
        c8PackedCoefficient<Phase,Transpose,3,Zeta>()>();
}

// Four simultaneous XOR-circulants: two injections and two emissions.
// Exactly five GFNI512 sites, three reused constant vectors, twelve XOR/XOR3.
// Inputs are copied before any result is exposed, so the packet hooks alias.
template<unsigned Phase,bool Transpose>
inline void c8PackedMaps(__m512i input0,__m512i input1,__m512i state,
                          __m512i& injection0,__m512i& injection1,
                          __m512i& emission0,__m512i& emission1) noexcept {
    __m512i z0,z1,z2,z3;
    c8TransposeRecords(input0,input1,state,state,z0,z1,z2,z3);
    c8Zeta(z0,z1,z2,z3);
    const auto g1=c8PackedConstants<Phase,Transpose,1>();
    const auto g2=c8PackedConstants<Phase,Transpose,2>();
    const auto g3=c8PackedConstants<Phase,Transpose,3>();
    const auto w1=_mm512_xor_si512(z1,_mm512_gf2p8mul_epi8(z0,g1));
    const auto w2=_mm512_xor_si512(z2,_mm512_gf2p8mul_epi8(z0,g2));
    auto w3=_mm512_ternarylogic_epi32(z3,_mm512_gf2p8mul_epi8(z2,g1),
                                    _mm512_gf2p8mul_epi8(z1,g2),0x96);
    w3=_mm512_xor_si512(w3,_mm512_gf2p8mul_epi8(z0,g3));
    z1=w1;z2=w2;z3=w3;c8Zeta(z0,z1,z2,z3);
    c8TransposeRecords(z0,z1,z2,z3,injection0,injection1,emission0,emission1);
}

// Input diagonal is external, exactly as in the existing packed kernels.
// Adjoint: x*=v+C^T*h; s*=D_s*(h+2*A^T*v). Seven GFNI including D_s/eta
// plus the two external input diagonals means NINE total per eight records.
// The adjoint has sixteen XOR/XOR3 (forward has fifteen), before routing.
template<unsigned Phase>
inline void c8InnerTransposePacked(const __m512i* input2,__m512i* output2,
                                    __m512i& state,__m512i diagonal) noexcept {
    const auto v0=input2[0],v1=input2[1],h=state;
    __m512i i0,i1,e0,e1;c8PackedMaps<Phase,true>(v0,v1,h,i0,i1,e0,e1);
    const auto moment=_mm512_gf2p8mul_epi8(_mm512_xor_si512(i0,i1),_mm512_set1_epi8(2));
    output2[0]=_mm512_xor_si512(v0,e0);output2[1]=_mm512_xor_si512(v1,e1);
    state=_mm512_gf2p8mul_epi8(_mm512_xor_si512(h,moment),diagonal);
}
template<unsigned Phase>
inline void c8InnerForwardPacked(const __m512i* input2,__m512i* output2,
                                  __m512i& state,__m512i diagonal) noexcept {
    const auto x0=input2[0],x1=input2[1];
    const auto r=_mm512_gf2p8mul_epi8(state,diagonal);
    const auto etaState=_mm512_gf2p8mul_epi8(r,_mm512_set1_epi8(2));
    __m512i i0,i1,e0,e1;c8PackedMaps<Phase,false>(x0,x1,etaState,i0,i1,e0,e1);
    output2[0]=_mm512_xor_si512(x0,e0);output2[1]=_mm512_xor_si512(x1,e1);
    state=_mm512_ternarylogic_epi32(r,i0,i1,0x96);
}

// Paired A/B schedule: all four lanes belong to INPUT maps or all four to
// STATE maps. Unlike c8PackedMaps, input moments have no register dependency
// on the incoming state. The map and its random setup remain identical.
template<unsigned FirstPhase,unsigned SecondPhase,bool Injection,unsigned Lane,unsigned Zeta>
constexpr std::uint8_t c8PairCoefficient() noexcept {
    static_assert(FirstPhase<8&&SecondPhase<8&&Lane<4&&Zeta>=1&&Zeta<=3);
    constexpr unsigned phase=Lane<2?FirstPhase:SecondPhase,first=4*(Lane%2);
    constexpr auto g3=c8Coefficient(phase,first+3,0,Injection);
    if constexpr(Zeta==3)return g3;
    else return c8Coefficient(phase,first+Zeta,0,Injection)^g3;
}
template<unsigned FirstPhase,unsigned SecondPhase,bool Injection,unsigned Zeta>
inline __m512i c8PairConstants() noexcept {
    return c8FourConstants<c8PairCoefficient<FirstPhase,SecondPhase,Injection,0,Zeta>(),
        c8PairCoefficient<FirstPhase,SecondPhase,Injection,1,Zeta>(),
        c8PairCoefficient<FirstPhase,SecondPhase,Injection,2,Zeta>(),
        c8PairCoefficient<FirstPhase,SecondPhase,Injection,3,Zeta>()>();
}
template<unsigned FirstPhase,unsigned SecondPhase,bool Injection>
inline void c8PairCirculants(__m512i a,__m512i b,__m512i c,__m512i d,
                              __m512i& out0,__m512i& out1,__m512i& out2,__m512i& out3) noexcept {
    __m512i z0,z1,z2,z3;c8TransposeRecords(a,b,c,d,z0,z1,z2,z3);
    c8Zeta(z0,z1,z2,z3);
    const auto g1=c8PairConstants<FirstPhase,SecondPhase,Injection,1>();
    const auto g2=c8PairConstants<FirstPhase,SecondPhase,Injection,2>();
    const auto g3=c8PairConstants<FirstPhase,SecondPhase,Injection,3>();
    const auto w1=_mm512_xor_si512(z1,_mm512_gf2p8mul_epi8(z0,g1));
    const auto w2=_mm512_xor_si512(z2,_mm512_gf2p8mul_epi8(z0,g2));
    auto w3=_mm512_ternarylogic_epi32(z3,_mm512_gf2p8mul_epi8(z2,g1),
                                    _mm512_gf2p8mul_epi8(z1,g2),0x96);
    w3=_mm512_xor_si512(w3,_mm512_gf2p8mul_epi8(z0,g3));
    z1=w1;z2=w2;z3=w3;c8Zeta(z0,z1,z2,z3);
    c8TransposeRecords(z0,z1,z2,z3,out0,out1,out2,out3);
}

// input4/output4 are ordered by reverse traversal: FirstPhase's eight
// records, then SecondPhase's eight. Exact hook alias remains supported.
// Per sixteen outputs:18GFNI (including four external input diagonals),
// 32XOR/XOR3,32record shuffles. Only two diagonal updates carry state.
template<unsigned FirstPhase,unsigned SecondPhase>
inline void c8InnerTransposePairedPacked(const __m512i* input4,__m512i* output4,
                                         __m512i& state,const __m512i* diagonal2) noexcept {
    static_assert(FirstPhase==SecondPhase+1&&FirstPhase%2==1);
    const auto v0=input4[0],v1=input4[1],v2=input4[2],v3=input4[3];
    __m512i i0,i1,i2,i3;
    // A^T moments for BOTH packets are ready independently of incoming h.
    c8PairCirculants<FirstPhase,SecondPhase,false>(v0,v1,v2,v3,i0,i1,i2,i3);
    const auto m0=_mm512_gf2p8mul_epi8(_mm512_xor_si512(i0,i1),_mm512_set1_epi8(2));
    const auto m1=_mm512_gf2p8mul_epi8(_mm512_xor_si512(i2,i3),_mm512_set1_epi8(2));
    const auto h0=state;
    const auto h1=_mm512_gf2p8mul_epi8(_mm512_xor_si512(h0,m0),diagonal2[0]);
    state=_mm512_gf2p8mul_epi8(_mm512_xor_si512(h1,m1),diagonal2[1]);
    __m512i e0,e1,e2,e3;
    c8PairCirculants<FirstPhase,SecondPhase,true>(h0,h0,h1,h1,e0,e1,e2,e3);
    output4[0]=_mm512_xor_si512(v0,e0);output4[1]=_mm512_xor_si512(v1,e1);
    output4[2]=_mm512_xor_si512(v2,e2);output4[3]=_mm512_xor_si512(v3,e3);
}
} // namespace spin::research::char2
