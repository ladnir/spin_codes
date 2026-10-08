#pragma once
#include "HdBlockwiseSmallInner.h"

namespace spin::research::char2::hdbwshared {
// Exact beta64 shared map, t64/s4; no distance certificate or timing claim.
// Q holds global novel-basis coefficients. All tiles use incoming Q:
// y_alpha=FFT8_(alpha^64)(IFFT8_alpha(x_alpha)+(Q,0));
// q_next=(I+u*v^T)Q+sum_alpha high4(IFFT8_alpha(x_alpha)).
using State=hdbwsmall::State<4>;
using Mixer=hdbwsmall::Mixer<4>;
using Scalar=hdbwsmall::Scalar;
inline constexpr unsigned packetWidth=64,stateSize=4,beta=64;
inline Mixer prepareMixer(const Scalar* u,const Scalar* v){return hdbwsmall::prepareMixer<4>(u,v);}
inline void loadState(State& s,const Block* p) noexcept {hdbwsmall::loadState(s,p);}
inline void storeState(Block* p,const State& s) noexcept {hdbwsmall::storeState(p,s);}
#ifndef HDBSH_FFT_CONTROL
#define HDBSH_FFT_CONTROL 0
#endif
#if defined(_MSC_VER)
#define HDBSH_INLINE __forceinline
#else
#define HDBSH_INLINE inline __attribute__((always_inline))
#endif

namespace detail {
template<unsigned I> HDBSH_INLINE __m128i lane(__m512i x) noexcept {
    if constexpr(I==0)return _mm512_castsi512_si128(x);
    else return _mm512_extracti64x2_epi64(x,I);
}
HDBSH_INLINE __m512i pack(__m128i a,__m128i b,__m128i c,__m128i d) noexcept {
    auto x=_mm512_castsi128_si512(a);
    x=_mm512_inserti64x2(x,b,1);x=_mm512_inserti64x2(x,c,2);
    return _mm512_inserti64x2(x,d,3);
}
template<unsigned A,unsigned B,unsigned C,unsigned D>
HDBSH_INLINE __m512i permute(__m512i x) noexcept {
    return _mm512_shuffle_i32x4(x,x,_MM_SHUFFLE(D,C,B,A));
}
template<unsigned A,unsigned B,unsigned C,unsigned D>
HDBSH_INLINE __m512i select(__m512i x,__m512i y) noexcept {
    return _mm512_permutex2var_epi64(x,_mm512_setr_epi64(
        2*A,2*A+1,2*B,2*B+1,2*C,2*C+1,2*D,2*D+1),y);
}
template<unsigned A,unsigned B,unsigned C,unsigned D>
HDBSH_INLINE __m512i product(__m512i x) noexcept {
    constexpr auto r=std::uint64_t(0x0101010101010101ULL);
    return _mm512_gf2p8mul_epi8(x,_mm512_setr_epi64(
        r*A,r*A,r*B,r*B,r*C,r*C,r*D,r*D));
}
struct Moments {
    // mask order: 0,8,16,32 /24,40,48,- /1,9,17,33 /
    //             2,10,18,34 /4,12,20,36 /3,5,6,-
    __m512i v0{},v1{},v2{},v3{},v4{},v5{};
};

HDBSH_INLINE __m512i feedback(const Moments& m) noexcept {
    const auto p0=product<2,4,8,16>(pack(lane<0>(m.v3),lane<0>(m.v4),lane<1>(m.v0),lane<2>(m.v0)));
    const auto p1=product<32,6,28,120>(pack(lane<3>(m.v0),lane<0>(m.v4),lane<1>(m.v0),lane<2>(m.v0)));
    const auto p2=product<235,2,24,6>(pack(lane<3>(m.v0),lane<0>(m.v3),lane<0>(m.v4),lane<1>(m.v5)));
    const auto p3=product<8,224,28,48>(pack(lane<2>(m.v5),lane<1>(m.v0),lane<1>(m.v2),lane<1>(m.v3)));
    const auto p4=product<64,193,120,224>(pack(lane<1>(m.v4),lane<2>(m.v0),lane<2>(m.v2),lane<2>(m.v3)));
    const auto p5=product<155,54,100,235>(pack(lane<2>(m.v4),lane<0>(m.v1),lane<3>(m.v0),lane<3>(m.v2)));
    const auto p6=product<237,65,180,171>(pack(lane<3>(m.v3),lane<3>(m.v4),lane<1>(m.v1),lane<2>(m.v1)));
    const auto z0=lane<0>(m.v0);
    const auto z1=_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(lane<0>(m.v2),lane<0>(p0)),lane<1>(p0)),_mm_xor_si128(_mm_xor_si128(lane<2>(p0),lane<3>(p0)),lane<0>(p1)));
    const auto z2=_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(lane<0>(m.v3),lane<1>(p1)),lane<2>(p1)),_mm_xor_si128(lane<3>(p1),lane<0>(p2)));
    const auto z3=_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(lane<0>(m.v5),lane<1>(p2)),lane<2>(p2)),_mm_xor_si128(lane<3>(p2),lane<0>(p3))),_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(lane<1>(p3),lane<2>(p3)),lane<3>(p3)),_mm_xor_si128(lane<0>(p4),lane<1>(p4)))),_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(lane<2>(p4),lane<3>(p4)),lane<0>(p5)),_mm_xor_si128(lane<1>(p5),lane<2>(p5))),_mm_xor_si128(_mm_xor_si128(_mm_xor_si128(lane<3>(p5),lane<0>(p6)),lane<1>(p6)),_mm_xor_si128(lane<2>(p6),lane<3>(p6)))));
    return pack(_mm_xor_si128(_mm_xor_si128(z0,z1),_mm_xor_si128(z2,z3)),
        _mm_xor_si128(z0,z2),_mm_xor_si128(z0,z1),z0);
}

// Eight GFNI calls prepare all state-emission products. The bounded 1KiB
// output is reused by the eight input tiles. This is not a heap allocation.
HDBSH_INLINE void emission(__m512i q,Block* output) noexcept {
    const auto p0=product<2,4,8,16>(permute<1,1,1,1>(q));
    const auto p1=product<32,64,6,28>(permute<1,1,2,2>(q));
    const auto p2=product<120,161,235,6>(permute<2,2,2,3>(q));
    const auto p3=product<8,22,27,28>(permute<3,3,3,3>(q));
    const auto p4=product<48,49,54,64>(permute<3,3,3,3>(q));
    const auto p5=product<65,75,120,155>(permute<3,3,3,3>(q));
    const auto p6=product<161,171,180,213>(permute<3,3,3,3>(q));
    const auto p7=product<222,224,235,237>(permute<3,3,3,3>(q));
    const auto a0=_mm_xor_si128(_mm_xor_si128(lane<0>(q),lane<1>(p1)),_mm_xor_si128(lane<1>(p2),lane<3>(p6)));
    const auto a1=_mm_xor_si128(lane<1>(q),lane<0>(p6));
    const auto a2=_mm_xor_si128(_mm_xor_si128(lane<0>(p0),lane<2>(q)),lane<2>(p3));
    const auto a3=lane<3>(q);
    const auto a4=_mm_xor_si128(_mm_xor_si128(lane<1>(p0),lane<2>(p1)),lane<1>(p4));
    const auto a5=lane<3>(p2);
    const auto a6=lane<0>(p3);
    const auto a8=_mm_xor_si128(_mm_xor_si128(lane<2>(p0),lane<3>(p1)),lane<0>(p7));
    const auto a9=lane<3>(p3);
    const auto a10=lane<0>(p4);
    const auto a12=lane<3>(p4);
    const auto a16=_mm_xor_si128(_mm_xor_si128(lane<3>(p0),lane<0>(p2)),lane<1>(p3));
    const auto a17=lane<2>(p5);
    const auto a18=lane<1>(p7);
    const auto a20=lane<3>(p5);
    const auto a24=lane<2>(p4);
    const auto a32=_mm_xor_si128(_mm_xor_si128(lane<0>(p1),lane<2>(p2)),lane<1>(p5));
    const auto a33=lane<2>(p7);
    const auto a34=lane<3>(p7);
    const auto a36=lane<0>(p5);
    const auto a40=lane<2>(p6);
    const auto a48=lane<1>(p6);
    const auto z0_1=_mm_xor_si128(a1,a0);
    const auto z0_3=_mm_xor_si128(a3,a2);
    const auto z0_5=_mm_xor_si128(a5,a4);
    const auto z0_9=_mm_xor_si128(a9,a8);
    const auto z0_17=_mm_xor_si128(a17,a16);
    const auto z0_33=_mm_xor_si128(a33,a32);
    const auto z1_2=_mm_xor_si128(a2,a0);
    const auto z1_3=_mm_xor_si128(z0_3,z0_1);
    const auto z1_6=_mm_xor_si128(a6,a4);
    const auto z1_7=_mm_xor_si128(a6,z0_5);
    const auto z1_10=_mm_xor_si128(a10,a8);
    const auto z1_11=_mm_xor_si128(a10,z0_9);
    const auto z1_18=_mm_xor_si128(a18,a16);
    const auto z1_19=_mm_xor_si128(a18,z0_17);
    const auto z1_34=_mm_xor_si128(a34,a32);
    const auto z1_35=_mm_xor_si128(a34,z0_33);
    const auto z2_4=_mm_xor_si128(a4,a0);
    const auto z2_5=_mm_xor_si128(z0_5,z0_1);
    const auto z2_6=_mm_xor_si128(z1_6,z1_2);
    const auto z2_7=_mm_xor_si128(z1_7,z1_3);
    const auto z2_12=_mm_xor_si128(a12,a8);
    const auto z2_13=_mm_xor_si128(a12,z0_9);
    const auto z2_14=_mm_xor_si128(a12,z1_10);
    const auto z2_15=_mm_xor_si128(a12,z1_11);
    const auto z2_20=_mm_xor_si128(a20,a16);
    const auto z2_21=_mm_xor_si128(a20,z0_17);
    const auto z2_22=_mm_xor_si128(a20,z1_18);
    const auto z2_23=_mm_xor_si128(a20,z1_19);
    const auto z2_36=_mm_xor_si128(a36,a32);
    const auto z2_37=_mm_xor_si128(a36,z0_33);
    const auto z2_38=_mm_xor_si128(a36,z1_34);
    const auto z2_39=_mm_xor_si128(a36,z1_35);
    const auto z3_8=_mm_xor_si128(a8,a0);
    const auto z3_9=_mm_xor_si128(z0_9,z0_1);
    const auto z3_10=_mm_xor_si128(z1_10,z1_2);
    const auto z3_11=_mm_xor_si128(z1_11,z1_3);
    const auto z3_12=_mm_xor_si128(z2_12,z2_4);
    const auto z3_13=_mm_xor_si128(z2_13,z2_5);
    const auto z3_14=_mm_xor_si128(z2_14,z2_6);
    const auto z3_15=_mm_xor_si128(z2_15,z2_7);
    const auto z3_24=_mm_xor_si128(a24,a16);
    const auto z3_25=_mm_xor_si128(a24,z0_17);
    const auto z3_26=_mm_xor_si128(a24,z1_18);
    const auto z3_27=_mm_xor_si128(a24,z1_19);
    const auto z3_28=_mm_xor_si128(a24,z2_20);
    const auto z3_29=_mm_xor_si128(a24,z2_21);
    const auto z3_30=_mm_xor_si128(a24,z2_22);
    const auto z3_31=_mm_xor_si128(a24,z2_23);
    const auto z3_40=_mm_xor_si128(a40,a32);
    const auto z3_41=_mm_xor_si128(a40,z0_33);
    const auto z3_42=_mm_xor_si128(a40,z1_34);
    const auto z3_43=_mm_xor_si128(a40,z1_35);
    const auto z3_44=_mm_xor_si128(a40,z2_36);
    const auto z3_45=_mm_xor_si128(a40,z2_37);
    const auto z3_46=_mm_xor_si128(a40,z2_38);
    const auto z3_47=_mm_xor_si128(a40,z2_39);
    const auto z4_16=_mm_xor_si128(a16,a0);
    const auto z4_17=_mm_xor_si128(z0_17,z0_1);
    const auto z4_18=_mm_xor_si128(z1_18,z1_2);
    const auto z4_19=_mm_xor_si128(z1_19,z1_3);
    const auto z4_20=_mm_xor_si128(z2_20,z2_4);
    const auto z4_21=_mm_xor_si128(z2_21,z2_5);
    const auto z4_22=_mm_xor_si128(z2_22,z2_6);
    const auto z4_23=_mm_xor_si128(z2_23,z2_7);
    const auto z4_24=_mm_xor_si128(z3_24,z3_8);
    const auto z4_25=_mm_xor_si128(z3_25,z3_9);
    const auto z4_26=_mm_xor_si128(z3_26,z3_10);
    const auto z4_27=_mm_xor_si128(z3_27,z3_11);
    const auto z4_28=_mm_xor_si128(z3_28,z3_12);
    const auto z4_29=_mm_xor_si128(z3_29,z3_13);
    const auto z4_30=_mm_xor_si128(z3_30,z3_14);
    const auto z4_31=_mm_xor_si128(z3_31,z3_15);
    const auto z4_48=_mm_xor_si128(a48,a32);
    const auto z4_49=_mm_xor_si128(a48,z0_33);
    const auto z4_50=_mm_xor_si128(a48,z1_34);
    const auto z4_51=_mm_xor_si128(a48,z1_35);
    const auto z4_52=_mm_xor_si128(a48,z2_36);
    const auto z4_53=_mm_xor_si128(a48,z2_37);
    const auto z4_54=_mm_xor_si128(a48,z2_38);
    const auto z4_55=_mm_xor_si128(a48,z2_39);
    const auto z4_56=_mm_xor_si128(a48,z3_40);
    const auto z4_57=_mm_xor_si128(a48,z3_41);
    const auto z4_58=_mm_xor_si128(a48,z3_42);
    const auto z4_59=_mm_xor_si128(a48,z3_43);
    const auto z4_60=_mm_xor_si128(a48,z3_44);
    const auto z4_61=_mm_xor_si128(a48,z3_45);
    const auto z4_62=_mm_xor_si128(a48,z3_46);
    const auto z4_63=_mm_xor_si128(a48,z3_47);
    const auto z5_32=_mm_xor_si128(a32,a0);
    const auto z5_33=_mm_xor_si128(z0_33,z0_1);
    const auto z5_34=_mm_xor_si128(z1_34,z1_2);
    const auto z5_35=_mm_xor_si128(z1_35,z1_3);
    const auto z5_36=_mm_xor_si128(z2_36,z2_4);
    const auto z5_37=_mm_xor_si128(z2_37,z2_5);
    const auto z5_38=_mm_xor_si128(z2_38,z2_6);
    const auto z5_39=_mm_xor_si128(z2_39,z2_7);
    const auto z5_40=_mm_xor_si128(z3_40,z3_8);
    const auto z5_41=_mm_xor_si128(z3_41,z3_9);
    const auto z5_42=_mm_xor_si128(z3_42,z3_10);
    const auto z5_43=_mm_xor_si128(z3_43,z3_11);
    const auto z5_44=_mm_xor_si128(z3_44,z3_12);
    const auto z5_45=_mm_xor_si128(z3_45,z3_13);
    const auto z5_46=_mm_xor_si128(z3_46,z3_14);
    const auto z5_47=_mm_xor_si128(z3_47,z3_15);
    const auto z5_48=_mm_xor_si128(z4_48,z4_16);
    const auto z5_49=_mm_xor_si128(z4_49,z4_17);
    const auto z5_50=_mm_xor_si128(z4_50,z4_18);
    const auto z5_51=_mm_xor_si128(z4_51,z4_19);
    const auto z5_52=_mm_xor_si128(z4_52,z4_20);
    const auto z5_53=_mm_xor_si128(z4_53,z4_21);
    const auto z5_54=_mm_xor_si128(z4_54,z4_22);
    const auto z5_55=_mm_xor_si128(z4_55,z4_23);
    const auto z5_56=_mm_xor_si128(z4_56,z4_24);
    const auto z5_57=_mm_xor_si128(z4_57,z4_25);
    const auto z5_58=_mm_xor_si128(z4_58,z4_26);
    const auto z5_59=_mm_xor_si128(z4_59,z4_27);
    const auto z5_60=_mm_xor_si128(z4_60,z4_28);
    const auto z5_61=_mm_xor_si128(z4_61,z4_29);
    const auto z5_62=_mm_xor_si128(z4_62,z4_30);
    const auto z5_63=_mm_xor_si128(z4_63,z4_31);
    _mm512_store_si512(output+0,pack(a0,z0_1,z1_2,z1_3));
    _mm512_store_si512(output+4,pack(z2_4,z2_5,z2_6,z2_7));
    _mm512_store_si512(output+8,pack(z3_8,z3_9,z3_10,z3_11));
    _mm512_store_si512(output+12,pack(z3_12,z3_13,z3_14,z3_15));
    _mm512_store_si512(output+16,pack(z4_16,z4_17,z4_18,z4_19));
    _mm512_store_si512(output+20,pack(z4_20,z4_21,z4_22,z4_23));
    _mm512_store_si512(output+24,pack(z4_24,z4_25,z4_26,z4_27));
    _mm512_store_si512(output+28,pack(z4_28,z4_29,z4_30,z4_31));
    _mm512_store_si512(output+32,pack(z5_32,z5_33,z5_34,z5_35));
    _mm512_store_si512(output+36,pack(z5_36,z5_37,z5_38,z5_39));
    _mm512_store_si512(output+40,pack(z5_40,z5_41,z5_42,z5_43));
    _mm512_store_si512(output+44,pack(z5_44,z5_45,z5_46,z5_47));
    _mm512_store_si512(output+48,pack(z5_48,z5_49,z5_50,z5_51));
    _mm512_store_si512(output+52,pack(z5_52,z5_53,z5_54,z5_55));
    _mm512_store_si512(output+56,pack(z5_56,z5_57,z5_58,z5_59));
    _mm512_store_si512(output+60,pack(z5_60,z5_61,z5_62,z5_63));
}

template<unsigned Tile,bool Routed,class Index>
HDBSH_INLINE void cubeTile(const Block* input,const Index* route,Block* output,
                          const Block* emitted,Moments& m) noexcept {
    constexpr unsigned begin=8*Tile;
    const auto x0=hdbw::detail::input4<Routed>(input,route,begin);
    const auto x1=hdbw::detail::input4<Routed>(input,route,begin+4);
    // Three four-record faces, then six edges and the whole-cube sum.
    const auto d1=_mm512_xor_si512(select<0,2,4,6>(x0,x1),select<1,3,5,7>(x0,x1));
    const auto d2=_mm512_xor_si512(select<0,1,4,5>(x0,x1),select<2,3,6,7>(x0,x1));
    const auto d4=_mm512_xor_si512(x0,x1);
    const auto e3=_mm512_xor_si512(d1,permute<1,0,3,2>(d1));
    const auto e5=_mm512_xor_si512(d1,permute<2,3,0,1>(d1));
    const auto e6=_mm512_xor_si512(d2,permute<2,3,0,1>(d2));
    const auto e7=hdbw::detail::repeatRecordSum(d1);
    constexpr unsigned active=3|((Tile&1)?12:0)|((Tile&2)?48:0)|((Tile&4)?192:0);
    constexpr unsigned pairs=((Tile&3)==3?3:0)|((Tile&5)==5?12:0)|((Tile&6)==6?48:0);
    m.v0=_mm512_mask_xor_epi64(m.v0,active,m.v0,e7);
    if constexpr(pairs)m.v1=_mm512_mask_xor_epi64(m.v1,pairs,m.v1,e7);
    m.v2=_mm512_mask_xor_epi64(m.v2,active,m.v2,permute<1,1,1,1>(e6));
    m.v3=_mm512_mask_xor_epi64(m.v3,active,m.v3,permute<1,1,1,1>(e5));
    m.v4=_mm512_mask_xor_epi64(m.v4,active,m.v4,permute<2,2,2,2>(e3));
    m.v5=_mm512_xor_si512(m.v5,select<0,1,7,7>(select<3,7,3,7>(d4,d2),d1));

    const auto p1=product<142,142,142,142>(d1);
    const auto p2=product<4,4,4,4>(d2);
    const auto p4=product<188,188,188,188>(d4);
    const auto p35=product<62,62,137,137>(select<0,2,4,5>(e3,e5));
    // The last lane is unobserved; multiplying its duplicate by zero is free
    // within this existing packed GFNI instruction, not an extra product.
    const auto p67=product<192,192,80,0>(select<0,1,4,4>(e6,e7));
    const auto common=_mm512_xor_si512(_mm512_xor_si512(p4,permute<2,2,3,3>(p35)),
        _mm512_xor_si512(permute<0,1,0,1>(p67),permute<2,2,2,2>(p67)));
    const auto y0=_mm512_xor_si512(_mm512_xor_si512(x0,common),
        _mm512_xor_si512(permute<0,0,1,1>(p1),
        _mm512_xor_si512(permute<0,1,0,1>(p2),permute<0,0,0,0>(p35))));
    const auto y1=_mm512_xor_si512(_mm512_xor_si512(x1,common),
        _mm512_xor_si512(permute<2,2,3,3>(p1),
        _mm512_xor_si512(permute<2,3,2,3>(p2),permute<1,1,1,1>(p35))));
    _mm512_storeu_si512(output+begin,_mm512_xor_si512(y0,_mm512_load_si512(emitted+begin)));
    _mm512_storeu_si512(output+begin+4,_mm512_xor_si512(y1,_mm512_load_si512(emitted+begin+4)));
}

template<unsigned Tile,bool Routed,class Index>
HDBSH_INLINE void fftTile(const Block* input,const Index* route,Block* output,
                         __m512i q,__m512i& h) noexcept {
    constexpr unsigned begin=8*Tile;
    auto a=hdbw::detail::input4<Routed>(input,route,begin);
    auto b=hdbw::detail::input4<Routed>(input,route,begin+4);
    hdbwsmall::detail::transform<8,begin,hdbw::detail::Kind::Inverse>(a,b);
    h=_mm512_xor_si512(h,b);a=_mm512_xor_si512(a,q);
    hdbwsmall::detail::transform<8,begin^64,hdbw::detail::Kind::Forward>(a,b);
    _mm512_storeu_si512(output+begin,a);_mm512_storeu_si512(output+begin+4,b);
}

template<bool Fft,bool Routed,class Index>
HDBSH_INLINE void packet(const Block* input,const Index* route,Block* output,
                         State& state,const Mixer& mixer) noexcept {
    const auto q=state.values[0];
    hdbwsmall::detail::mix(state,mixer);
    if constexpr(Fft){
        auto h=_mm512_setzero_si512();
#define HDBSH_TILE(I) fftTile<I,Routed>(input,route,output,q,h)
        HDBSH_TILE(0);HDBSH_TILE(1);HDBSH_TILE(2);HDBSH_TILE(3);
        HDBSH_TILE(4);HDBSH_TILE(5);HDBSH_TILE(6);HDBSH_TILE(7);
#undef HDBSH_TILE
        state.values[0]=_mm512_xor_si512(state.values[0],h);
    }else{
        alignas(64) Block emitted[64];
        emission(q,emitted);
        Moments m{};
#define HDBSH_TILE(I) cubeTile<I,Routed>(input,route,output,emitted,m)
        HDBSH_TILE(0);HDBSH_TILE(1);HDBSH_TILE(2);HDBSH_TILE(3);
        HDBSH_TILE(4);HDBSH_TILE(5);HDBSH_TILE(6);HDBSH_TILE(7);
#undef HDBSH_TILE
        state.values[0]=_mm512_xor_si512(state.values[0],feedback(m));
    }
}
} // namespace detail

// Raw exact alias is supported. Routed output must be disjoint from the
// readable region. State and mixer are always disjoint from payload.
template<bool Fft=(HDBSH_FFT_CONTROL!=0)>
HDBSH_INLINE void packetForwardRaw(const Block* input,Block* output,
                                   State& state,const Mixer& mixer) noexcept {
    detail::packet<Fft,false,std::uint16_t>(input,nullptr,output,state,mixer);
}
template<bool Fft=(HDBSH_FFT_CONTROL!=0),class Index>
HDBSH_INLINE void packetForwardRouted(const Block* input,const Index* route,Block* output,
                                      State& state,const Mixer& mixer) noexcept {
    detail::packet<Fft,true>(input,route,output,state,mixer);
}
#undef HDBSH_INLINE
} // namespace spin::research::char2::hdbwshared
