#pragma once
// Retained fused parity/output circuit, shared by native-packing experiments.
#include "K16CodeSign.h"
namespace spin::research::k16codesign::nativeouter {
struct PlaneValues { __m512i v[8]; };
template<unsigned Plane>
static SPIN_FORCEINLINE PlaneValues parityValues(const __m512i* packed) {
    static_assert(Plane<4);
    const auto m2=_mm512_set1_epi64(0x0204080320408030ULL);
    const auto m5=_mm512_set1_epi64(0x050a070e50a070e0ULL);
    const auto m8=_mm512_set1_epi64(0x0803060c803060c0ULL);
    const auto m11=_mm512_set1_epi64(0x0b050a07b050a070ULL);
    const auto m12=_mm512_set1_epi64(0x0c0b050ac0b050a0ULL);
    const auto m13=_mm512_set1_epi64(0x0d090102d0901020ULL);
    auto x0=_mm512_load_si512(packed+32+Plane);
    auto x1=_mm512_load_si512(packed+36+Plane);
    auto x2=_mm512_load_si512(packed+40+Plane);
    auto x3=_mm512_load_si512(packed+44+Plane);
    auto x4=_mm512_load_si512(packed+48+Plane);
    auto x5=_mm512_load_si512(packed+52+Plane);
    auto x6=_mm512_load_si512(packed+56+Plane);
    auto x7=_mm512_load_si512(packed+60+Plane);
    x0=_mm512_xor_si512(x0,x1);x2=_mm512_xor_si512(x2,x3);
    x4=_mm512_xor_si512(x4,x5);x6=_mm512_xor_si512(x6,x7);
    x0=_mm512_xor_si512(x0,x2);x1=_mm512_xor_si512(x1,x3);
    x4=_mm512_xor_si512(x4,x6);x5=_mm512_xor_si512(x5,x7);
    x0=_mm512_xor_si512(x0,x4);x1=_mm512_xor_si512(x1,x5);
    x2=_mm512_xor_si512(x2,x6);x3=_mm512_xor_si512(x3,x7);
    auto y0=x0;
    auto y1=_mm512_xor_si512(x1,_mm512_gf2p8affine_epi64_epi8(x0,m11,0));
    auto y2=_mm512_xor_si512(x2,_mm512_gf2p8affine_epi64_epi8(x0,m2,0));
    auto y3=_mm512_xor_si512(x3,_mm512_gf2p8affine_epi64_epi8(x0,m13,0));
    y3=_mm512_xor_si512(y3,_mm512_gf2p8affine_epi64_epi8(x1,m2,0));
    y3=_mm512_xor_si512(y3,_mm512_gf2p8affine_epi64_epi8(x2,m11,0));
    auto y4=_mm512_xor_si512(x4,_mm512_gf2p8affine_epi64_epi8(x0,m5,0));
    auto y5=_mm512_xor_si512(x5,_mm512_gf2p8affine_epi64_epi8(x0,m12,0));
    y5=_mm512_xor_si512(y5,_mm512_gf2p8affine_epi64_epi8(x1,m5,0));
    y5=_mm512_xor_si512(y5,_mm512_gf2p8affine_epi64_epi8(x4,m11,0));
    auto y6=_mm512_xor_si512(x6,_mm512_gf2p8affine_epi64_epi8(x0,m11,0));
    y6=_mm512_xor_si512(y6,_mm512_gf2p8affine_epi64_epi8(x2,m5,0));
    y6=_mm512_xor_si512(y6,_mm512_gf2p8affine_epi64_epi8(x4,m2,0));
    auto y7=_mm512_xor_si512(x7,_mm512_gf2p8affine_epi64_epi8(x0,m8,0));
    y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x1,m11,0));
    y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x2,m12,0));
    y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x3,m5,0));
    y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x4,m13,0));
    y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x5,m2,0));
    y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x6,m11,0));
    y0=_mm512_xor_si512(y0,y1);y2=_mm512_xor_si512(y2,y3);
    y4=_mm512_xor_si512(y4,y5);y6=_mm512_xor_si512(y6,y7);
    y0=_mm512_xor_si512(y0,y2);y1=_mm512_xor_si512(y1,y3);
    y4=_mm512_xor_si512(y4,y6);y5=_mm512_xor_si512(y5,y7);
    y0=_mm512_xor_si512(y0,y4);y1=_mm512_xor_si512(y1,y5);
    y2=_mm512_xor_si512(y2,y6);y3=_mm512_xor_si512(y3,y7);
    return {{
        _mm512_xor_si512(_mm512_load_si512(packed+0+Plane),y0),
        _mm512_xor_si512(_mm512_load_si512(packed+4+Plane),y1),
        _mm512_xor_si512(_mm512_load_si512(packed+8+Plane),y2),
        _mm512_xor_si512(_mm512_load_si512(packed+12+Plane),y3),
        _mm512_xor_si512(_mm512_load_si512(packed+16+Plane),y4),
        _mm512_xor_si512(_mm512_load_si512(packed+20+Plane),y5),
        _mm512_xor_si512(_mm512_load_si512(packed+24+Plane),y6),
        _mm512_xor_si512(_mm512_load_si512(packed+28+Plane),y7)}};
}
template<unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storePair(const PlaneValues& a,const PlaneValues& b,Block* out) {
    static_assert(Symbol<8 && (Plane==0 || Plane==2));
    const auto i0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
    const auto i1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x0=_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0);
    const auto x1=_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0);
    _mm512_storeu_si512(out+32*Plane+4*Symbol,_mm512_permutex2var_epi8(x0,i0,x1));
    _mm512_storeu_si512(out+32*(Plane+1)+4*Symbol,_mm512_permutex2var_epi8(x0,i1,x1));
}
template<unsigned Plane>
static SPIN_NOINLINE void finishPair(const __m512i* packed,Block* out) {
    const auto a=parityValues<Plane>(packed);
    const auto b=parityValues<Plane+1>(packed);
    storePair<0,Plane>(a,b,out);storePair<1,Plane>(a,b,out);
    storePair<2,Plane>(a,b,out);storePair<3,Plane>(a,b,out);
    storePair<4,Plane>(a,b,out);storePair<5,Plane>(a,b,out);
    storePair<6,Plane>(a,b,out);storePair<7,Plane>(a,b,out);
}
}

