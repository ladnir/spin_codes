// Selected native-field RS circuit, with the measured loop/unrolled policies.
#include "Plan.h"
#include "NativeOuterFinish.h"
#include "../packet/PacketLargeInner.h"
#include <cstdint>

namespace spin::detail::certified22 {
namespace {
using Coeff=std::array<std::uint8_t,3>;
using nativeouter::PlaneValues;
struct Pair { __m512i lo,hi; };
static SPIN_FORCEINLINE Pair quadraticMultiply(__m512i x0,__m512i x1,const std::uint8_t* c) {
    const auto p0=_mm512_gf2p8mul_epi8(x0,_mm512_set1_epi8(char(c[0])));
    const auto p1=_mm512_gf2p8mul_epi8(x1,_mm512_set1_epi8(char(c[1])));
    const auto p2=_mm512_gf2p8mul_epi8(_mm512_xor_si512(x0,x1),_mm512_set1_epi8(char(c[2])));
    return {_mm512_xor_si512(p0,p1),_mm512_xor_si512(p0,p2)};
}
static SPIN_FORCEINLINE void mix(const Block* in, __m512i* out, const Coeff& coeff) {
    using detail::packet::large::packIndex;
    const auto a = _mm512_loadu_si512(in);
    const auto b = _mm512_loadu_si512(in+4);
    const auto c = _mm512_loadu_si512(in+8);
    const auto d = _mm512_loadu_si512(in+12);
    const auto basis = _mm512_set1_epi64(0x0102040810204080ULL);
    const auto v0 = _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0);
    const auto v1 = _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0);
    const auto v2 = _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(c,packIndex<0>(),d),0);
    const auto v3 = _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(c,packIndex<8>(),d),0);
    const auto lo = quadraticMultiply(v0,v2,coeff.data());
    const auto hi = quadraticMultiply(v1,v3,coeff.data());
    out[0]=lo.lo;out[1]=hi.lo;out[2]=lo.hi;out[3]=hi.hi;
}

static SPIN_FORCEINLINE __m512i xor3(__m512i a,__m512i b,__m512i c) {
    return _mm512_ternarylogic_epi64(a,b,c,0x96);
}

// Same nilpotent-factor RS16 circuit as parityValues. The six distinct
// multiplications of x0 become only 2*x0, 4*x0, 8*x0; the other products are
// their binary sums. This removes three GFNI affine instructions per plane
// after ordinary common-subexpression elimination (18 -> 15), at the cost
// of four XOR/ternary-logic instructions. All four planes are independent.
template<unsigned Plane>
static SPIN_FORCEINLINE PlaneValues sharedParity(const __m512i* packed) {
    static_assert(Plane<4);
    const auto m2=_mm512_set1_epi64(0x0204080320408030ULL);
    const auto m4=_mm512_set1_epi64(0x0408030640803060ULL);
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
    const auto t2=_mm512_gf2p8affine_epi64_epi8(x0,m2,0);
    const auto t4=_mm512_gf2p8affine_epi64_epi8(x0,m4,0);
    const auto t8=_mm512_gf2p8affine_epi64_epi8(x0,m8,0);
    const auto t5=_mm512_xor_si512(t4,x0);
    const auto t11=xor3(t8,t2,x0);
    const auto t12=_mm512_xor_si512(t8,t4);
    const auto t13=_mm512_xor_si512(t12,x0);
    auto y0=x0;
    auto y1=_mm512_xor_si512(x1,t11);
    auto y2=_mm512_xor_si512(x2,t2);
    auto y3=xor3(_mm512_xor_si512(x3,t13),
        _mm512_gf2p8affine_epi64_epi8(x1,m2,0),_mm512_gf2p8affine_epi64_epi8(x2,m11,0));
    auto y4=_mm512_xor_si512(x4,t5);
    auto y5=xor3(_mm512_xor_si512(x5,t12),
        _mm512_gf2p8affine_epi64_epi8(x1,m5,0),_mm512_gf2p8affine_epi64_epi8(x4,m11,0));
    auto y6=xor3(_mm512_xor_si512(x6,t11),
        _mm512_gf2p8affine_epi64_epi8(x2,m5,0),_mm512_gf2p8affine_epi64_epi8(x4,m2,0));
    auto y7=xor3(_mm512_xor_si512(x7,t8),
        _mm512_gf2p8affine_epi64_epi8(x1,m11,0),_mm512_gf2p8affine_epi64_epi8(x2,m12,0));
    y7=xor3(y7,_mm512_gf2p8affine_epi64_epi8(x3,m5,0),_mm512_gf2p8affine_epi64_epi8(x4,m13,0));
    y7=xor3(y7,_mm512_gf2p8affine_epi64_epi8(x5,m2,0),_mm512_gf2p8affine_epi64_epi8(x6,m11,0));
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


template<bool Streaming,unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storePair(const PlaneValues& a,const PlaneValues& b,Block* out) {
    static_assert(Symbol<8 && (Plane==0 || Plane==2));
    if constexpr(!Streaming) {
        nativeouter::storePair<Symbol,Plane>(a,b,out);
    } else {
    const auto i0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
    const auto i1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x0=_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0);
    const auto x1=_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0);
    _mm512_stream_si512(reinterpret_cast<__m512i*>(out+32*Plane+4*Symbol),
        _mm512_permutex2var_epi8(x0,i0,x1));
    _mm512_stream_si512(reinterpret_cast<__m512i*>(out+32*(Plane+1)+4*Symbol),
        _mm512_permutex2var_epi8(x0,i1,x1));
    }
}

template<bool Streaming,unsigned Plane>
static SPIN_FORCEINLINE void finish(const __m512i* packed,Block* out) {
    const auto a=sharedParity<Plane>(packed);
    const auto b=sharedParity<Plane+1>(packed);
    storePair<Streaming,0,Plane>(a,b,out);storePair<Streaming,1,Plane>(a,b,out);
    storePair<Streaming,2,Plane>(a,b,out);storePair<Streaming,3,Plane>(a,b,out);
    storePair<Streaming,4,Plane>(a,b,out);storePair<Streaming,5,Plane>(a,b,out);
    storePair<Streaming,6,Plane>(a,b,out);storePair<Streaming,7,Plane>(a,b,out);
}

template<bool Streaming>
static SPIN_NOINLINE void groupLoop(const Block* __restrict in,Block* __restrict out,
    const Coeff* __restrict coeff) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol)
        mix(in+16*symbol,packed+4*symbol,coeff[symbol]);
    finish<Streaming,0>(packed,out);finish<Streaming,2>(packed,out);
}

template<bool Streaming>
static SPIN_NOINLINE void groupUnrolled(const Block* __restrict in,Block* __restrict out,
    const Coeff* __restrict coeff) {
    alignas(64) __m512i packed[64];
    mix(in+16*0,packed+4*0,coeff[0]);mix(in+16*1,packed+4*1,coeff[1]);
    mix(in+16*2,packed+4*2,coeff[2]);mix(in+16*3,packed+4*3,coeff[3]);
    mix(in+16*4,packed+4*4,coeff[4]);mix(in+16*5,packed+4*5,coeff[5]);
    mix(in+16*6,packed+4*6,coeff[6]);mix(in+16*7,packed+4*7,coeff[7]);
    mix(in+16*8,packed+4*8,coeff[8]);mix(in+16*9,packed+4*9,coeff[9]);
    mix(in+16*10,packed+4*10,coeff[10]);mix(in+16*11,packed+4*11,coeff[11]);
    mix(in+16*12,packed+4*12,coeff[12]);mix(in+16*13,packed+4*13,coeff[13]);
    mix(in+16*14,packed+4*14,coeff[14]);mix(in+16*15,packed+4*15,coeff[15]);
    finish<Streaming,0>(packed,out);finish<Streaming,2>(packed,out);
}
}
void outerFastCached(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t g=0;g<plan.groups;++g)
        groupLoop<false>(scratch+groupStride*g,output+128*g,plan.outerField.data()+16*g);
}
void outerFastCachedUnrolled(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t g=0;g<plan.groups;++g)
        groupUnrolled<false>(scratch+groupStride*g,output+128*g,plan.outerField.data()+16*g);
}
void outerFastStreaming(const Block* scratch,Block* output,const Plan& plan) {
    static_assert(sizeof(Block)==16 && (128*sizeof(Block))%64==0);
    if((reinterpret_cast<std::uintptr_t>(output)&63U)!=0) {
        outerFastCached(scratch,output,plan);
        return;
    }
    for(std::size_t g=0;g<plan.groups;++g)
        groupUnrolled<true>(scratch+groupStride*g,output+128*g,plan.outerField.data()+16*g);
    _mm_sfence();
}
}
