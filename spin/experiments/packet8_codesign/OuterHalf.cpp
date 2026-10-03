// Same shared-parity circuit as the frozen mode52 outer, with compact
// two-plane indexing. Process the two 64-bit payload halves separately.
#include "Packet8.h"
#include "../k16_codesign_100us/kernel/NativeOuterFinish.h"
#include "../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"

namespace spin::research::packet8 {
namespace {
using k16codesign::nativeouter::PlaneValues;
static SPIN_FORCEINLINE __m512i xor3(__m512i a,__m512i b,__m512i c) {
    return _mm512_ternarylogic_epi64(a,b,c,0x96);
}
template<unsigned Plane>
static SPIN_FORCEINLINE PlaneValues sharedParity(const __m512i* packed) {
    static_assert(Plane<2);
    const auto m2=_mm512_set1_epi64(0x0204080320408030ULL);
    const auto m4=_mm512_set1_epi64(0x0408030640803060ULL);
    const auto m5=_mm512_set1_epi64(0x050a070e50a070e0ULL);
    const auto m8=_mm512_set1_epi64(0x0803060c803060c0ULL);
    const auto m11=_mm512_set1_epi64(0x0b050a07b050a070ULL);
    const auto m12=_mm512_set1_epi64(0x0c0b050ac0b050a0ULL);
    const auto m13=_mm512_set1_epi64(0x0d090102d0901020ULL);
    auto x0=_mm512_load_si512(packed+16+Plane);
    auto x1=_mm512_load_si512(packed+18+Plane);
    auto x2=_mm512_load_si512(packed+20+Plane);
    auto x3=_mm512_load_si512(packed+22+Plane);
    auto x4=_mm512_load_si512(packed+24+Plane);
    auto x5=_mm512_load_si512(packed+26+Plane);
    auto x6=_mm512_load_si512(packed+28+Plane);
    auto x7=_mm512_load_si512(packed+30+Plane);
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
        _mm512_xor_si512(_mm512_load_si512(packed+2+Plane),y1),
        _mm512_xor_si512(_mm512_load_si512(packed+4+Plane),y2),
        _mm512_xor_si512(_mm512_load_si512(packed+6+Plane),y3),
        _mm512_xor_si512(_mm512_load_si512(packed+8+Plane),y4),
        _mm512_xor_si512(_mm512_load_si512(packed+10+Plane),y5),
        _mm512_xor_si512(_mm512_load_si512(packed+12+Plane),y6),
        _mm512_xor_si512(_mm512_load_si512(packed+14+Plane),y7)}};
}


template<unsigned Symbol,unsigned Plane,unsigned Half>
static SPIN_FORCEINLINE void storeHalf(const PlaneValues& values,Block* output) {
    static_assert(Symbol<8&&Plane<2&&Half<2);
    const auto i0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
    const auto i1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x=_mm512_gf2p8affine_epi64_epi8(basis,values.v[Symbol],0);
    constexpr __mmask8 mask=Half?0xaa:0x55;
    _mm512_mask_storeu_epi64(output+64*Plane+4*Symbol,mask,_mm512_permutexvar_epi8(i0,x));
    _mm512_mask_storeu_epi64(output+64*Plane+32+4*Symbol,mask,_mm512_permutexvar_epi8(i1,x));
}
template<unsigned Plane,unsigned Half>
static SPIN_FORCEINLINE void finishHalf(const __m512i* packed,Block* output) {
    const auto values=sharedParity<Plane>(packed);
    storeHalf<0,Plane,Half>(values,output);storeHalf<1,Plane,Half>(values,output);
    storeHalf<2,Plane,Half>(values,output);storeHalf<3,Plane,Half>(values,output);
    storeHalf<4,Plane,Half>(values,output);storeHalf<5,Plane,Half>(values,output);
    storeHalf<6,Plane,Half>(values,output);storeHalf<7,Plane,Half>(values,output);
}
template<unsigned Half>
static SPIN_FORCEINLINE void half(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,3>* __restrict coefficients,__m512i* packed) {
    for(unsigned symbol=0;symbol<16;++symbol) {
        const auto a=_mm512_load_si512(input+16*symbol+4*Half);
        const auto b=_mm512_load_si512(input+16*symbol+8+4*Half);
        const auto mixed=rs::tower32byte::quadraticMultiply(a,b,coefficients[symbol].data());
        packed[2*symbol]=mixed.lo;packed[2*symbol+1]=mixed.hi;
    }
    finishHalf<0,Half>(packed,output);
    finishHalf<1,Half>(packed,output);
}
static SPIN_NOINLINE void groupHalf(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,3>* __restrict coefficients) {
    alignas(64) __m512i packed[32];
    half<0>(input,output,coefficients,packed);
    half<1>(input,output,coefficients,packed);
}
}
void outerFastHalf(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t g=0;g<plan.outer.groups;++g)
        groupHalf(scratch+rs::groupStride*g,output+128*g,plan.native.field.data()+16*g);
}
}
