#include "OuterTraffic.h"
#include "../k16_codesign_100us/kernel/NativeOuterFinish.h"
#include "../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"
#include <stdexcept>

namespace spin::research::packet8wide24opt {
namespace {
using Coefficients=std::array<std::uint8_t,9>;
using k16codesign::nativeouter::PlaneValues;

static SPIN_FORCEINLINE __m512i vx(__m512i a,__m512i b) {
    return _mm512_xor_si512(a,b);
}
static SPIN_FORCEINLINE __m512i vx3(__m512i a,__m512i b,__m512i c) {
    return _mm512_ternarylogic_epi64(a,b,c,0x96);
}

// Exactly the retained 15-GFNI parity circuit, excluding the systematic
// inputs. Stride and Plane are compile-time layout choices only.
template<unsigned Stride,unsigned Plane>
static SPIN_FORCEINLINE PlaneValues parity(const __m512i* packed) {
    static_assert(Plane<Stride);
    const auto m2=_mm512_set1_epi64(0x0204080320408030ULL);
    const auto m4=_mm512_set1_epi64(0x0408030640803060ULL);
    const auto m5=_mm512_set1_epi64(0x050a070e50a070e0ULL);
    const auto m8=_mm512_set1_epi64(0x0803060c803060c0ULL);
    const auto m11=_mm512_set1_epi64(0x0b050a07b050a070ULL);
    const auto m12=_mm512_set1_epi64(0x0c0b050ac0b050a0ULL);
    const auto m13=_mm512_set1_epi64(0x0d090102d0901020ULL);
    auto x0=_mm512_load_si512(packed+0*Stride+Plane);
    auto x1=_mm512_load_si512(packed+1*Stride+Plane);
    auto x2=_mm512_load_si512(packed+2*Stride+Plane);
    auto x3=_mm512_load_si512(packed+3*Stride+Plane);
    auto x4=_mm512_load_si512(packed+4*Stride+Plane);
    auto x5=_mm512_load_si512(packed+5*Stride+Plane);
    auto x6=_mm512_load_si512(packed+6*Stride+Plane);
    auto x7=_mm512_load_si512(packed+7*Stride+Plane);
    x0=vx(x0,x1);x2=vx(x2,x3);x4=vx(x4,x5);x6=vx(x6,x7);
    x0=vx(x0,x2);x1=vx(x1,x3);x4=vx(x4,x6);x5=vx(x5,x7);
    x0=vx(x0,x4);x1=vx(x1,x5);x2=vx(x2,x6);x3=vx(x3,x7);
    const auto t2=_mm512_gf2p8affine_epi64_epi8(x0,m2,0);
    const auto t4=_mm512_gf2p8affine_epi64_epi8(x0,m4,0);
    const auto t8=_mm512_gf2p8affine_epi64_epi8(x0,m8,0);
    const auto t5=vx(t4,x0),t11=vx3(t8,t2,x0),t12=vx(t8,t4),t13=vx(t12,x0);
    auto y0=x0,y1=vx(x1,t11),y2=vx(x2,t2);
    auto y3=vx3(vx(x3,t13),_mm512_gf2p8affine_epi64_epi8(x1,m2,0),_mm512_gf2p8affine_epi64_epi8(x2,m11,0));
    auto y4=vx(x4,t5);
    auto y5=vx3(vx(x5,t12),_mm512_gf2p8affine_epi64_epi8(x1,m5,0),_mm512_gf2p8affine_epi64_epi8(x4,m11,0));
    auto y6=vx3(vx(x6,t11),_mm512_gf2p8affine_epi64_epi8(x2,m5,0),_mm512_gf2p8affine_epi64_epi8(x4,m2,0));
    auto y7=vx3(vx(x7,t8),_mm512_gf2p8affine_epi64_epi8(x1,m11,0),_mm512_gf2p8affine_epi64_epi8(x2,m12,0));
    y7=vx3(y7,_mm512_gf2p8affine_epi64_epi8(x3,m5,0),_mm512_gf2p8affine_epi64_epi8(x4,m13,0));
    y7=vx3(y7,_mm512_gf2p8affine_epi64_epi8(x5,m2,0),_mm512_gf2p8affine_epi64_epi8(x6,m11,0));
    y0=vx(y0,y1);y2=vx(y2,y3);y4=vx(y4,y5);y6=vx(y6,y7);
    y0=vx(y0,y2);y1=vx(y1,y3);y4=vx(y4,y6);y5=vx(y5,y7);
    y0=vx(y0,y4);y1=vx(y1,y5);y2=vx(y2,y6);y3=vx(y3,y7);
    return {{y0,y1,y2,y3,y4,y5,y6,y7}};
}

template<unsigned Stride,unsigned Plane>
static SPIN_FORCEINLINE void transformParity(__m512i* packed) {
    const auto p=parity<Stride,Plane>(packed);
    _mm512_store_si512(packed+0*Stride+Plane,p.v[0]);
    _mm512_store_si512(packed+1*Stride+Plane,p.v[1]);
    _mm512_store_si512(packed+2*Stride+Plane,p.v[2]);
    _mm512_store_si512(packed+3*Stride+Plane,p.v[3]);
    _mm512_store_si512(packed+4*Stride+Plane,p.v[4]);
    _mm512_store_si512(packed+5*Stride+Plane,p.v[5]);
    _mm512_store_si512(packed+6*Stride+Plane,p.v[6]);
    _mm512_store_si512(packed+7*Stride+Plane,p.v[7]);
}

struct Four { __m512i a,b,c,d; };
struct Eight { __m512i a,b,c,d,e,f,g,h; };

template<unsigned Half>
static SPIN_FORCEINLINE Four mixHalf(const Block* input,const std::uint8_t* coeff) {
    static_assert(Half<2);
    auto a=_mm512_load_si512(input+4*Half);
    auto b=_mm512_load_si512(input+8+4*Half);
    auto c=_mm512_load_si512(input+16+4*Half);
    auto d=_mm512_load_si512(input+24+4*Half);
    rs::tower32byte::applyMultiply(a,b,c,d,coeff);
    return {a,b,c,d};
}
static SPIN_FORCEINLINE Eight mixFull(const Block* input,const std::uint8_t* coeff) {
    auto a=_mm512_load_si512(input),b=_mm512_load_si512(input+4);
    auto c=_mm512_load_si512(input+8),d=_mm512_load_si512(input+12);
    auto e=_mm512_load_si512(input+16),f=_mm512_load_si512(input+20);
    auto g=_mm512_load_si512(input+24),h=_mm512_load_si512(input+28);
    rs::tower32byte::applyMultiply(a,c,e,g,coeff);
    rs::tower32byte::applyMultiply(b,d,f,h,coeff);
    return {a,b,c,d,e,f,g,h};
}

// Exact flat9 MUL circuit from RandomizerVariants.cpp, privately instantiated
// here so the measured original variants and the other source stay unchanged.
template<unsigned Index>
static SPIN_FORCEINLINE __m512i flatProduct(__m512i value,const std::uint8_t* c) {
    return _mm512_gf2p8mul_epi8(value,_mm512_set1_epi8(char(c[Index])));
}
static SPIN_FORCEINLINE void flat9(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,
    const std::uint8_t* c) {
    const auto x01=vx(x0,x1),x23=vx(x2,x3);
    const auto p0=flatProduct<0>(x0,c),p1=flatProduct<1>(x1,c),p2=flatProduct<2>(x01,c);
    const auto q0=flatProduct<3>(x2,c),q1=flatProduct<4>(x3,c),q2=flatProduct<5>(x23,c);
    const auto r0=flatProduct<6>(vx(x0,x2),c),r1=flatProduct<7>(vx(x1,x3),c);
    const auto r2=flatProduct<8>(vx(x01,x23),c);
    const auto pq=vx(p0,q0),pr=vx(p0,r0);
    x0=vx3(pq,p1,q1);x1=vx3(pq,p2,q2);
    x2=vx3(pr,p1,r1);x3=vx3(pr,p2,r2);
}
template<bool Identity=false>
static SPIN_FORCEINLINE Eight mixFullFlat(const Block* input,const std::uint8_t* c) {
    auto a=_mm512_load_si512(input),b=_mm512_load_si512(input+4);
    auto d0=_mm512_load_si512(input+8),d1=_mm512_load_si512(input+12);
    auto e=_mm512_load_si512(input+16),f=_mm512_load_si512(input+20);
    auto g=_mm512_load_si512(input+24),h=_mm512_load_si512(input+28);
    if constexpr(!Identity) {
        flat9(a,d0,e,g,c);flat9(b,d1,f,h,c);
    }
    return {a,b,d0,d1,e,f,g,h};
}
static SPIN_FORCEINLINE void storeFour(__m512i* out,Four x) {
    _mm512_store_si512(out,x.a);_mm512_store_si512(out+1,x.b);
    _mm512_store_si512(out+2,x.c);_mm512_store_si512(out+3,x.d);
}
static SPIN_FORCEINLINE void storeEight(__m512i* out,Eight x) {
    _mm512_store_si512(out,x.a);_mm512_store_si512(out+1,x.b);
    _mm512_store_si512(out+2,x.c);_mm512_store_si512(out+3,x.d);
    _mm512_store_si512(out+4,x.e);_mm512_store_si512(out+5,x.f);
    _mm512_store_si512(out+6,x.g);_mm512_store_si512(out+7,x.h);
}

static SPIN_FORCEINLINE __m512i unpackIndex0() {
    return _mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
}
static SPIN_FORCEINLINE __m512i unpackIndex1() {
    return _mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
}
template<unsigned Symbol,unsigned Byte>
static SPIN_FORCEINLINE void unpackPair(__m512i low,__m512i high,Block* out) {
    static_assert(Symbol<8&&Byte<4);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto a=_mm512_gf2p8affine_epi64_epi8(basis,low,0);
    const auto b=_mm512_gf2p8affine_epi64_epi8(basis,high,0);
    _mm512_storeu_si512(out+64*Byte+4*Symbol,_mm512_permutex2var_epi8(a,unpackIndex0(),b));
    _mm512_storeu_si512(out+64*Byte+32+4*Symbol,_mm512_permutex2var_epi8(a,unpackIndex1(),b));
}
template<unsigned Symbol,unsigned Byte,unsigned Half>
static SPIN_FORCEINLINE void unpackHalf(__m512i value,Block* out) {
    static_assert(Symbol<8&&Byte<4&&Half<2);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x=_mm512_gf2p8affine_epi64_epi8(basis,value,0);
    constexpr __mmask8 mask=Half?0xaa:0x55;
    _mm512_mask_storeu_epi64(out+64*Byte+4*Symbol,mask,_mm512_permutexvar_epi8(unpackIndex0(),x));
    _mm512_mask_storeu_epi64(out+64*Byte+32+4*Symbol,mask,_mm512_permutexvar_epi8(unpackIndex1(),x));
}

template<unsigned Symbol,unsigned Byte,unsigned Half>
static SPIN_FORCEINLINE void emitHalf(const PlaneValues& p,const __m512i* packed,Block* output) {
    unpackHalf<Symbol,Byte,Half>(vx(p.v[Symbol],_mm512_load_si512(packed+4*Symbol+Byte)),output);
}
template<unsigned Byte,unsigned Half>
static SPIN_FORCEINLINE void finishHalf(const __m512i* packed,Block* output) {
    const auto p=parity<4,Byte>(packed+32);
    emitHalf<0,Byte,Half>(p,packed,output);emitHalf<1,Byte,Half>(p,packed,output);
    emitHalf<2,Byte,Half>(p,packed,output);emitHalf<3,Byte,Half>(p,packed,output);
    emitHalf<4,Byte,Half>(p,packed,output);emitHalf<5,Byte,Half>(p,packed,output);
    emitHalf<6,Byte,Half>(p,packed,output);emitHalf<7,Byte,Half>(p,packed,output);
}
template<unsigned Half>
static SPIN_FORCEINLINE void allSymbolsHalf(const Block* input,Block* output,
    const Coefficients* coefficients,__m512i* packed) {
    for(unsigned symbol=0;symbol<16;++symbol)
        storeFour(packed+4*symbol,mixHalf<Half>(input+32*symbol,coefficients[symbol].data()));
    finishHalf<0,Half>(packed,output);finishHalf<1,Half>(packed,output);
    finishHalf<2,Half>(packed,output);finishHalf<3,Half>(packed,output);
}
static SPIN_NOINLINE void groupHalf(const Block* __restrict input,Block* __restrict output,
    const Coefficients* __restrict coefficients) {
    alignas(64) __m512i packed[64];
    allSymbolsHalf<0>(input,output,coefficients,packed);
    allSymbolsHalf<1>(input,output,coefficients,packed);
}

template<unsigned Symbol>
static SPIN_FORCEINLINE void systematicFull(const Block* input,Block* output,
    const Coefficients* coefficients,const __m512i* packed) {
    const auto x=mixFull(input+32*Symbol,coefficients[Symbol].data());
    const auto* p=packed+8*Symbol;
    unpackPair<Symbol,0>(vx(x.a,_mm512_load_si512(p)),vx(x.b,_mm512_load_si512(p+1)),output);
    unpackPair<Symbol,1>(vx(x.c,_mm512_load_si512(p+2)),vx(x.d,_mm512_load_si512(p+3)),output);
    unpackPair<Symbol,2>(vx(x.e,_mm512_load_si512(p+4)),vx(x.f,_mm512_load_si512(p+5)),output);
    unpackPair<Symbol,3>(vx(x.g,_mm512_load_si512(p+6)),vx(x.h,_mm512_load_si512(p+7)),output);
}
static SPIN_NOINLINE void groupParityFirst(const Block* __restrict input,Block* __restrict output,
    const Coefficients* __restrict coefficients) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=8;symbol<16;++symbol)
        storeEight(packed+8*(symbol-8),mixFull(input+32*symbol,coefficients[symbol].data()));
    transformParity<8,0>(packed);transformParity<8,1>(packed);
    transformParity<8,2>(packed);transformParity<8,3>(packed);
    transformParity<8,4>(packed);transformParity<8,5>(packed);
    transformParity<8,6>(packed);transformParity<8,7>(packed);
    systematicFull<0>(input,output,coefficients,packed);systematicFull<1>(input,output,coefficients,packed);
    systematicFull<2>(input,output,coefficients,packed);systematicFull<3>(input,output,coefficients,packed);
    systematicFull<4>(input,output,coefficients,packed);systematicFull<5>(input,output,coefficients,packed);
    systematicFull<6>(input,output,coefficients,packed);systematicFull<7>(input,output,coefficients,packed);
}

template<unsigned Symbol,bool OmitZero>
static SPIN_FORCEINLINE void systematicFullFlat(const Block* input,Block* output,
    const Coefficients* coefficients,const __m512i* packed) {
    const auto x=mixFullFlat<OmitZero&&Symbol==0>(input+32*Symbol,coefficients[Symbol].data());
    const auto* p=packed+8*Symbol;
    unpackPair<Symbol,0>(vx(x.a,_mm512_load_si512(p)),vx(x.b,_mm512_load_si512(p+1)),output);
    unpackPair<Symbol,1>(vx(x.c,_mm512_load_si512(p+2)),vx(x.d,_mm512_load_si512(p+3)),output);
    unpackPair<Symbol,2>(vx(x.e,_mm512_load_si512(p+4)),vx(x.f,_mm512_load_si512(p+5)),output);
    unpackPair<Symbol,3>(vx(x.g,_mm512_load_si512(p+6)),vx(x.h,_mm512_load_si512(p+7)),output);
}
template<bool OmitZero>
static SPIN_NOINLINE void groupParityFirstFlat(const Block* __restrict input,Block* __restrict output,
    const Coefficients* __restrict coefficients) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=8;symbol<16;++symbol)
        storeEight(packed+8*(symbol-8),mixFullFlat(input+32*symbol,coefficients[symbol].data()));
    transformParity<8,0>(packed);transformParity<8,1>(packed);
    transformParity<8,2>(packed);transformParity<8,3>(packed);
    transformParity<8,4>(packed);transformParity<8,5>(packed);
    transformParity<8,6>(packed);transformParity<8,7>(packed);
    systematicFullFlat<0,OmitZero>(input,output,coefficients,packed);systematicFullFlat<1,OmitZero>(input,output,coefficients,packed);
    systematicFullFlat<2,OmitZero>(input,output,coefficients,packed);systematicFullFlat<3,OmitZero>(input,output,coefficients,packed);
    systematicFullFlat<4,OmitZero>(input,output,coefficients,packed);systematicFullFlat<5,OmitZero>(input,output,coefficients,packed);
    systematicFullFlat<6,OmitZero>(input,output,coefficients,packed);systematicFullFlat<7,OmitZero>(input,output,coefficients,packed);
}

template<unsigned Symbol,unsigned Half>
static SPIN_FORCEINLINE void systematicHalf(const Block* input,Block* output,
    const Coefficients* coefficients,const __m512i* packed) {
    const auto x=mixHalf<Half>(input+32*Symbol,coefficients[Symbol].data());
    const auto* p=packed+4*Symbol;
    unpackHalf<Symbol,0,Half>(vx(x.a,_mm512_load_si512(p)),output);
    unpackHalf<Symbol,1,Half>(vx(x.b,_mm512_load_si512(p+1)),output);
    unpackHalf<Symbol,2,Half>(vx(x.c,_mm512_load_si512(p+2)),output);
    unpackHalf<Symbol,3,Half>(vx(x.d,_mm512_load_si512(p+3)),output);
}
template<unsigned Half>
static SPIN_FORCEINLINE void parityFirstHalf(const Block* input,Block* output,
    const Coefficients* coefficients,__m512i* packed) {
    for(unsigned symbol=8;symbol<16;++symbol)
        storeFour(packed+4*(symbol-8),mixHalf<Half>(input+32*symbol,coefficients[symbol].data()));
    transformParity<4,0>(packed);transformParity<4,1>(packed);
    transformParity<4,2>(packed);transformParity<4,3>(packed);
    systematicHalf<0,Half>(input,output,coefficients,packed);systematicHalf<1,Half>(input,output,coefficients,packed);
    systematicHalf<2,Half>(input,output,coefficients,packed);systematicHalf<3,Half>(input,output,coefficients,packed);
    systematicHalf<4,Half>(input,output,coefficients,packed);systematicHalf<5,Half>(input,output,coefficients,packed);
    systematicHalf<6,Half>(input,output,coefficients,packed);systematicHalf<7,Half>(input,output,coefficients,packed);
}
static SPIN_NOINLINE void groupParityFirstHalf(const Block* __restrict input,Block* __restrict output,
    const Coefficients* __restrict coefficients) {
    alignas(64) __m512i packed[32];
    parityFirstHalf<0>(input,output,coefficients,packed);
    parityFirstHalf<1>(input,output,coefficients,packed);
}

template<unsigned Variant>
static void run(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t group=0;group<plan.groups;++group) {
        const auto* input=scratch+packet8wide24::groupStride*group;
        auto* out=output+256*group;
        const auto* coefficients=plan.outerCoefficients.data()+16*group;
        if constexpr(Variant==1)groupHalf(input,out,coefficients);
        else if constexpr(Variant==2)groupParityFirst(input,out,coefficients);
        else if constexpr(Variant==3)groupParityFirstHalf(input,out,coefficients);
        else if constexpr(Variant==4)groupParityFirstFlat<false>(input,out,coefficients);
        else {static_assert(Variant==5);groupParityFirstFlat<true>(input,out,coefficients);}
    }
}
}

void outerTraffic(const Block* scratch,Block* output,const Plan& plan,unsigned variant) {
    switch(variant) {
    case 0:packet8wide24::outerFast(scratch,output,plan);return;
    case 1:run<1>(scratch,output,plan);return;
    case 2:run<2>(scratch,output,plan);return;
    case 3:run<3>(scratch,output,plan);return;
    case 4:run<4>(scratch,output,plan);return;
    case 5:run<5>(scratch,output,plan);return;
    default:throw std::invalid_argument("outerTraffic variant must be0..5");
    }
}
}
