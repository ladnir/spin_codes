// Isolated A/B controls: identical maps, packing, parity, and output layout.
#include "RsWide.h"
#include "Tower32Randomizer.h"
#include "Tower32ByteRandomizer.h"
#include "../Rs16OuterFast.h"

namespace spin::research::rswide {
namespace {
// Exact pack half from retained fast16::mixSymbol, with its GL16 factored out.
static SPIN_FORCEINLINE void packHalf(const Block* input,__m512i* v) {
    const auto i0=_mm512_setr_epi64(0x0141115100401050ULL,0x0343135302421252ULL,0x0545155504441454ULL,0x0747175706461656ULL,0x0949195908481858ULL,0x0b4b1b5b0a4a1a5aULL,0x0d4d1d5d0c4c1c5cULL,0x0f4f1f5f0e4e1e5eULL);
    const auto i1=_mm512_setr_epi64(0x2161317120603070ULL,0x2363337322623272ULL,0x2565357524643474ULL,0x2767377726663676ULL,0x2969397928683878ULL,0x2b6b3b7b2a6a3a7aULL,0x2d6d3d7d2c6c3c7cULL,0x2f6f3f7f2e6e3e7eULL);
    const auto i2=_mm512_setr_epi64(0x0302434201004140ULL,0x0706474605044544ULL,0x0b0a4b4a09084948ULL,0x0f0e4f4e0d0c4d4cULL,0x1312535211105150ULL,0x1716575615145554ULL,0x1b1a5b5a19185958ULL,0x1f1e5f5e1d1c5d5cULL);
    const auto i3=_mm512_setr_epi64(0x2322636221206160ULL,0x2726676625246564ULL,0x2b2a6b6a29286968ULL,0x2f2e6f6e2d2c6d6cULL,0x3332737231307170ULL,0x3736777635347574ULL,0x3b3a7b7a39387978ULL,0x3f3e7f7e3d3c7d7cULL);
    const auto r0=_mm512_loadu_si512(input),r1=_mm512_loadu_si512(input+4);
    const auto r2=_mm512_loadu_si512(input+8),r3=_mm512_loadu_si512(input+12);
    const auto a0=_mm512_permutex2var_epi8(r0,i0,r1),a1=_mm512_permutex2var_epi8(r0,i1,r1);
    const auto a2=_mm512_permutex2var_epi8(r2,i0,r3),a3=_mm512_permutex2var_epi8(r2,i1,r3);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    v[0]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a0,i2,a2),0);
    v[1]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a0,i3,a2),0);
    v[2]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a1,i2,a3),0);
    v[3]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a1,i3,a3),0);
}
template<unsigned Out>
static SPIN_FORCEINLINE void mixByte(const __m512i* v,__m512i* packed,const std::uint64_t* coeff) {
    const auto a=_mm512_set1_epi64(coeff[4*Out]);
    const auto b=_mm512_set1_epi64(coeff[4*Out+1]);
    const auto c=_mm512_set1_epi64(coeff[4*Out+2]);
    const auto d=_mm512_set1_epi64(coeff[4*Out+3]);
    packed[0]=_mm512_xor_si512(
        _mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[0],a,0),_mm512_gf2p8affine_epi64_epi8(v[2],b,0)),
        _mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[4],c,0),_mm512_gf2p8affine_epi64_epi8(v[6],d,0)));
    packed[1]=_mm512_xor_si512(
        _mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[1],a,0),_mm512_gf2p8affine_epi64_epi8(v[3],b,0)),
        _mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[5],c,0),_mm512_gf2p8affine_epi64_epi8(v[7],d,0)));
}
template<Randomizer Mode>
static SPIN_FORCEINLINE void mixSymbol(const Block* input,__m512i* low,__m512i* high,const std::uint8_t* bytes) {
    __m512i v[8];packHalf(input,v);packHalf(input+16,v+4);
    if constexpr(Mode==Randomizer::DenseGl32) {
        const auto* coeff=reinterpret_cast<const std::uint64_t*>(bytes);
        mixByte<0>(v,low,coeff);mixByte<1>(v,low+2,coeff);
        mixByte<2>(v,high,coeff);mixByte<3>(v,high+2,coeff);
    } else {
        if constexpr(Mode==Randomizer::TowerByte32) {
            rs::tower32byte::applyMultiply(v[0],v[2],v[4],v[6],bytes);
            rs::tower32byte::applyMultiply(v[1],v[3],v[5],v[7],bytes);
        } else {
            const auto* coeff=reinterpret_cast<const std::uint64_t*>(bytes);
            rs::tower32::applyTranspose(v[0],v[2],v[4],v[6],coeff);
            rs::tower32::applyTranspose(v[1],v[3],v[5],v[7],coeff);
        }
        low[0]=v[0];low[1]=v[1];low[2]=v[2];low[3]=v[3];
        high[0]=v[4];high[1]=v[5];high[2]=v[6];high[3]=v[7];
    }
}
static SPIN_FORCEINLINE void finishHalf(__m512i* packed,Block* out) {
    rs::fast16::parityPlane<0>(packed);rs::fast16::parityPlane<1>(packed);
    rs::fast16::parityPlane<2>(packed);rs::fast16::parityPlane<3>(packed);
    rs::fast16::storeSymbol<0>(packed,out);rs::fast16::storeSymbol<1>(packed,out);
    rs::fast16::storeSymbol<2>(packed,out);rs::fast16::storeSymbol<3>(packed,out);
    rs::fast16::storeSymbol<4>(packed,out);rs::fast16::storeSymbol<5>(packed,out);
    rs::fast16::storeSymbol<6>(packed,out);rs::fast16::storeSymbol<7>(packed,out);
}

static SPIN_NOINLINE void finishHalfIsolated(__m512i* packed,Block* out) {
    finishHalf(packed,out);
}
template<Randomizer Mode,bool Pair,bool IsolateFinish>
static SPIN_NOINLINE void outerGroupVariant(const Block* __restrict in,Block* __restrict out,
    const std::uint8_t* __restrict coeff) {
    alignas(64) __m512i low[64],high[64];
    constexpr unsigned stride=Mode==Randomizer::DenseGl32?128:(Mode==Randomizer::TowerByte32?9:72);
    if constexpr(Pair) {
        for(unsigned symbol=0;symbol<16;symbol+=2) {
            mixSymbol<Mode>(in+32*symbol,low+4*symbol,high+4*symbol,coeff+stride*symbol);
            mixSymbol<Mode>(in+32*(symbol+1),low+4*(symbol+1),high+4*(symbol+1),coeff+stride*(symbol+1));
        }
    } else {
        for(unsigned symbol=0;symbol<16;++symbol)
            mixSymbol<Mode>(in+32*symbol,low+4*symbol,high+4*symbol,coeff+stride*symbol);
    }
    if constexpr(IsolateFinish) {
        finishHalfIsolated(low,out);finishHalfIsolated(high,out+128);
    } else {
        finishHalf(low,out);finishHalf(high,out+128);
    }
}
template<bool Pair,bool IsolateFinish>
void outerVariant(const Block* scratch,Block* output,const Plan& plan) {
    const auto* coeff=reinterpret_cast<const std::uint8_t*>(plan.compactCoefficients());
    if(plan.randomizer==Randomizer::DenseGl32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupVariant<Randomizer::DenseGl32,Pair,IsolateFinish>(
                scratch+groupStride*group,output+256*group,coeff+2048*group);
    } else if(plan.randomizer==Randomizer::TowerField32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupVariant<Randomizer::TowerField32,Pair,IsolateFinish>(
                scratch+groupStride*group,output+256*group,coeff+1152*group);
    } else {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupVariant<Randomizer::TowerByte32,Pair,IsolateFinish>(
                scratch+groupStride*group,output+256*group,coeff+144*group);
    }
}

template<unsigned Plane>
static SPIN_NOINLINE void parityPlaneIsolated(__m512i* packed) {
    rs::fast16::parityPlane<Plane>(packed);
}
static SPIN_NOINLINE void finishHalfPlaneIsolated(__m512i* packed,Block* out) {
    parityPlaneIsolated<0>(packed);parityPlaneIsolated<1>(packed);
    parityPlaneIsolated<2>(packed);parityPlaneIsolated<3>(packed);
    rs::fast16::storeSymbol<0>(packed,out);rs::fast16::storeSymbol<1>(packed,out);
    rs::fast16::storeSymbol<2>(packed,out);rs::fast16::storeSymbol<3>(packed,out);
    rs::fast16::storeSymbol<4>(packed,out);rs::fast16::storeSymbol<5>(packed,out);
    rs::fast16::storeSymbol<6>(packed,out);rs::fast16::storeSymbol<7>(packed,out);
}
// Exact retained factorization, returning final values instead of storing
// them into packed. The paired finish can go directly to final row stores.
struct PlaneValues { __m512i v[8]; };
template<unsigned plane> static SPIN_FORCEINLINE PlaneValues parityValues(const __m512i* packed){
static_assert(plane<4);
const auto m2=_mm512_set1_epi64(0x0204080320408030ULL);
const auto m5=_mm512_set1_epi64(0x050a070e50a070e0ULL);
const auto m8=_mm512_set1_epi64(0x0803060c803060c0ULL);
const auto m11=_mm512_set1_epi64(0x0b050a07b050a070ULL);
const auto m12=_mm512_set1_epi64(0x0c0b050ac0b050a0ULL);
const auto m13=_mm512_set1_epi64(0x0d090102d0901020ULL);
auto x0=_mm512_load_si512(packed+32+plane);
auto x1=_mm512_load_si512(packed+36+plane);
auto x2=_mm512_load_si512(packed+40+plane);
auto x3=_mm512_load_si512(packed+44+plane);
auto x4=_mm512_load_si512(packed+48+plane);
auto x5=_mm512_load_si512(packed+52+plane);
auto x6=_mm512_load_si512(packed+56+plane);
auto x7=_mm512_load_si512(packed+60+plane);
x0=_mm512_xor_si512(x0,x1);
x2=_mm512_xor_si512(x2,x3);
x4=_mm512_xor_si512(x4,x5);
x6=_mm512_xor_si512(x6,x7);
x0=_mm512_xor_si512(x0,x2);
x1=_mm512_xor_si512(x1,x3);
x4=_mm512_xor_si512(x4,x6);
x5=_mm512_xor_si512(x5,x7);
x0=_mm512_xor_si512(x0,x4);
x1=_mm512_xor_si512(x1,x5);
x2=_mm512_xor_si512(x2,x6);
x3=_mm512_xor_si512(x3,x7);
auto y0=x0;
auto y1=x1;
y1=_mm512_xor_si512(y1,_mm512_gf2p8affine_epi64_epi8(x0,m11,0));
auto y2=x2;
y2=_mm512_xor_si512(y2,_mm512_gf2p8affine_epi64_epi8(x0,m2,0));
auto y3=x3;
y3=_mm512_xor_si512(y3,_mm512_gf2p8affine_epi64_epi8(x0,m13,0));
y3=_mm512_xor_si512(y3,_mm512_gf2p8affine_epi64_epi8(x1,m2,0));
y3=_mm512_xor_si512(y3,_mm512_gf2p8affine_epi64_epi8(x2,m11,0));
auto y4=x4;
y4=_mm512_xor_si512(y4,_mm512_gf2p8affine_epi64_epi8(x0,m5,0));
auto y5=x5;
y5=_mm512_xor_si512(y5,_mm512_gf2p8affine_epi64_epi8(x0,m12,0));
y5=_mm512_xor_si512(y5,_mm512_gf2p8affine_epi64_epi8(x1,m5,0));
y5=_mm512_xor_si512(y5,_mm512_gf2p8affine_epi64_epi8(x4,m11,0));
auto y6=x6;
y6=_mm512_xor_si512(y6,_mm512_gf2p8affine_epi64_epi8(x0,m11,0));
y6=_mm512_xor_si512(y6,_mm512_gf2p8affine_epi64_epi8(x2,m5,0));
y6=_mm512_xor_si512(y6,_mm512_gf2p8affine_epi64_epi8(x4,m2,0));
auto y7=x7;
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x0,m8,0));
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x1,m11,0));
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x2,m12,0));
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x3,m5,0));
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x4,m13,0));
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x5,m2,0));
y7=_mm512_xor_si512(y7,_mm512_gf2p8affine_epi64_epi8(x6,m11,0));
y0=_mm512_xor_si512(y0,y1);
y2=_mm512_xor_si512(y2,y3);
y4=_mm512_xor_si512(y4,y5);
y6=_mm512_xor_si512(y6,y7);
y0=_mm512_xor_si512(y0,y2);
y1=_mm512_xor_si512(y1,y3);
y4=_mm512_xor_si512(y4,y6);
y5=_mm512_xor_si512(y5,y7);
y0=_mm512_xor_si512(y0,y4);
y1=_mm512_xor_si512(y1,y5);
y2=_mm512_xor_si512(y2,y6);
y3=_mm512_xor_si512(y3,y7);
return {{
    _mm512_xor_si512(_mm512_load_si512(packed+0+plane),y0),
    _mm512_xor_si512(_mm512_load_si512(packed+4+plane),y1),
    _mm512_xor_si512(_mm512_load_si512(packed+8+plane),y2),
    _mm512_xor_si512(_mm512_load_si512(packed+12+plane),y3),
    _mm512_xor_si512(_mm512_load_si512(packed+16+plane),y4),
    _mm512_xor_si512(_mm512_load_si512(packed+20+plane),y5),
    _mm512_xor_si512(_mm512_load_si512(packed+24+plane),y6),
    _mm512_xor_si512(_mm512_load_si512(packed+28+plane),y7)
}};
}

template<unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storePairSymbol(const PlaneValues& a,const PlaneValues& b,Block* out) {
    static_assert(Symbol<8 && (Plane==0 || Plane==2));
const auto index0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
const auto index1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);

    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x0=_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0);
    const auto x1=_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0);
    _mm512_storeu_si512(out+32*Plane+4*Symbol,_mm512_permutex2var_epi8(x0,index0,x1));
    _mm512_storeu_si512(out+32*(Plane+1)+4*Symbol,_mm512_permutex2var_epi8(x0,index1,x1));
}
template<unsigned Plane>
static SPIN_NOINLINE void finishPairFused(const __m512i* packed,Block* out) {
    static_assert(Plane==0 || Plane==2);
    const auto a=parityValues<Plane>(packed);
    const auto b=parityValues<Plane+1>(packed);
    storePairSymbol<0,Plane>(a,b,out);storePairSymbol<1,Plane>(a,b,out);
    storePairSymbol<2,Plane>(a,b,out);storePairSymbol<3,Plane>(a,b,out);
    storePairSymbol<4,Plane>(a,b,out);storePairSymbol<5,Plane>(a,b,out);
    storePairSymbol<6,Plane>(a,b,out);storePairSymbol<7,Plane>(a,b,out);
}
template<Randomizer Mode,bool FusedPair>
static SPIN_NOINLINE void outerGroupFurther(const Block* __restrict in,Block* __restrict out,
    const std::uint8_t* __restrict coeff) {
    alignas(64) __m512i low[64],high[64];
    constexpr unsigned stride=Mode==Randomizer::DenseGl32?128:(Mode==Randomizer::TowerByte32?9:72);
    for(unsigned symbol=0;symbol<16;++symbol)
        mixSymbol<Mode>(in+32*symbol,low+4*symbol,high+4*symbol,coeff+stride*symbol);
    if constexpr(FusedPair) {
        finishPairFused<0>(low,out);finishPairFused<2>(low,out);
        finishPairFused<0>(high,out+128);finishPairFused<2>(high,out+128);
    } else {
        finishHalfPlaneIsolated(low,out);finishHalfPlaneIsolated(high,out+128);
    }
}
template<bool FusedPair>
void outerFurther(const Block* scratch,Block* output,const Plan& plan) {
    const auto* coeff=reinterpret_cast<const std::uint8_t*>(plan.compactCoefficients());
    if(plan.randomizer==Randomizer::DenseGl32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupFurther<Randomizer::DenseGl32,FusedPair>(
                scratch+groupStride*group,output+256*group,coeff+2048*group);
    } else if(plan.randomizer==Randomizer::TowerField32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupFurther<Randomizer::TowerField32,FusedPair>(
                scratch+groupStride*group,output+256*group,coeff+1152*group);
    } else {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupFurther<Randomizer::TowerByte32,FusedPair>(
                scratch+groupStride*group,output+256*group,coeff+144*group);
    }
}

// Separate helpers keep the ordinary-store A/B control unchanged. Every store
// covers a complete cache line; the outer entry checks alignment once.
template<unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storePairSymbolNt(const PlaneValues& a,const PlaneValues& b,Block* out) {
    static_assert(Symbol<8 && (Plane==0 || Plane==2));
    const auto index0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
    const auto index1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x0=_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0);
    const auto x1=_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0);
    _mm512_stream_si512(reinterpret_cast<__m512i*>(out+32*Plane+4*Symbol),
        _mm512_permutex2var_epi8(x0,index0,x1));
    _mm512_stream_si512(reinterpret_cast<__m512i*>(out+32*(Plane+1)+4*Symbol),
        _mm512_permutex2var_epi8(x0,index1,x1));
}

template<unsigned Plane>
static SPIN_NOINLINE void finishPairFusedNt(const __m512i* packed,Block* out) {
    static_assert(Plane==0 || Plane==2);
    const auto a=parityValues<Plane>(packed);
    const auto b=parityValues<Plane+1>(packed);
    storePairSymbolNt<0,Plane>(a,b,out);storePairSymbolNt<1,Plane>(a,b,out);
    storePairSymbolNt<2,Plane>(a,b,out);storePairSymbolNt<3,Plane>(a,b,out);
    storePairSymbolNt<4,Plane>(a,b,out);storePairSymbolNt<5,Plane>(a,b,out);
    storePairSymbolNt<6,Plane>(a,b,out);storePairSymbolNt<7,Plane>(a,b,out);
}

template<Randomizer Mode>
static SPIN_NOINLINE void outerGroupFusedPairNt(const Block* __restrict in,Block* __restrict out,
    const std::uint8_t* __restrict coeff) {
    alignas(64) __m512i low[64],high[64];
    constexpr unsigned stride=Mode==Randomizer::DenseGl32?128:(Mode==Randomizer::TowerByte32?9:72);
    for(unsigned symbol=0;symbol<16;++symbol)
        mixSymbol<Mode>(in+32*symbol,low+4*symbol,high+4*symbol,coeff+stride*symbol);
    finishPairFusedNt<0>(low,out);finishPairFusedNt<2>(low,out);
    finishPairFusedNt<0>(high,out+128);finishPairFusedNt<2>(high,out+128);
}

void outerFusedPairNtAligned(const Block* scratch,Block* output,const Plan& plan) {
    const auto* coeff=reinterpret_cast<const std::uint8_t*>(plan.compactCoefficients());
    if(plan.randomizer==Randomizer::DenseGl32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupFusedPairNt<Randomizer::DenseGl32>(
                scratch+groupStride*group,output+256*group,coeff+2048*group);
    } else if(plan.randomizer==Randomizer::TowerField32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupFusedPairNt<Randomizer::TowerField32>(
                scratch+groupStride*group,output+256*group,coeff+1152*group);
    } else {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroupFusedPairNt<Randomizer::TowerByte32>(
                scratch+groupStride*group,output+256*group,coeff+144*group);
    }
    // Complete the streaming stores before the caller reads or publishes output.
    _mm_sfence();
}
}
void outerFastPair(const Block* scratch,Block* output,const Plan& plan) {
    outerVariant<true,false>(scratch,output,plan);
}
void outerFastNoInline(const Block* scratch,Block* output,const Plan& plan) {
    outerVariant<false,true>(scratch,output,plan);
}
void outerFastPairNoInline(const Block* scratch,Block* output,const Plan& plan) {
    outerVariant<true,true>(scratch,output,plan);
}
void transposeFastPair(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFastPair(scratch,output,plan);
}
void transposeFastNoInline(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFastNoInline(scratch,output,plan);
}
void transposeFastPairNoInline(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFastPairNoInline(scratch,output,plan);
}
void outerFastPlaneNoInline(const Block* scratch,Block* output,const Plan& plan) {
    outerFurther<false>(scratch,output,plan);
}
void outerFastFusedPair(const Block* scratch,Block* output,const Plan& plan) {
    outerFurther<true>(scratch,output,plan);
}
void transposeFastPlaneNoInline(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFastPlaneNoInline(scratch,output,plan);
}
void transposeFastFusedPair(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFastFusedPair(scratch,output,plan);
}
void outerFastFusedPairNt(const Block* scratch,Block* output,const Plan& plan) {
    if((reinterpret_cast<std::uintptr_t>(output)&63u)!=0) {
        outerFastFusedPair(scratch,output,plan);
        return;
    }
    outerFusedPairNtAligned(scratch,output,plan);
}
void transposeFastFusedPairNt(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFastFusedPairNt(scratch,output,plan);
}

}
