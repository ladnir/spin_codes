#include "RsWide.h"
#include "Tower32Randomizer.h"
#include "Tower32ByteRandomizer.h"
#include "InnerField16.h"
#include "../../../../../spin/src/packet/PacketLargeInner.h"
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
template<Randomizer Mode>
static SPIN_NOINLINE void outerGroup(const Block* __restrict in,Block* __restrict out,const std::uint8_t* __restrict coeff) {
    alignas(64) __m512i low[64],high[64];
    constexpr unsigned stride=Mode==Randomizer::DenseGl32?128:(Mode==Randomizer::TowerByte32?9:72);
    for(unsigned symbol=0;symbol<16;++symbol)
        mixSymbol<Mode>(in+32*symbol,low+4*symbol,high+4*symbol,coeff+stride*symbol);
    finishHalf(low,out);finishHalf(high,out+128);
}
struct StreamingRoute {
    Block* scratch;const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+bases[packet]),value);
    }
};
}
void reverseRoute(const Block* input,Block* scratch,const Plan& plan) {
    StreamingRoute emit{scratch,plan.route.data()};
    if(plan.innerRandomizer==InnerRandomizer::TowerByte16)
        reverseField16(input,plan.n,plan.fieldUpdates.data(),emit);
    else detail::packet::large::reverse(input,plan.n,plan.denseUpdates.data(),emit);
    _mm_sfence();
}
void outerFast(const Block* scratch,Block* output,const Plan& plan) {
    const auto* coeff=reinterpret_cast<const std::uint8_t*>(plan.compactCoefficients());
    if(plan.randomizer==Randomizer::DenseGl32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroup<Randomizer::DenseGl32>(scratch+groupStride*group,output+256*group,coeff+2048*group);
    } else if(plan.randomizer==Randomizer::TowerField32) {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroup<Randomizer::TowerField32>(scratch+groupStride*group,output+256*group,coeff+1152*group);
    } else {
        for(std::size_t group=0;group<plan.groups;++group)
            outerGroup<Randomizer::TowerByte32>(scratch+groupStride*group,output+256*group,coeff+144*group);
    }
}
void transposeFast(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRoute(input,scratch,plan);outerFast(scratch,output,plan);
}
}
