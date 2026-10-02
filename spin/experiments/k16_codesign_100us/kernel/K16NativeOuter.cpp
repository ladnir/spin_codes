#include "K16NativeOuter.h"
#include "NativeOuterFinish.h"
#include "../../../src/packet/PacketLargeInner.h"
#include "../../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"

namespace spin::research::k16codesign {
namespace {
static SPIN_FORCEINLINE void packNative(const Block* input,__m512i* v) {
    using detail::packet::large::packIndex;
    const auto a=_mm512_loadu_si512(input);
    const auto b=_mm512_loadu_si512(input+4);
    const auto c=_mm512_loadu_si512(input+8);
    const auto d=_mm512_loadu_si512(input+12);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    v[0]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0);
    v[1]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0);
    v[2]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(c,packIndex<0>(),d),0);
    v[3]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(c,packIndex<8>(),d),0);
}
static SPIN_FORCEINLINE void mixGl(const Block* input,__m512i* packed,const std::uint64_t* coeff) {
    __m512i v[4];packNative(input,v);
    const auto m0=_mm512_set1_epi64(coeff[0]),m1=_mm512_set1_epi64(coeff[1]);
    const auto m2=_mm512_set1_epi64(coeff[2]),m3=_mm512_set1_epi64(coeff[3]);
    packed[0]=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[0],m0,0),_mm512_gf2p8affine_epi64_epi8(v[2],m1,0));
    packed[1]=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[1],m0,0),_mm512_gf2p8affine_epi64_epi8(v[3],m1,0));
    packed[2]=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[0],m2,0),_mm512_gf2p8affine_epi64_epi8(v[2],m3,0));
    packed[3]=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(v[1],m2,0),_mm512_gf2p8affine_epi64_epi8(v[3],m3,0));
}
static SPIN_FORCEINLINE void mixField(const Block* input,__m512i* packed,const std::uint8_t* coeff) {
    __m512i v[4];packNative(input,v);
    const auto lo=rs::tower32byte::quadraticMultiply(v[0],v[2],coeff);
    const auto hi=rs::tower32byte::quadraticMultiply(v[1],v[3],coeff);
    packed[0]=lo.lo;packed[1]=hi.lo;packed[2]=lo.hi;packed[3]=hi.hi;
}
static SPIN_NOINLINE void nativeGlGroup(const Block* __restrict in,Block* __restrict out,
    const std::uint64_t* __restrict coeff) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol)mixGl(in+16*symbol,packed+4*symbol,coeff+4*symbol);
    nativeouter::finishPair<0>(packed,out);nativeouter::finishPair<2>(packed,out);
}
static SPIN_NOINLINE void nativeFieldGroup(const Block* __restrict in,Block* __restrict out,
    const std::array<std::uint8_t,3>* __restrict coeff) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol)mixField(in+16*symbol,packed+4*symbol,coeff[symbol].data());
    nativeouter::finishPair<0>(packed,out);nativeouter::finishPair<2>(packed,out);
}
}
void outerNativeGl16(const Block* scratch,Block* output,const rs::Plan& plan,const NativeOuterTables& tables) {
    for(std::size_t group=0;group<plan.groups;++group)
        nativeGlGroup(scratch+rs::groupStride*group,output+128*group,tables.gl.data()+64*group);
}
void outerNativeField16(const Block* scratch,Block* output,const rs::Plan& plan,const NativeOuterTables& tables) {
    for(std::size_t group=0;group<plan.groups;++group)
        nativeFieldGroup(scratch+rs::groupStride*group,output+128*group,tables.field.data()+16*group);
}
}
