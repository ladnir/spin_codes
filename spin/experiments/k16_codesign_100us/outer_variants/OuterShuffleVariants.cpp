#include "OuterShuffleVariants.h"
#include "../kernel/NativeOuterFinish.h"
#include "../../../src/packet/PacketLargeInner.h"
#include "../../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"

namespace spin::research::k16codesign::outervariants {
namespace {
static SPIN_FORCEINLINE __m512i packUnaryIndex() {
    alignas(64) static constexpr auto index=[] {
        std::array<std::uint8_t,64> result{};
        for(unsigned byte=0;byte<8;++byte)
            for(unsigned coordinate=0;coordinate<8;++coordinate) {
                const auto c=7-coordinate;
                result[8*byte+coordinate]=std::uint8_t(8*(2*(c%4)+c/4)+byte);
            }
        return result;
    }();
    return _mm512_load_si512(index.data());
}
static SPIN_FORCEINLINE __m512i unpackUnaryIndex() {
    alignas(64) static constexpr auto index=[] {
        std::array<std::uint8_t,64> result{};
        for(unsigned lane=0;lane<8;++lane)
            for(unsigned byte=0;byte<8;++byte)
                result[8*lane+byte]=std::uint8_t(8*byte+lane/2+4*(lane%2));
        return result;
    }();
    return _mm512_load_si512(index.data());
}

template<bool Unary>
static SPIN_FORCEINLINE void mixShuffle(const Block* input,__m512i* out,const std::uint8_t* coeff) {
    using detail::packet::large::packIndex;
    const auto a=_mm512_loadu_si512(input),b=_mm512_loadu_si512(input+4);
    const auto c=_mm512_loadu_si512(input+8),d=_mm512_loadu_si512(input+12);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    __m512i w0,w1,w2,w3;
    if constexpr(Unary) {
        const auto index=packUnaryIndex();
        w0=_mm512_permutexvar_epi8(index,_mm512_unpacklo_epi64(a,b));
        w1=_mm512_permutexvar_epi8(index,_mm512_unpackhi_epi64(a,b));
        w2=_mm512_permutexvar_epi8(index,_mm512_unpacklo_epi64(c,d));
        w3=_mm512_permutexvar_epi8(index,_mm512_unpackhi_epi64(c,d));
    } else {
        w0=_mm512_permutex2var_epi8(a,packIndex<0>(),b);
        w1=_mm512_permutex2var_epi8(a,packIndex<8>(),b);
        w2=_mm512_permutex2var_epi8(c,packIndex<0>(),d);
        w3=_mm512_permutex2var_epi8(c,packIndex<8>(),d);
    }
    const auto v0=_mm512_gf2p8affine_epi64_epi8(basis,w0,0);
    const auto v1=_mm512_gf2p8affine_epi64_epi8(basis,w1,0);
    const auto v2=_mm512_gf2p8affine_epi64_epi8(basis,w2,0);
    const auto v3=_mm512_gf2p8affine_epi64_epi8(basis,w3,0);
    const auto lo=rs::tower32byte::quadraticMultiply(v0,v2,coeff);
    const auto hi=rs::tower32byte::quadraticMultiply(v1,v3,coeff);
    out[0]=lo.lo;out[1]=hi.lo;out[2]=lo.hi;out[3]=hi.hi;
}

template<unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storeUnary(const nativeouter::PlaneValues& a,
    const nativeouter::PlaneValues& b,Block* out) {
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto index=unpackUnaryIndex();
    const auto x0=_mm512_permutexvar_epi8(index,_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0));
    const auto x1=_mm512_permutexvar_epi8(index,_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0));
    _mm512_storeu_si512(out+32*Plane+4*Symbol,_mm512_unpacklo_epi64(x0,x1));
    _mm512_storeu_si512(out+32*(Plane+1)+4*Symbol,_mm512_unpackhi_epi64(x0,x1));
}

template<unsigned Plane>
static SPIN_NOINLINE void finishUnary(const __m512i* packed,Block* out) {
    const auto a=nativeouter::parityValues<Plane>(packed);
    const auto b=nativeouter::parityValues<Plane+1>(packed);
    storeUnary<0,Plane>(a,b,out);storeUnary<1,Plane>(a,b,out);
    storeUnary<2,Plane>(a,b,out);storeUnary<3,Plane>(a,b,out);
    storeUnary<4,Plane>(a,b,out);storeUnary<5,Plane>(a,b,out);
    storeUnary<6,Plane>(a,b,out);storeUnary<7,Plane>(a,b,out);
}

template<bool Pack,bool Unpack>
static SPIN_NOINLINE void shuffleGroup(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,3>* __restrict coeff) {
    alignas(64) __m512i packed[64];
    // Retain the compact fixed-width symbol loop: the earlier fully unrolled
    // alternatives increased register pressure and instruction footprint.
    for(unsigned symbol=0;symbol<16;++symbol)
        mixShuffle<Pack>(input+16*symbol,packed+4*symbol,coeff[symbol].data());
    if constexpr(Unpack) {finishUnary<0>(packed,output);finishUnary<2>(packed,output);}
    else {nativeouter::finishPair<0>(packed,output);nativeouter::finishPair<2>(packed,output);}
}
template<bool Pack,bool Unpack>
static SPIN_FORCEINLINE void shuffleRun(const Block* in,Block* out,const rs::Plan& p,
    const NativeOuterTables& t) {
    for(std::size_t g=0;g<p.groups;++g)
        shuffleGroup<Pack,Unpack>(in+rs::groupStride*g,out+128*g,t.field.data()+16*g);
}
}
void fieldUnaryPack(const Block* in,Block* out,const rs::Plan& p,const NativeOuterTables& t) {
    shuffleRun<true,false>(in,out,p,t);
}
void fieldUnaryUnpack(const Block* in,Block* out,const rs::Plan& p,const NativeOuterTables& t) {
    shuffleRun<false,true>(in,out,p,t);
}
void fieldUnaryBoth(const Block* in,Block* out,const rs::Plan& p,const NativeOuterTables& t) {
    shuffleRun<true,true>(in,out,p,t);
}
}
