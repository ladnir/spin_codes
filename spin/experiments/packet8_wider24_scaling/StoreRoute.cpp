#include "StoreVariants.h"
#include "../../src/packet/PacketLargeInner.h"

namespace spin::research::packet8wide24scaling {
namespace {
using packet8wide24::adjointMatrix;
struct Pair { __m512i lo,hi; };

static SPIN_FORCEINLINE Pair vx(Pair a,Pair b) {
    return {_mm512_xor_si512(a.lo,b.lo),_mm512_xor_si512(a.hi,b.hi)};
}
static SPIN_FORCEINLINE Pair vx3(Pair a,Pair b,Pair c) {
    return {_mm512_ternarylogic_epi64(a.lo,b.lo,c.lo,0x96),
            _mm512_ternarylogic_epi64(a.hi,b.hi,c.hi,0x96)};
}
static SPIN_FORCEINLINE Pair affine(Pair x,std::uint64_t matrix) {
    const auto m=_mm512_set1_epi64(matrix);
    return {_mm512_gf2p8affine_epi64_epi8(x.lo,m,0),
            _mm512_gf2p8affine_epi64_epi8(x.hi,m,0)};
}
static SPIN_FORCEINLINE Pair multiply(Pair x,std::uint8_t scalar) {
    const auto m=_mm512_set1_epi8(static_cast<char>(scalar));
    return {_mm512_gf2p8mul_epi8(x.lo,m),_mm512_gf2p8mul_epi8(x.hi,m)};
}
static SPIN_FORCEINLINE Pair pack(const Block* input) {
    using detail::packet::large::packIndex;
    const auto a=_mm512_loadu_si512(input),b=_mm512_loadu_si512(input+4);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    return {_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0),
            _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0)};
}
static SPIN_FORCEINLINE void store(Block* output,Pair x) {
    _mm512_stream_si512(reinterpret_cast<__m512i*>(output),x.lo);
    _mm512_stream_si512(reinterpret_cast<__m512i*>(output+4),x.hi);
}

template<unsigned P>
static SPIN_FORCEINLINE Pair pairStep(const Block* input,Block* scratch,
    const std::uint32_t* route,Pair even,Pair oddDirection,Pair& odd) {
    const auto x=pack(input+16*P),y=pack(input+16*P+8);
    store(scratch+route[2*P],vx(x,even));
    store(scratch+route[2*P+1],vx3(y,even,oddDirection));
    if constexpr(P==0)odd=y;else odd=vx(odd,y);
    return vx(x,y);
}

// Arithmetic and schedule copied from the frozen reverseRouteFast. Only
// the two aligned stores per packet change, followed by one final fence.
static void reverseRouteStreaming(const Block* input,Block* scratch,const Plan& plan) {
    Pair a{_mm512_setzero_si512(),_mm512_setzero_si512()},b=a,c=a;
    for(std::size_t epoch=plan.n/64;epoch-->0;) {
        const auto* in=input+64*epoch;
        const auto* route=plan.route.data()+8*epoch;
        Pair odd,p01,p23,p45,p67;
        {
            const auto one=vx(b,c);
            const auto two=vx(affine(b,adjointMatrix(2)),affine(c,adjointMatrix(4)));
            const auto four=vx(affine(b,adjointMatrix(4)),affine(c,adjointMatrix(16)));
            p01=pairStep<0>(in,scratch,route,a,one,odd);
            p23=pairStep<1>(in,scratch,route,vx(a,two),one,odd);
            p45=pairStep<2>(in,scratch,route,vx(a,four),one,odd);
            p67=pairStep<3>(in,scratch,route,vx3(a,two,four),one,odd);
        }
        if(epoch) {
            const auto middle=vx(p23,p67),high=vx(p45,p67);
            const auto d0=vx3(p01,p23,high);
            const auto d1=vx3(odd,affine(middle,adjointMatrix(2)),affine(high,adjointMatrix(4)));
            const auto d2=vx3(odd,affine(middle,adjointMatrix(4)),affine(high,adjointMatrix(16)));
            const auto& r=plan.updates[epoch].coefficients;
            const auto m0=multiply(a,r[0]),m1=multiply(b,r[1]),m2=multiply(c,r[2]);
            const auto m01=multiply(vx(a,b),r[3]);
            const auto m02=multiply(vx(a,c),r[4]);
            const auto m12=multiply(vx(b,c),r[5]);
            const auto nextA=vx3(vx3(m0,m12,m1),m2,d0);
            const auto nextB=vx3(m01,m0,vx(m12,d1));
            const auto nextC=vx3(m02,m0,vx(m1,d2));
            a=nextA;b=nextB;c=nextC;
        }
    }
    _mm_sfence();
}
}

void reverseRoute(const Block* input,Block* scratch,const Plan& plan,bool stream) {
    if(stream)reverseRouteStreaming(input,scratch,plan);
    else packet8wide24::reverseRouteFast(input,scratch,plan);
}
}
