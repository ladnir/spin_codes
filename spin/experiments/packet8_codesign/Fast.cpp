#include "Packet8.h"
#include "../../src/packet/PacketLargeInner.h"
#include "../k16_codesign_100us/kernel/NativeOuterFinish.h"
#include "../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"

namespace spin::research::packet8 {
namespace {
struct Pair {__m512i lo,hi;};
static SPIN_FORCEINLINE Pair vx(Pair a,Pair b) {
    return {_mm512_xor_si512(a.lo,b.lo),_mm512_xor_si512(a.hi,b.hi)};
}
static SPIN_FORCEINLINE Pair vx3(Pair a,Pair b,Pair c) {
    return {_mm512_ternarylogic_epi64(a.lo,b.lo,c.lo,0x96),_mm512_ternarylogic_epi64(a.hi,b.hi,c.hi,0x96)};
}
static SPIN_FORCEINLINE Pair affine(Pair x,std::uint64_t matrix) {
    const auto m=_mm512_set1_epi64(matrix);
    return {_mm512_gf2p8affine_epi64_epi8(x.lo,m,0),_mm512_gf2p8affine_epi64_epi8(x.hi,m,0)};
}
static SPIN_FORCEINLINE Pair pack(const Block* input) {
    using detail::packet::large::packIndex;
    const auto a=_mm512_loadu_si512(input),b=_mm512_loadu_si512(input+4);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    return {_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0),
            _mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0)};
}
static SPIN_FORCEINLINE void store(Block* out,Pair v) {
    _mm512_store_si512(out,v.lo);_mm512_store_si512(out+4,v.hi);
}
template<unsigned P>
static SPIN_FORCEINLINE Pair pairStep(const Block* in,Block* out,const std::uint32_t* route,
    Pair even,Pair stateB,Pair& odd) {
    const auto x=pack(in+16*P),y=pack(in+16*P+8);
    store(out+route[2*P],vx(x,even));
    store(out+route[2*P+1],vx3(y,even,stateB));
    if constexpr(P==0)odd=y;else odd=vx(odd,y);
    return vx(x,y);
}
static SPIN_NOINLINE void outerGroup(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,3>* __restrict coefficients) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol) {
        const auto a=_mm512_load_si512(input+16*symbol);
        const auto b=_mm512_load_si512(input+16*symbol+4);
        const auto c=_mm512_load_si512(input+16*symbol+8);
        const auto d=_mm512_load_si512(input+16*symbol+12);
        const auto lo=rs::tower32byte::quadraticMultiply(a,c,coefficients[symbol].data());
        const auto hi=rs::tower32byte::quadraticMultiply(b,d,coefficients[symbol].data());
        packed[4*symbol]=lo.lo;packed[4*symbol+1]=hi.lo;
        packed[4*symbol+2]=lo.hi;packed[4*symbol+3]=hi.hi;
    }
    k16codesign::nativeouter::finishPair<0>(packed,output);
    k16codesign::nativeouter::finishPair<2>(packed,output);
}
}
void reverseRouteFast(const Block* input,Block* scratch,const Plan& plan) {
    Pair a{_mm512_setzero_si512(),_mm512_setzero_si512()},b=a;
    for(std::size_t epoch=plan.n()/64;epoch-->0;) {
        const auto* in=input+64*epoch;
        const auto* route=plan.route.data()+8*epoch;
        const auto two=affine(b,adjointMatrix(2)),four=affine(b,adjointMatrix(4));
        Pair odd;
        const auto p01=pairStep<0>(in,scratch,route,a,b,odd);
        const auto p23=pairStep<1>(in,scratch,route,vx(a,two),b,odd);
        const auto p45=pairStep<2>(in,scratch,route,vx(a,four),b,odd);
        const auto p67=pairStep<3>(in,scratch,route,vx3(a,two,four),b,odd);
        if(epoch) {
            const auto high=vx(p45,p67);
            const auto d0=vx3(p01,p23,high);
            const auto d1=vx3(odd,affine(vx(p23,p67),adjointMatrix(2)),affine(high,adjointMatrix(4)));
            const auto& m=plan.updates[epoch].adjoint;
            const auto nextA=vx3(affine(a,m[0]),affine(b,m[2]),d0);
            b=vx3(affine(a,m[1]),affine(b,m[3]),d1);
            a=nextA;
        }
    }
}
void routeOnlyFast(const Block* input,Block* scratch,const Plan& plan) {
    for(std::size_t packet=plan.n()/8;packet-->0;)
        store(scratch+plan.route[packet],pack(input+8*packet));
}
namespace {
static SPIN_FORCEINLINE void prefetchPair(Block* destination) {
#ifdef _MSC_VER
    _m_prefetchw(destination);_m_prefetchw(destination+4);
#else
    __builtin_prefetch(destination,1,3);__builtin_prefetch(destination+4,1,3);
#endif
}
template<std::size_t... I>
static SPIN_FORCEINLINE void prefetchEpoch(Block* scratch,const std::uint32_t* route,
    std::index_sequence<I...>) {
    (prefetchPair(scratch+route[I]),...);
}
// Keep the retained baseline above byte-for-byte in place. These separately
// instantiated variants add only a compile-time lookahead to the same steps.
static SPIN_FORCEINLINE void prefetchStep(const Block* input,Block* scratch,const Plan& plan,
    std::size_t epoch,Pair& a,Pair& b) {
        const auto* in=input+64*epoch;
        const auto* route=plan.route.data()+8*epoch;
        const auto two=affine(b,adjointMatrix(2)),four=affine(b,adjointMatrix(4));
        Pair odd;
        const auto p01=pairStep<0>(in,scratch,route,a,b,odd);
        const auto p23=pairStep<1>(in,scratch,route,vx(a,two),b,odd);
        const auto p45=pairStep<2>(in,scratch,route,vx(a,four),b,odd);
        const auto p67=pairStep<3>(in,scratch,route,vx3(a,two,four),b,odd);
        if(epoch) {
            const auto high=vx(p45,p67);
            const auto d0=vx3(p01,p23,high);
            const auto d1=vx3(odd,affine(vx(p23,p67),adjointMatrix(2)),affine(high,adjointMatrix(4)));
            const auto& m=plan.updates[epoch].adjoint;
            const auto nextA=vx3(affine(a,m[0]),affine(b,m[2]),d0);
            b=vx3(affine(a,m[1]),affine(b,m[3]),d1);
            a=nextA;
        }
}
template<unsigned Packets>
void reversePrefetch(const Block* input,Block* scratch,const Plan& plan) {
    static_assert(Packets==8||Packets==16||Packets==32);
    Pair a{_mm512_setzero_si512(),_mm512_setzero_si512()},b=a;
    auto epoch=plan.n()/64;
    for(;epoch>Packets/8;) {
        --epoch;
        prefetchEpoch(scratch,plan.route.data()+8*epoch-Packets,std::make_index_sequence<8>{});
        prefetchStep(input,scratch,plan,epoch,a,b);
    }
    for(;epoch-->0;)prefetchStep(input,scratch,plan,epoch,a,b);
}
template<unsigned Packets>
void routePrefetch(const Block* input,Block* scratch,const Plan& plan) {
    static_assert(Packets==8||Packets==16||Packets==32);
    auto packet=plan.n()/8;
    for(;packet>Packets;) {
        --packet;prefetchPair(scratch+plan.route[packet-Packets]);
        store(scratch+plan.route[packet],pack(input+8*packet));
    }
    for(;packet-->0;)store(scratch+plan.route[packet],pack(input+8*packet));
}
}
void reverseRoutePrefetch8(const Block* in,Block* scratch,const Plan& p){reversePrefetch<8>(in,scratch,p);}
void reverseRoutePrefetch16(const Block* in,Block* scratch,const Plan& p){reversePrefetch<16>(in,scratch,p);}
void reverseRoutePrefetch32(const Block* in,Block* scratch,const Plan& p){reversePrefetch<32>(in,scratch,p);}
void routeOnlyPrefetch8(const Block* in,Block* scratch,const Plan& p){routePrefetch<8>(in,scratch,p);}
void routeOnlyPrefetch16(const Block* in,Block* scratch,const Plan& p){routePrefetch<16>(in,scratch,p);}
void routeOnlyPrefetch32(const Block* in,Block* scratch,const Plan& p){routePrefetch<32>(in,scratch,p);}
void outerFast(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t g=0;g<plan.outer.groups;++g)
        outerGroup(scratch+rs::groupStride*g,output+128*g,plan.native.field.data()+16*g);
}
void transposeFast(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRouteFast(input,scratch,plan);outerFast(scratch,output,plan);
}
}
