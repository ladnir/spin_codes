#include "K16CodeSign.h"
#include "../../../src/packet/PacketLargeInner.h"
#include "../../../../research/workstreams/k16_design/implementation/Rs16OuterFast.h"
#include "../../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"

namespace spin::research::k16codesign {
namespace {
using namespace detail::packet::large;
struct CachedRoute {
    Block* scratch;
    const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+bases[packet]),value);
    }
};

template<class Emit>
static SPIN_NOINLINE void reverseField16(const Block* input,std::size_t n,
    const std::array<std::uint8_t,3>* rows,Emit& emit) {
    alignas(64) __m128i words[16],moments[64],syndrome[16];
    alignas(64) __m512i packets[16];
    PackedState state{},feedback;
    for(std::size_t epoch=n/64;epoch-->0;) {
        const auto* raw=input+64*epoch;
        const bool first=epoch+1==n/64;
        if(first) {
            loadPackets(raw,packets,std::make_index_sequence<16>{});
            for(unsigned h=16;h-->0;)emit(16*epoch+h,packets[h]);
        } else {
            wideUnpackVbmi(state,words);
            streamStep(words,raw,16*epoch,emit,moments);
        }
        if(!epoch)break;
        if(first)packetMoments(packets,moments);
        finish(moments,syndrome);
        widePackVbmi(syndrome,feedback);
        if(first)state=feedback;
        else {
            const auto lo=rs::tower32byte::quadraticMultiply(state.v[0],state.v[2],rows[epoch].data());
            const auto hi=rs::tower32byte::quadraticMultiply(state.v[1],state.v[3],rows[epoch].data());
            state.v[0]=_mm512_xor_si512(lo.lo,feedback.v[0]);
            state.v[1]=_mm512_xor_si512(hi.lo,feedback.v[1]);
            state.v[2]=_mm512_xor_si512(lo.hi,feedback.v[2]);
            state.v[3]=_mm512_xor_si512(hi.hi,feedback.v[3]);
        }
    }
}

// Same fixed packing network as retained Rs16OuterFast::mixSymbol. Keeping
// the 16 symbol bits in byteplanes makes both GL16 and field maps lane-local.
static SPIN_FORCEINLINE void packSymbol(const Block* input,__m512i* v) {
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
static SPIN_FORCEINLINE void mixField(const Block* input,__m512i* packed,const std::uint8_t* coeff) {
    __m512i v[4];packSymbol(input,v);
    const auto lo=rs::tower32byte::quadraticMultiply(v[0],v[2],coeff);
    const auto hi=rs::tower32byte::quadraticMultiply(v[1],v[3],coeff);
    packed[0]=lo.lo;packed[1]=hi.lo;packed[2]=lo.hi;packed[3]=hi.hi;
}

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
static SPIN_NOINLINE void groupGl16(const Block* __restrict in,Block* __restrict out,
    const std::uint64_t* __restrict coeff) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol)
        rs::fast16::mixSymbol(in+16*symbol,packed+4*symbol,coeff+4*symbol);
    finishPair<0>(packed,out);finishPair<2>(packed,out);
}
static SPIN_NOINLINE void groupField16(const Block* __restrict in,Block* __restrict out,
    const std::array<std::uint8_t,3>* __restrict coeff) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol)
        mixField(in+16*symbol,packed+4*symbol,coeff[symbol].data());
    finishPair<0>(packed,out);finishPair<2>(packed,out);
}
}

void reverseRouteCachedGl16(const Block* input,Block* scratch,const rs::Plan& plan) {
    CachedRoute emit{scratch,plan.route.data()};
    detail::packet::large::reverse(input,plan.n,plan.denseUpdates.data(),emit);
}
void reverseRouteCachedField16(const Block* input,Block* scratch,const rs::Plan& plan,const Tables& tables) {
    CachedRoute emit{scratch,plan.route.data()};
    reverseField16(input,plan.n,tables.innerField.data(),emit);
}
void outerFusedGl16(const Block* scratch,Block* output,const rs::Plan& plan) {
    const auto* coefficients=plan.compactCoefficients();
    for(std::size_t group=0;group<plan.groups;++group)
        groupGl16(scratch+rs::groupStride*group,output+128*group,coefficients+64*group);
}
void outerFusedField16(const Block* scratch,Block* output,const rs::Plan& plan,const Tables& tables) {
    for(std::size_t group=0;group<plan.groups;++group)
        groupField16(scratch+rs::groupStride*group,output+128*group,tables.outerField.data()+16*group);
}
}
