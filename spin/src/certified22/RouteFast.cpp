// Exact packed-carry K20Route.cpp circuit; cached/streaming specialization only.
#include "Plan.h"
#include "productionMap.h"
#include "../packet/PacketInnerFast.h"
#include "LaneMajor22.h"
#include "StepCircuit.h"
#include <algorithm>
#include <cstddef>
#include <array>
#include <cstdint>
#include <type_traits>
#include <utility>

namespace spin::detail::certified22 {
namespace {
using namespace step;
using detail::packet::fast::packIndex;
using detail::packet::fast::unpackIndex;
using detail::packet::fast::firstStage;
struct alignas(64) Words22 { __m512i v[6]; };
struct alignas(64) Packed22 { __m512i v[6]; };

template<unsigned I>
static SPIN_FORCEINLINE __m512i broadcast(const Words22& state) {
    static_assert(I<22);
    return _mm512_shuffle_i32x4(state.v[I/4],state.v[I/4],0x55*(I%4));
}
template<unsigned G>
static SPIN_FORCEINLINE void pack(const Words22& words,Packed22& packed) {
    static_assert(G<3);
    const auto lo=_mm512_permutex2var_epi8(words.v[2*G],packIndex<0>(),words.v[2*G+1]);
    const auto hi=_mm512_permutex2var_epi8(words.v[2*G],packIndex<8>(),words.v[2*G+1]);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    packed.v[2*G]=_mm512_gf2p8affine_epi64_epi8(basis,lo,0);
    packed.v[2*G+1]=_mm512_gf2p8affine_epi64_epi8(basis,hi,0);
}
struct Four { __m512i a,b,c,d; };
static SPIN_FORCEINLINE Four transpose4(__m512i a,__m512i b,__m512i c,__m512i d) {
    const auto ab0=_mm512_shuffle_i32x4(a,b,0x44),ab1=_mm512_shuffle_i32x4(a,b,0xee);
    const auto cd0=_mm512_shuffle_i32x4(c,d,0x44),cd1=_mm512_shuffle_i32x4(c,d,0xee);
    return {_mm512_shuffle_i32x4(ab0,cd0,0x88),_mm512_shuffle_i32x4(ab0,cd0,0xdd),
            _mm512_shuffle_i32x4(ab1,cd1,0x88),_mm512_shuffle_i32x4(ab1,cd1,0xdd)};
}
static SPIN_FORCEINLINE __m512i totals4(__m512i a,__m512i b,__m512i c,__m512i d) {
    const auto ab=_mm512_xor_si512(_mm512_shuffle_i32x4(a,b,0x44),_mm512_shuffle_i32x4(a,b,0xee));
    const auto cd=_mm512_xor_si512(_mm512_shuffle_i32x4(c,d,0x44),_mm512_shuffle_i32x4(c,d,0xee));
    return _mm512_xor_si512(_mm512_shuffle_i32x4(ab,cd,0x88),_mm512_shuffle_i32x4(ab,cd,0xdd));
}
static SPIN_FORCEINLINE Words22 finishMonomials(const __m512i* high) {
    const auto l0=_mm512_mask_xor_epi64(high[0],0x33,high[0],_mm512_shuffle_i32x4(high[0],high[0],0xf5));
    const auto low=_mm512_mask_xor_epi64(l0,0x0f,l0,_mm512_shuffle_i32x4(l0,l0,0xee));
    const auto middle=transpose4(high[1],high[2],high[4],high[8]);
    const auto odd=_mm512_xor_si512(middle.b,middle.d);
    const auto upper=_mm512_xor_si512(middle.c,middle.d);
    const auto total=_mm512_ternarylogic_epi64(middle.a,middle.c,odd,0x96);
    const auto quadratic=totals4(high[3],high[5],high[6],high[9]);
    const auto pair=_mm512_xor_si512(_mm512_shuffle_i32x4(high[10],high[12],0x44),
                                    _mm512_shuffle_i32x4(high[10],high[12],0xee));
    const auto both=_mm512_xor_si512(pair,_mm512_shuffle_i32x4(pair,pair,0xb1));
    // (0,1,2,3);(4,8,16,32);(5,9,17,33);(6,10,18,34);
    // (12,20,24,36);(40,48,0,0).
    return {{low,total,odd,upper,quadratic,
        _mm512_maskz_mov_epi64(0x0f,_mm512_shuffle_i32x4(both,both,0x08))}};
}
static SPIN_FORCEINLINE void highFromPackets(const __m512i* packets,__m512i* high) {
    const auto a=firstStage<0>(packets),b=firstStage<4>(packets),c=firstStage<8>(packets),d=firstStage<12>(packets);
    const auto z4=_mm512_xor_si512(b.z0,d.z0),z8=_mm512_xor_si512(c.z0,d.z0);
    high[0]=_mm512_ternarylogic_epi64(a.z0,c.z0,z4,0x96);high[4]=z4;high[8]=z8;high[12]=d.z0;
    const auto z5=_mm512_xor_si512(b.z1,d.z1),z9=_mm512_xor_si512(c.z1,d.z1);
    high[1]=_mm512_ternarylogic_epi64(a.z1,c.z1,z5,0x96);high[5]=z5;high[9]=z9;
    const auto z6=_mm512_xor_si512(b.z2,d.z2),z10=_mm512_xor_si512(c.z2,d.z2);
    high[2]=_mm512_ternarylogic_epi64(a.z2,c.z2,z6,0x96);high[6]=z6;high[10]=z10;
    high[3]=_mm512_ternarylogic_epi64(a.z3,b.z3,_mm512_xor_si512(c.z3,d.z3),0x96);
}
// This order takes fourteen coordinates directly from four existing vectors.
// The remaining coordinate3 occupies an otherwise-zero padding slot.
constexpr unsigned quadraticCoordinate(unsigned i) { return i<14 ? i+8 : 3; }
static_assert([] {
    constexpr unsigned affine[]={0,1,2,4,5,6,7};
    for(auto i:affine) if(map::feedback[i]!=(1u<<i)) return false;
    for(unsigned j=0;j<15;++j)
        for(auto i:affine) if(map::feedback[quadraticCoordinate(j)]&(1u<<i)) return false;
    return true;
}(),"The fast H path requires affine identity and a separate quadratic block");
// Produce H feedback in the same packed output lanes used by the random
// update. The hot loop combines both inputs before the six final unpacks.
struct Feedback22 { __m512i word0,word1; __m512i packed[4]; };
static_assert(map::feedback[3]==0x128008U);
template<unsigned Group,unsigned Input>
static constexpr std::uint64_t feedbackBlock() {
    static_assert(Group==1 || Group==2);
    std::uint64_t result=0;
    for(unsigned r=0;r<8;++r) for(unsigned c=0;c<8;++c) {
        const unsigned output=8*Group+r,input=8*Input+c;
        if(output<22 && input<15)
            result|=std::uint64_t((map::feedback[output]>>quadraticCoordinate(input))&1u)<<(8*r+c);
    }
    return result;
}
template<unsigned G>
static SPIN_FORCEINLINE void feedbackGroup(const Packed22& before,Feedback22& after) {
    const auto m0=_mm512_set1_epi64(feedbackBlock<G,0>());
    const auto m1=_mm512_set1_epi64(feedbackBlock<G,1>());
    after.packed[2*(G-1)]=_mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0),
        _mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0));
    after.packed[2*(G-1)+1]=_mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0),
        _mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0));
}
static SPIN_FORCEINLINE Feedback22 feedbackPacked(const Words22& mono) {
    // H row3 is z3+z34+z20+z40. The other seven affine coordinates
    // pass through unchanged; they need no packed feedback calculation.
    Feedback22 feedback;
    // Only lane 3 contributes here, and mono.v[3] already holds z34 in that lane.
    const auto extra=_mm512_ternarylogic_epi64(
        mono.v[3],broadcast<17>(mono),broadcast<20>(mono),0x96);
    feedback.word0=_mm512_mask_xor_epi64(mono.v[0],0xc0,mono.v[0],extra);
    feedback.word1=mono.v[1];
    // Pack8,...,21,3,0: no cross-vector gathers are needed.
    Words22 quadratic;
    quadratic.v[0]=mono.v[2];quadratic.v[1]=mono.v[3];quadratic.v[2]=mono.v[4];
    quadratic.v[3]=_mm512_mask_mov_epi64(mono.v[5],0x30,broadcast<3>(mono));
    Packed22 before;
    pack<0>(quadratic,before);pack<1>(quadratic,before);
    feedbackGroup<1>(before,feedback);feedbackGroup<2>(before,feedback);
    return feedback;
}
static SPIN_FORCEINLINE Words22 initialFeedback(const Feedback22& feedback) {
    return {{feedback.word0,feedback.word1,
        _mm512_permutex2var_epi8(feedback.packed[0],unpackIndex<0>(),feedback.packed[1]),
        _mm512_permutex2var_epi8(feedback.packed[0],unpackIndex<4>(),feedback.packed[1]),
        _mm512_permutex2var_epi8(feedback.packed[2],unpackIndex<0>(),feedback.packed[3]),
        _mm512_permutex2var_epi8(feedback.packed[2],unpackIndex<4>(),feedback.packed[3])}};
}
template<unsigned G>
static SPIN_FORCEINLINE void updateGroup(const Packed22& before,const std::uint64_t* row,
    const Feedback22& feedback,Words22& words) {
    const auto m0=_mm512_set1_epi64(row[3*G]),m1=_mm512_set1_epi64(row[3*G+1]),m2=_mm512_set1_epi64(row[3*G+2]);
    auto lo=_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0),
        _mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0),
        _mm512_gf2p8affine_epi64_epi8(m2,before.v[4],0),0x96);
    auto hi=_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0),
        _mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0),
        _mm512_gf2p8affine_epi64_epi8(m2,before.v[5],0),0x96);
    if constexpr(G==0) {
        words.v[0]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi),feedback.word0);
        words.v[1]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo,unpackIndex<4>(),hi),feedback.word1);
    } else {
        lo=_mm512_xor_si512(lo,feedback.packed[2*(G-1)]);
        hi=_mm512_xor_si512(hi,feedback.packed[2*(G-1)+1]);
        words.v[2*G]=_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi);
        words.v[2*G+1]=_mm512_permutex2var_epi8(lo,unpackIndex<4>(),hi);
    }
}
// The update's raw lo/hi bytes are payload bytes interleaved by coordinate.
// They are not Packed22: pack(unpack(lo,hi)) reverses each eight-byte row
// and transposes its bits. Keeping them avoids four cross-vector repacks.
struct alignas(64) RawQuadratic22 { __m512i v[4]; };
static SPIN_FORCEINLINE __m512i packInterleaved(__m512i value) {
    const auto reverse=_mm512_setr_epi64(
        0x0001020304050607ULL,0x08090a0b0c0d0e0fULL,
        0x0001020304050607ULL,0x08090a0b0c0d0e0fULL,
        0x0001020304050607ULL,0x08090a0b0c0d0e0fULL,
        0x0001020304050607ULL,0x08090a0b0c0d0e0fULL);
    return _mm512_gf2p8affine_epi64_epi8(
        _mm512_set1_epi64(0x0102040810204080ULL),
        _mm512_shuffle_epi8(value,reverse),0);
}
template<unsigned G>
static SPIN_FORCEINLINE void packCarried(const RawQuadratic22& raw,Packed22& packed) {
    static_assert(G==1 || G==2);
    packed.v[2*G]=packInterleaved(raw.v[2*(G-1)]);
    packed.v[2*G+1]=packInterleaved(raw.v[2*(G-1)+1]);
}
template<unsigned G>
static SPIN_FORCEINLINE void updateGroupCarried(const Packed22& before,const std::uint64_t* row,
    const Feedback22& feedback,Words22& words,RawQuadratic22& raw) {
    static_assert(G==1 || G==2);
    const auto m0=_mm512_set1_epi64(row[3*G]),m1=_mm512_set1_epi64(row[3*G+1]),m2=_mm512_set1_epi64(row[3*G+2]);
    auto lo=_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0),
        _mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0),
        _mm512_gf2p8affine_epi64_epi8(m2,before.v[4],0),0x96);
    auto hi=_mm512_ternarylogic_epi64(
        _mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0),
        _mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0),
        _mm512_gf2p8affine_epi64_epi8(m2,before.v[5],0),0x96);
    lo=_mm512_xor_si512(lo,feedback.packed[2*(G-1)]);
    hi=_mm512_xor_si512(hi,feedback.packed[2*(G-1)+1]);
    raw.v[2*(G-1)]=lo;
    raw.v[2*(G-1)+1]=hi;
    words.v[2*G]=_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi);
    words.v[2*G+1]=_mm512_permutex2var_epi8(lo,unpackIndex<4>(),hi);
}

// Software prefetch distance, in epochs, of the routed packets gathered by
// the forward inner (measured on Zen 4 at K = 2^16): two for one lane, one
// for two lanes, none for the lane-major kernel or the reverse scatter.
#ifndef SPIN_RS8_PF_NARROW
#define SPIN_RS8_PF_NARROW 2
#endif
#ifndef SPIN_RS8_PF_WIDE
#define SPIN_RS8_PF_WIDE 1
#endif
#ifndef SPIN_RS8_PF_LANE_MAJOR
#define SPIN_RS8_PF_LANE_MAJOR 0
#endif
template<bool Write>
static SPIN_FORCEINLINE void prefetchPackets(const void* base,const std::uint32_t* routeBytes) {
    const auto* b=static_cast<const char*>(base);
    for(unsigned i=0; i<16; ++i) SPIN_PREFETCH(b+routeBytes[i],Write?1:0);
}
template<bool Streaming>
struct ByteRoute {
    std::byte* scratch; const std::uint32_t* offsets;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) const {
        auto* destination=reinterpret_cast<__m512i*>(scratch+offsets[packet]);
        if constexpr(Streaming)_mm512_stream_si512(destination,value);
        else _mm512_store_si512(destination,value);
    }
};
template<bool RawFeedback=false,bool Gather=false,class Emit> static SPIN_FORCEINLINE void streamHigh(const Words22& state,
const Block* raw,const std::uint32_t* route,std::size_t packetBase,Emit& emit,__m512i* high) {
// Each of these registers holds one coefficient for each of four affine
// terms. Evaluate at 00, 01, 10, 11 before transposing to one term per register.
const auto row0=state.v[1];
const auto row1=_mm512_xor_si512(row0,state.v[2]);
const auto row2=_mm512_xor_si512(row0,state.v[3]);
const auto row3=_mm512_xor_si512(row1,state.v[3]);
const auto terms=transpose4(row0,row1,row2,row3);
const auto d1=terms.a,d2=terms.b,d4=terms.c,d8=terms.d;
const auto c12=broadcast<16>(state),c20=broadcast<17>(state);
const auto c24=broadcast<18>(state),c36=broadcast<19>(state);
const auto c40=broadcast<20>(state),c48=broadcast<21>(state);
const auto d0a=_mm512_mask_xor_epi64(state.v[0],0xcc,state.v[0],_mm512_shuffle_i32x4(state.v[0],state.v[0],0xa0));
const auto d0=_mm512_mask_xor_epi64(d0a,0xf0,d0a,_mm512_shuffle_i32x4(d0a,d0a,0x44));
auto d3=c12;
auto d5=c20;
auto d6=c24;
auto d9=c36;
auto d10=c40;
auto d12=c48;
__m512i g30,g31,g32,g33;
{
const auto b0=_mm512_xor_si512(_mm512_xor_si512(_mm512_xor_si512(d0,d4),d8),d12);
const auto b1=_mm512_xor_si512(_mm512_xor_si512(d1,d5),d9);
const auto b2=_mm512_xor_si512(_mm512_xor_si512(d2,d6),d10);
const auto e1=_mm512_xor_si512(b0,b1);
const auto e2=_mm512_xor_si512(b0,b2);
const auto e3=_mm512_xor_si512(_mm512_xor_si512(e1,b2),d3);
const auto x3=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[15]):_mm512_loadu_si512(raw+60));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+15,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[14]):_mm512_loadu_si512(raw+56));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+14,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[13]):_mm512_loadu_si512(raw+52));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+13,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[12]):_mm512_loadu_si512(raw+48));
const auto y0=RawFeedback?x0:_mm512_xor_si512(b0,x0);
emit(packetBase+12,RawFeedback?_mm512_xor_si512(b0,x0):y0);
g31=_mm512_xor_si512(y1,y3);
g32=_mm512_xor_si512(y2,y3);
g30=_mm512_ternarylogic_epi64(y0,y2,g31,0x96);
g33=y3;
}
__m512i g20,g21,g22,g23;
{
const auto b0=_mm512_xor_si512(d0,d8);
const auto b1=_mm512_xor_si512(d1,d9);
const auto b2=_mm512_xor_si512(d2,d10);
const auto e1=_mm512_xor_si512(b0,b1);
const auto e2=_mm512_xor_si512(b0,b2);
const auto e3=_mm512_xor_si512(_mm512_xor_si512(e1,b2),d3);
const auto x3=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[11]):_mm512_loadu_si512(raw+44));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+11,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[10]):_mm512_loadu_si512(raw+40));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+10,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[9]):_mm512_loadu_si512(raw+36));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+9,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[8]):_mm512_loadu_si512(raw+32));
const auto y0=RawFeedback?x0:_mm512_xor_si512(b0,x0);
emit(packetBase+8,RawFeedback?_mm512_xor_si512(b0,x0):y0);
g21=_mm512_xor_si512(y1,y3);
g22=_mm512_xor_si512(y2,y3);
g20=_mm512_ternarylogic_epi64(y0,y2,g21,0x96);
g23=y3;
}
const auto cd0=_mm512_xor_si512(g20,g30);
const auto cd1=_mm512_xor_si512(g21,g31);
const auto cd2=_mm512_xor_si512(g22,g32);
const auto cd3=_mm512_xor_si512(g23,g33);
__m512i g10,g11,g12,g13;
{
const auto b0=_mm512_xor_si512(d0,d4);
const auto b1=_mm512_xor_si512(d1,d5);
const auto b2=_mm512_xor_si512(d2,d6);
const auto e1=_mm512_xor_si512(b0,b1);
const auto e2=_mm512_xor_si512(b0,b2);
const auto e3=_mm512_xor_si512(_mm512_xor_si512(e1,b2),d3);
const auto x3=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[7]):_mm512_loadu_si512(raw+28));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+7,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[6]):_mm512_loadu_si512(raw+24));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+6,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[5]):_mm512_loadu_si512(raw+20));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+5,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[4]):_mm512_loadu_si512(raw+16));
const auto y0=RawFeedback?x0:_mm512_xor_si512(b0,x0);
emit(packetBase+4,RawFeedback?_mm512_xor_si512(b0,x0):y0);
g11=_mm512_xor_si512(y1,y3);
g12=_mm512_xor_si512(y2,y3);
g10=_mm512_ternarylogic_epi64(y0,y2,g11,0x96);
g13=y3;
}
__m512i g00,g01,g02,g03;
{
const auto b0=d0;
const auto b1=d1;
const auto b2=d2;
const auto e1=_mm512_xor_si512(b0,b1);
const auto e2=_mm512_xor_si512(b0,b2);
const auto e3=_mm512_xor_si512(_mm512_xor_si512(e1,b2),d3);
const auto x3=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[3]):_mm512_loadu_si512(raw+12));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+3,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[2]):_mm512_loadu_si512(raw+8));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+2,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[1]):_mm512_loadu_si512(raw+4));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+1,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=(Gather?_mm512_loadu_si512(reinterpret_cast<const std::byte*>(raw)+route[0]):_mm512_loadu_si512(raw+0));
const auto y0=RawFeedback?x0:_mm512_xor_si512(b0,x0);
emit(packetBase+0,RawFeedback?_mm512_xor_si512(b0,x0):y0);
g01=_mm512_xor_si512(y1,y3);
g02=_mm512_xor_si512(y2,y3);
g00=_mm512_ternarylogic_epi64(y0,y2,g01,0x96);
g03=y3;
}
high[0]=_mm512_ternarylogic_epi64(g00,g10,cd0,0x96);
high[1]=_mm512_ternarylogic_epi64(g01,g11,cd1,0x96);
high[2]=_mm512_ternarylogic_epi64(g02,g12,cd2,0x96);
high[3]=_mm512_ternarylogic_epi64(g03,g13,cd3,0x96);
high[4]=_mm512_xor_si512(g10,g30);
high[5]=_mm512_xor_si512(g11,g31);
high[6]=_mm512_xor_si512(g12,g32);
high[8]=cd0;
high[9]=cd1;
high[10]=cd2;
high[12]=g30;
}
template<bool Streaming,std::size_t... I>
static SPIN_FORCEINLINE void firstStep(const Block* raw,std::size_t packetBase,
    const ByteRoute<Streaming>& emit,__m512i* packets,std::index_sequence<I...>) {
    ((packets[I]=_mm512_loadu_si512(raw+4*I)),...);
    (emit(packetBase+15-I,packets[15-I]),...);
}
template<bool Streaming>
static SPIN_NOINLINE void reversePackedPeeled(const Block* __restrict input,
    Block* __restrict scratch,const Plan& plan) {
    const auto epochs=plan.n/64;
    if(!epochs)return;
    ByteRoute<Streaming> emit{reinterpret_cast<std::byte*>(scratch),plan.routeBytes.data()};
    alignas(64) __m512i packets[16],high[16];
    const auto first=epochs-1;
    firstStep(input+64*first,16*first,emit,packets,std::make_index_sequence<16>{});
    if(!first)return;
    highFromPackets(packets,high);
    const auto initial=feedbackPacked(finishMonomials(high));
    Words22 state=initialFeedback(initial);
    RawQuadratic22 raw{{initial.packed[0],initial.packed[1],initial.packed[2],initial.packed[3]}};
    for(std::size_t epoch=first;--epoch;) {
        Packed22 before;
        pack<0>(state,before);packCarried<1>(raw,before);packCarried<2>(raw,before);
        streamHigh<true>(state,input+64*epoch,nullptr,16*epoch,emit,high);
        const auto feedback=feedbackPacked(finishMonomials(high));
        updateGroup<0>(before,plan.updates.data()+9*epoch,feedback,state);
        updateGroupCarried<1>(before,plan.updates.data()+9*epoch,feedback,state,raw);
        updateGroupCarried<2>(before,plan.updates.data()+9*epoch,feedback,state,raw);
    }
    streamHigh<true>(state,input,nullptr,0,emit,high);
}

// Forward inner. The circuit basis w = basis*s turns the fixed map into the
// same packed-carry schedule as the reverse pass: outputs y = x + A s,
// feedback basis*A^T*x from the raw packets (the H circuit), and
// w' = basis*R^T*inverseBasis*w + feedback. Epochs run forward, routed
// packets are gathered from scratch, and output is written sequentially.
#ifndef SPIN_RS8_W2_LANE_MAJOR
#define SPIN_RS8_W2_LANE_MAJOR 0
#endif
#ifndef SPIN_RS8_LANE_MAJOR
#define SPIN_RS8_LANE_MAJOR 1
#endif
#ifndef SPIN_RS8_L4_SHUF
#define SPIN_RS8_L4_SHUF 0
#endif
#ifndef SPIN_RS8_WIDE_EPOCHS
#define SPIN_RS8_WIDE_EPOCHS 4
#endif
struct ForwardLane {
    Words22 state;
    RawQuadratic22 raw;
};
template<bool Stream>
static SPIN_FORCEINLINE void storeOutput(Block* output,__m512i value) {
    if constexpr(Stream) _mm512_stream_si512(reinterpret_cast<__m512i*>(output),value);
    else _mm512_storeu_si512(output,value);
}
template<bool Stream> struct DirectOutput {
    Block* output;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) const {
        storeOutput<Stream>(output+4*packet,value);
    }
};
struct EpochOutput {
    __m512i* output;
    SPIN_FORCEINLINE void operator()(std::size_t packet,__m512i value) const {
        output[packet&15]=value;
    }
};
template<std::size_t... I>
static SPIN_FORCEINLINE void gatherPackets(const Block* input,const std::uint32_t* route,
    __m512i* packets,std::index_sequence<I...>) {
    ((packets[I]=_mm512_loadu_si512(reinterpret_cast<const std::byte*>(input)+route[I])),...);
}
template<class Emit>
static SPIN_FORCEINLINE void forwardFirst(const Block* input,const Plan& plan,
    ForwardLane& lane,Emit& emit) {
    alignas(64) __m512i packets[16],high[16];
    gatherPackets(input,plan.routeBytes.data(),packets,std::make_index_sequence<16>{});
    // The initial state is zero, so the first epoch passes through unchanged.
    for(unsigned packet=0; packet<16; ++packet) emit(packet,packets[packet]);
    highFromPackets(packets,high);
    const auto initial=feedbackPacked(finishMonomials(high));
    lane.state=initialFeedback(initial);
    lane.raw={{initial.packed[0],initial.packed[1],initial.packed[2],initial.packed[3]}};
}
template<bool Last,class Emit>
static SPIN_FORCEINLINE void forwardStep(const Block* input,const Plan& plan,std::size_t epoch,
    ForwardLane& lane,Emit& emit) {
    alignas(64) __m512i high[16];
    if constexpr(Last) {
        streamHigh<true,true>(lane.state,input,plan.routeBytes.data()+16*epoch,16*epoch,emit,high);
    } else {
        Packed22 before;
        pack<0>(lane.state,before);packCarried<1>(lane.raw,before);packCarried<2>(lane.raw,before);
        streamHigh<true,true>(lane.state,input,plan.routeBytes.data()+16*epoch,16*epoch,emit,high);
        const auto feedback=feedbackPacked(finishMonomials(high));
        const auto* row=plan.forwardUpdates.data()+9*epoch;
        updateGroup<0>(before,row,feedback,lane.state);
        updateGroupCarried<1>(before,row,feedback,lane.state,lane.raw);
        updateGroupCarried<2>(before,row,feedback,lane.state,lane.raw);
    }
}
template<bool Stream>
static SPIN_NOINLINE void forwardNarrow(const Block* __restrict input,Block* __restrict output,
    const Plan& plan) {
    const auto epochs=plan.n/64;
    DirectOutput<Stream> emit{output};
    ForwardLane lane;
    forwardFirst(input,plan,lane,emit);
    if(epochs==1) return;
    for(std::size_t epoch=1; epoch+1<epochs; ++epoch) {
        if constexpr(SPIN_RS8_PF_NARROW>0) if(epoch+SPIN_RS8_PF_NARROW<epochs)
            prefetchPackets<false>(input,plan.routeBytes.data()+16*(epoch+SPIN_RS8_PF_NARROW));
        forwardStep<false>(input,plan,epoch,lane,emit);
    }
    forwardStep<true>(input,plan,epochs-1,lane,emit);
}
// Interleave staged epochs of L lane planes into L-wide records: record i
// of lane l is output block i*L+l. Every store writes a complete cache line.
template<bool Stream,unsigned L,unsigned E>
static SPIN_FORCEINLINE void storeEpoch(Block* output,const __m512i (&v)[L][16*E],unsigned e) {
    const auto low=_mm512_setr_epi64(0,1,8,9,2,3,10,11);
    const auto high=_mm512_setr_epi64(4,5,12,13,6,7,14,15);
    for(unsigned packet=0; packet<16; ++packet) {
        const auto i=16*e+packet;
#if SPIN_RS8_L4_SHUF
        if constexpr(L==4) {
            const auto t=transpose4(v[0][i],v[1][i],v[2][i],v[3][i]);
            storeOutput<Stream>(output+16*packet,t.a);
            storeOutput<Stream>(output+16*packet+4,t.b);
            storeOutput<Stream>(output+16*packet+8,t.c);
            storeOutput<Stream>(output+16*packet+12,t.d);
            continue;
        }
#endif
        const auto a=_mm512_permutex2var_epi64(v[0][i],low,v[1][i]);
        const auto b=_mm512_permutex2var_epi64(v[0][i],high,v[1][i]);
        if constexpr(L==2) {
            storeOutput<Stream>(output+8*packet,a);
            storeOutput<Stream>(output+8*packet+4,b);
        } else {
            const auto c=_mm512_permutex2var_epi64(v[2][i],low,v[3][i]);
            const auto d=_mm512_permutex2var_epi64(v[2][i],high,v[3][i]);
            storeOutput<Stream>(output+16*packet,_mm512_shuffle_i32x4(a,c,0x44));
            storeOutput<Stream>(output+16*packet+4,_mm512_shuffle_i32x4(a,c,0xee));
            storeOutput<Stream>(output+16*packet+8,_mm512_shuffle_i32x4(b,d,0x44));
            storeOutput<Stream>(output+16*packet+12,_mm512_shuffle_i32x4(b,d,0xee));
        }
    }
}
// Each lane runs up to E consecutive middle epochs with its state in
// registers, staging outputs in L1; the chunk is then interleaved. The first
// and last epochs are peeled so the hot loop contains one step variant.
template<bool Stream,unsigned L,unsigned E>
static SPIN_FORCEINLINE void storeChunk(Block* output,const __m512i (&values)[L][16*E],
    std::size_t e0,unsigned count) {
    for(unsigned e=0; e<count; ++e) storeEpoch<Stream,L,E>(output+64*L*(e0+e),values,e);
}
template<bool Stream,unsigned L,unsigned E>
static SPIN_NOINLINE void forwardWide(const Block* __restrict input,Block* __restrict output,
    const Plan& plan,std::size_t plane) {
    const auto epochs=plan.n/64;
    ForwardLane lanes[L];
    alignas(64) __m512i values[L][16*E];
    for(unsigned lane=0; lane<L; ++lane) {
        EpochOutput emit{values[lane]};
        forwardFirst(input+lane*plane,plan,lanes[lane],emit);
    }
    storeChunk<Stream,L,E>(output,values,0,1);
    if(epochs==1) return;
    for(std::size_t e0=1; e0+1<epochs; e0+=E) {
        const auto count=unsigned(std::min<std::size_t>(E,epochs-1-e0));
        for(unsigned lane=0; lane<L; ++lane) {
            const auto* in=input+lane*plane;
            auto state=lanes[lane];
            for(unsigned e=0; e<count; ++e) {
                if constexpr(SPIN_RS8_PF_WIDE>0) if(e0+e+SPIN_RS8_PF_WIDE<epochs)
                    prefetchPackets<false>(in,plan.routeBytes.data()+16*(e0+e+SPIN_RS8_PF_WIDE));
                EpochOutput emit{values[lane]+16*e};
                forwardStep<false>(in,plan,e0+e,state,emit);
            }
            lanes[lane]=state;
        }
        storeChunk<Stream,L,E>(output,values,e0,count);
    }
    for(unsigned lane=0; lane<L; ++lane) {
        EpochOutput emit{values[lane]};
        forwardStep<true>(input+lane*plane,plan,epochs-1,lanes[lane],emit);
    }
    storeChunk<Stream,L,E>(output,values,epochs-1,1);
}
// Lane-major forward inner for 2 and 4 lanes. One vector holds one step
// coordinate for every lane (256-bit for two lanes, 512-bit for four), so
// the output records are already interleaved and the per-step evaluation and
// feedback are the pure XOR networks in LaneMajor22.h. Only the GL22 update
// converts the 22 state coordinates to the narrow packed layout and back.
struct X512 {
    static SPIN_FORCEINLINE __m512i x2(__m512i a,__m512i b) { return _mm512_xor_si512(a,b); }
    static SPIN_FORCEINLINE __m512i x3(__m512i a,__m512i b,__m512i c) { return _mm512_ternarylogic_epi64(a,b,c,0x96); }
};
struct X256 {
    static SPIN_FORCEINLINE __m256i x2(__m256i a,__m256i b) { return _mm256_xor_si256(a,b); }
    static SPIN_FORCEINLINE __m256i x3(__m256i a,__m256i b,__m256i c) { return _mm256_ternarylogic_epi64(a,b,c,0x96); }
};
template<unsigned L> struct LaneVector;
template<> struct LaneVector<4> { using V=__m512i; using X=X512; };
template<> struct LaneVector<2> { using V=__m256i; using X=X256; };

// The GL22 forward update in the narrow layout, without feedback.
static SPIN_FORCEINLINE Words22 updateNarrow(const Words22& state,const std::uint64_t* row) {
    Packed22 before;
    pack<0>(state,before);pack<1>(state,before);pack<2>(state,before);
    Words22 result;
    const auto group=[&](auto g) {
        constexpr unsigned G=decltype(g)::value;
        const auto m0=_mm512_set1_epi64(row[3*G]),m1=_mm512_set1_epi64(row[3*G+1]),m2=_mm512_set1_epi64(row[3*G+2]);
        const auto lo=_mm512_ternarylogic_epi64(
            _mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0),
            _mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0),
            _mm512_gf2p8affine_epi64_epi8(m2,before.v[4],0),0x96);
        const auto hi=_mm512_ternarylogic_epi64(
            _mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0),
            _mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0),
            _mm512_gf2p8affine_epi64_epi8(m2,before.v[5],0),0x96);
        result.v[2*G]=_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi);
        result.v[2*G+1]=_mm512_permutex2var_epi8(lo,unpackIndex<4>(),hi);
    };
    group(std::integral_constant<unsigned,0>{});
    group(std::integral_constant<unsigned,1>{});
    group(std::integral_constant<unsigned,2>{});
    return result;
}
// Four step coordinates of one routed packet, one vector per coordinate.
template<unsigned L>
static SPIN_FORCEINLINE void gatherPacket(const std::byte* at,std::size_t stride,
    typename LaneVector<L>::V* x) {
    if constexpr(L==4) {
        const auto t=transpose4(_mm512_loadu_si512(at),_mm512_loadu_si512(at+stride),
                                _mm512_loadu_si512(at+2*stride),_mm512_loadu_si512(at+3*stride));
        x[0]=t.a;x[1]=t.b;x[2]=t.c;x[3]=t.d;
    } else {
        const auto q0=_mm512_loadu_si512(at),q1=_mm512_loadu_si512(at+stride);
        const auto t0=_mm512_permutex2var_epi64(q0,_mm512_setr_epi64(0,1,8,9,2,3,10,11),q1);
        const auto t1=_mm512_permutex2var_epi64(q0,_mm512_setr_epi64(4,5,12,13,6,7,14,15),q1);
        x[0]=_mm512_castsi512_si256(t0);x[1]=_mm512_extracti64x4_epi64(t0,1);
        x[2]=_mm512_castsi512_si256(t1);x[3]=_mm512_extracti64x4_epi64(t1,1);
    }
}
template<unsigned L>
static SPIN_FORCEINLINE void toNarrow(const typename LaneVector<L>::V* w,Words22 (&n)[L]) {
    for(unsigned k=0; k<6; ++k) {
        if constexpr(L==4) {
            const auto zero=_mm512_setzero_si512();
            const auto t=transpose4(w[4*k],w[4*k+1],4*k+2<22?w[4*k+2]:zero,4*k+3<22?w[4*k+3]:zero);
            n[0].v[k]=t.a;n[1].v[k]=t.b;n[2].v[k]=t.c;n[3].v[k]=t.d;
        } else {
            const auto zero=_mm256_setzero_si256();
            const auto a=_mm512_inserti64x4(_mm512_castsi256_si512(w[4*k]),w[4*k+1],1);
            const auto b=_mm512_inserti64x4(_mm512_castsi256_si512(4*k+2<22?w[4*k+2]:zero),4*k+3<22?w[4*k+3]:zero,1);
            n[0].v[k]=_mm512_shuffle_i32x4(a,b,0x88);
            n[1].v[k]=_mm512_shuffle_i32x4(a,b,0xdd);
        }
    }
}
template<unsigned L>
static SPIN_FORCEINLINE void fromNarrow(const Words22 (&n)[L],const typename LaneVector<L>::V* f,
    typename LaneVector<L>::V* w) {
    using X=typename LaneVector<L>::X;
    for(unsigned k=0; k<6; ++k) {
        if constexpr(L==4) {
            const auto t=transpose4(n[0].v[k],n[1].v[k],n[2].v[k],n[3].v[k]);
            w[4*k]=X::x2(t.a,f[4*k]);w[4*k+1]=X::x2(t.b,f[4*k+1]);
            if(4*k+2<22) { w[4*k+2]=X::x2(t.c,f[4*k+2]);w[4*k+3]=X::x2(t.d,f[4*k+3]); }
        } else {
            const auto a=_mm512_permutex2var_epi64(n[0].v[k],_mm512_setr_epi64(0,1,8,9,2,3,10,11),n[1].v[k]);
            w[4*k]=X::x2(_mm512_castsi512_si256(a),f[4*k]);
            w[4*k+1]=X::x2(_mm512_extracti64x4_epi64(a,1),f[4*k+1]);
            if(4*k+2<22) {
                const auto b=_mm512_permutex2var_epi64(n[0].v[k],_mm512_setr_epi64(4,5,12,13,6,7,14,15),n[1].v[k]);
                w[4*k+2]=X::x2(_mm512_castsi512_si256(b),f[4*k+2]);
                w[4*k+3]=X::x2(_mm512_extracti64x4_epi64(b,1),f[4*k+3]);
            }
        }
    }
}
template<bool Stream,unsigned L> struct LaneStore {
    Block* output;
    template<class V> SPIN_FORCEINLINE void operator()(unsigned j,V value) const {
        auto* at=output+L*j;
        if constexpr(L==4) storeOutput<Stream>(at,value);
        else if constexpr(Stream) _mm256_stream_si256(reinterpret_cast<__m256i*>(at),value);
        else _mm256_storeu_si256(reinterpret_cast<__m256i*>(at),value);
    }
};
template<bool Stream,unsigned L>
static SPIN_NOINLINE void forwardLaneMajor(const Block* __restrict input,Block* __restrict output,
    const Plan& plan,std::size_t plane) {
    using V=typename LaneVector<L>::V;
    using X=typename LaneVector<L>::X;
    const auto epochs=plan.n/64;
    const auto* route=plan.routeBytes.data();
    const auto* base=reinterpret_cast<const std::byte*>(input);
    const auto stride=plane*sizeof(Block);
    alignas(64) V w[22],f[22],mono[22];
    {
        MonoAccumulator<V,X> acc;
        LaneStore<Stream,L> st{output};
SPIN_UNROLL8
        for(unsigned h=0; h<8; ++h) {
            V x[8];
            gatherPacket<L>(base+route[2*h],stride,x);gatherPacket<L>(base+route[2*h+1],stride,x+4);
            for(unsigned l=0; l<8; ++l) st(8*h+l,x[l]);
            acc.block(h,x);
        }
        acc.finish(mono);
        lanemajor::feedbackH<V,X>(mono,w);
    }
    if(epochs==1) return;
    for(std::size_t epoch=1; epoch+1<epochs; ++epoch) {
        if constexpr(SPIN_RS8_PF_LANE_MAJOR>0) if(epoch+SPIN_RS8_PF_LANE_MAJOR<epochs)
            for(unsigned lane=0; lane<L; ++lane)
                prefetchPackets<false>(base+lane*stride,route+16*(epoch+SPIN_RS8_PF_LANE_MAJOR));
        const auto* r=route+16*epoch;
        StepWeights<V,X> weights;
        weights.prepare(w);
        MonoAccumulator<V,X> acc;
        LaneStore<Stream,L> st{output+64*L*epoch};
SPIN_UNROLL8
        for(unsigned h=0; h<8; ++h) {
            V x[8];
            gatherPacket<L>(base+r[2*h],stride,x);gatherPacket<L>(base+r[2*h+1],stride,x+4);
            evalBlock(weights,h,x,st,8*h);
            acc.block(h,x);
        }
        acc.finish(mono);
        lanemajor::feedbackH<V,X>(mono,f);
        Words22 n[L];
        toNarrow<L>(w,n);
        for(unsigned lane=0; lane<L; ++lane) n[lane]=updateNarrow(n[lane],plan.forwardUpdates.data()+9*epoch);
        fromNarrow<L>(n,f,w);
    }
    const auto* r=route+16*(epochs-1);
    StepWeights<V,X> weights;
    weights.prepare(w);
    LaneStore<Stream,L> st{output+64*L*(epochs-1)};
SPIN_UNROLL8
    for(unsigned h=0; h<8; ++h) {
        V x[8];
        gatherPacket<L>(base+r[2*h],stride,x);gatherPacket<L>(base+r[2*h+1],stride,x+4);
        evalBlock(weights,h,x,st,8*h);
    }
}
template<bool Stream>
static void forwardInner(const Block* scratch,Block* output,const Plan& plan,unsigned lanes,std::size_t plane) {
    if(lanes==1) forwardNarrow<Stream>(scratch,output,plan);
    // Measured selection: packet layout with staged interleave for two
    // lanes, lane-major XOR networks for four lanes.
#if SPIN_RS8_W2_LANE_MAJOR
    else if(lanes==2) forwardLaneMajor<Stream,2>(scratch,output,plan,plane);
#else
    else if(lanes==2) forwardWide<Stream,2,SPIN_RS8_WIDE_EPOCHS>(scratch,output,plan,plane);
#endif
#if SPIN_RS8_LANE_MAJOR
    else forwardLaneMajor<Stream,4>(scratch,output,plan,plane);
#else
    else forwardWide<Stream,4,SPIN_RS8_WIDE_EPOCHS>(scratch,output,plan,plane);
#endif
    if constexpr(Stream) _mm_sfence();
}

}
void forwardInnerFast(const Block* scratch,Block* output,const Plan& plan,unsigned lanes,bool stream,
                      std::size_t plane) {
    if(!plane) plane=plan.scratchBlocks();
    if(stream && !(reinterpret_cast<std::uintptr_t>(output)&63U))
        forwardInner<true>(scratch,output,plan,lanes,plane);
    else
        forwardInner<false>(scratch,output,plan,lanes,plane);
}
void reverseRouteFast(const Block* input,Block* scratch,const Plan& plan) {
    if(plan.k<=262144) {
        reversePackedPeeled<false>(input,scratch,plan);
    } else {
        reversePackedPeeled<true>(input,scratch,plan);
        // Routing must complete before the outer reads the scratch buffer.
        _mm_sfence();
    }
}
}
