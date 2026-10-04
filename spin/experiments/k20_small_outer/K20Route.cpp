// Exact certified t64/s22 circuit copied from ../k18_certified22/K18Fast.cpp.
// Routing stores and optional state repacking are specialized here; A, C,
// H=P P^T, and the fixed emission/feedback schedule remain unchanged.
// reverseRouteNt retains the store-only control. The frozen source is retained.
#include "K20Route.h"
#include "K18Map.h"
#include "../../src/packet/PacketInnerFast.h"
#include <array>
#include <cstdint>
#include <utility>

namespace spin::research::k20small {
using Block = rs::Block;
using k18cert22::Tables;
namespace map = k18cert22::map;
namespace {
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
    const auto extra=_mm512_ternarylogic_epi64(
        broadcast<15>(mono),broadcast<17>(mono),broadcast<20>(mono),0x96);
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

static bool checkCarryPair(__m512i lo,__m512i hi) noexcept {
    Words22 words{};
    words.v[0]=_mm512_permutex2var_epi8(lo,unpackIndex<0>(),hi);
    words.v[1]=_mm512_permutex2var_epi8(lo,unpackIndex<4>(),hi);
    Packed22 reference;
    pack<0>(words,reference);
    return _mm512_cmpeq_epi64_mask(reference.v[0],packInterleaved(lo))==0xff
        && _mm512_cmpeq_epi64_mask(reference.v[1],packInterleaved(hi))==0xff;
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
template<bool RawFeedback=false,class Emit> static SPIN_FORCEINLINE void streamHigh(const Words22& state,
const Block* raw,std::size_t packetBase,Emit& emit,__m512i* high) {
const auto c4=broadcast<4>(state),c8=broadcast<5>(state);
const auto c16=broadcast<6>(state),c32=broadcast<7>(state);
const auto c5=broadcast<8>(state),c9=broadcast<9>(state);
const auto c17=broadcast<10>(state),c33=broadcast<11>(state);
const auto c6=broadcast<12>(state),c10=broadcast<13>(state);
const auto c18=broadcast<14>(state),c34=broadcast<15>(state);
const auto c12=broadcast<16>(state),c20=broadcast<17>(state);
const auto c24=broadcast<18>(state),c36=broadcast<19>(state);
const auto c40=broadcast<20>(state),c48=broadcast<21>(state);
const auto d0a=_mm512_mask_xor_epi64(state.v[0],0xcc,state.v[0],_mm512_shuffle_i32x4(state.v[0],state.v[0],0xa0));
const auto d0=_mm512_mask_xor_epi64(d0a,0xf0,d0a,_mm512_shuffle_i32x4(d0a,d0a,0x44));
auto d1=c4;
d1=_mm512_mask_xor_epi64(d1,0xcc,d1,c5);
d1=_mm512_mask_xor_epi64(d1,0xf0,d1,c6);
auto d2=c8;
d2=_mm512_mask_xor_epi64(d2,0xcc,d2,c9);
d2=_mm512_mask_xor_epi64(d2,0xf0,d2,c10);
auto d3=c12;
auto d4=c16;
d4=_mm512_mask_xor_epi64(d4,0xf0,d4,c18);
d4=_mm512_mask_xor_epi64(d4,0xcc,d4,c17);
auto d5=c20;
auto d6=c24;
auto d8=c32;
d8=_mm512_mask_xor_epi64(d8,0xcc,d8,c33);
d8=_mm512_mask_xor_epi64(d8,0xf0,d8,c34);
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
const auto x3=_mm512_loadu_si512(raw+60);
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+15,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+56);
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+14,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+52);
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+13,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+48);
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
const auto x3=_mm512_loadu_si512(raw+44);
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+11,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+40);
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+10,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+36);
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+9,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+32);
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
const auto x3=_mm512_loadu_si512(raw+28);
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+7,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+24);
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+6,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+20);
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+5,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+16);
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
const auto x3=_mm512_loadu_si512(raw+12);
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+3,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+8);
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+2,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+4);
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+1,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+0);
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
static SPIN_NOINLINE void reversePeeled(const Block* __restrict input,
    Block* __restrict scratch,const rs::Plan& plan,const Tables& tables) {
    const auto epochs=plan.n/64;
    if(!epochs)return;
    ByteRoute<Streaming> emit{reinterpret_cast<std::byte*>(scratch),tables.routeBytes.data()};
    alignas(64) __m512i packets[16],high[16];
    const auto first=epochs-1;
    firstStep(input+64*first,16*first,emit,packets,std::make_index_sequence<16>{});
    if(!first)return;
    highFromPackets(packets,high);
    Words22 state=initialFeedback(feedbackPacked(finishMonomials(high)));
    for(std::size_t epoch=first;--epoch;) {
        Packed22 before;
        pack<0>(state,before);pack<1>(state,before);pack<2>(state,before);
        streamHigh<true>(state,input+64*epoch,16*epoch,emit,high);
        const auto feedback=feedbackPacked(finishMonomials(high));
        updateGroup<0>(before,tables.updates.data()+9*epoch,feedback,state);
        updateGroup<1>(before,tables.updates.data()+9*epoch,feedback,state);
        updateGroup<2>(before,tables.updates.data()+9*epoch,feedback,state);
    }
    streamHigh<true>(state,input,0,emit,high);
}
}
template<bool Streaming>
static SPIN_NOINLINE void reversePackedPeeled(const Block* __restrict input,
    Block* __restrict scratch,const rs::Plan& plan,const Tables& tables) {
    const auto epochs=plan.n/64;
    if(!epochs)return;
    ByteRoute<Streaming> emit{reinterpret_cast<std::byte*>(scratch),tables.routeBytes.data()};
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
        streamHigh<true>(state,input+64*epoch,16*epoch,emit,high);
        const auto feedback=feedbackPacked(finishMonomials(high));
        updateGroup<0>(before,tables.updates.data()+9*epoch,feedback,state);
        updateGroupCarried<1>(before,tables.updates.data()+9*epoch,feedback,state,raw);
        updateGroupCarried<2>(before,tables.updates.data()+9*epoch,feedback,state,raw);
    }
    streamHigh<true>(state,input,0,emit,high);
}

void reverseRouteNt(const Block* input,Block* scratch,const rs::Plan& plan,const Tables& tables) {
    reversePeeled<true>(input,scratch,plan,tables);
    // Complete all non-temporal routing stores before the outer reads scratch.
    _mm_sfence();
}
void reverseRoutePackedNt(const Block* input,Block* scratch,const rs::Plan& plan,const Tables& tables) {
    reversePackedPeeled<true>(input,scratch,plan,tables);
    _mm_sfence();
}
void reverseRoutePackedCached(const Block* input,Block* scratch,const rs::Plan& plan,const Tables& tables) {
    reversePackedPeeled<false>(input,scratch,plan,tables);
}
bool packedCarrySelfTest() noexcept {
    alignas(64) std::array<std::uint64_t,16> raw{};
    for(unsigned bit=0;bit<1024;++bit) {
        raw.fill(0);
        raw[bit/64]=std::uint64_t{1}<<(bit%64);
        if(!checkCarryPair(_mm512_load_si512(raw.data()),_mm512_load_si512(raw.data()+8)))return false;
    }
    std::uint64_t seed=0x528be7367bcb24d5ULL;
    for(unsigned trial=0;trial<256;++trial) {
        for(auto& value:raw) {
            seed+=0x9e3779b97f4a7c15ULL;
            auto word=seed;
            word=(word^(word>>30))*0xbf58476d1ce4e5b9ULL;
            word=(word^(word>>27))*0x94d049bb133111ebULL;
            value=word^(word>>31);
        }
        if(!checkCarryPair(_mm512_load_si512(raw.data()),_mm512_load_si512(raw.data()+8)))return false;
        // Group 2 has six true coordinates and two zero SIMD padding slots.
        for(auto& value:raw)value&=0x0000ffffffffffffULL;
        const auto lo=_mm512_load_si512(raw.data()),hi=_mm512_load_si512(raw.data()+8);
        if(!checkCarryPair(lo,hi))return false;
        const auto padding=_mm512_set1_epi64(0xc0c0c0c0c0c0c0c0ULL);
        if(_mm512_test_epi64_mask(packInterleaved(lo),padding)
            || _mm512_test_epi64_mask(packInterleaved(hi),padding))return false;
    }
    return true;
}
}
