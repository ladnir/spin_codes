// Imported exact mode-52 paired-s15 state/feedback circuits.
// Regenerate with src/paired15/import_frozen.py.
// K16Paired15Fold.cpp SHA256: 4a1437d996b4e3b20668b0c05b9c7e780befd195d7fd8b56b3ff3c359590042c
// T64Paired15ShuffleMap.h SHA256: 076d4b38792ad0f0c4b9533440d86109a17124c49c06a8121086d1b2becb379f

#pragma once
#include "Paired15.h"
#include "../packet/PacketInnerFast.h"
namespace spin::detail::paired15::fast {
using namespace detail::packet::fast;
struct Four {__m512i a,b,c,d;};
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
static SPIN_FORCEINLINE void finishWide(const __m512i* high,WordPackets16& feedback) {
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
    // low=(z0,z1,z2,z3); total=(z4,z8,z16,z32)
    feedback.v[0]=low;
    feedback.v[1]=total;
    // odd=(z5,z9,z17,z33), upper=(z6,z10,z18,z34)
    // quadratic=(z12,z20,z24,z36), both=(z40,z40,z48,z48)
    // v2=(z9,0,z33+z12,z17+z36), v3=(z10+z20,z34+z24,z6+z40,z48).
    const auto left2=_mm512_maskz_shuffle_i32x4(0xff0f,odd,odd,0xb1);
    const auto right2=_mm512_shuffle_i32x4(quadratic,quadratic,0xc0);
    feedback.v[2]=_mm512_mask_xor_epi64(left2,0xf0,left2,right2);
    const auto left3=_mm512_shuffle_i32x4(upper,both,0x8d);
    const auto right3=_mm512_shuffle_i32x4(quadratic,upper,0x09);
    feedback.v[3]=_mm512_mask_xor_epi64(left3,0x3f,left3,right3);
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
template<bool Field>
static SPIN_FORCEINLINE void updateWide(const PackedState& before,const std::uint64_t* row,
    const WordPackets16& feedback,WordPackets16& words) {
    __m512i lo0,hi0,lo1,hi1;
    if constexpr(Field) {
        const auto m0=_mm512_set1_epi64(row[0]),m1=_mm512_set1_epi64(row[1]),m2=_mm512_set1_epi64(row[2]);
        const auto p00=_mm512_gf2p8affine_epi64_epi8(m0,before.v[0],0),p01=_mm512_gf2p8affine_epi64_epi8(m0,before.v[1],0);
        const auto p10=_mm512_gf2p8affine_epi64_epi8(m1,before.v[2],0),p11=_mm512_gf2p8affine_epi64_epi8(m1,before.v[3],0);
        const auto p20=_mm512_gf2p8affine_epi64_epi8(m2,_mm512_xor_si512(before.v[0],before.v[2]),0);
        const auto p21=_mm512_gf2p8affine_epi64_epi8(m2,_mm512_xor_si512(before.v[1],before.v[3]),0);
        lo0=_mm512_xor_si512(p00,p10);hi0=_mm512_xor_si512(p01,p11);
        lo1=_mm512_xor_si512(p00,p20);hi1=_mm512_xor_si512(p01,p21);
    } else {
        const auto m00=_mm512_set1_epi64(row[0]),m01=_mm512_set1_epi64(row[1]);
        const auto m10=_mm512_set1_epi64(row[2]),m11=_mm512_set1_epi64(row[3]);
        lo0=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m00,before.v[0],0),_mm512_gf2p8affine_epi64_epi8(m01,before.v[2],0));
        hi0=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m00,before.v[1],0),_mm512_gf2p8affine_epi64_epi8(m01,before.v[3],0));
        lo1=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m10,before.v[0],0),_mm512_gf2p8affine_epi64_epi8(m11,before.v[2],0));
        hi1=_mm512_xor_si512(_mm512_gf2p8affine_epi64_epi8(m10,before.v[1],0),_mm512_gf2p8affine_epi64_epi8(m11,before.v[3],0));
    }
    words.v[0]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo0,unpackIndex<0>(),hi0),feedback.v[0]);
    words.v[1]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo0,unpackIndex<4>(),hi0),feedback.v[1]);
    words.v[2]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo1,unpackIndex<0>(),hi1),feedback.v[2]);
    words.v[3]=_mm512_xor_si512(_mm512_permutex2var_epi8(lo1,unpackIndex<4>(),hi1),feedback.v[3]);
}
template<bool Gather,bool RawFeedback=false,class Emit> static SPIN_FORCEINLINE void streamHigh(const WordPackets16& state,
const Block* raw,const std::uint32_t* route,std::size_t packetBase,Emit& emit,__m512i* high) {
const auto q3=_mm512_shuffle_i32x4(state.v[1],state.v[1],0x00);
const auto c4=q3;
const auto q4=_mm512_shuffle_i32x4(state.v[1],state.v[1],0x55);
const auto c8=q4;
const auto q5=_mm512_shuffle_i32x4(state.v[1],state.v[1],0xaa);
const auto c16=q5;
const auto q6=_mm512_shuffle_i32x4(state.v[1],state.v[1],0xff);
const auto c32=q6;
const auto q8=_mm512_shuffle_i32x4(state.v[2],state.v[2],0x00);
const auto c9=q8;
const auto q9=_mm512_shuffle_i32x4(state.v[3],state.v[3],0xff);
const auto c48=q9;
const auto q11=_mm512_shuffle_i32x4(state.v[2],state.v[2],0xaa);
const auto c33=q11;
const auto c12=q11;
const auto q12=_mm512_shuffle_i32x4(state.v[3],state.v[3],0x00);
const auto c10=q12;
const auto c20=q12;
const auto q13=_mm512_shuffle_i32x4(state.v[2],state.v[2],0xff);
const auto c17=q13;
const auto c36=q13;
const auto q14=_mm512_shuffle_i32x4(state.v[3],state.v[3],0x55);
const auto c34=q14;
const auto c24=q14;
const auto q15=_mm512_shuffle_i32x4(state.v[3],state.v[3],0xaa);
const auto c6=q15;
const auto c40=q15;
const auto d0a=_mm512_mask_xor_epi64(state.v[0],0xcc,state.v[0],_mm512_shuffle_i32x4(state.v[0],state.v[0],0xa0));
const auto d0=_mm512_mask_xor_epi64(d0a,0xf0,d0a,_mm512_shuffle_i32x4(d0a,d0a,0x44));
auto d1=c4;
d1=_mm512_mask_xor_epi64(d1,0xf0,d1,c6);
auto d2=c8;
d2=_mm512_mask_xor_epi64(d2,0xcc,d2,c9);
d2=_mm512_mask_xor_epi64(d2,0xf0,d2,c10);
auto d3=c12;
auto d4=c16;
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
const auto x3=_mm512_loadu_si512(raw+(Gather?route[15]:60));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+15,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+(Gather?route[14]:56));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+14,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+(Gather?route[13]:52));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+13,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+(Gather?route[12]:48));
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
const auto x3=_mm512_loadu_si512(raw+(Gather?route[11]:44));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+11,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+(Gather?route[10]:40));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+10,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+(Gather?route[9]:36));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+9,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+(Gather?route[8]:32));
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
const auto x3=_mm512_loadu_si512(raw+(Gather?route[7]:28));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+7,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+(Gather?route[6]:24));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+6,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+(Gather?route[5]:20));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+5,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+(Gather?route[4]:16));
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
const auto x3=_mm512_loadu_si512(raw+(Gather?route[3]:12));
const auto y3=RawFeedback?x3:_mm512_xor_si512(e3,x3);
emit(packetBase+3,RawFeedback?_mm512_xor_si512(e3,x3):y3);
const auto x2=_mm512_loadu_si512(raw+(Gather?route[2]:8));
const auto y2=RawFeedback?x2:_mm512_xor_si512(e2,x2);
emit(packetBase+2,RawFeedback?_mm512_xor_si512(e2,x2):y2);
const auto x1=_mm512_loadu_si512(raw+(Gather?route[1]:4));
const auto y1=RawFeedback?x1:_mm512_xor_si512(e1,x1);
emit(packetBase+1,RawFeedback?_mm512_xor_si512(e1,x1):y1);
const auto x0=_mm512_loadu_si512(raw+(Gather?route[0]:0));
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

}
