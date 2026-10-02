#pragma once
// Paired t64/s15: remove old literal row 10 and compact remaining coordinates.
// SIMD state additionally swaps compact coordinates 3 and 7, as PairedBasis.
#include "T64PairedMap.h"
namespace spin::research::k16codesign::paired15 {
using namespace detail::packet::fast;
inline constexpr std::uint16_t columns[64]={0x1,0x3,0x5,0x87,0x9,0xb,0x400d,0x408f,0x11,0x113,0x815,0x997,0x419,0x51b,0x4c1d,0x4d9f,0x21,0x1023,0x25,0x10a7,0x829,0x182b,0x482d,0x58af,0x2031,0x3133,0x2835,0x39b7,0x2c39,0x3d3b,0x643d,0x75bf,0x41,0x443,0x2045,0x24c7,0x1049,0x144b,0x704d,0x74cf,0x4051,0x4553,0x6855,0x6dd7,0x5459,0x515b,0x3c5d,0x39df,0x261,0x1663,0x2265,0x36e7,0x1a69,0xe6b,0x7a6d,0x6eef,0x6271,0x7773,0x4a75,0x5ff7,0x7e79,0x6b7b,0x167d,0x3ff};
template<bool RawFeedback=false,class Emit> static SPIN_FORCEINLINE void streamHigh(const WordPackets16& state,
const Block* raw,std::size_t packetBase,Emit& emit,__m512i* high) {
const auto q3=_mm512_shuffle_i32x4(state.v[1],state.v[1],0xff);
const auto c4=q3;
const auto q4=_mm512_shuffle_i32x4(state.v[1],state.v[1],0x00);
const auto c8=q4;
const auto q5=_mm512_shuffle_i32x4(state.v[1],state.v[1],0x55);
const auto c16=q5;
const auto q6=_mm512_shuffle_i32x4(state.v[1],state.v[1],0xaa);
const auto c32=q6;
const auto q8=_mm512_shuffle_i32x4(state.v[2],state.v[2],0x00);
const auto c9=q8;
const auto q9=_mm512_shuffle_i32x4(state.v[2],state.v[2],0x55);
const auto c48=q9;
const auto q11=_mm512_shuffle_i32x4(state.v[2],state.v[2],0xaa);
const auto c33=q11;
const auto c12=q11;
const auto q12=_mm512_shuffle_i32x4(state.v[2],state.v[2],0xff);
const auto c10=q12;
const auto c20=q12;
const auto q13=_mm512_shuffle_i32x4(state.v[3],state.v[3],0x00);
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
}
