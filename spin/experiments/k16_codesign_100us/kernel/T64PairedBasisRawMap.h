#pragma once
// Exact paired expansion; reduce the raw input moments in the same load pass.
// C*A=0 makes the final feedback identical to reducing emitted packets.
// Intermediate moments need not agree, so only the complete selected feedback is used.
#include "T64PairedMap.h"
namespace spin::research::k16codesign::pairedbasisraw {
using namespace detail::packet::fast;
template<class Emit> static SPIN_FORCEINLINE void streamHigh(const WordPackets16& state,
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
const auto q10=_mm512_shuffle_i32x4(state.v[2],state.v[2],0xaa);
const auto c5=q10;
const auto c18=q10;
const auto q11=_mm512_shuffle_i32x4(state.v[2],state.v[2],0xff);
const auto c33=q11;
const auto c12=q11;
const auto q12=_mm512_shuffle_i32x4(state.v[3],state.v[3],0x00);
const auto c10=q12;
const auto c20=q12;
const auto q13=_mm512_shuffle_i32x4(state.v[3],state.v[3],0x55);
const auto c17=q13;
const auto c36=q13;
const auto q14=_mm512_shuffle_i32x4(state.v[3],state.v[3],0xaa);
const auto c34=q14;
const auto c24=q14;
const auto q15=_mm512_shuffle_i32x4(state.v[3],state.v[3],0xff);
const auto c6=q15;
const auto c40=q15;
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
const auto y3=_mm512_loadu_si512(raw+60);
emit(packetBase+15,_mm512_xor_si512(e3,y3));
const auto y2=_mm512_loadu_si512(raw+56);
emit(packetBase+14,_mm512_xor_si512(e2,y2));
const auto y1=_mm512_loadu_si512(raw+52);
emit(packetBase+13,_mm512_xor_si512(e1,y1));
const auto y0=_mm512_loadu_si512(raw+48);
emit(packetBase+12,_mm512_xor_si512(b0,y0));
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
const auto y3=_mm512_loadu_si512(raw+44);
emit(packetBase+11,_mm512_xor_si512(e3,y3));
const auto y2=_mm512_loadu_si512(raw+40);
emit(packetBase+10,_mm512_xor_si512(e2,y2));
const auto y1=_mm512_loadu_si512(raw+36);
emit(packetBase+9,_mm512_xor_si512(e1,y1));
const auto y0=_mm512_loadu_si512(raw+32);
emit(packetBase+8,_mm512_xor_si512(b0,y0));
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
const auto y3=_mm512_loadu_si512(raw+28);
emit(packetBase+7,_mm512_xor_si512(e3,y3));
const auto y2=_mm512_loadu_si512(raw+24);
emit(packetBase+6,_mm512_xor_si512(e2,y2));
const auto y1=_mm512_loadu_si512(raw+20);
emit(packetBase+5,_mm512_xor_si512(e1,y1));
const auto y0=_mm512_loadu_si512(raw+16);
emit(packetBase+4,_mm512_xor_si512(b0,y0));
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
const auto y3=_mm512_loadu_si512(raw+12);
emit(packetBase+3,_mm512_xor_si512(e3,y3));
const auto y2=_mm512_loadu_si512(raw+8);
emit(packetBase+2,_mm512_xor_si512(e2,y2));
const auto y1=_mm512_loadu_si512(raw+4);
emit(packetBase+1,_mm512_xor_si512(e1,y1));
const auto y0=_mm512_loadu_si512(raw+0);
emit(packetBase+0,_mm512_xor_si512(b0,y0));
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
