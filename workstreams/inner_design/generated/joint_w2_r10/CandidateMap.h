#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x800,0x92,0x801,0x5d5a0,0x10800,0x14f00,0x10a21,0x49812,0x10082,0x10974,0x343c7,0x69f02,0x2428a,0x204ee,0x3ef,0x590b8,0x50080,0xd0bb,0x2309,0x2601,0x1ecb3,0x4731a,0x4cd1a,0x48780,0x1ebe7,0x43ab8,0x68b2a,0x68f46,0x745dc,0x2db11,0x2731,0x6ccf,0x804,0x5e273,0x2f040,0x2cf04,0x18b05,0x42ee0,0x37161,0x301b7,0x37550,0x69e43,0x3ce50,0x3f070,0xb459,0x510d8,0xd79,0x7ccb,0x16be,0x2460,0x7cd72,0x22a9f,0x4798c,0x404c0,0x3a060,0x6081f,0x6880f,0x6bbb5,0x31087,0x6f60e,0xa535,0xd91d,0x53f9d,0x9686,0x10080,0x1cb5f,0x49a4c,0x184a0,0xce57,0x4a1a,0x556bb,0x7c5,0xa66,0xc0dd,0x7d3ee,0x2cc66,0x386b9,0x30390,0x45d11,0x10d0b,0x44920,0x15a56,0x4f064,0x43621,0x6bc4,0x53720,0xd0a0,0x5977,0xa023,0x5b231,0x25a23,0x29d02,0x6c0cf,0x39d4f,0x438ef,0x4b05c,0x62e9e,0x307a4,0x14c17,0x1b01e,0x76348,0x205e0,0x3e1,0xb07a,0x551ae,0x79f0,0x7063,0x8d0e,0x65e70,0x339bc,0x37d9d,0x3cf62,0x67104,0x68097,0x43005,0x114a5,0x2d0e1,0x26ee0,0x93c0,0x5f8f2,0xedd1,0x1d26,0xef94,0x5ca50,0x60e3c,0x6b159,0x60e59,0x3640f};
static constexpr auto groupedColumns=columns;
static constexpr std::array<unsigned,S> groupOrder{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
static inline void finish(const __m128i* z,__m128i* out) {
const auto v0=z[0];
const auto v1=z[1];
const auto v2=z[2];
const auto v3=z[3];
const auto v4=z[4];
const auto v5=z[5];
const auto v6=z[6];
const auto v7=z[8];
const auto v8=z[9];
const auto v9=z[10];
const auto v10=z[12];
const auto v11=z[16];
const auto v12=z[17];
const auto v13=z[18];
const auto v14=z[20];
const auto v15=z[24];
const auto v16=z[32];
const auto v17=z[33];
const auto v18=z[34];
const auto v19=z[36];
const auto v20=z[40];
const auto v21=z[48];
const auto v22=z[64];
const auto v23=z[65];
const auto v24=z[66];
const auto v25=z[68];
const auto v26=z[72];
const auto v27=z[80];
const auto v28=z[96];
const auto v29=vx(v10,v18);
const auto v30=vx(v13,v19);
const auto v31=vx(v18,v21);
const auto v32=vx(v11,v14);
const auto v33=vx(v9,v29);
const auto v34=vx(v6,v30);
const auto v35=vx(v19,v27);
const auto v36=vx(v26,v31);
const auto v37=vx(v23,v25);
out[0]=v2;
out[1]=vx(v7, v22);
out[2]=v16;
out[3]=vx(vx(v3, v13), vx(v23, v36));
out[4]=vx(vx(v1, v4), v7);
out[5]=vx(vx(vx(v5, v15), vx(v17, v26)), v37);
out[6]=vx(vx(vx(v16, v20), vx(v24, v31)), v34);
out[7]=vx(vx(v0, v4), v22);
out[8]=vx(vx(v20, v28), v33);
out[9]=vx(vx(v8, v29), v35);
out[10]=vx(vx(v1, v9), v36);
out[11]=v0;
out[12]=vx(vx(v12, v25), vx(v32, v35));
out[13]=vx(vx(vx(v0, v5), vx(v10, v19)), vx(v28, v32));
out[14]=vx(vx(v3, v12), v34);
out[15]=vx(vx(vx(v2, v6), vx(v24, v33)), v37);
out[16]=v4;
out[17]=vx(v17, v30);
out[18]=vx(v11, v22);
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x4,0x90,0x40,0x27700,0x1a,0x59e00,0x7eb40,0x89,0x63c00,0x43900,0x3c502,0x1,0x45020,0x3d321,0x5df00,0x68e04,0x8,0x7a100,0xa0};
constexpr u32 inverse[19]={0x800,0x892,0x1,0x10000,0x10882,0x50880,0x4,0x10880,0x522c4,0x6c,0x27592,0x762cc,0x759ec,0x2386c,0x562e4,0x5e0ed,0x8d9b,0x75456,0x72388};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
