#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x100,0x108,0x2100,0x7e4d9,0x102,0x6a092,0x51b25,0x47f64,0x140,0x2196c,0x70703,0x2dafe,0x3d9ec,0x76058,0x1e588,0x299ed,0x300,0x1b2d,0x7a4fc,0x7900,0x33d43,0x584f6,0x1a098,0xdcfc,0x1f00,0x21f09,0x9ebf,0x55b67,0xf9ed,0x4587c,0x54275,0x62635,0x2180,0x544de,0x25ab1,0xfa3e,0x4de76,0x71ab0,0x39f60,0x79e77,0x7fd09,0x8073,0x2a07b,0x218d0,0xda51,0x106b3,0xbd04,0x6a437,0x71dc7,0x260bc,0x2e10a,0x59a0,0xdc70,0x30093,0x1a9a,0x403a8,0xdd0e,0x7b851,0x20780,0x2a70e,0x4c417,0x500d0,0x324be,0x525a8,0x900,0x7d185,0x889b,0x95cf,0xc7b,0x17566,0x5b7c7,0x30b0b,0x2107c,0x7d0dd,0x5b7a4,0x7b2d4,0x1cda9,0x2ac90,0x35056,0x7f4be,0x790bd,0x501d,0x96da,0x93ab,0x4ab87,0x5cabf,0x697c7,0x332e,0x59581,0x4d05,0x5b5a5,0x7a8f0,0x57615,0x60f09,0x6c16,0x4d0db,0x473aa,0x6ce79,0x6a900,0x3d102,0x8925,0x4956e,0x769a8,0x4b032,0x1b61f,0x113e8,0x44af6,0x32ad0,0x6943e,0x9051,0x652f0,0x7934e,0x4d450,0x671a6,0x18906,0x4e921,0x3109e,0x714f0,0x377ef,0xb650,0x10da5,0x1b077,0x376b0,0x40eb3,0x511c5,0x30d8f,0x250f7,0x3896c};
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
const auto v29=vx(v14,v27);
const auto v30=vx(v9,v19);
const auto v31=vx(v5,v28);
const auto v32=vx(v3,v25);
const auto v33=vx(v3,v17);
const auto v34=vx(v8,v24);
const auto v35=vx(v21,v26);
out[0]=vx(vx(v19, v26), vx(v28, v32));
out[1]=v4;
out[2]=vx(vx(v10, v31), v34);
out[3]=v1;
out[4]=vx(vx(vx(v0, v1), vx(v8, v20)), vx(vx(v27, v28), v30));
out[5]=vx(vx(v23, v30), v33);
out[6]=v7;
out[7]=vx(v2, v16);
out[8]=v0;
out[9]=v11;
out[10]=vx(vx(v6, v7), vx(v12, v33));
out[11]=v22;
out[12]=vx(vx(vx(v8, v15), vx(v22, v29)), v35);
out[13]=v2;
out[14]=vx(vx(vx(v1, v9), vx(v11, v15)), v31);
out[15]=vx(vx(vx(v4, v12), vx(v16, v19)), vx(v25, v29));
out[16]=vx(vx(vx(v2, v11), vx(v14, v19)), v34);
out[17]=vx(vx(v6, v13), vx(v29, v32));
out[18]=vx(vx(v4, v9), vx(v18, v35));
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x52e00,0x8,0x15900,0x2,0x20503,0x61300,0x10,0x44,0x1,0x20,0x27e10,0x80,0x28c80,0x4,0x6a322,0x45048,0x28a24,0x18f00,0x3c508};
constexpr u32 inverse[19]={0x100,0x8,0x2000,0x2,0x40,0x200,0x2080,0x800,0x37911,0x1c5d2,0xffd2,0x5d9b1,0x41806,0x2e396,0x19bba,0x6b9b8,0x32318,0x387db,0x5233e};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
