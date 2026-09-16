#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x800,0x808,0x900,0x3193d,0x802,0xdc0f,0x45936,0x79d0e,0x840,0x1c849,0x60d51,0x4dd6d,0xc63,0x1186f,0x25946,0x5d7f,0x2800,0x5a828,0x66930,0xf92d,0x7ec03,0x2b82e,0x5fd07,0x3b91f,0x53c60,0x17c49,0x57941,0x2295d,0x2fc42,0x6686e,0x6e957,0x16d4e,0xa00,0x79a0d,0x48f01,0xf39,0x4da03,0x39e0b,0x40f36,0x5b0b,0x7da75,0x18a79,0x55b65,0x1b5c,0x30e57,0x58a5e,0x5df73,0x4b4f,0x16a05,0x37a28,0x3af34,0x2af2c,0x27e07,0xba2f,0x4eb02,0x53f1f,0x3ae50,0x7e7c,0x76f70,0x7af69,0xbe73,0x3ba5a,0x2f67,0x3b7b,0x880,0x58889,0x3ddb0,0x54d8c,0x31896,0x64c9a,0x49d92,0x2d9ab,0x388c5,0x7c8cd,0x659e4,0x109d9,0x9cf2,0x408ff,0x11de7,0x699df,0x3eca5,0x3ec8c,0x679a5,0x569b9,0x738b2,0x7ec9e,0x6fd86,0x5399f,0x578c0,0x4b8e8,0x6e9d1,0x439cc,0x1a8f6,0xbcdb,0x669d3,0x46dcb,0x5ceb0,0x7debc,0x29f81,0x39fb8,0x20ea7,0xcaae,0x10fa2,0xdb9e,0x19ec0,0x24ecd,0xcbe0,0xbd8,0x65af6,0x55efe,0x35fe2,0x34bdf,0x76a90,0xfabc,0x67b91,0x2fb88,0x76e86,0x2aaf,0x22fb3,0x67baf,0x62ec0,0x7eed,0x13bd0,0x47bc8,0x62ef7,0xaadf,0x56bd3,0xffce};
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
const auto v29=vx(v20,v24);
const auto v30=vx(v14,v27);
const auto v31=vx(v28,v29);
const auto v32=vx(v3,v25);
const auto v33=vx(v17,v26);
const auto v34=vx(v6,v13);
const auto v35=vx(v8,v30);
const auto v36=vx(v5,v19);
const auto v37=vx(v35,v36);
const auto v38=vx(v21,v33);
const auto v39=vx(v18,v23);
const auto v40=vx(v12,v31);
const auto v41=vx(v15,v34);
const auto v42=vx(v3,v10);
const auto v43=vx(v20,v39);
const auto v44=vx(v9,v32);
const auto v45=vx(v9,v43);
const auto v46=vx(v15,v32);
const auto v47=vx(v31,v34);
out[0]=vx(vx(v37, v38), vx(v42, v45));
out[1]=v4;
out[2]=vx(vx(vx(v5, v6), vx(v20, v27)), vx(v32, v38));
out[3]=v1;
out[4]=vx(v44, v47);
out[5]=vx(vx(v27, v40), vx(v41, v42));
out[6]=v7;
out[7]=v22;
out[8]=v2;
out[9]=v16;
out[10]=vx(vx(vx(v5, v9), vx(v10, v15)), vx(vx(v18, v24), vx(v28, v30)));
out[11]=v0;
out[12]=vx(vx(vx(v6, v17), vx(v29, v36)), v46);
out[13]=v11;
out[14]=vx(vx(v21, v37), v47);
out[15]=vx(vx(v33, v37), vx(v39, v40));
out[16]=vx(vx(vx(v23, v35), vx(v38, v40)), v46);
out[17]=vx(vx(vx(v13, v29), vx(v30, v33)), v44);
out[18]=vx(vx(vx(v12, v14), vx(v17, v19)), vx(vx(v28, v41), v45));
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x40000,0x8,0x400,0x2,0x10000,0x2000,0x10,0x80,0x4,0x40,0x1000,0x1,0x200,0x20,0x100,0x800,0x4000,0x20000,0x8000};
constexpr u32 inverse[19]={0x800,0x8,0x100,0x2,0x40,0x2000,0x200,0x80,0x4000,0x1000,0x4,0x8000,0x400,0x20,0x10000,0x40000,0x10,0x20000,0x1};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
