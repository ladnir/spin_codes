#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x7f,0x80,0x100,0x200,0x38f,0x400,0x5b3,0x6d5,0x769,0x800,0x1000,0x2000,0x380f,0x4000,0x5833,0x6855,0x7069,0x8000,0x9983,0xaa85,0xb309,0xcc91,0xd521,0xe641,0xfffe,0x10000,0x20000,0x27b24,0x17b2b,0x1b929,0x2b91a,0x2c258,0x1c264,0x155d4,0x25457,0x22c75,0x12df9,0x1e86c,0x2e9dc,0x29198,0x19027,0x2ff03,0x1e700,0x1ac22,0x2b42e,0x20e3b,0x1160b,0x15d4f,0x24570,0x22256,0x13bd6,0x173f2,0x26a7d,0x2d7ff,0x1ce4c,0x1860e,0x29fb2,0x40000,0x72850,0x4d8f0,0x7f0af,0x7304e,0x4182d,0x7e8eb,0x4c087,0x4df8a,0x7f659,0x405ff,0x72c23,0x7eb55,0x4c2b5,0x73175,0x4189a,0x4c921,0x7f972,0x439d4,0x70988,0x7b17e,0x4811e,0x741de,0x471b1,0x49e2a,0x7affa,0x46c5a,0x75d85,0x7e2e4,0x4d307,0x710c1,0x4212d,0x67d6e,0x6553d,0x5debf,0x5f6e3,0x5f418,0x5dc78,0x6579c,0x67ff3,0x6f7b1,0x6de61,0x556e5,0x57f3a,0x57a56,0x553b5,0x6db57,0x6f2bb,0x5434d,0x5731d,0x6c899,0x6f8c6,0x6822a,0x6b249,0x509ab,0x539c7,0x54113,0x570c0,0x6c842,0x6f99e,0x684e5,0x6b505,0x50de1,0x53c0e};
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
const auto v29=vx(v0,v21);
const auto v30=vx(v27,v29);
const auto v31=vx(v15,v17);
const auto v32=vx(v3,v30);
const auto v33=vx(v8,v26);
const auto v34=vx(v18,v32);
const auto v35=vx(v10,v23);
const auto v36=vx(v14,v31);
const auto v37=vx(v13,v28);
const auto v38=vx(v1,v9);
const auto v39=vx(v20,v34);
const auto v40=vx(v5,v25);
const auto v41=vx(v6,v33);
const auto v42=vx(v24,v37);
const auto v43=vx(v12,v41);
const auto v44=vx(v38,v40);
const auto v45=vx(v19,v35);
const auto v46=vx(v11,v36);
const auto v47=vx(v0,v4);
const auto v48=vx(v7,v43);
const auto v49=vx(v36,v39);
const auto v50=vx(v42,v44);
const auto v51=vx(v31,v45);
const auto v52=vx(v12,v42);
const auto v53=vx(v2,v9);
out[0]=v0;
out[1]=vx(v0, v1);
out[2]=vx(v0, v2);
out[3]=vx(vx(v2, v49), v50);
out[4]=v47;
out[5]=vx(vx(vx(v1, v4), vx(v37, v39)), v43);
out[6]=vx(vx(vx(v2, v18), vx(v27, v33)), vx(vx(v40, v47), v51));
out[7]=vx(v0, v7);
out[8]=vx(vx(vx(v6, v7), vx(v26, v29)), vx(v50, v51));
out[9]=vx(vx(vx(v7, v8), vx(v17, v34)), vx(v52, v53));
out[10]=vx(vx(vx(v4, v13), vx(v17, v30)), vx(v45, v48));
out[11]=vx(v0, v11);
out[12]=vx(vx(vx(v24, v30), vx(v35, v44)), v46);
out[13]=vx(vx(vx(v5, v11), vx(v13, v14)), vx(vx(v33, v35), vx(v39, v53)));
out[14]=vx(vx(vx(v3, v4), vx(v29, v46)), v52);
out[15]=vx(vx(v32, v46), v48);
out[16]=vx(v0, v16);
out[17]=vx(vx(vx(v10, v16), vx(v28, v38)), vx(v43, v49));
out[18]=vx(v0, v22);
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x1,0x3,0x5,0x76607,0x9,0x41f0b,0x1a70d,0x11,0x44913,0x79015,0x43019,0x21,0xe823,0x3a125,0x6d929,0xf231,0x41,0x6ce43,0x81};
constexpr u32 inverse[19]={0x1,0x3,0x5,0x11,0x81,0x801,0x10001,0x40001,0x2e11,0xd8c9,0x8a0b,0x37c8a,0x334c1,0x304fb,0x367b3,0x307c1,0x7fbc,0x7a92,0x34ab};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
