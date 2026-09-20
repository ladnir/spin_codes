#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x800,0x2000,0x808,0x30409,0x908,0x39950,0x48a09,0x43e50,0x40,0x6f1ca,0x6508b,0x38500,0x209d2,0x74000,0xda10,0x6b7c3,0x2,0x2d30,0x76434,0x46d07,0x33444,0xa12e,0xd37b,0x6210,0x5e4bc,0x31004,0x4d049,0x100f0,0x4d860,0x19480,0x16f9c,0x7077d,0xa00,0xe50d,0x7a856,0x4635a,0x4c329,0x7947c,0x7e276,0x79122,0x36025,0x556a2,0x292b0,0x78036,0x5a196,0x2f49,0xd00a,0x67ad4,0x3fc60,0x3165f,0x33a08,0xf436,0x40007,0x75260,0x4566,0x3300,0x57abb,0x3490e,0x3ec10,0x6fba4,0x8e46,0x505ab,0x29be4,0x43408,0x80,0x308ba,0x1e75b,0x1cb60,0x81,0xb0e3,0x56453,0x6f030,0x77190,0x2a020,0xc688,0x63339,0x5790b,0x310e3,0x64d1a,0x300f3,0x3b0e0,0xbde8,0x53305,0x51a0c,0x85af,0x30ff,0x28543,0x11412,0x12d0e,0x4f98c,0x1fe28,0x70eab,0x10db,0x67c01,0x440f4,0x1082f,0x45f8f,0x790b8,0x21a0a,0x2f13c,0x97af,0xe0c0,0x25123,0x1024d,0x4cfa,0x55a47,0x59bc,0x66b00,0x68c40,0x22a5,0x21a0f,0x790eb,0x4118d,0x7db88,0x53036,0x5de32,0x3ece3,0x39ebe,0x64e51,0x5180d,0x5ee06,0xfd89,0x29f7e,0x4a8f0,0x1bf2,0x6b025,0x3e983,0x66655};
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
const auto v29=vx(v6,v23);
const auto v30=vx(v5,v15);
const auto v31=vx(v8,v19);
const auto v32=vx(v21,v26);
const auto v33=vx(v20,v25);
const auto v34=vx(v4,v24);
const auto v35=vx(v18,v30);
const auto v36=vx(v12,v17);
const auto v37=vx(v17,v23);
const auto v38=vx(v6,v7);
const auto v39=vx(v29,v33);
const auto v40=vx(v11,v32);
const auto v41=vx(v27,v31);
const auto v42=vx(v0,v22);
out[0]=vx(vx(vx(v3, v4), vx(v12, v26)), vx(v39, v41));
out[1]=vx(v0, v11);
out[2]=vx(vx(vx(v5, v8), vx(v9, v17)), vx(v21, v38));
out[3]=v2;
out[4]=vx(vx(vx(v2, v9), vx(v10, v40)), v41);
out[5]=vx(vx(vx(v0, v25), vx(v26, v30)), v37);
out[6]=vx(v0, v7);
out[7]=v42;
out[8]=vx(v2, v4);
out[9]=v16;
out[10]=vx(vx(v14, v31), v34);
out[11]=v0;
out[12]=vx(vx(vx(v3, v18), vx(v19, v37)), v40);
out[13]=vx(v0, v1);
out[14]=vx(vx(v16, v28), vx(v34, v39));
out[15]=vx(vx(v7, v10), vx(v20, v35));
out[16]=vx(vx(vx(v1, v22), vx(v28, v32)), v35);
out[17]=vx(vx(vx(v9, v13), vx(v29, v36)), v42);
out[18]=vx(vx(v3, v16), vx(v36, v38));
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x34608,0x21,0x72f10,0x4,0x7fc24,0x59e01,0x11,0x81,0xc,0x40,0x28a08,0x1,0x5d620,0x3,0xb248,0x9f10,0x56682,0x3cc81,0x27e50};
constexpr u32 inverse[19]={0x800,0x2800,0x8,0x108,0x840,0x802,0x200,0x880,0x33c67,0x17b8a,0x983c,0x5bfd4,0x12d92,0x6d0c7,0x5aff6,0x6c5d7,0x649c8,0x20481,0x5f5cf};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
