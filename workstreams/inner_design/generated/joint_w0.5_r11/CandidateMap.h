#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x8,0x30000,0x20a,0x2d756,0x2000c,0x50407,0x360ef,0x5b1b0,0x800,0x2dc9c,0xee8,0x30f20,0x619c4,0xc95b,0x77dcd,0x7806,0x802,0x6c095,0x20bd8,0x5161b,0x8c6b,0x240ff,0x3ed50,0xf490,0x70c49,0x104a,0x50b79,0x3c22e,0x399e0,0x81e0,0xfc31,0x23165,0x10000,0x3bdb,0x3a086,0x34e09,0x7c0c3,0x2ff1b,0x402a4,0xe828,0x2e34f,0x20c00,0x4723,0x17d38,0x324c,0x4d900,0x3f4c1,0x6cad9,0x4b30e,0x404a,0x41250,0x13440,0xf7a0,0xe7,0x1341f,0x160c,0x5c02,0x57bd2,0xf9b6,0x40b32,0x96c,0x12abf,0x1ce39,0x138be,0x700c3,0x55821,0x22492,0x1a924,0x45e6,0x61907,0x40356,0x38ae3,0x35c90,0xd0e6,0x67c2b,0x42509,0x875,0x78000,0x44a2f,0x2170e,0x75070,0xc00d,0x75f9,0x630d0,0x29138,0x10546,0x4d650,0x6977a,0x40060,0x24489,0x32103,0x4b0be,0x5d0e8,0x79002,0x3936a,0x6d4,0xd0,0x63e1,0x78605,0x63060,0x38532,0x7e200,0x56106,0xd360,0x7b7c4,0x60061,0x35fb,0x570a,0x23e6,0x59040,0x6c338,0x2a5ca,0x5eb67,0x40c9,0x6c6a,0x41290,0x4eae8,0x54545,0xf04,0x75fd,0x55030,0x12f0a,0xd3d7,0x579b9,0x407f,0x3b46,0x4a179,0x50f14};
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
const auto v29=vx(v17,v20);
const auto v30=vx(v19,v21);
const auto v31=vx(v6,v23);
const auto v32=vx(v7,v11);
const auto v33=vx(v1,v13);
const auto v34=vx(v26,v29);
const auto v35=vx(v14,v18);
const auto v36=vx(v3,v30);
const auto v37=vx(v35,v36);
const auto v38=vx(v5,v8);
const auto v39=vx(v8,v20);
const auto v40=vx(v9,v22);
const auto v41=vx(v15,v29);
const auto v42=vx(v9,v16);
const auto v43=vx(v28,v33);
const auto v44=vx(v27,v28);
const auto v45=vx(v0,v6);
const auto v46=vx(v1,v16);
out[0]=vx(vx(v11, v31), vx(v37, v44));
out[1]=v32;
out[2]=vx(v4, v46);
out[3]=v0;
out[4]=vx(vx(vx(v30, v33), vx(v38, v42)), v45);
out[5]=vx(vx(vx(v1, v24), vx(v32, v34)), vx(v40, v45));
out[6]=vx(vx(v7, v35), vx(v41, v43));
out[7]=vx(vx(vx(v3, v7), vx(v15, v20)), vx(vx(v21, v22), vx(v27, v31)));
out[8]=vx(vx(vx(v10, v22), vx(v30, v34)), v43);
out[9]=vx(v2, v32);
out[10]=vx(vx(vx(v8, v12), vx(v14, v23)), v34);
out[11]=vx(v0, v7);
out[12]=vx(vx(v7, v12), vx(v37, v39));
out[13]=vx(vx(vx(v2, v11), vx(v13, v27)), vx(v38, v41));
out[14]=vx(vx(vx(v10, v18), vx(v24, v25)), vx(v31, v42));
out[15]=vx(vx(vx(v4, v19), vx(v39, v40)), v44);
out[16]=vx(v0, v16);
out[17]=v46;
out[18]=vx(vx(v7, v13), vx(v17, v19));
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x6d620,0x30,0x4a,0x1,0x8e43,0x38fb3,0x4ed12,0x59a90,0x25882,0x34,0x56800,0x11,0xbe10,0x30724,0x68e40,0x20588,0x41,0x42,0x7a110};
constexpr u32 inverse[19]={0x8,0x30008,0x202,0x20004,0x808,0x80a,0x10008,0x700cb,0x614f2,0xf704,0x6ddfd,0x2b26d,0x81b5,0x299a5,0x37bdb,0x6988c,0x55dc3,0x549c0,0x609d0};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
