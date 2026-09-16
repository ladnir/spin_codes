#pragma once
#include "Inner.h"
namespace bare_spin {
struct CandidateMap : Map128S19 {
static constexpr std::array<std::uint32_t,T> columns{0x1000,0x280,0x1142,0x7c0d4,0x1010,0x70f8b,0xee57,0x32da,0x200,0x7cd1e,0x58c01,0x58009,0x9a09,0x580c,0x5eb0d,0x2ea1e,0x40,0xc655,0x205ef,0x500ec,0x7c342,0x84c,0x539e8,0x531f0,0x4cb1b,0x3d090,0x341f7,0x3996a,0x39000,0x38690,0x4e5e9,0x3306f,0x100,0x1f209,0x1ad67,0x79d78,0x11ba0,0x7e5b2,0x48c2,0x175c6,0xbe05,0x69092,0x49d21,0x570a0,0x13cbc,0x1f30,0x5e09d,0x30007,0x2f03c,0x3d7a0,0x158b6,0x7bc3c,0x4298e,0x20309,0x77e01,0x69790,0x69662,0x6c60,0xb1ab,0x188bf,0xd7c9,0x120d0,0x60f05,0x3b0a,0x81,0x30860,0x43d01,0xf6f6,0x40d73,0x889,0xcff6,0x3091a,0x6d58f,0x200f0,0x7674c,0x47125,0x24064,0x19800,0x30da2,0x716d0,0x3da0d,0x679,0x5e360,0x1fc02,0x14ed,0x4c582,0x6d285,0x5c0fc,0x1d658,0x5d7b2,0x26076,0x1a28a,0x280a1,0x18c50,0x1c98a,0x5066d,0x397ed,0x17e85,0x60748,0x32d36,0x680af,0x364dc,0x3ef0f,0x1c86a,0x5efe6,0xdb10,0x5f000,0x707e0,0x60bd,0x25950,0x805e,0x57aa5,0x2ac1d,0x91e0,0x53855,0xc6be,0x784d,0x548ab,0x71300,0x5e0f0,0xd4d,0x5ed2e,0x21646,0x3533,0x24104,0xac7c,0xa50a,0x58b64};
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
const auto v29=vx(v10,v24);
const auto v30=vx(v12,v20);
const auto v31=vx(v8,v19);
const auto v32=vx(v6,v23);
const auto v33=vx(v4,v30);
const auto v34=vx(v9,v32);
const auto v35=vx(v7,v22);
const auto v36=vx(v14,v31);
const auto v37=vx(v9,v27);
const auto v38=vx(v10,v26);
const auto v39=vx(v18,v29);
const auto v40=vx(v13,v23);
const auto v41=vx(v21,v38);
out[0]=vx(vx(v0, v1), v35);
out[1]=vx(vx(v2, v11), v16);
out[2]=vx(vx(vx(v5, v18), vx(v25, v30)), v37);
out[3]=vx(vx(v12, v19), v34);
out[4]=v4;
out[5]=vx(vx(v25, v34), v39);
out[6]=vx(v0, v11);
out[7]=vx(v1, v7);
out[8]=vx(v0, v16);
out[9]=vx(v0, v7);
out[10]=vx(vx(vx(v3, v18), vx(v21, v33)), v36);
out[11]=vx(vx(vx(v1, v4), vx(v31, v37)), v41);
out[12]=v0;
out[13]=vx(vx(v14, v16), vx(v27, v39));
out[14]=vx(vx(vx(v3, v7), vx(v8, v15)), vx(v29, v40));
out[15]=vx(vx(v13, v17), vx(v19, v35));
out[16]=vx(vx(v20, v28), v41);
out[17]=vx(vx(vx(v11, v22), vx(v24, v26)), vx(v33, v40));
out[18]=vx(vx(v2, v7), vx(v24, v36));
}
static void conjugate(u32* rows) {
constexpr u32 V[19]={0x93,0x64,0x1bf00,0x46d00,0x8,0x68e00,0x21,0x12,0x41,0x11,0xbe08,0x7fc0a,0x1,0x6b340,0x64810,0x7a190,0x5f900,0x43a8,0x28a14};
constexpr u32 inverse[19]={0x1000,0x1280,0x142,0x10,0x1200,0x1040,0x1100,0x1081,0x713b7,0x4f2cf,0x6950e,0x3a9ec,0x356ed,0x14898,0x1e1a9,0x3d448,0x717a3,0x9c29,0x2866c};
u32 left[19]{},result[19]{};
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((V[i]>>j)&1) left[i]^=rows[j];
for(unsigned i=0;i<19;++i) for(unsigned j=0;j<19;++j) if((left[i]>>j)&1) result[i]^=inverse[j];
std::memcpy(rows,result,sizeof(result));
}
};
}
