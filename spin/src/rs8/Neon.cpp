// NEON Turbo (PacketRs8T64) kernels for AArch64 (NEON is baseline there).
// Same linear map as the oracle and the other kernels (spin_rs8_test neon).
//
// Inner: the shared generic inner (InnerGeneric.h) on 128-bit NEON vectors,
// with EOR3 three-way XOR where the SHA3 extension is available.
//
// Outer: the byte-plane circuit of rs8/Avx2.cpp (itself the AVX-512 circuit
// of rs8/Fast.cpp and rs8/ForwardFast.cpp), one group per 128-bit vector:
// eight coordinates are bit-transposed across registers so that each byte
// holds those coordinates of one OT, and GF(256) constant and fixed nibble
// maps are two TBL lookups. NEON's byte shifts make the transposes and the
// nibble split one instruction shorter than on AVX2.
#include "Rs8.h"
#if SPIN_AARCH64
#include "InnerGeneric.h"
#include "OuterMaps.h"
#include <arm_neon.h>
#include <cstddef>

namespace spin::detail::rs8 {
namespace {
using namespace outer_maps;
using V=uint8x16_t;

static SPIN_FORCEINLINE V xor2(V a,V b) { return veorq_u8(a,b); }
static SPIN_FORCEINLINE V xor3(V a,V b,V c) {
#if defined(__ARM_FEATURE_SHA3)
    return veor3q_u8(a,b,c);
#else
    return veorq_u8(veorq_u8(a,b),c);
#endif
}
struct Ops {
    using V=uint8x16_t;
    static SPIN_FORCEINLINE V load(const void* p) { return vld1q_u8(static_cast<const std::uint8_t*>(p)); }
    static SPIN_FORCEINLINE void store(void* p,V v) { vst1q_u8(static_cast<std::uint8_t*>(p),v); }
    static SPIN_FORCEINLINE V x2(V a,V b) { return xor2(a,b); }
    static SPIN_FORCEINLINE V x3(V a,V b,V c) { return xor3(a,b,c); }
    static SPIN_FORCEINLINE V zero() { return vdupq_n_u8(0); }
};

// ------------------------------------------------------------ primitives --
struct Table { V lo,hi; };
static SPIN_FORCEINLINE Table table(const NibbleMap& m) { return {vld1q_u8(m.lo),vld1q_u8(m.hi)}; }
struct Nibbles { V lo,hi; };
static SPIN_FORCEINLINE Nibbles nibbles(V x) { return {vandq_u8(x,vdupq_n_u8(0x0f)),vshrq_n_u8(x,4)}; }
static SPIN_FORCEINLINE V apply(const Nibbles& x,const Table& t) {
    return veorq_u8(vqtbl1q_u8(t.lo,x.lo),vqtbl1q_u8(t.hi,x.hi));
}
static SPIN_FORCEINLINE V apply(V x,const Table& t) { return apply(nibbles(x),t); }
static SPIN_FORCEINLINE Nibbles xorN(const Nibbles& a,const Nibbles& b) { return {xor2(a.lo,b.lo),xor2(a.hi,b.hi)}; }

// Bit k+D of a <-> bit k of b in every byte, for the bits k selected by Mask.
template<int D,int Mask>
static SPIN_FORCEINLINE void swapBits(V& a,V& b) {
    const auto t=vandq_u8(veorq_u8(vshrq_n_u8(a,D),b),vdupq_n_u8(Mask));
    b=veorq_u8(b,t);
    a=veorq_u8(a,vshlq_n_u8(t,D));
}
static SPIN_FORCEINLINE void inBlockRounds(V* x) {
    swapBits<2,0x33>(x[0],x[2]);swapBits<2,0x33>(x[1],x[3]);swapBits<2,0x33>(x[4],x[6]);swapBits<2,0x33>(x[5],x[7]);
    swapBits<1,0x55>(x[0],x[1]);swapBits<1,0x55>(x[2],x[3]);swapBits<1,0x55>(x[4],x[5]);swapBits<1,0x55>(x[6],x[7]);
}
// In place 8x8 bit transpose across registers, independently in every byte.
static SPIN_FORCEINLINE void transpose8(V* x) {
    swapBits<4,0x0f>(x[0],x[4]);swapBits<4,0x0f>(x[1],x[5]);swapBits<4,0x0f>(x[2],x[6]);swapBits<4,0x0f>(x[3],x[7]);
    inBlockRounds(x);
}
// Coordinates b[0..7] to byte planes: bit j of byte u of chunk r is bit
// 8u+r of coordinate j (a fixed OT order shared by every plane and symbol).
static SPIN_FORCEINLINE void toPlanes(const V* b,V* chunk) {
    for(unsigned m=0;m<8;++m) chunk[m]=b[m];
    transpose8(chunk);
}
static SPIN_FORCEINLINE void fromPlanes(const V* chunk,V* b) {
    for(unsigned m=0;m<8;++m) b[m]=chunk[m];
    transpose8(b);
}
// toPlanes then nibbles(), without the last round: the low (high) nibbles of
// x[r] are coordinates 0..3 of OT sets r (r+4), those of x[r+4] 4..7.
static SPIN_FORCEINLINE void toNibblePlanes(V* x,Nibbles* n) {
    inBlockRounds(x);
    const auto mask=vdupq_n_u8(0x0f);
    for(unsigned r=0;r<4;++r) {
        n[r]={vandq_u8(x[r],mask),vandq_u8(x[r+4],mask)};
        n[r+4]={vshrq_n_u8(x[r],4),vshrq_n_u8(x[r+4],4)};
    }
}

// ------------------------------------------------------------------ outer --
struct alignas(64) Planes { V v[8][4][8]; };
using Coeff=std::array<std::uint8_t,9>;
static SPIN_FORCEINLINE void flat9(const Nibbles& x0,const Nibbles& x1,const Nibbles& x2,const Nibbles& x3,
    const Table* c,V& y0,V& y1,V& y2,V& y3) {
    const auto x01=xorN(x0,x1),x23=xorN(x2,x3);
    const auto p0=apply(x0,c[0]),p1=apply(x1,c[1]),p2=apply(x01,c[2]);
    const auto q0=apply(x2,c[3]),q1=apply(x3,c[4]),q2=apply(x23,c[5]);
    const auto r0=apply(xorN(x0,x2),c[6]),r1=apply(xorN(x1,x3),c[7]);
    const auto r2=apply(xorN(x01,x23),c[8]);
    const auto pq=xor2(p0,q0),pr=xor2(p0,r0);
    y0=xor3(pq,p1,q1);y1=xor3(pq,p2,q2);y2=xor3(pr,p1,r1);y3=xor3(pr,p2,r2);
}
static SPIN_FORCEINLINE void flat9Adjoint(V& x0,V& x1,V& x2,V& x3,const Table* m) {
    const auto y0=x0,y1=x1,y2=x2,y3=x3;
    const auto y01=xor2(y0,y1),y23=xor2(y2,y3);
    const auto a0=apply(xor2(y01,y23),m[0]),a1=apply(xor2(y0,y2),m[1]),a2=apply(xor2(y1,y3),m[2]);
    const auto a3=apply(y01,m[3]),a4=apply(y0,m[4]),a5=apply(y1,m[5]);
    const auto a6=apply(y23,m[6]),a7=apply(y2,m[7]),a8=apply(y3,m[8]);
    const auto s68=xor2(a6,a8),s78=xor2(a7,a8);
    x0=xor3(a0,a2,s68);x1=xor3(a1,a2,s78);x2=xor3(a3,a5,s68);x3=xor3(a4,a5,s78);
}
template<bool Adjoint>
static SPIN_FORCEINLINE void coefficients(const Coeff& c,Table* t) {
    const auto* m=Adjoint?maps.adjoint:maps.multiply;
    for(unsigned i=0;i<9;++i) t[i]=table(m[c[i]]);
}
constexpr unsigned messageBlock(unsigned q,unsigned i,unsigned j) {
    return 16*(2*i+(j>>2))+4*q+(j&3);
}
static SPIN_NOINLINE void groupTranspose(const Block* in,Block* out,const Coeff* c) {
    Planes p;
    for(unsigned s=0;s<8;++s) {
        Nibbles n[4][8];
        for(unsigned i=0;i<4;++i) {
            V b[8];
            for(unsigned j=0;j<8;++j) b[j]=Ops::load(in+32*s+8*i+j);
            toNibblePlanes(b,n[i]);
        }
        Table t[9];
        coefficients<false>(c[s],t);
        for(unsigned m=0;m<8;++m)
            flat9(n[0][m],n[1][m],n[2][m],n[3][m],t,p.v[s][0][m],p.v[s][1][m],p.v[s][2][m],p.v[s][3][m]);
    }
    const auto t6=table(a6),t8=table(a8);
    for(unsigned i=0;i<4;++i)
        for(unsigned m=0;m<8;++m) {
            auto x0=p.v[4][i][m],x1=p.v[5][i][m],x2=p.v[6][i][m],x3=p.v[7][i][m];
            x0=xor2(x0,x1);x2=xor2(x2,x3);x0=xor2(x0,x2);x1=xor2(x1,x3);
            const auto n0=nibbles(x0);
            const auto u6=apply(n0,t6),u8=apply(n0,t8);
            auto y0=x0,y1=xor2(x1,u8),y2=xor2(x2,u6);
            auto y3=xor3(xor3(x3,u6,x0),apply(x1,t6),apply(x2,t8));
            y0=xor2(y0,y1);y2=xor2(y2,y3);y0=xor2(y0,y2);y1=xor2(y1,y3);
            p.v[0][i][m]=xor2(p.v[0][i][m],y0);p.v[1][i][m]=xor2(p.v[1][i][m],y1);
            p.v[2][i][m]=xor2(p.v[2][i][m],y2);p.v[3][i][m]=xor2(p.v[3][i][m],y3);
        }
    for(unsigned q=0;q<4;++q)
        for(unsigned i=0;i<4;++i) {
            V b[8];
            fromPlanes(p.v[q][i],b);
            for(unsigned j=0;j<8;++j) Ops::store(out+messageBlock(q,i,j),b[j]);
        }
}
static SPIN_NOINLINE void groupForward(const Block* message,unsigned lanes,unsigned lane,Block* out,const Coeff* c) {
    Planes p;
    for(unsigned q=0;q<4;++q)
        for(unsigned i=0;i<4;++i) {
            V b[8];
            for(unsigned j=0;j<8;++j) b[j]=Ops::load(message+messageBlock(q,i,j)*lanes+lane);
            toPlanes(b,p.v[q][i]);
        }
    const auto t6=table(a6t),t8=table(a8t);
    for(unsigned i=0;i<4;++i)
        for(unsigned m=0;m<8;++m) {
            auto y0=p.v[0][i][m],y1=p.v[1][i][m],y2=p.v[2][i][m],y3=p.v[3][i][m];
            y3=xor2(y3,y1);y2=xor2(y2,y0);y3=xor2(y3,y2);y1=xor2(y1,y0);
            const auto n3=nibbles(y3);
            auto x0=xor3(xor3(y0,y3,apply(y1,t8)),apply(y2,t6),apply(n3,t6));
            auto x1=xor2(y1,apply(n3,t6)),x2=xor2(y2,apply(n3,t8)),x3=y3;
            x3=xor2(x3,x1);x2=xor2(x2,x0);x3=xor2(x3,x2);x1=xor2(x1,x0);
            p.v[4][i][m]=x0;p.v[5][i][m]=x1;p.v[6][i][m]=x2;p.v[7][i][m]=x3;
        }
    for(unsigned s=0;s<8;++s) {
        Table t[9];
        coefficients<true>(c[s],t);
        for(unsigned m=0;m<8;++m) flat9Adjoint(p.v[s][0][m],p.v[s][1][m],p.v[s][2][m],p.v[s][3][m],t);
        for(unsigned i=0;i<4;++i) {
            V b[8];
            fromPlanes(p.v[s][i],b);
            for(unsigned j=0;j<8;++j) Ops::store(out+32*s+8*i+j,b[j]);
        }
    }
}
}

void transposeNeon(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    inner::reverse<Ops>(input,scratch,*plan.certified22Plan);
    for(std::size_t g=0;g<plan.groups;++g)
        groupTranspose(scratch+groupStride*g,output+128*g,plan.field.data()+8*g);
}
void forwardNeon(const Block* input,Block* output,Block* scratch,const Plan& plan,unsigned lanes) {
    const auto plane=plan.scratchBlocks();
    for(std::size_t g=0;g<plan.groups;++g)
        for(unsigned lane=0;lane<lanes;++lane)
            groupForward(input+128*lanes*g,lanes,lane,scratch+lane*plane+groupStride*g,plan.field.data()+8*g);
    inner::forwardLanes<Ops>(scratch,output,*plan.certified22Plan,lanes,plane);
}
}
#endif
