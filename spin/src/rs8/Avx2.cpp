// AVX2 Turbo (PacketRs8T64) kernels: transpose and forward, no AVX-512,
// VBMI or GFNI. Same linear map as the scalar oracle and the AVX-512 kernels.
//
// Inner: the fixed s22 step circuit is the generic lane-major XOR network of
// certified22/StepCircuit.h in the fast basis w = basis*s. The per-epoch GL22
// update w' = U w + f reads the packed 3x3-block tables prepared at setup
// (updates for the transpose, forwardUpdates for forward) and applies them by
// the method of four Russians: six data tables of nibble combinations of w,
// six lookups per row.
//
// Outer: the AVX-512 circuit (rs8/Fast.cpp, rs8/ForwardFast.cpp) on byte
// planes. Eight 128-bit coordinates are bit-transposed across registers so
// that each byte holds those eight coordinates of one OT; GF(256) constant and
// fixed affine byte maps are two nibble lookups (pshufb). Each 256-bit vector
// carries two independent jobs in its 128-bit halves (two groups, or two
// record lanes of one group); pshufb is lane-local, so each half has its own
// tables.
#include "Rs8.h"
#include "../certified22/productionMap.h"
#include "../certified22/LaneMajor22.h"
#include "../certified22/StepCircuit.h"
#include <bit>
#include <cstddef>

namespace spin::detail::rs8 {
namespace {
using namespace certified22::step;
using V=__m256i;

// ---------------------------------------------------------------- tables --
constexpr std::uint8_t mul8(std::uint8_t a,std::uint8_t b) {
    std::uint8_t result=0;
    for(unsigned i=0;i<8;++i) {
        if(b&1U) result^=a;
        a=std::uint8_t((unsigned(a)<<1)^((a&128U)?0x11bU:0U));
        b>>=1;
    }
    return result;
}
// GFNI affine semantics: bit i of the result is parity(byte 7-i of m & x).
constexpr std::uint8_t affine(std::uint64_t m,std::uint8_t x) {
    std::uint8_t result=0;
    for(unsigned i=0;i<8;++i)
        result|=std::uint8_t((std::popcount(unsigned((m>>(8*(7-i)))&x)&255U)&1U)<<i);
    return result;
}
constexpr std::uint64_t transposeAffine(std::uint64_t m) {
    std::uint64_t result=0;
    for(unsigned i=0;i<8;++i) for(unsigned j=0;j<8;++j)
        result|=((m>>(8*(7-i)+j))&1ULL)<<(8*(7-j)+i);
    return result;
}
// A GF(2)-linear byte map f as two nibble tables: f(x)=lo[x&15]^hi[x>>4].
struct alignas(32) NibbleMap { std::uint8_t lo[16]{},hi[16]{}; };
template<class F> constexpr NibbleMap nibbleMap(F f) {
    NibbleMap result{};
    for(unsigned v=0;v<16;++v) {
        result.lo[v]=f(std::uint8_t(v));
        result.hi[v]=f(std::uint8_t(v<<4));
    }
    return result;
}
struct Maps { NibbleMap multiply[256],adjoint[256]; };
// Built from the eight columns c*e_b of each multiplication (cheap enough for
// default constexpr step limits).
constexpr Maps makeMaps() {
    Maps maps{};
    for(unsigned c=0;c<256;++c) {
        std::uint8_t column[8]{};
        for(unsigned b=0;b<8;++b) column[b]=mul8(std::uint8_t(c),std::uint8_t(1U<<b));
        const auto times=[&](unsigned x) {
            std::uint8_t r=0;
            for(unsigned b=0;b<8;++b) if((x>>b)&1U) r^=column[b];
            return r;
        };
        // Adjoint of multiplication by c: bit b is parity(c*e_b & y).
        const auto adjoint=[&](unsigned y) {
            std::uint8_t r=0;
            for(unsigned b=0;b<8;++b) r|=std::uint8_t((std::popcount(unsigned(column[b]&y))&1U)<<b);
            return r;
        };
        for(unsigned v=0;v<16;++v) {
            maps.multiply[c].lo[v]=times(v);maps.multiply[c].hi[v]=times(v<<4);
            maps.adjoint[c].lo[v]=adjoint(v);maps.adjoint[c].hi[v]=adjoint(v<<4);
        }
    }
    return maps;
}
alignas(64) constexpr Maps maps=makeMaps();
// The parity network's adjoint nibble multiplications (rs8/Fast.cpp) and
// their transposes (rs8/ForwardFast.cpp).
constexpr std::uint64_t m6=0x060c0b0560c0b050ULL,m8=0x0803060c803060c0ULL;
alignas(32) constexpr NibbleMap a6=nibbleMap([](std::uint8_t x){ return affine(m6,x); });
alignas(32) constexpr NibbleMap a8=nibbleMap([](std::uint8_t x){ return affine(m8,x); });
alignas(32) constexpr NibbleMap a6t=nibbleMap([](std::uint8_t x){ return affine(transposeAffine(m6),x); });
alignas(32) constexpr NibbleMap a8t=nibbleMap([](std::uint8_t x){ return affine(transposeAffine(m8),x); });

// ------------------------------------------------------------ primitives --
static SPIN_FORCEINLINE __m128i load(const void* p) { return _mm_loadu_si128(static_cast<const __m128i*>(p)); }
static SPIN_FORCEINLINE V join(__m128i a,__m128i b) { return _mm256_inserti128_si256(_mm256_castsi128_si256(a),b,1); }
struct Table { V lo,hi; };
static SPIN_FORCEINLINE Table table(const NibbleMap& a,const NibbleMap& b) {
    return {join(load(a.lo),load(b.lo)),join(load(a.hi),load(b.hi))};
}
static SPIN_FORCEINLINE Table table(const NibbleMap& a) {
    return {_mm256_broadcastsi128_si256(load(a.lo)),_mm256_broadcastsi128_si256(load(a.hi))};
}
struct Nibbles { V lo,hi; };
static SPIN_FORCEINLINE Nibbles nibbles(V x) {
    const auto mask=_mm256_set1_epi8(0x0f);
    return {_mm256_and_si256(x,mask),_mm256_and_si256(_mm256_srli_epi16(x,4),mask)};
}
static SPIN_FORCEINLINE V apply(const Nibbles& x,const Table& t) {
    return _mm256_xor_si256(_mm256_shuffle_epi8(t.lo,x.lo),_mm256_shuffle_epi8(t.hi,x.hi));
}
static SPIN_FORCEINLINE V apply(V x,const Table& t) { return apply(nibbles(x),t); }
static SPIN_FORCEINLINE V xor2(V a,V b) { return _mm256_xor_si256(a,b); }
static SPIN_FORCEINLINE V xor3(V a,V b,V c) { return _mm256_xor_si256(_mm256_xor_si256(a,b),c); }

// Bit k+D of a <-> bit k of b in every byte, for the bits k selected by Mask.
template<int D,int Mask>
static SPIN_FORCEINLINE void swapBits(V& a,V& b) {
    const auto t=_mm256_and_si256(_mm256_xor_si256(_mm256_srli_epi16(a,D),b),_mm256_set1_epi8(char(Mask)));
    b=_mm256_xor_si256(b,t);
    a=_mm256_xor_si256(a,_mm256_slli_epi16(t,D));
}
// In place 8x8 bit transpose across registers, independently in every byte:
// bit k of byte u of x[j] <-> bit j of byte u of x[k].
static SPIN_FORCEINLINE void transpose8(V* x) {
    swapBits<4,0x0f>(x[0],x[4]);swapBits<4,0x0f>(x[1],x[5]);swapBits<4,0x0f>(x[2],x[6]);swapBits<4,0x0f>(x[3],x[7]);
    swapBits<2,0x33>(x[0],x[2]);swapBits<2,0x33>(x[1],x[3]);swapBits<2,0x33>(x[4],x[6]);swapBits<2,0x33>(x[5],x[7]);
    swapBits<1,0x55>(x[0],x[1]);swapBits<1,0x55>(x[2],x[3]);swapBits<1,0x55>(x[4],x[5]);swapBits<1,0x55>(x[6],x[7]);
}
// Coordinates b[0..7] to byte planes: bit j of byte u of chunk r is bit
// 8u+r of coordinate j. The byte maps act on every byte alike, so this OT
// order (a fixed permutation, shared by every plane and symbol) is free.
static SPIN_FORCEINLINE void toPlanes(const V* b,V* chunk) {
    for(unsigned m=0;m<8;++m) chunk[m]=b[m];
    transpose8(chunk);
}
// Exact inverse of toPlanes (the transpose is an involution).
static SPIN_FORCEINLINE void fromPlanes(const V* chunk,V* b) {
    for(unsigned m=0;m<8;++m) b[m]=chunk[m];
    transpose8(b);
}
// toPlanes followed by nibbles(), without the last transpose round: after the
// in-block rounds, the low (high) nibbles of x[r] are coordinates 0..3 of OT
// sets r (r+4), and those of x[r+4] coordinates 4..7.
static SPIN_FORCEINLINE void toNibblePlanes(V* x,Nibbles* n) {
    swapBits<2,0x33>(x[0],x[2]);swapBits<2,0x33>(x[1],x[3]);swapBits<2,0x33>(x[4],x[6]);swapBits<2,0x33>(x[5],x[7]);
    swapBits<1,0x55>(x[0],x[1]);swapBits<1,0x55>(x[2],x[3]);swapBits<1,0x55>(x[4],x[5]);swapBits<1,0x55>(x[6],x[7]);
    const auto mask=_mm256_set1_epi8(0x0f);
    for(unsigned r=0;r<4;++r) {
        n[r]={_mm256_and_si256(x[r],mask),_mm256_and_si256(x[r+4],mask)};
        n[r+4]={_mm256_and_si256(_mm256_srli_epi16(x[r],4),mask),_mm256_and_si256(_mm256_srli_epi16(x[r+4],4),mask)};
    }
}
static SPIN_FORCEINLINE Nibbles xorN(const Nibbles& a,const Nibbles& b) {
    return {_mm256_xor_si256(a.lo,b.lo),_mm256_xor_si256(a.hi,b.hi)};
}

// ------------------------------------------------------------------ outer --
// Byte planes of one group per half: [symbol][byte][chunk].
struct alignas(64) Planes { V v[8][4][8]; };
using Coeff=std::array<std::uint8_t,9>;

// The exact nine-product tower map of rs8/Fast.cpp (multiplication by the
// sampled GF(2^32) scalar), flattened.
// Inputs are nibble-split byte planes (linear combinations are taken on the
// split nibbles), outputs are whole bytes.
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
// Transposed flat9 of rs8/ForwardFast.cpp with adjoint byte maps.
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
static SPIN_FORCEINLINE void coefficients(const Coeff& a,const Coeff& b,Table* t) {
    const auto* m=Adjoint?maps.adjoint:maps.multiply;
    for(unsigned i=0;i<9;++i) t[i]=table(m[a[i]],m[b[i]]);
}

// Systematic message bit j of byte plane i of column q (one RS row nibble
// per half byte): the message block written by storePair in rs8/Fast.cpp.
constexpr unsigned messageBlock(unsigned q,unsigned i,unsigned j) {
    return 16*(2*i+(j>>2))+4*q+(j&3);
}

// Transpose outer of one group per half: in(k) loads scratch coordinate k
// (0..255) of both jobs, out(k,v) stores message block k (0..127).
template<class In,class Out>
static SPIN_FORCEINLINE void outerTranspose(In in,Out out,const Coeff* ca,const Coeff* cb) {
    Planes p;
    for(unsigned s=0;s<8;++s) {
        alignas(64) Nibbles n[4][8];
        for(unsigned i=0;i<4;++i) {
            V b[8];
            for(unsigned j=0;j<8;++j) b[j]=in(32*s+8*i+j);
            toNibblePlanes(b,n[i]);
        }
        Table t[9];
        coefficients<false>(ca[s],cb[s],t);
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
            for(unsigned j=0;j<8;++j) out(messageBlock(q,i,j),b[j]);
        }
}
// Forward outer, the exact adjoint: in(k) loads message block k, out(k,v)
// stores scratch coordinate k.
template<class In,class Out>
static SPIN_FORCEINLINE void outerForward(In in,Out out,const Coeff* ca,const Coeff* cb) {
    Planes p;
    for(unsigned q=0;q<4;++q)
        for(unsigned i=0;i<4;++i) {
            V b[8];
            for(unsigned j=0;j<8;++j) b[j]=in(messageBlock(q,i,j));
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
        coefficients<true>(ca[s],cb[s],t);
        for(unsigned m=0;m<8;++m) flat9Adjoint(p.v[s][0][m],p.v[s][1][m],p.v[s][2][m],p.v[s][3][m],t);
        for(unsigned i=0;i<4;++i) {
            V b[8];
            fromPlanes(p.v[s][i],b);
            for(unsigned j=0;j<8;++j) out(32*s+8*i+j,b[j]);
        }
    }
}

// Two groups (or one group twice) of contiguous blocks.
struct GroupsIn {
    const Block* a; const Block* b;
    SPIN_FORCEINLINE V operator()(unsigned k) const { return join(load(a+k),load(b+k)); }
};
// Streaming stores are used only with 64-byte aligned output: each (q,i)
// plane of the transpose outer writes two complete cache lines.
template<bool Both,bool Stream=false> struct GroupsOut {
    Block* a; Block* b;
    SPIN_FORCEINLINE void operator()(unsigned k,V v) const {
        if constexpr(Stream) {
            _mm_stream_si128(&a[k].mData,_mm256_castsi256_si128(v));
            if constexpr(Both) _mm_stream_si128(&b[k].mData,_mm256_extracti128_si256(v,1));
        } else {
            _mm_storeu_si128(&a[k].mData,_mm256_castsi256_si128(v));
            if constexpr(Both) _mm_storeu_si128(&b[k].mData,_mm256_extracti128_si256(v,1));
        }
    }
};
// Two adjacent record lanes of an L-lane interleaved message.
template<unsigned L> struct RecordsIn {
    const Block* base;
    SPIN_FORCEINLINE V operator()(unsigned k) const { return _mm256_loadu_si256(reinterpret_cast<const V*>(base+L*k)); }
};

template<bool Stream>
static SPIN_NOINLINE void transposeOuter(const Block* scratch,Block* output,const Plan& plan) {
    const auto groups=plan.groups;
    const auto* c=plan.field.data();
    std::size_t g=0;
    for(;g+1<groups;g+=2)
        outerTranspose(GroupsIn{scratch+groupStride*g,scratch+groupStride*(g+1)},
                       GroupsOut<true,Stream>{output+128*g,output+128*(g+1)},c+8*g,c+8*(g+1));
    if(g<groups)
        outerTranspose(GroupsIn{scratch+groupStride*g,scratch+groupStride*g},
                       GroupsOut<false,Stream>{output+128*g,nullptr},c+8*g,c+8*g);
    if constexpr(Stream) _mm_sfence();
}
// Large outputs stream, as the AVX-512 outer does above K18.
static void transposeOuter(const Block* scratch,Block* output,const Plan& plan) {
    if(plan.k>262144 && !(reinterpret_cast<std::uintptr_t>(output)&63U)) transposeOuter<true>(scratch,output,plan);
    else transposeOuter<false>(scratch,output,plan);
}
static SPIN_NOINLINE void forwardOuter(const Block* message,Block* scratch,const Plan& plan,unsigned lanes) {
    const auto groups=plan.groups,plane=plan.scratchBlocks();
    const auto* c=plan.field.data();
    if(lanes==1) {
        std::size_t g=0;
        for(;g+1<groups;g+=2)
            outerForward(GroupsIn{message+128*g,message+128*(g+1)},
                         GroupsOut<true>{scratch+groupStride*g,scratch+groupStride*(g+1)},c+8*g,c+8*(g+1));
        if(g<groups)
            outerForward(GroupsIn{message+128*g,message+128*g},
                         GroupsOut<false>{scratch+groupStride*g,nullptr},c+8*g,c+8*g);
    } else if(lanes==2) {
        for(std::size_t g=0;g<groups;++g)
            outerForward(RecordsIn<2>{message+2*128*g},
                         GroupsOut<true>{scratch+groupStride*g,scratch+plane+groupStride*g},c+8*g,c+8*g);
    } else {
        for(std::size_t g=0;g<groups;++g)
            for(unsigned pair=0;pair<2;++pair)
                outerForward(RecordsIn<4>{message+4*128*g+2*pair},
                             GroupsOut<true>{scratch+2*pair*plane+groupStride*g,scratch+(2*pair+1)*plane+groupStride*g},
                             c+8*g,c+8*g);
    }
}

// ------------------------------------------------------------------ inner --
struct X128 {
    static SPIN_FORCEINLINE __m128i x2(__m128i a,__m128i b) { return _mm_xor_si128(a,b); }
    static SPIN_FORCEINLINE __m128i x3(__m128i a,__m128i b,__m128i c) { return _mm_xor_si128(_mm_xor_si128(a,b),c); }
    static SPIN_FORCEINLINE __m128i zero() { return _mm_setzero_si128(); }
};
struct X256 {
    static SPIN_FORCEINLINE __m256i x2(__m256i a,__m256i b) { return _mm256_xor_si256(a,b); }
    static SPIN_FORCEINLINE __m256i x3(__m256i a,__m256i b,__m256i c) { return _mm256_xor_si256(_mm256_xor_si256(a,b),c); }
    static SPIN_FORCEINLINE __m256i zero() { return _mm256_setzero_si256(); }
};
struct Y2 { __m256i a,b; };
struct X512 {
    static SPIN_FORCEINLINE Y2 x2(Y2 p,Y2 q) { return {_mm256_xor_si256(p.a,q.a),_mm256_xor_si256(p.b,q.b)}; }
    static SPIN_FORCEINLINE Y2 x3(Y2 p,Y2 q,Y2 r) { return x2(x2(p,q),r); }
    static SPIN_FORCEINLINE Y2 zero() { return {_mm256_setzero_si256(),_mm256_setzero_si256()}; }
};

// w = U w + f for the packed 3x3-block GL22 table of one epoch (qword 3o+i,
// byte j is byte i of row 8o+j; bit c of row r is U[r][c]).
// The tables of nibble chunks 2i and 2i+1 (columns 8i..8i+3 and 8i+4..8i+7)
// are built together, one pair per entry, at twice the width.
struct Y4 { Y2 a,b; };
struct X1024 {
    static SPIN_FORCEINLINE Y4 x2(Y4 p,Y4 q) { return {X512::x2(p.a,q.a),X512::x2(p.b,q.b)}; }
    static SPIN_FORCEINLINE Y4 zero() { return {X512::zero(),X512::zero()}; }
};
template<class W> struct PairOf;
template<> struct PairOf<__m128i> {
    using E=__m256i; using XE=X256;
    static SPIN_FORCEINLINE E make(__m128i a,__m128i b) { return join(a,b); }
};
template<> struct PairOf<__m256i> {
    using E=Y2; using XE=X512;
    static SPIN_FORCEINLINE E make(__m256i a,__m256i b) { return {a,b}; }
};
template<> struct PairOf<Y2> {
    using E=Y4; using XE=X1024;
    static SPIN_FORCEINLINE E make(Y2 a,Y2 b) { return {a,b}; }
};
// t[v] = t[v&(v-1)] + p[ctz v]: every nonzero nibble combination, one XOR each.
template<class E,class XE,std::size_t... V>
static SPIN_FORCEINLINE void gray(E* t,const E* p,std::index_sequence<V...>) {
    ((t[V+1]=XE::x2(t[(V+1)&V],p[std::countr_zero(unsigned(V+1))])),...);
}
// Scaled table offsets of the 72 packed bytes: lo = (b&15)<<Shift and
// hi = (b>>4)<<Shift, so that 8*offset is the byte offset of the entry.
template<unsigned Shift>
static SPIN_FORCEINLINE void nibbleOffsets(const std::uint64_t* packed,std::uint8_t* lo,std::uint8_t* hi) {
    const auto* p=reinterpret_cast<const char*>(packed);
    const auto m=_mm256_set1_epi8(0x0f);
    for(unsigned c=0;c<64;c+=32) {
        const auto x=_mm256_loadu_si256(reinterpret_cast<const V*>(p+c));
        _mm256_store_si256(reinterpret_cast<V*>(lo+c),_mm256_slli_epi16(_mm256_and_si256(x,m),Shift));
        _mm256_store_si256(reinterpret_cast<V*>(hi+c),_mm256_srli_epi16(_mm256_andnot_si256(m,x),4-Shift));
    }
    const auto x=_mm_loadl_epi64(reinterpret_cast<const __m128i*>(p+64));
    const auto m128=_mm_set1_epi8(0x0f);
    _mm_storel_epi64(reinterpret_cast<__m128i*>(lo+64),_mm_slli_epi16(_mm_and_si128(x,m128),Shift));
    _mm_storel_epi64(reinterpret_cast<__m128i*>(hi+64),_mm_srli_epi16(_mm_andnot_si128(m128,x),4-Shift));
}
template<class W,class X>
static SPIN_FORCEINLINE void update(W* w,const W* f,const std::uint64_t* packed) {
    using P=PairOf<W>;
    using E=typename P::E;
    using XE=typename P::XE;
    constexpr unsigned S=sizeof(W);
    static_assert(sizeof(E)==2*S && (S==16 || S==32 || S==64));
    alignas(64) E t[3][16];
    alignas(32) std::uint8_t lo[96],hi[96];
    nibbleOffsets<std::countr_zero(S/4)>(packed,lo,hi);
    for(unsigned i=0;i<3;++i) {
        E p[4];
        for(unsigned b=0;b<4;++b) p[b]=P::make(w[8*i+b],8*i+4+b<22?w[8*i+4+b]:X::zero());
        t[i][0]=XE::zero();
        gray<E,XE>(t[i],p,std::make_index_sequence<15>{});
    }
    const auto* base=reinterpret_cast<const char*>(t);
    for(unsigned r=0;r<22;++r) {
        const unsigned o=r/8,j=r%8;
        auto acc=f[r];
        for(unsigned i=0;i<3;++i) {
            const unsigned q=8*(3*o+i)+j;
            const auto* entry=base+sizeof(t[0])*i;
            acc=X::x3(acc,*reinterpret_cast<const W*>(entry+8*std::size_t(lo[q])),
                          *reinterpret_cast<const W*>(entry+S+8*std::size_t(hi[q])));
        }
        w[r]=acc;
    }
}

// Transpose inner: epochs run backwards, input is read sequentially and every
// four-coordinate packet is scattered to its routed scratch offset.
// Streaming stores write each 64-byte packet completely, four in a row.
template<bool Stream> struct Scatter {
    std::byte* base; const std::uint32_t* route;
    SPIN_FORCEINLINE void operator()(unsigned j,__m128i v) const {
        auto* at=reinterpret_cast<__m128i*>(base+route[j>>2]+16*(j&3));
        if constexpr(Stream) _mm_stream_si128(at,v);
        else _mm_store_si128(at,v);
    }
};
template<bool Stream>
static SPIN_NOINLINE void reverseInner(const Block* __restrict input,Block* __restrict scratch,
    const certified22::Plan& plan) {
    using W=__m128i;
    using Scatter=rs8::Scatter<Stream>;
    const auto epochs=plan.n/64;
    auto* base=reinterpret_cast<std::byte*>(scratch);
    const auto* route=plan.routeBytes.data();
    alignas(64) W w[22],f[22],mono[22];
    {
        const auto e=epochs-1;
        const Scatter st{base,route+16*e};
        MonoAccumulator<W,X128> acc;
#pragma GCC unroll 8
        for(unsigned h=0;h<8;++h) {
            W x[8];
            for(unsigned l=0;l<8;++l) { x[l]=load(input+64*e+8*h+l);st(8*h+l,x[l]); }
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<W,X128>(mono,w);
    }
    if(epochs==1) return;
    for(std::size_t e=epochs-2;e>0;--e) {
        const Scatter st{base,route+16*e};
        StepWeights<W,X128> weights;
        weights.prepare(w);
        MonoAccumulator<W,X128> acc;
#pragma GCC unroll 8
        for(unsigned h=0;h<8;++h) {
            W x[8];
            for(unsigned l=0;l<8;++l) x[l]=load(input+64*e+8*h+l);
            evalBlock(weights,h,x,st,8*h);
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<W,X128>(mono,f);
        update<W,X128>(w,f,plan.updates.data()+9*e);
    }
    const Scatter st{base,route};
    StepWeights<W,X128> weights;
    weights.prepare(w);
#pragma GCC unroll 8
    for(unsigned h=0;h<8;++h) {
        W x[8];
        for(unsigned l=0;l<8;++l) x[l]=load(input+8*h+l);
        evalBlock(weights,h,x,st,8*h);
    }
}
// Large scratch is written with streaming stores, as in the AVX-512 reverse
// route above K18; routing completes before the outer reads the scratch.
static void reverseInner(const Block* input,Block* scratch,const certified22::Plan& plan) {
    if(plan.k>262144) { reverseInner<true>(input,scratch,plan); _mm_sfence(); }
    else reverseInner<false>(input,scratch,plan);
}

// Forward inner: lane-major, one vector holds one step coordinate of every
// lane, so output records are already interleaved.
template<unsigned L> struct Lanes;
// Streaming output requires 16-byte (one lane) or 32-byte (wider) alignment.
template<> struct Lanes<1> {
    using W=__m128i; using X=X128;
    static constexpr unsigned prefetch=2;
    static SPIN_FORCEINLINE void gather(const std::byte* at,std::size_t,W* x) {
        for(unsigned k=0;k<4;++k) x[k]=load(at+16*k);
    }
    template<bool Stream> static SPIN_FORCEINLINE void store(Block* at,W v) {
        if constexpr(Stream) _mm_stream_si128(&at->mData,v);
        else _mm_storeu_si128(&at->mData,v);
    }
};
template<> struct Lanes<2> {
    using W=__m256i; using X=X256;
    static constexpr unsigned prefetch=1;
    static SPIN_FORCEINLINE void gather(const std::byte* at,std::size_t stride,W* x) {
        for(unsigned k=0;k<4;++k) x[k]=join(load(at+16*k),load(at+stride+16*k));
    }
    template<bool Stream> static SPIN_FORCEINLINE void store(Block* at,W v) {
        if constexpr(Stream) _mm256_stream_si256(reinterpret_cast<V*>(at),v);
        else _mm256_storeu_si256(reinterpret_cast<V*>(at),v);
    }
};
template<> struct Lanes<4> {
    using W=Y2; using X=X512;
    static constexpr unsigned prefetch=0;
    static SPIN_FORCEINLINE void gather(const std::byte* at,std::size_t stride,W* x) {
        for(unsigned k=0;k<4;++k)
            x[k]={join(load(at+16*k),load(at+stride+16*k)),join(load(at+2*stride+16*k),load(at+3*stride+16*k))};
    }
    template<bool Stream> static SPIN_FORCEINLINE void store(Block* at,W v) {
        Lanes<2>::store<Stream>(at,v.a);Lanes<2>::store<Stream>(at+2,v.b);
    }
};
template<unsigned L,bool Stream> struct LaneStore {
    Block* output;
    template<class W> SPIN_FORCEINLINE void operator()(unsigned j,W v) const {
        Lanes<L>::template store<Stream>(output+L*j,v);
    }
};
// Prefetch (software, measured on the AVX-512 kernels) of the routed packets
// gathered Distance epochs ahead, for every lane plane.
template<unsigned L>
static SPIN_FORCEINLINE void prefetchEpoch(const std::byte* base,std::size_t stride,const std::uint32_t* route) {
    for(unsigned lane=0;lane<L;++lane)
        for(unsigned i=0;i<16;++i)
            _mm_prefetch(reinterpret_cast<const char*>(base+lane*stride+route[i]),_MM_HINT_T0);
}
template<unsigned L,bool Stream>
static SPIN_NOINLINE void forwardInner(const Block* __restrict input,Block* __restrict output,
    const certified22::Plan& plan,std::size_t plane) {
    using W=typename Lanes<L>::W;
    using X=typename Lanes<L>::X;
    using LaneStore=rs8::LaneStore<L,Stream>;
    constexpr unsigned distance=Lanes<L>::prefetch;
    const auto epochs=plan.n/64;
    const auto* route=plan.routeBytes.data();
    const auto* base=reinterpret_cast<const std::byte*>(input);
    const auto stride=plane*sizeof(Block);
    alignas(64) W w[22],f[22],mono[22];
    {
        MonoAccumulator<W,X> acc;
        const LaneStore st{output};
#pragma GCC unroll 8
        for(unsigned h=0;h<8;++h) {
            W x[8];
            Lanes<L>::gather(base+route[2*h],stride,x);Lanes<L>::gather(base+route[2*h+1],stride,x+4);
            for(unsigned l=0;l<8;++l) st(8*h+l,x[l]);
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<W,X>(mono,w);
    }
    if(epochs==1) return;
    for(std::size_t e=1;e+1<epochs;++e) {
        if constexpr(distance>0) if(e+distance<epochs) prefetchEpoch<L>(base,stride,route+16*(e+distance));
        const auto* r=route+16*e;
        StepWeights<W,X> weights;
        weights.prepare(w);
        MonoAccumulator<W,X> acc;
        const LaneStore st{output+64*L*e};
#pragma GCC unroll 8
        for(unsigned h=0;h<8;++h) {
            W x[8];
            Lanes<L>::gather(base+r[2*h],stride,x);Lanes<L>::gather(base+r[2*h+1],stride,x+4);
            evalBlock(weights,h,x,st,8*h);
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<W,X>(mono,f);
        update<W,X>(w,f,plan.forwardUpdates.data()+9*e);
    }
    const auto* r=route+16*(epochs-1);
    StepWeights<W,X> weights;
    weights.prepare(w);
    const LaneStore st{output+64*L*(epochs-1)};
#pragma GCC unroll 8
    for(unsigned h=0;h<8;++h) {
        W x[8];
        Lanes<L>::gather(base+r[2*h],stride,x);Lanes<L>::gather(base+r[2*h+1],stride,x+4);
        evalBlock(weights,h,x,st,8*h);
    }
}
}

// Streaming output needs 16-byte alignment for one lane and 32 for wider
// records (forward_bytes guarantees 16).
static void forwardInnerAvx2(const Block* scratch,Block* output,const certified22::Plan& plan,
    unsigned lanes,bool stream,std::size_t plane) {
    if(stream && (lanes==1 || !(reinterpret_cast<std::uintptr_t>(output)&31U))) {
        if(lanes==1) forwardInner<1,true>(scratch,output,plan,plane);
        else if(lanes==2) forwardInner<2,true>(scratch,output,plan,plane);
        else forwardInner<4,true>(scratch,output,plan,plane);
        _mm_sfence();
    } else {
        if(lanes==1) forwardInner<1,false>(scratch,output,plan,plane);
        else if(lanes==2) forwardInner<2,false>(scratch,output,plan,plane);
        else forwardInner<4,false>(scratch,output,plan,plane);
    }
}
void transposeAvx2(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseInner(input,scratch,*plan.certified22Plan);
    transposeOuter(scratch,output,plan);
}
void forwardAvx2(const Block* input,Block* output,Block* scratch,const Plan& plan,unsigned lanes,bool stream) {
    forwardOuter(input,scratch,plan,lanes);
    forwardInnerAvx2(scratch,output,*plan.certified22Plan,lanes,stream,plan.scratchBlocks());
}
// Stage entry points for differential checks and stage timing.
void reverseRouteAvx2(const Block* input,Block* scratch,const Plan& plan) {
    reverseInner(input,scratch,*plan.certified22Plan);
}
void outerAvx2(const Block* scratch,Block* output,const Plan& plan) { transposeOuter(scratch,output,plan); }
void outerForwardAvx2(const Block* input,Block* scratch,const Plan& plan,unsigned lanes) {
    forwardOuter(input,scratch,plan,lanes);
}
void forwardInnerAvx2(const Block* scratch,Block* output,const Plan& plan,unsigned lanes) {
    forwardInnerAvx2(scratch,output,*plan.certified22Plan,lanes,false,plan.scratchBlocks());
}
}
