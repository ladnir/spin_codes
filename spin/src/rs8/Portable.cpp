// Portable Turbo (PacketRs8T64) kernels: plain C++ on the common 128-bit
// block, for any 64-bit CPU. Same linear map as the oracle and the ISA
// kernels (spin_rs8_test portable checks it against the oracle).
//
// Inner: the fixed s22 step XOR network of certified22/StepCircuit.h in the
// fast basis, with the per-epoch GL22 update applied by four-Russians tables
// (as in rs8/Avx2.cpp). Forward records are lane-major: one vector of L
// blocks per step coordinate.
//
// Outer, bitsliced: multiplication of the 32 coordinates of a symbol by its
// GF(2^32) scalar is a 32x32 binary matrix, applied by four-Russians tables
// over nibbles of the input coordinates; its rows (transpose) or columns
// (forward) are linear in the scalar and come from nibble-indexed tables.
// The RS[8,4] parity uses the factored circuit of rs8/Fast.cpp, its nibble
// maps applied as 4x4 binary matrices to four coordinates.
#include "Rs8.h"
#include "Field.h"
#include "../certified22/productionMap.h"
#include "../certified22/LaneMajor22.h"
#include "../certified22/StepCircuit.h"
#include <bit>
#include <cstddef>
#include <utility>

namespace spin::detail::rs8 {
namespace {
using namespace certified22::step;
using B=Block;

// ------------------------------------------------------------ vectors --
struct XB {
    static SPIN_FORCEINLINE B x2(B a,B b) { return a^b; }
    static SPIN_FORCEINLINE B x3(B a,B b,B c) { return a^b^c; }
    static SPIN_FORCEINLINE B zero() { return B::zero(); }
};
// One step coordinate of L record lanes.
template<unsigned L> struct alignas(16) R { B v[L]; };
template<unsigned L> struct XR {
    static SPIN_FORCEINLINE R<L> x2(const R<L>& a,const R<L>& b) {
        R<L> r; for(unsigned i=0;i<L;++i) r.v[i]=a.v[i]^b.v[i]; return r;
    }
    static SPIN_FORCEINLINE R<L> x3(const R<L>& a,const R<L>& b,const R<L>& c) {
        R<L> r; for(unsigned i=0;i<L;++i) r.v[i]=a.v[i]^b.v[i]^c.v[i]; return r;
    }
    static SPIN_FORCEINLINE R<L> zero() { R<L> r; for(auto& x:r.v) x=B::zero(); return r; }
};
template<unsigned L> struct LaneTypes { using W=R<L>; using X=XR<L>; };
template<> struct LaneTypes<1> { using W=B; using X=XB; };

static SPIN_FORCEINLINE B load(const void* p) { B b; std::memcpy(&b,p,16); return b; }
static SPIN_FORCEINLINE void store(void* p,B b) { std::memcpy(p,&b,16); }

// t[v] = t[v&(v-1)] + p[ctz v] for every nonzero nibble value v < N.
template<class W,class X,std::size_t... V>
static SPIN_FORCEINLINE void gray(W* t,const W* p,std::index_sequence<V...>) {
    ((t[V+1]=X::x2(t[(V+1)&V],p[std::countr_zero(unsigned(V+1))])),...);
}

// ------------------------------------------------------------- inner --
// w = U w + f for the packed 3x3-block GL22 table of one epoch (qword 3o+i,
// byte j is byte i of row 8o+j; bit c of row r is U[r][c]).
template<class W,class X>
static SPIN_FORCEINLINE void update(W* w,const W* f,const std::uint64_t* packed) {
    W t[6][16];
    for(unsigned k=0;k<6;++k) {
        t[k][0]=X::zero();
        if(k<5) gray<W,X>(t[k],w+4*k,std::make_index_sequence<15>{});
        else gray<W,X>(t[k],w+4*k,std::make_index_sequence<3>{});
    }
    const auto* bytes=reinterpret_cast<const std::uint8_t*>(packed);
    for(unsigned r=0;r<22;++r) {
        const auto* row=bytes+24*(r/8)+(r%8);
        const unsigned b0=row[0],b1=row[8],b2=row[16];
        w[r]=X::x3(X::x3(f[r],t[0][b0&15],t[1][b0>>4]),X::x2(t[2][b1&15],t[3][b1>>4]),
                   X::x2(t[4][b2&15],t[5][b2>>4]));
    }
}
// Transpose inner: epochs backwards, sequential input, packets scattered.
struct Scatter {
    std::byte* base; const std::uint32_t* route;
    SPIN_FORCEINLINE void operator()(unsigned j,B v) const { store(base+route[j>>2]+16*(j&3),v); }
};
static SPIN_NOINLINE void reverseInner(const B* __restrict input,B* __restrict scratch,
    const certified22::Plan& plan) {
    const auto epochs=plan.n/64;
    auto* base=reinterpret_cast<std::byte*>(scratch);
    const auto* route=plan.routeBytes.data();
    B w[22],f[22],mono[22];
    {
        const auto e=epochs-1;
        const Scatter st{base,route+16*e};
        MonoAccumulator<B,XB> acc;
        for(unsigned h=0;h<8;++h) {
            B x[8];
            for(unsigned l=0;l<8;++l) { x[l]=input[64*e+8*h+l];st(8*h+l,x[l]); }
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<B,XB>(mono,w);
    }
    if(epochs==1) return;
    for(std::size_t e=epochs-2;e>0;--e) {
        const Scatter st{base,route+16*e};
        StepWeights<B,XB> weights;
        weights.prepare(w);
        MonoAccumulator<B,XB> acc;
        for(unsigned h=0;h<8;++h) {
            B x[8];
            for(unsigned l=0;l<8;++l) x[l]=input[64*e+8*h+l];
            evalBlock(weights,h,x,st,8*h);
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<B,XB>(mono,f);
        update<B,XB>(w,f,plan.updates.data()+9*e);
    }
    const Scatter st{base,route};
    StepWeights<B,XB> weights;
    weights.prepare(w);
    for(unsigned h=0;h<8;++h) {
        B x[8];
        for(unsigned l=0;l<8;++l) x[l]=input[8*h+l];
        evalBlock(weights,h,x,st,8*h);
    }
}
// Forward inner, lane-major: gathers routed packets from L scratch planes.
template<unsigned L>
static SPIN_FORCEINLINE void gather(const std::byte* at,std::size_t stride,typename LaneTypes<L>::W* x) {
    for(unsigned k=0;k<4;++k) {
        if constexpr(L==1) x[k]=load(at+16*k);
        else for(unsigned lane=0;lane<L;++lane) x[k].v[lane]=load(at+lane*stride+16*k);
    }
}
template<unsigned L> struct LaneStore {
    B* output;
    template<class W> SPIN_FORCEINLINE void operator()(unsigned j,const W& v) const {
        if constexpr(L==1) output[j]=v;
        else for(unsigned lane=0;lane<L;++lane) output[L*j+lane]=v.v[lane];
    }
};
template<unsigned L>
static SPIN_NOINLINE void forwardInner(const B* __restrict input,B* __restrict output,
    const certified22::Plan& plan,std::size_t plane) {
    using W=typename LaneTypes<L>::W;
    using X=typename LaneTypes<L>::X;
    const auto epochs=plan.n/64;
    const auto* route=plan.routeBytes.data();
    const auto* base=reinterpret_cast<const std::byte*>(input);
    const auto stride=plane*sizeof(B);
    W w[22],f[22],mono[22];
    {
        MonoAccumulator<W,X> acc;
        const LaneStore<L> st{output};
        for(unsigned h=0;h<8;++h) {
            W x[8];
            gather<L>(base+route[2*h],stride,x);gather<L>(base+route[2*h+1],stride,x+4);
            for(unsigned l=0;l<8;++l) st(8*h+l,x[l]);
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<W,X>(mono,w);
    }
    if(epochs==1) return;
    for(std::size_t e=1;e+1<epochs;++e) {
        const auto* r=route+16*e;
        StepWeights<W,X> weights;
        weights.prepare(w);
        MonoAccumulator<W,X> acc;
        const LaneStore<L> st{output+64*L*e};
        for(unsigned h=0;h<8;++h) {
            W x[8];
            gather<L>(base+r[2*h],stride,x);gather<L>(base+r[2*h+1],stride,x+4);
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
    const LaneStore<L> st{output+64*L*(epochs-1)};
    for(unsigned h=0;h<8;++h) {
        W x[8];
        gather<L>(base+r[2*h],stride,x);gather<L>(base+r[2*h+1],stride,x+4);
        evalBlock(weights,h,x,st,8*h);
    }
}

// ------------------------------------------------------------- outer --
// Rows (Mult(c)[r] as a mask over input coordinates) and columns (c*e_col)
// of multiplication by c, as XORs of nibble-indexed tables: both are linear
// in c. Built once from the literal field arithmetic.
struct ScalarTables { std::uint32_t rows[8][16][32],columns[8][16][32]; };
static const ScalarTables& scalarTables() {
    static const ScalarTables tables=[] {
        ScalarTables t{};
        for(unsigned p=0;p<8;++p)
            for(unsigned v=1;v<16;++v)
                for(unsigned col=0;col<32;++col) {
                    const auto image=field::multiply32(std::uint32_t(v)<<(4*p),std::uint32_t{1}<<col);
                    t.columns[p][v][col]=image;
                    for(unsigned r=0;r<32;++r) t.rows[p][v][r]|=((image>>r)&1U)<<col;
                }
        return t;
    }();
    return tables;
}
// The 32 rows (or columns) of multiplication by c.
static SPIN_FORCEINLINE void matrixOf(const std::uint32_t (&table)[8][16][32],std::uint32_t c,std::uint32_t* out) {
    for(unsigned i=0;i<32;++i) out[i]=0;
    for(unsigned p=0;p<8;++p) {
        const auto* t=table[p][(c>>(4*p))&15U];
        for(unsigned i=0;i<32;++i) out[i]^=t[i];
    }
}
// out[i] = XOR of in[c] over the bits c of masks[i], by four-Russians
// tables over nibbles of in (32 inputs and outputs).
static SPIN_FORCEINLINE void apply32(const B* in,const std::uint32_t* masks,B* out) {
    B t[8][16];
    for(unsigned k=0;k<8;++k) {
        t[k][0]=B::zero();
        gray<B,XB>(t[k],in+4*k,std::make_index_sequence<15>{});
    }
    for(unsigned i=0;i<32;++i) {
        const auto m=masks[i];
        out[i]=t[0][m&15U]^t[1][(m>>4)&15U]^t[2][(m>>8)&15U]^t[3][(m>>12)&15U]
              ^t[4][(m>>16)&15U]^t[5][(m>>20)&15U]^t[6][(m>>24)&15U]^t[7][m>>28];
    }
}

// The parity network's nibble maps as 4x4 binary matrices: column b is the
// image of e_b (rs8/Fast.cpp m6/m8, their transposes in rs8/ForwardFast.cpp).
constexpr std::uint8_t affine(std::uint64_t m,unsigned x) {
    unsigned result=0;
    for(unsigned i=0;i<8;++i) result|=(std::popcount(unsigned((m>>(8*(7-i)))&x)&255U)&1U)<<i;
    return std::uint8_t(result);
}
constexpr std::uint64_t transposeAffine(std::uint64_t m) {
    std::uint64_t result=0;
    for(unsigned i=0;i<8;++i) for(unsigned j=0;j<8;++j)
        result|=((m>>(8*(7-i)+j))&1ULL)<<(8*(7-j)+i);
    return result;
}
using Nibble=std::array<std::uint8_t,4>;
constexpr Nibble nibbleMatrix(std::uint64_t m) {
    Nibble columns{};
    for(unsigned b=0;b<4;++b) columns[b]=std::uint8_t(affine(m,1U<<b)&15U);
    return columns;
}
constexpr std::uint64_t m6=0x060c0b0560c0b050ULL,m8=0x0803060c803060c0ULL;
constexpr Nibble a6=nibbleMatrix(m6),a8=nibbleMatrix(m8);
constexpr Nibble a6t=nibbleMatrix(transposeAffine(m6)),a8t=nibbleMatrix(transposeAffine(m8));
struct Quad { B v[4]; };
static SPIN_FORCEINLINE Quad operator^(const Quad& a,const Quad& b) {
    return {{a.v[0]^b.v[0],a.v[1]^b.v[1],a.v[2]^b.v[2],a.v[3]^b.v[3]}};
}
template<const Nibble& M>
static SPIN_FORCEINLINE Quad map(const Quad& x) {
    Quad y{{B::zero(),B::zero(),B::zero(),B::zero()}};
    for(unsigned b=0;b<4;++b)
        for(unsigned t=0;t<4;++t)
            if((M[b]>>t)&1U) y.v[t]^=x.v[b];
    return y;
}

// Transpose outer of one group: scratch coordinates 32s+c to message blocks
// 16 row + 4q + t.
static SPIN_FORCEINLINE void groupTranspose(const B* in,B* out,const std::uint32_t* scalars) {
    const auto& tables=scalarTables();
    B mixed[256];
    for(unsigned s=0;s<8;++s) {
        std::uint32_t rows[32];
        matrixOf(tables.rows,scalars[s],rows);
        apply32(in+32*s,rows,mixed+32*s);
    }
    const auto nibble=[&](unsigned s,unsigned row) {
        const auto* p=mixed+32*s+4*row; return Quad{{p[0],p[1],p[2],p[3]}};
    };
    for(unsigned row=0;row<8;++row) {
        auto x0=nibble(4,row),x1=nibble(5,row),x2=nibble(6,row),x3=nibble(7,row);
        x0=x0^x1;x2=x2^x3;x0=x0^x2;x1=x1^x3;
        const auto u6=map<a6>(x0),u8=map<a8>(x0);
        auto y0=x0,y1=x1^u8,y2=x2^u6,y3=x3^u6^x0^map<a6>(x1)^map<a8>(x2);
        y0=y0^y1;y2=y2^y3;y0=y0^y2;y1=y1^y3;
        const Quad y[4]{y0,y1,y2,y3};
        for(unsigned q=0;q<4;++q) {
            const auto v=nibble(q,row)^y[q];
            for(unsigned t=0;t<4;++t) out[16*row+4*q+t]=v.v[t];
        }
    }
}
// Forward outer (exact adjoint) of one group and one record lane.
static SPIN_FORCEINLINE void groupForward(const B* message,unsigned lanes,unsigned lane,B* out,
    const std::uint32_t* scalars) {
    const auto& tables=scalarTables();
    B coded[256];
    const auto nibble=[&](unsigned q,unsigned row) {
        Quad v;
        for(unsigned t=0;t<4;++t) v.v[t]=message[(16*row+4*q+t)*lanes+lane];
        return v;
    };
    for(unsigned row=0;row<8;++row) {
        Quad y[4];
        for(unsigned q=0;q<4;++q) {
            y[q]=nibble(q,row);
            for(unsigned t=0;t<4;++t) coded[32*q+4*row+t]=y[q].v[t];
        }
        auto y0=y[0],y1=y[1],y2=y[2],y3=y[3];
        y3=y3^y1;y2=y2^y0;y3=y3^y2;y1=y1^y0;
        auto x0=y0^y3^map<a8t>(y1)^map<a6t>(y2)^map<a6t>(y3);
        auto x1=y1^map<a6t>(y3),x2=y2^map<a8t>(y3),x3=y3;
        x3=x3^x1;x2=x2^x0;x3=x3^x2;x1=x1^x0;
        const Quad x[4]{x0,x1,x2,x3};
        for(unsigned q=0;q<4;++q)
            for(unsigned t=0;t<4;++t) coded[32*(4+q)+4*row+t]=x[q].v[t];
    }
    for(unsigned s=0;s<8;++s) {
        std::uint32_t columns[32];
        matrixOf(tables.columns,scalars[s],columns);
        apply32(coded+32*s,columns,out+32*s);
    }
}
}

void transposePortable(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseInner(input,scratch,*plan.certified22Plan);
    for(std::size_t g=0;g<plan.groups;++g)
        groupTranspose(scratch+groupStride*g,output+128*g,plan.scalars.data()+8*g);
}
void forwardPortable(const Block* input,Block* output,Block* scratch,const Plan& plan,unsigned lanes) {
    const auto plane=plan.scratchBlocks();
    for(std::size_t g=0;g<plan.groups;++g)
        for(unsigned lane=0;lane<lanes;++lane)
            groupForward(input+128*lanes*g,lanes,lane,scratch+lane*plane+groupStride*g,plan.scalars.data()+8*g);
    if(lanes==1) forwardInner<1>(scratch,output,*plan.certified22Plan,plane);
    else if(lanes==2) forwardInner<2>(scratch,output,*plan.certified22Plan,plane);
    else forwardInner<4>(scratch,output,*plan.certified22Plan,plane);
}
}
