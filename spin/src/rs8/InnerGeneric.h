// The s22 inner (transpose and forward) over a 128-bit vector type, shared by
// the Portable and NEON Turbo kernels. Ops supplies V (one 128-bit block),
// load/store from block memory, and x2/x3/zero; every instantiating TU must
// use a TU-local Ops so that no two ISA builds of one instantiation merge.
//
// The fixed step circuit is the XOR network of certified22/StepCircuit.h in
// the fast basis; the per-epoch GL22 update w' = U w + f applies the packed
// 3x3-block table (updates or forwardUpdates) by four-Russians tables.
// Forward records are lane-major: one vector of L blocks per coordinate.
#pragma once
#include "Rs8.h"
#include "../certified22/LaneMajor22.h"
#include "../certified22/StepCircuit.h"
#include <bit>
#include <cstddef>
#include <utility>

namespace spin::detail::rs8::inner {
using namespace certified22::step;

template<class Ops> struct X1 {
    using V=typename Ops::V;
    static SPIN_FORCEINLINE V x2(V a,V b) { return Ops::x2(a,b); }
    static SPIN_FORCEINLINE V x3(V a,V b,V c) { return Ops::x3(a,b,c); }
    static SPIN_FORCEINLINE V zero() { return Ops::zero(); }
};
template<class V,unsigned L> struct R { V v[L]; };
template<class Ops,unsigned L> struct XR {
    using W=R<typename Ops::V,L>;
    static SPIN_FORCEINLINE W x2(const W& a,const W& b) {
        W r; for(unsigned i=0;i<L;++i) r.v[i]=Ops::x2(a.v[i],b.v[i]); return r;
    }
    static SPIN_FORCEINLINE W x3(const W& a,const W& b,const W& c) {
        W r; for(unsigned i=0;i<L;++i) r.v[i]=Ops::x3(a.v[i],b.v[i],c.v[i]); return r;
    }
    static SPIN_FORCEINLINE W zero() { W r; for(auto& x:r.v) x=Ops::zero(); return r; }
};
template<class Ops,unsigned L> struct Lane { using W=R<typename Ops::V,L>; using X=XR<Ops,L>; };
template<class Ops> struct Lane<Ops,1> { using W=typename Ops::V; using X=X1<Ops>; };

// t[v] = t[v&(v-1)] + p[ctz v] for every nonzero nibble value v.
template<class W,class X,std::size_t... I>
SPIN_FORCEINLINE void gray(W* t,const W* p,std::index_sequence<I...>) {
    ((t[I+1]=X::x2(t[(I+1)&I],p[std::countr_zero(unsigned(I+1))])),...);
}
// w = U w + f (qword 3o+i, byte j is byte i of row 8o+j; bit c of row r is U[r][c]).
template<class W,class X>
SPIN_FORCEINLINE void update(W* w,const W* f,const std::uint64_t* packed) {
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

// Transpose inner: epochs backwards, sequential input, packets scattered to
// their routed scratch offsets.
template<class Ops> struct Scatter {
    std::byte* base; const std::uint32_t* route;
    SPIN_FORCEINLINE void operator()(unsigned j,typename Ops::V v) const { Ops::store(base+route[j>>2]+16*(j&3),v); }
};
template<class Ops>
SPIN_NOINLINE void reverse(const Block* __restrict input,Block* __restrict scratch,const certified22::Plan& plan) {
    using V=typename Ops::V;
    using X=X1<Ops>;
    const auto epochs=plan.n/64;
    auto* base=reinterpret_cast<std::byte*>(scratch);
    const auto* route=plan.routeBytes.data();
    V w[22],f[22],mono[22];
    {
        const auto e=epochs-1;
        const Scatter<Ops> st{base,route+16*e};
        MonoAccumulator<V,X> acc;
        for(unsigned h=0;h<8;++h) {
            V x[8];
            for(unsigned l=0;l<8;++l) { x[l]=Ops::load(input+64*e+8*h+l);st(8*h+l,x[l]); }
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<V,X>(mono,w);
    }
    if(epochs==1) return;
    for(std::size_t e=epochs-2;e>0;--e) {
        const Scatter<Ops> st{base,route+16*e};
        StepWeights<V,X> weights;
        weights.prepare(w);
        MonoAccumulator<V,X> acc;
        for(unsigned h=0;h<8;++h) {
            V x[8];
            for(unsigned l=0;l<8;++l) x[l]=Ops::load(input+64*e+8*h+l);
            evalBlock(weights,h,x,st,8*h);
            acc.block(h,x);
        }
        acc.finish(mono);
        certified22::lanemajor::feedbackH<V,X>(mono,f);
        update<V,X>(w,f,plan.updates.data()+9*e);
    }
    const Scatter<Ops> st{base,route};
    StepWeights<V,X> weights;
    weights.prepare(w);
    for(unsigned h=0;h<8;++h) {
        V x[8];
        for(unsigned l=0;l<8;++l) x[l]=Ops::load(input+8*h+l);
        evalBlock(weights,h,x,st,8*h);
    }
}

// Forward inner, lane-major: gathers routed packets from L scratch planes.
template<class Ops,unsigned L>
SPIN_FORCEINLINE void gather(const std::byte* at,std::size_t stride,typename Lane<Ops,L>::W* x) {
    for(unsigned k=0;k<4;++k) {
        if constexpr(L==1) x[k]=Ops::load(at+16*k);
        else for(unsigned lane=0;lane<L;++lane) x[k].v[lane]=Ops::load(at+lane*stride+16*k);
    }
}
template<class Ops,unsigned L> struct LaneStore {
    Block* output;
    template<class W> SPIN_FORCEINLINE void operator()(unsigned j,const W& v) const {
        if constexpr(L==1) Ops::store(output+j,v);
        else for(unsigned lane=0;lane<L;++lane) Ops::store(output+L*j+lane,v.v[lane]);
    }
};
template<class Ops,unsigned L>
SPIN_NOINLINE void forward(const Block* __restrict input,Block* __restrict output,
    const certified22::Plan& plan,std::size_t plane) {
    using W=typename Lane<Ops,L>::W;
    using X=typename Lane<Ops,L>::X;
    const auto epochs=plan.n/64;
    const auto* route=plan.routeBytes.data();
    const auto* base=reinterpret_cast<const std::byte*>(input);
    const auto stride=plane*sizeof(Block);
    W w[22],f[22],mono[22];
    {
        MonoAccumulator<W,X> acc;
        const LaneStore<Ops,L> st{output};
        for(unsigned h=0;h<8;++h) {
            W x[8];
            gather<Ops,L>(base+route[2*h],stride,x);gather<Ops,L>(base+route[2*h+1],stride,x+4);
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
        const LaneStore<Ops,L> st{output+64*L*e};
        for(unsigned h=0;h<8;++h) {
            W x[8];
            gather<Ops,L>(base+r[2*h],stride,x);gather<Ops,L>(base+r[2*h+1],stride,x+4);
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
    const LaneStore<Ops,L> st{output+64*L*(epochs-1)};
    for(unsigned h=0;h<8;++h) {
        W x[8];
        gather<Ops,L>(base+r[2*h],stride,x);gather<Ops,L>(base+r[2*h+1],stride,x+4);
        evalBlock(weights,h,x,st,8*h);
    }
}
template<class Ops>
void forwardLanes(const Block* scratch,Block* output,const certified22::Plan& plan,unsigned lanes,std::size_t plane) {
    if(lanes==1) forward<Ops,1>(scratch,output,plan,plane);
    else if(lanes==2) forward<Ops,2>(scratch,output,plan,plane);
    else forward<Ops,4>(scratch,output,plan,plane);
}
}
