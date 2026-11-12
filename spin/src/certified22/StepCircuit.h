// The fixed s22 step circuit in the fast basis w = basis*s, generic over the
// vector type V and its XOR helpers X (x2, x3). Pure XOR networks: shared by
// the AVX-512 lane-major forward (RouteFast.cpp) and the AVX2 kernels
// (rs8/Avx2.cpp). Every instantiating TU must use TU-local V/X types so no
// two ISA builds of one instantiation can be merged by the linker.
#pragma once
#include "../kernels/Block.h"
#include "productionMap.h"

namespace spin::detail::certified22::step {
// Mask-indexed access to the fast-circuit state coordinates.
constexpr unsigned stateIndex(unsigned mask) {
    for(unsigned i=0; i<22; ++i) if(map::masks[i]==mask) return i;
    return 22;
}
template<unsigned Mask> static constexpr unsigned SI=[] {
    static_assert(stateIndex(Mask)<22,"monomial outside the s22 state");
    return stateIndex(Mask);
}();
// One step coordinate j = 8h + l splits into three high bits h (masks bits
// 3..5) and three low bits l (bits 0..2). Per step, the evaluation weights are
// base A(h) = sum over monomials inside h and q_b(h) = w_b + sum_{a in h} w_ab,
// so each 8-coordinate block needs a handful of live values; mono is its exact
// transpose, accumulated per block. Checked against the definitions in the
// generator's docstring and by the forward/transpose tests.
template<class V,class X>
struct StepWeights {
    V A[8],q[3][8],p01,p02,p12,p0212;
    SPIN_FORCEINLINE void prepare(const V* w) {
        A[0]=w[SI<0>];
        A[1]=X::x2(w[SI<0>],w[SI<8>]);A[2]=X::x2(w[SI<0>],w[SI<16>]);A[4]=X::x2(w[SI<0>],w[SI<32>]);
        A[3]=X::x3(A[1],w[SI<16>],w[SI<24>]);A[5]=X::x3(A[1],w[SI<32>],w[SI<40>]);
        A[6]=X::x3(A[2],w[SI<32>],w[SI<48>]);A[7]=X::x2(X::x3(A[3],w[SI<32>],w[SI<40>]),w[SI<48>]);
        cross<0>(w);cross<1>(w);cross<2>(w);
        p01=w[SI<3>];p02=w[SI<5>];p12=w[SI<6>];p0212=X::x2(p02,p12);
    }
    template<unsigned B> SPIN_FORCEINLINE void cross(const V* w) {
        auto* r=q[B];
        r[0]=w[SI<(1u<<B)>];
        r[1]=X::x2(r[0],w[SI<(8u|(1u<<B))>]);r[2]=X::x2(r[0],w[SI<(16u|(1u<<B))>]);r[4]=X::x2(r[0],w[SI<(32u|(1u<<B))>]);
        r[3]=X::x2(r[1],w[SI<(16u|(1u<<B))>]);r[5]=X::x2(r[1],w[SI<(32u|(1u<<B))>]);
        r[6]=X::x2(r[2],w[SI<(32u|(1u<<B))>]);r[7]=X::x2(r[3],w[SI<(32u|(1u<<B))>]);
    }
};
template<class V,class X>
struct MonoAccumulator {
    V B[8],Q[3][8],p01,p02,p12,p0212;
    // Adjoint of one block's evaluation for raw inputs x0..x7 at block h.
    SPIN_FORCEINLINE void block(unsigned h,const V* x) {
        const auto a37=X::x2(x[3],x[7]),t2=X::x2(x[2],x[6]);
        const auto t1=X::x3(x[1],x[5],a37);
        B[h]=X::x2(X::x3(x[0],x[4],t2),t1);
        Q[0][h]=t1;Q[1][h]=X::x2(a37,t2);Q[2][h]=X::x2(X::x3(x[4],x[5],x[6]),x[7]);
        if(h==0) { p01=a37;p02=x[5];p12=x[6];p0212=x[7]; }
        else { p01=X::x2(p01,a37);p02=X::x2(p02,x[5]);p12=X::x2(p12,x[6]);p0212=X::x2(p0212,x[7]); }
    }
    // mono[i] for every state monomial: superset sums over the high bits.
    SPIN_FORCEINLINE void finish(V* mono) const {
        const V* b=B;
        const auto s1=X::x3(b[1],b[3],b[5]),s2=X::x3(b[2],b[3],b[6]),s4=X::x3(b[4],b[5],b[6]);
        mono[SI<8>]=X::x2(s1,b[7]);mono[SI<16>]=X::x2(s2,b[7]);mono[SI<32>]=X::x2(s4,b[7]);
        mono[SI<24>]=X::x2(b[3],b[7]);mono[SI<40>]=X::x2(b[5],b[7]);mono[SI<48>]=X::x2(b[6],b[7]);
        mono[SI<0>]=X::x3(X::x3(b[0],b[1],b[2]),X::x3(b[3],b[4],b[5]),X::x2(b[6],b[7]));
        lowCross<0>(mono);lowCross<1>(mono);lowCross<2>(mono);
        mono[SI<3>]=p01;mono[SI<5>]=X::x2(p02,p0212);mono[SI<6>]=X::x2(p12,p0212);
    }
    template<unsigned Bit> SPIN_FORCEINLINE void lowCross(V* mono) const {
        const V* r=Q[Bit];
        const auto r13=X::x2(r[1],r[3]),r57=X::x2(r[5],r[7]),r26=X::x2(r[2],r[6]);
        mono[SI<(8u|(1u<<Bit))>]=X::x2(r13,r57);
        mono[SI<(16u|(1u<<Bit))>]=X::x3(r26,r[3],r[7]);
        mono[SI<(32u|(1u<<Bit))>]=X::x2(X::x3(r[4],r[5],r[6]),r[7]);
        mono[SI<(1u<<Bit)>]=X::x3(X::x3(r[0],r13,r57),r26,r[4]);
    }
};
template<class V,class X,class Store>
static SPIN_FORCEINLINE void evalBlock(const StepWeights<V,X>& s,unsigned h,const V* x,Store& st,unsigned j) {
    const auto base=s.A[h],q0=s.q[0][h],q1=s.q[1][h],q2=s.q[2][h];
    const auto t1=X::x2(base,q0),t2=X::x2(base,q1);
    const auto v3=X::x3(t1,q1,s.p01);
    st(j+0,X::x2(x[0],base));st(j+1,X::x2(x[1],t1));st(j+2,X::x2(x[2],t2));st(j+3,X::x2(x[3],v3));
    st(j+4,X::x3(x[4],base,q2));st(j+5,X::x3(x[5],t1,X::x2(q2,s.p02)));
    st(j+6,X::x3(x[6],t2,X::x2(q2,s.p12)));st(j+7,X::x3(x[7],v3,X::x2(q2,s.p0212)));
}
}
