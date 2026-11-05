#pragma once
#include "SparseByteField128Simd.h"

// Isolated arithmetic alternative. No setup sampler, key expansion, basis
// conversion, or change to a measured encoder is introduced here.
namespace spin::research::char2::field128_lift::sparse_byte::deferred {
#if defined(_MSC_VER)
#define SPARSE128_DEFER_INLINE __forceinline
#else
#define SPARSE128_DEFER_INLINE inline __attribute__((always_inline))
#endif

// Four independent polynomial products over GF2[X], with AES coefficients
// spaced by16 bits. Vectors hold Y degrees0..7,8..15,16..23,24..31.
// Each slot has degree<=14 in X. XORing products cannot increase that degree.
struct Product { __m512i p0,p1,p2,p3; };

SPARSE128_DEFER_INLINE Product multiplyRaw4(__m512i a,__m512i b) noexcept {
    const auto zero=_mm512_setzero_si512();
    const auto a0=_mm512_unpacklo_epi8(a,zero),a1=_mm512_unpackhi_epi8(a,zero);
    const auto b0=_mm512_unpacklo_epi8(b,zero),b1=_mm512_unpackhi_epi8(b,zero);
    // Cross first: absorb each diagonal product before computing the next.
    // This avoids keeping three separate two-vector products plus the cross.
    auto cross=carryless128(_mm512_xor_si512(a0,a1),_mm512_xor_si512(b0,b1));
    const auto low=carryless128(a0,b0);
    cross.low=_mm512_xor_si512(cross.low,low.low);
    cross.high=_mm512_xor_si512(cross.high,low.high);
    const auto high=carryless128(a1,b1);
    cross.low=_mm512_xor_si512(cross.low,high.low);
    cross.high=_mm512_xor_si512(cross.high,high.high);
    return {low.low,_mm512_xor_si512(low.high,cross.low),
            _mm512_xor_si512(high.low,cross.high),high.high};
}

SPARSE128_DEFER_INLINE Product sum(Product a,Product b) noexcept {
    return {_mm512_xor_si512(a.p0,b.p0),_mm512_xor_si512(a.p1,b.p1),
            _mm512_xor_si512(a.p2,b.p2),_mm512_xor_si512(a.p3,b.p3)};
}

SPARSE128_DEFER_INLINE __m512i reduce(Product a) noexcept {
    return reduceExtension(reduceAndPack(a.p0,a.p1),reduceAndPack(a.p2,a.p3));
}

SPARSE128_DEFER_INLINE __m512i reduceSum3(Product a,Product b,Product c) noexcept {
    const auto low=reduceAndPack(_mm512_ternarylogic_epi64(a.p0,b.p0,c.p0,0x96),
                                _mm512_ternarylogic_epi64(a.p1,b.p1,c.p1,0x96));
    const auto high=reduceAndPack(_mm512_ternarylogic_epi64(a.p2,b.p2,c.p2,0x96),
                                 _mm512_ternarylogic_epi64(a.p3,b.p3,c.p3,0x96));
    return reduceExtension(low,high);
}

// Same nine-product Toeplitz4 circuit and seven compact keys as the frozen
// LiteralEncoder::toeplitz4. Layout: seven vectors, four independent field
// coefficients each: A0,A1,A2,V0,V2,W0,W1. V1=V0+A1+A2, W2=W0+A1+A2.
// The represented matrix is T[i][j]=t[3+i-j], with compact keys obtained from
// Plan::nineKeys(t) at positions0,1,2,3,5,6,7. Any t (singular or not) is valid
// arithmetic input. Sampling/rejection and stored key bytes remain unchanged.
//
// AES reduction and sparse extension reduction are GF2-linear. Thus each
// output can XOR its raw products before one reduction. Static arithmetic per
// four records:9*9=81 VPCLMUL,4*6=24 GFNI (old81/54), no key/data conversions.
// U0/U1 are the only raw sums retained across the V/W phases. Each shared
// V/W product dies at its phase boundary; outputs replace dead input vectors.
// There is no array of nine raw products. Compiler spill behavior is unmeasured.
// Stores are last: x is in-place; key storage may overlap x. Unaligned key
// bytes are accepted. x itself addresses four properly aligned __m512i values.
SPARSE128_DEFER_INLINE void toeplitz4(__m512i* x,const void* packedKeys) noexcept {
    const auto* p=static_cast<const unsigned char*>(packedKeys);
    const auto x0=x[0],x1=x[1],x2=x[2],x3=x[3];
    Product u0,u1;
    __m512i keyXor;
    {
        const auto a0=_mm512_loadu_si512(p),a1=_mm512_loadu_si512(p+64),a2=_mm512_loadu_si512(p+128);
        keyXor=_mm512_xor_si512(a1,a2);
        const auto low=_mm512_xor_si512(x0,x2),high=_mm512_xor_si512(x1,x3);
        const auto shared=multiplyRaw4(_mm512_xor_si512(low,high),a0);
        u0=sum(shared,multiplyRaw4(high,a1));
        u1=sum(shared,multiplyRaw4(low,a2));
    }
    __m512i out2,out3;
    {
        const auto v0=_mm512_loadu_si512(p+192);
        const auto shared=multiplyRaw4(_mm512_xor_si512(x0,x1),v0);
        out2=reduceSum3(u0,shared,multiplyRaw4(x1,_mm512_xor_si512(v0,keyXor)));
        out3=reduceSum3(u1,shared,multiplyRaw4(x0,_mm512_loadu_si512(p+256)));
    }
    __m512i out0,out1;
    {
        const auto w0=_mm512_loadu_si512(p+320);
        const auto shared=multiplyRaw4(_mm512_xor_si512(x2,x3),w0);
        out0=reduceSum3(u0,shared,multiplyRaw4(x3,_mm512_loadu_si512(p+384)));
        out1=reduceSum3(u1,shared,multiplyRaw4(x2,_mm512_xor_si512(w0,keyXor)));
    }
    x[0]=out0;x[1]=out1;x[2]=out2;x[3]=out3;
}

// Byte-addressable boundary:256 input/output bytes and448 packed-key bytes.
// All reads precede output stores; arbitrary input/output/key overlaps work.
SPARSE128_DEFER_INLINE void toeplitz4Records(const void* input,const void* packedKeys,void* output) noexcept {
    const auto* a=static_cast<const unsigned char*>(input);
    auto* c=static_cast<unsigned char*>(output);
    alignas(64) __m512i x[4]={_mm512_loadu_si512(a),_mm512_loadu_si512(a+64),
                            _mm512_loadu_si512(a+128),_mm512_loadu_si512(a+192)};
    toeplitz4(x,packedKeys);
    _mm512_storeu_si512(c,x[0]);_mm512_storeu_si512(c+64,x[1]);
    _mm512_storeu_si512(c+128,x[2]);_mm512_storeu_si512(c+192,x[3]);
}

#undef SPARSE128_DEFER_INLINE
} // namespace spin::research::char2::field128_lift::sparse_byte::deferred
