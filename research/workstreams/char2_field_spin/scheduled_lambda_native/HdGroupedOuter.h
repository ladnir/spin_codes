#pragma once
#include "HdSumrankInner.h"
#include "Hd35Outer.h"
#ifndef HDGROUP_FIELD32
#define HDGROUP_FIELD32 0
#endif
#ifndef HDGROUP_NORMALIZE
#define HDGROUP_NORMALIZE 1
#endif
#ifndef HDGROUP_FUSE_LAST
#define HDGROUP_FUSE_LAST 0
#endif
#ifndef HDGROUP_SOURCE_PREFIX
#define HDGROUP_SOURCE_PREFIX 0
#endif
#ifndef HDGROUP_SUBFIELD
#define HDGROUP_SUBFIELD 0
#endif
#ifndef HDGROUP_TOEPLITZ
#define HDGROUP_TOEPLITZ 0
#endif
#ifndef HDGROUP_BYTE_KEYS
#define HDGROUP_BYTE_KEYS 0
#endif
#ifndef HDGROUP_TOEPLITZ_COMPACT7
#define HDGROUP_TOEPLITZ_COMPACT7 0
#endif
#ifndef HDGROUP_XOR_MIXER
#define HDGROUP_XOR_MIXER 0
#endif
#ifndef HDGROUP_BIT_SYNDROMES
#define HDGROUP_BIT_SYNDROMES 0
#endif
#ifndef HDGROUP_SCHEDULED_LAMBDA
#define HDGROUP_SCHEDULED_LAMBDA 0
#endif
#ifndef HDGROUP_CROSS_LAYERS
#define HDGROUP_CROSS_LAYERS 0
#endif
#ifndef HDGROUP_CROSS_VARYING
#define HDGROUP_CROSS_VARYING 1
#endif
#ifndef HDGROUP_LONG_RS
#define HDGROUP_LONG_RS 0
#endif
#ifndef HDGROUP_COSET_PIPELINE
#define HDGROUP_COSET_PIPELINE 0
#endif
#ifndef HDGROUP_COSET_PIPELINE_STREAM
#define HDGROUP_COSET_PIPELINE_STREAM 1
#endif

namespace spin::research::char2::hdgroup {
#if defined(_MSC_VER)
#define HDGROUP_INLINE __forceinline
#else
#define HDGROUP_INLINE inline __attribute__((always_inline))
#endif
using hdfast::outer_detail::normalized;
using Kind=hdbw::detail::Kind;
static_assert(HDGROUP_SOURCE_PREFIX>=0&&HDGROUP_SOURCE_PREFIX<=2);
static_assert(!HDGROUP_SOURCE_PREFIX||(HDGROUP_TOEPLITZ&&HDGROUP_SUBFIELD&&HDGROUP_BYTE_KEYS
    &&HDGROUP_FUSE_LAST&&!HDGROUP_LONG_RS),
    "Source-prefix fusion requires the certified GF256 byte-key Toeplitz/fused-last path");
static_assert(HDGROUP_COSET_PIPELINE==0||HDGROUP_COSET_PIPELINE==1);
static_assert(HDGROUP_COSET_PIPELINE_STREAM==0||HDGROUP_COSET_PIPELINE_STREAM==1);
static_assert(!(HDGROUP_SUBFIELD&&HDGROUP_FIELD32),"Subfield DHD and extension-field mixers are distinct ensembles");
static_assert(!(HDGROUP_TOEPLITZ&&HDGROUP_FIELD32),"Toeplitz and extension-field mixers are distinct ensembles");
static_assert(!HDGROUP_BYTE_KEYS||HDGROUP_SUBFIELD,"Byte setup keys require the GF256 ensemble");
#if HDGROUP_TOEPLITZ_COMPACT7
static_assert(HDGROUP_TOEPLITZ_COMPACT7==1);
static_assert(HDGROUP_TOEPLITZ&&HDGROUP_BYTE_KEYS&&HDGROUP_SUBFIELD,
    "Compact7 requires the GF256 byte-key Toeplitz4 path");
#endif
static_assert(!HDGROUP_LONG_RS||(HDGROUP_SUBFIELD&&HDGROUP_BYTE_KEYS&&!HDGROUP_TOEPLITZ
    &&!HDGROUP_XOR_MIXER&&!HDGROUP_FIELD32),
    "Long RS is a distinct GF256 byte-key outer ensemble, without a column mixer");
static_assert(HDGROUP_XOR_MIXER>=0&&HDGROUP_XOR_MIXER<=5);
static_assert(HDGROUP_BIT_SYNDROMES==0||HDGROUP_BIT_SYNDROMES==1);
static_assert(HDGROUP_SCHEDULED_LAMBDA==0||HDGROUP_SCHEDULED_LAMBDA==1);
static_assert(!HDGROUP_SCHEDULED_LAMBDA||(HDGROUP_XOR_MIXER==5&&HDGROUP_SUBFIELD&&HDGROUP_BYTE_KEYS
    &&!HDGROUP_TOEPLITZ&&!HDGROUP_FIELD32&&!HDGROUP_CROSS_LAYERS&&!HDGROUP_LONG_RS
    &&!HDGROUP_COSET_PIPELINE&&!HDGROUP_SOURCE_PREFIX),
    "Scheduled lambda requires the isolated no-cross GF256 G16 path");
static_assert(!HDGROUP_BIT_SYNDROMES||HDGROUP_XOR_MIXER==5,
    "Bit syndromes require the experimental G16 GF256 mixer");
static_assert(HDGROUP_CROSS_LAYERS>=0&&HDGROUP_CROSS_LAYERS<=2);
static_assert(HDGROUP_CROSS_VARYING==0||HDGROUP_CROSS_VARYING==1);
static_assert(!HDGROUP_CROSS_LAYERS||(HDGROUP_XOR_MIXER==5&&HDGROUP_BIT_SYNDROMES
    &&HDGROUP_SUBFIELD&&HDGROUP_BYTE_KEYS&&!HDGROUP_LONG_RS&&!HDGROUP_COSET_PIPELINE),
    "Cross-column pairs require the isolated byte-key H-star G16 path");
static_assert(!HDGROUP_XOR_MIXER||(HDGROUP_SUBFIELD&&!HDGROUP_TOEPLITZ&&!HDGROUP_FIELD32
    &&(HDGROUP_XOR_MIXER>=2||!HDGROUP_NORMALIZE)),
    "The XOR-column candidates require GF256; only D2 H D1 permits gauge normalization");
#if HDGROUP_BYTE_KEYS
using Coefficient=std::uint8_t;
inline Coefficient prepareCoefficient(std::uint16_t value){
    if(value>255)throw std::invalid_argument("Non-byte outer coefficient");
    return static_cast<std::uint8_t>(value);
}
#else
using Coefficient=wide16::Coefficient;
inline Coefficient prepareCoefficient(std::uint16_t value){return wide16::prepareCoefficient(value);}
#endif
template<unsigned Group>
inline constexpr unsigned toeplitzProducts=Group==2?3:9;
#if HDGROUP_TOEPLITZ_COMPACT7
// Stored coefficients differ from the unchanged nine GFNI products.
template<unsigned Group>
inline constexpr unsigned toeplitzStoredCoefficients=7;
template<unsigned Group>
inline constexpr unsigned beforeRegions=256/Group*toeplitzStoredCoefficients<Group>;
#else
template<unsigned Group>
inline constexpr unsigned beforeBaseRegions=HDGROUP_TOEPLITZ?256/Group*toeplitzProducts<Group>:
    (HDGROUP_NORMALIZE&&!(HDGROUP_FIELD32&&Group==2)?256-256/Group:256);
template<unsigned Group>
inline constexpr unsigned beforeRegions=beforeBaseRegions<Group>+256*HDGROUP_CROSS_LAYERS;
#endif

HDGROUP_INLINE __m512i diagonalProduct(__m512i value,const Coefficient* coefficients){
#if HDGROUP_BYTE_KEYS
    std::uint32_t four;std::memcpy(&four,coefficients,sizeof(four));
    const auto packed=_mm512_set1_epi32(static_cast<int>(four));
    const auto select=_mm512_setr_epi64(0,0,0x0101010101010101LL,0x0101010101010101LL,
        0x0202020202020202LL,0x0202020202020202LL,0x0303030303030303LL,0x0303030303030303LL);
    return _mm512_gf2p8mul_epi8(value,_mm512_shuffle_epi8(packed,select));
#else
    if constexpr(HDGROUP_SUBFIELD){
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(coefficients)));
        return _mm512_gf2p8mul_epi8(value,wide16::expandPairs<0>(packed));
    }else return wide16::multiply512(value,coefficients);
#endif
}

#if HDGROUP_LONG_RS
// Deterministic RS[256,128], followed by independent nonzero GF256 output
// diagonals. Retain the established explicit sixteen-register FFT leaves.
// Four lane coefficients for each output are contiguous in the hot key table.
template<unsigned Offset,bool Stream>
HD35_NOINLINE void longLeaf16(const __m512i* source,Block* planes,std::size_t stride,
                             const Coefficient* diagonal){
#define HDGROUP_RS_LOAD(I) auto x##I=_mm512_load_si512(source+I)
    HDGROUP_RS_LOAD(0);HDGROUP_RS_LOAD(1);HDGROUP_RS_LOAD(2);HDGROUP_RS_LOAD(3);
    HDGROUP_RS_LOAD(4);HDGROUP_RS_LOAD(5);HDGROUP_RS_LOAD(6);HDGROUP_RS_LOAD(7);
    HDGROUP_RS_LOAD(8);HDGROUP_RS_LOAD(9);HDGROUP_RS_LOAD(10);HDGROUP_RS_LOAD(11);
    HDGROUP_RS_LOAD(12);HDGROUP_RS_LOAD(13);HDGROUP_RS_LOAD(14);HDGROUP_RS_LOAD(15);
#undef HDGROUP_RS_LOAD
#define HDGROUP_RS_PAIR(A,B,LEVEL,DELTA) hdfast::outer_detail::pair<true,normalized(LEVEL,Offset^DELTA)>(x##A,x##B)
    HDGROUP_RS_PAIR(0,8,3,0);HDGROUP_RS_PAIR(1,9,3,0);
    HDGROUP_RS_PAIR(2,10,3,0);HDGROUP_RS_PAIR(3,11,3,0);
    HDGROUP_RS_PAIR(4,12,3,0);HDGROUP_RS_PAIR(5,13,3,0);
    HDGROUP_RS_PAIR(6,14,3,0);HDGROUP_RS_PAIR(7,15,3,0);
    HDGROUP_RS_PAIR(0,4,2,0);HDGROUP_RS_PAIR(1,5,2,0);
    HDGROUP_RS_PAIR(2,6,2,0);HDGROUP_RS_PAIR(3,7,2,0);
    HDGROUP_RS_PAIR(8,12,2,8);HDGROUP_RS_PAIR(9,13,2,8);
    HDGROUP_RS_PAIR(10,14,2,8);HDGROUP_RS_PAIR(11,15,2,8);
    HDGROUP_RS_PAIR(0,2,1,0);HDGROUP_RS_PAIR(1,3,1,0);
    HDGROUP_RS_PAIR(4,6,1,4);HDGROUP_RS_PAIR(5,7,1,4);
    HDGROUP_RS_PAIR(8,10,1,8);HDGROUP_RS_PAIR(9,11,1,8);
    HDGROUP_RS_PAIR(12,14,1,12);HDGROUP_RS_PAIR(13,15,1,12);
    HDGROUP_RS_PAIR(0,1,0,0);HDGROUP_RS_PAIR(2,3,0,2);
    HDGROUP_RS_PAIR(4,5,0,4);HDGROUP_RS_PAIR(6,7,0,6);
    HDGROUP_RS_PAIR(8,9,0,8);HDGROUP_RS_PAIR(10,11,0,10);
    HDGROUP_RS_PAIR(12,13,0,12);HDGROUP_RS_PAIR(14,15,0,14);
#undef HDGROUP_RS_PAIR
#define HDGROUP_RS_STORE(I) do { \
    const auto scaled=diagonalProduct(x##I,diagonal+4*I); \
    if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(planes+I*stride),scaled); \
    else _mm512_storeu_si512(planes+I*stride,scaled); \
} while(false)
    HDGROUP_RS_STORE(0);HDGROUP_RS_STORE(1);HDGROUP_RS_STORE(2);HDGROUP_RS_STORE(3);
    HDGROUP_RS_STORE(4);HDGROUP_RS_STORE(5);HDGROUP_RS_STORE(6);HDGROUP_RS_STORE(7);
    HDGROUP_RS_STORE(8);HDGROUP_RS_STORE(9);HDGROUP_RS_STORE(10);HDGROUP_RS_STORE(11);
    HDGROUP_RS_STORE(12);HDGROUP_RS_STORE(13);HDGROUP_RS_STORE(14);HDGROUP_RS_STORE(15);
#undef HDGROUP_RS_STORE
}

template<unsigned Offset,bool Stream>
HD35_NOINLINE void longCoset128(const Block* source,Block* planes,std::size_t stride,
                               const Coefficient* diagonal){
    static_assert(Offset==0||Offset==128);
    alignas(64) __m512i values[128];
    for(unsigned base=0;base!=16;base+=2){
        c8Unroll(std::make_index_sequence<2>{},[&](auto index){
            const auto j=base+unsigned(index);
            auto a=_mm512_loadu_si512(source+4*j);
            auto b=_mm512_loadu_si512(source+4*(16+j));
            auto c=_mm512_loadu_si512(source+4*(32+j));
            auto d=_mm512_loadu_si512(source+4*(48+j));
            auto e=_mm512_loadu_si512(source+4*(64+j));
            auto f=_mm512_loadu_si512(source+4*(80+j));
            auto g=_mm512_loadu_si512(source+4*(96+j));
            auto h=_mm512_loadu_si512(source+4*(112+j));
            hdfast::outer_detail::pair<true,normalized(6,Offset)>(a,e);
            hdfast::outer_detail::pair<true,normalized(6,Offset)>(b,f);
            hdfast::outer_detail::pair<true,normalized(6,Offset)>(c,g);
            hdfast::outer_detail::pair<true,normalized(6,Offset)>(d,h);
            hdfast::outer_detail::pair<true,normalized(5,Offset)>(a,c);
            hdfast::outer_detail::pair<true,normalized(5,Offset)>(b,d);
            hdfast::outer_detail::pair<true,normalized(5,Offset+64)>(e,g);
            hdfast::outer_detail::pair<true,normalized(5,Offset+64)>(f,h);
            hdfast::outer_detail::pair<true,normalized(4,Offset)>(a,b);
            hdfast::outer_detail::pair<true,normalized(4,Offset+32)>(c,d);
            hdfast::outer_detail::pair<true,normalized(4,Offset+64)>(e,f);
            hdfast::outer_detail::pair<true,normalized(4,Offset+96)>(g,h);
            _mm512_store_si512(values+j,a);_mm512_store_si512(values+16+j,b);
            _mm512_store_si512(values+32+j,c);_mm512_store_si512(values+48+j,d);
            _mm512_store_si512(values+64+j,e);_mm512_store_si512(values+80+j,f);
            _mm512_store_si512(values+96+j,g);_mm512_store_si512(values+112+j,h);
        });
    }
    longLeaf16<Offset,Stream>(values,planes,stride,diagonal);
    longLeaf16<Offset+16,Stream>(values+16,planes+16*stride,stride,diagonal+64);
    longLeaf16<Offset+32,Stream>(values+32,planes+32*stride,stride,diagonal+128);
    longLeaf16<Offset+48,Stream>(values+48,planes+48*stride,stride,diagonal+192);
    longLeaf16<Offset+64,Stream>(values+64,planes+64*stride,stride,diagonal+256);
    longLeaf16<Offset+80,Stream>(values+80,planes+80*stride,stride,diagonal+320);
    longLeaf16<Offset+96,Stream>(values+96,planes+96*stride,stride,diagonal+384);
    longLeaf16<Offset+112,Stream>(values+112,planes+112*stride,stride,diagonal+448);
}
#endif

// [[a,b],[c,a]] times (x,y), with (a,a+b,a+c) prepared in setup.
// Each product applies four independently sampled coefficients to four records.
HDGROUP_INLINE void toeplitz2(__m512i x,__m512i y,const Coefficient* p,
                             __m512i& out0,__m512i& out1){
    const auto shared=diagonalProduct(_mm512_xor_si512(x,y),p);
    out0=_mm512_xor_si512(shared,diagonalProduct(y,p+4));
    out1=_mm512_xor_si512(shared,diagonalProduct(x,p+8));
}
#if HDGROUP_TOEPLITZ_COMPACT7
// Each four-byte key contains one scalar for each 128-bit record. This is
// the unchanged byte-key expansion, separated from its GFNI product.
HDGROUP_INLINE __m512i compactToeplitzKey(const Coefficient* coefficients){
    std::uint32_t four;std::memcpy(&four,coefficients,sizeof(four));
    const auto packed=_mm512_set1_epi32(static_cast<int>(four));
    const auto select=_mm512_setr_epi64(0,0,0x0101010101010101LL,0x0101010101010101LL,
        0x0202020202020202LL,0x0202020202020202LL,0x0303030303030303LL,0x0303030303030303LL);
    return _mm512_shuffle_epi8(packed,select);
}
// Keys are (A0,A1,A2,V0,V2,W0,W1). For the old nine-key circuit,
// V1=V0^(A1^A2) and W2=W0^(A1^A2). Reconstruct expanded coefficients
// without changing the nine products, data layout, or sampled matrix.
HDGROUP_INLINE void compactToeplitz4(__m512i* x,const Coefficient* p){
    const auto a0=compactToeplitzKey(p);
    const auto a1=compactToeplitzKey(p+4),a2=compactToeplitzKey(p+8);
    const auto keyXor=_mm512_xor_si512(a1,a2);
    const auto low=_mm512_xor_si512(x[0],x[2]),high=_mm512_xor_si512(x[1],x[3]);
    const auto sharedU=_mm512_gf2p8mul_epi8(_mm512_xor_si512(low,high),a0);
    const auto u0=_mm512_xor_si512(sharedU,_mm512_gf2p8mul_epi8(high,a1));
    const auto u1=_mm512_xor_si512(sharedU,_mm512_gf2p8mul_epi8(low,a2));
    const auto v0=compactToeplitzKey(p+12),v2=compactToeplitzKey(p+16);
    const auto sharedV=_mm512_gf2p8mul_epi8(_mm512_xor_si512(x[0],x[1]),v0);
    const auto outV0=_mm512_xor_si512(sharedV,_mm512_gf2p8mul_epi8(x[1],_mm512_xor_si512(v0,keyXor)));
    const auto outV1=_mm512_xor_si512(sharedV,_mm512_gf2p8mul_epi8(x[0],v2));
    const auto w0=compactToeplitzKey(p+20),w1=compactToeplitzKey(p+24);
    const auto sharedW=_mm512_gf2p8mul_epi8(_mm512_xor_si512(x[2],x[3]),w0);
    const auto outW0=_mm512_xor_si512(sharedW,_mm512_gf2p8mul_epi8(x[3],w1));
    const auto outW1=_mm512_xor_si512(sharedW,_mm512_gf2p8mul_epi8(x[2],_mm512_xor_si512(w0,keyXor)));
    x[0]=_mm512_xor_si512(u0,outW0);x[1]=_mm512_xor_si512(u1,outW1);
    x[2]=_mm512_xor_si512(u0,outV0);x[3]=_mm512_xor_si512(u1,outV1);
}
#endif
template<unsigned Group>
HDGROUP_INLINE void toeplitzColumn(__m512i* x,const Coefficient* p){
#if HDGROUP_TOEPLITZ_COMPACT7
    static_assert(Group==4,"Compact7 supports only four-symbol Toeplitz columns");
    compactToeplitz4(x,p);
#else
    static_assert(Group==2||Group==4,"Toeplitz prototype supports two or four symbols");
    if constexpr(Group==2)toeplitz2(x[0],x[1],p,x[0],x[1]);
    else{
        // T=[[A,B],[C,A]]. Three 2x2 products, nine scalar products total.
        __m512i u0,u1,v0,v1,w0,w1;
        toeplitz2(_mm512_xor_si512(x[0],x[2]),_mm512_xor_si512(x[1],x[3]),p,u0,u1);
        toeplitz2(x[0],x[1],p+12,v0,v1); // (C+A) x_low
        toeplitz2(x[2],x[3],p+24,w0,w1); // (B+A) x_high
        x[0]=_mm512_xor_si512(u0,w0);x[1]=_mm512_xor_si512(u1,w1);
        x[2]=_mm512_xor_si512(u0,v0);x[3]=_mm512_xor_si512(u1,v1);
    }
#endif
}

// Every index and butterfly coefficient is compile-time fixed. Small column
// transforms stay in registers; RS leaves retain the established explicit
// sixteen-register schedule rather than a runtime FFT loop.
template<unsigned Size,unsigned Offset,bool Inverse=false>
HDGROUP_INLINE void tinyFft(__m512i* x){
    if constexpr(Size>1){
        constexpr unsigned half=Size/2,level=std::countr_zero(Size)-1;
        if constexpr(Inverse){tinyFft<half,Offset,true>(x);tinyFft<half,Offset+half,true>(x+half);}
        c8Unroll(std::make_index_sequence<half>{},[&](auto j){
            hdsumrank::affinePair<normalized(level,Offset),Inverse?Kind::Inverse:Kind::Forward>(x[j],x[j+half]);
        });
        if constexpr(!Inverse){tinyFft<half,Offset>(x);tinyFft<half,Offset+half>(x+half);}
    }
}
template<unsigned Size,unsigned Offset>
HD35_NOINLINE void rsFft(__m512i* x){
    if constexpr(Size==16)hdfast::outer_detail::leaf16<Offset,true>(x,reinterpret_cast<Block*>(x),4);
    else if constexpr(Size<16)tinyFft<Size,Offset>(x);
    else{
        constexpr unsigned half=Size/2,level=std::countr_zero(Size)-1;
        c8Unroll(std::make_index_sequence<half>{},[&](auto j){
            hdfast::outer_detail::pair<true,normalized(level,Offset)>(x[j],x[j+half]);
        });
        rsFft<half,Offset>(x);rsFft<half,Offset+half>(x+half);
    }
}

// Leave only the final size-two butterflies for the column consumer. The
// bounded leaf is fully unrolled, like tinyFft; no runtime recursion or callback.
template<unsigned Size,unsigned Offset>
HDGROUP_INLINE void tinyPrefix(__m512i* x){
    if constexpr(Size>2){
        constexpr unsigned half=Size/2,level=std::countr_zero(Size)-1;
        c8Unroll(std::make_index_sequence<half>{},[&](auto j){
            hdfast::outer_detail::pair<true,normalized(level,Offset)>(x[j],x[j+half]);
        });
        tinyPrefix<half,Offset>(x);tinyPrefix<half,Offset+half>(x+half);
    }
}
template<unsigned Size,unsigned Offset>
HD35_NOINLINE void rsPrefix(__m512i* x){
    if constexpr(Size<=16)tinyPrefix<Size,Offset>(x);
    else{
        constexpr unsigned half=Size/2,level=std::countr_zero(Size)-1;
        c8Unroll(std::make_index_sequence<half>{},[&](auto j){
            hdfast::outer_detail::pair<true,normalized(level,Offset)>(x[j],x[j+half]);
        });
        rsPrefix<half,Offset>(x);rsPrefix<half,Offset+half>(x+half);
    }
}

#if HDGROUP_SOURCE_PREFIX
// Exact RS32 prefix: fuse the source loads with its first size-32 stage.
// Store only the transformed pair, then retain the two existing RS16 leaves.
// Mode 1 loads one pair at a time; mode 2 stages the low half in registers
// before consuming the high half. Source/output layouts are unchanged.
template<unsigned Offset>
HDGROUP_INLINE void sourcePrefix32(const Block* source,__m512i* values){
    static_assert(Offset==0||Offset==32,"Certified Group4 cosets have offsets zero and thirty-two");
#if HDGROUP_SOURCE_PREFIX==2
#define HDGROUP_SOURCE_LOW(J) auto a##J=_mm512_loadu_si512(source+4*(J))
    HDGROUP_SOURCE_LOW(0);HDGROUP_SOURCE_LOW(1);HDGROUP_SOURCE_LOW(2);HDGROUP_SOURCE_LOW(3);
    HDGROUP_SOURCE_LOW(4);HDGROUP_SOURCE_LOW(5);HDGROUP_SOURCE_LOW(6);HDGROUP_SOURCE_LOW(7);
    HDGROUP_SOURCE_LOW(8);HDGROUP_SOURCE_LOW(9);HDGROUP_SOURCE_LOW(10);HDGROUP_SOURCE_LOW(11);
    HDGROUP_SOURCE_LOW(12);HDGROUP_SOURCE_LOW(13);HDGROUP_SOURCE_LOW(14);HDGROUP_SOURCE_LOW(15);
#undef HDGROUP_SOURCE_LOW
#define HDGROUP_SOURCE_HIGH(J) do { \
    auto b=_mm512_loadu_si512(source+4*((J)+16)); \
    hdfast::outer_detail::pair<true,normalized(4,Offset)>(a##J,b); \
    _mm512_store_si512(values+(J),a##J); \
    _mm512_store_si512(values+(J)+16,b); \
} while(false)
    HDGROUP_SOURCE_HIGH(0);HDGROUP_SOURCE_HIGH(1);HDGROUP_SOURCE_HIGH(2);HDGROUP_SOURCE_HIGH(3);
    HDGROUP_SOURCE_HIGH(4);HDGROUP_SOURCE_HIGH(5);HDGROUP_SOURCE_HIGH(6);HDGROUP_SOURCE_HIGH(7);
    HDGROUP_SOURCE_HIGH(8);HDGROUP_SOURCE_HIGH(9);HDGROUP_SOURCE_HIGH(10);HDGROUP_SOURCE_HIGH(11);
    HDGROUP_SOURCE_HIGH(12);HDGROUP_SOURCE_HIGH(13);HDGROUP_SOURCE_HIGH(14);HDGROUP_SOURCE_HIGH(15);
#undef HDGROUP_SOURCE_HIGH
#else
#define HDGROUP_SOURCE_PAIR(J) do { \
    auto a=_mm512_loadu_si512(source+4*(J)); \
    auto b=_mm512_loadu_si512(source+4*((J)+16)); \
    hdfast::outer_detail::pair<true,normalized(4,Offset)>(a,b); \
    _mm512_store_si512(values+(J),a); \
    _mm512_store_si512(values+(J)+16,b); \
} while(false)
    HDGROUP_SOURCE_PAIR(0);HDGROUP_SOURCE_PAIR(1);HDGROUP_SOURCE_PAIR(2);HDGROUP_SOURCE_PAIR(3);
    HDGROUP_SOURCE_PAIR(4);HDGROUP_SOURCE_PAIR(5);HDGROUP_SOURCE_PAIR(6);HDGROUP_SOURCE_PAIR(7);
    HDGROUP_SOURCE_PAIR(8);HDGROUP_SOURCE_PAIR(9);HDGROUP_SOURCE_PAIR(10);HDGROUP_SOURCE_PAIR(11);
    HDGROUP_SOURCE_PAIR(12);HDGROUP_SOURCE_PAIR(13);HDGROUP_SOURCE_PAIR(14);HDGROUP_SOURCE_PAIR(15);
#undef HDGROUP_SOURCE_PAIR
#endif
    rsPrefix<16,Offset>(values);rsPrefix<16,Offset+16>(values+16);
}
#endif

#if HDGROUP_XOR_MIXER==3 || HDGROUP_XOR_MIXER==4
// H=I+11^T+uu^T in the AES polynomial basis. Mode 3 uses u=(0,...,7),
// four fixed byte products and 28 binary XORs per eight vectors. Mode 4
// uses u=(8,...,15), adding two fixed products and two XORs. Neither
// circuit samples keys or changes the diagonal layout.
HDGROUP_INLINE void twoSumColumn(__m512i* x){
    const auto p01=_mm512_xor_si512(x[0],x[1]);
    const auto p23=_mm512_xor_si512(x[2],x[3]);
    const auto p45=_mm512_xor_si512(x[4],x[5]);
    const auto p67=_mm512_xor_si512(x[6],x[7]);
    const auto b=_mm512_xor_si512(p23,p67);
    const auto s=_mm512_xor_si512(_mm512_xor_si512(p01,p45),b);
    const auto c=_mm512_xor_si512(p45,p67);
    const auto a=_mm512_xor_si512(_mm512_xor_si512(x[1],x[3]),_mm512_xor_si512(x[5],x[7]));
    const auto tLow=_mm512_xor_si512(a,_mm512_xor_si512(
        _mm512_gf2p8mul_epi8(b,_mm512_set1_epi8(2)),
        _mm512_gf2p8mul_epi8(c,_mm512_set1_epi8(4))));
#if HDGROUP_XOR_MIXER==4
    const auto t=_mm512_xor_si512(tLow,_mm512_gf2p8mul_epi8(s,_mm512_set1_epi8(8)));
    const auto common=_mm512_xor_si512(s,_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(8)));
#else
    const auto t=tLow,common=s;
#endif
    const auto t2=_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(2));
    const auto t4=_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(4));
    const auto c1=_mm512_xor_si512(common,t),c2=_mm512_xor_si512(common,t2);
    const auto c3=_mm512_xor_si512(c1,t2),c4=_mm512_xor_si512(common,t4);
    const auto c5=_mm512_xor_si512(c1,t4),c6=_mm512_xor_si512(c2,t4);
    const auto c7=_mm512_xor_si512(c3,t4);
    x[0]=_mm512_xor_si512(x[0],common);x[1]=_mm512_xor_si512(x[1],c1);
    x[2]=_mm512_xor_si512(x[2],c2);x[3]=_mm512_xor_si512(x[3],c3);
    x[4]=_mm512_xor_si512(x[4],c4);x[5]=_mm512_xor_si512(x[5],c5);
    x[6]=_mm512_xor_si512(x[6],c6);x[7]=_mm512_xor_si512(x[7],c7);
}
#endif

#if HDGROUP_XOR_MIXER==5
// Sixteen labels 06 + span_GF2(01,0c,50,b0) in AES GF256. The input
// bit-sums share 26 XORs; the complete mixer uses 8 fixed products/62 XORs.
// Bit syndromes reuse these sums, adding four XORs without field products.
HDGROUP_INLINE void twoSum16Column(__m512i* x){
    const auto p0=_mm512_xor_si512(x[0],x[1]),p1=_mm512_xor_si512(x[2],x[3]);
    const auto p2=_mm512_xor_si512(x[4],x[5]),p3=_mm512_xor_si512(x[6],x[7]);
    const auto p4=_mm512_xor_si512(x[8],x[9]),p5=_mm512_xor_si512(x[10],x[11]);
    const auto p6=_mm512_xor_si512(x[12],x[13]),p7=_mm512_xor_si512(x[14],x[15]);
    const auto a0=_mm512_xor_si512(
        _mm512_xor_si512(_mm512_xor_si512(x[1],x[3]),_mm512_xor_si512(x[5],x[7])),
        _mm512_xor_si512(_mm512_xor_si512(x[9],x[11]),_mm512_xor_si512(x[13],x[15])));
    const auto a1=_mm512_xor_si512(_mm512_xor_si512(p1,p3),_mm512_xor_si512(p5,p7));
    const auto q0=_mm512_xor_si512(p0,p1),q1=_mm512_xor_si512(p2,p3);
    const auto q2=_mm512_xor_si512(p4,p5),q3=_mm512_xor_si512(p6,p7);
    const auto a2=_mm512_xor_si512(q1,q3),a3=_mm512_xor_si512(q2,q3);
    const auto s=_mm512_xor_si512(_mm512_xor_si512(q0,q1),a3);
    const auto t=_mm512_xor_si512(a0,_mm512_xor_si512(
        _mm512_xor_si512(_mm512_gf2p8mul_epi8(s,_mm512_set1_epi8(0x06)),
                        _mm512_gf2p8mul_epi8(a1,_mm512_set1_epi8(0x0c))),
        _mm512_xor_si512(_mm512_gf2p8mul_epi8(a2,_mm512_set1_epi8(0x50)),
                        _mm512_gf2p8mul_epi8(a3,_mm512_set1_epi8(static_cast<char>(0xb0))))));
    const auto common=_mm512_xor_si512(s,_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(0x06)));
#if HDGROUP_BIT_SYNDROMES
    const auto t1=_mm512_xor_si512(t,a0);
    const auto t2=_mm512_xor_si512(_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(0x0c)),a1);
    const auto t4=_mm512_xor_si512(_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(0x50)),a2);
    const auto t8=_mm512_xor_si512(_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(static_cast<char>(0xb0))),a3);
    const auto c1=_mm512_xor_si512(common,t1),c2=_mm512_xor_si512(common,t2);
#else
    const auto t2=_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(0x0c));
    const auto t4=_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(0x50));
    const auto t8=_mm512_gf2p8mul_epi8(t,_mm512_set1_epi8(static_cast<char>(0xb0)));
    const auto c1=_mm512_xor_si512(common,t),c2=_mm512_xor_si512(common,t2);
#endif
    const auto c3=_mm512_xor_si512(c1,t2),c4=_mm512_xor_si512(common,t4);
    const auto c5=_mm512_xor_si512(c1,t4),c6=_mm512_xor_si512(c2,t4);
    const auto c7=_mm512_xor_si512(c3,t4),c8=_mm512_xor_si512(common,t8);
    const auto c9=_mm512_xor_si512(c1,t8),c10=_mm512_xor_si512(c2,t8);
    const auto c11=_mm512_xor_si512(c3,t8),c12=_mm512_xor_si512(c4,t8);
    const auto c13=_mm512_xor_si512(c5,t8),c14=_mm512_xor_si512(c6,t8);
    const auto c15=_mm512_xor_si512(c7,t8);
    x[0]=_mm512_xor_si512(x[0],common);x[1]=_mm512_xor_si512(x[1],c1);
    x[2]=_mm512_xor_si512(x[2],c2);x[3]=_mm512_xor_si512(x[3],c3);
    x[4]=_mm512_xor_si512(x[4],c4);x[5]=_mm512_xor_si512(x[5],c5);
    x[6]=_mm512_xor_si512(x[6],c6);x[7]=_mm512_xor_si512(x[7],c7);
    x[8]=_mm512_xor_si512(x[8],c8);x[9]=_mm512_xor_si512(x[9],c9);
    x[10]=_mm512_xor_si512(x[10],c10);x[11]=_mm512_xor_si512(x[11],c11);
    x[12]=_mm512_xor_si512(x[12],c12);x[13]=_mm512_xor_si512(x[13],c13);
    x[14]=_mm512_xor_si512(x[14],c14);x[15]=_mm512_xor_si512(x[15],c15);
}
#endif

#if HDGROUP_SCHEDULED_LAMBDA
// Label arithmetic only: GF16=GF2[x]/(x^4+x+1), not AES byte arithmetic.
constexpr unsigned scheduledLabelMultiply(unsigned a,unsigned b){
    unsigned value=0;
    while(b){if(b&1)value^=a;b>>=1;a<<=1;if(a&16)a^=0x13;}
    return value;
}
constexpr unsigned scheduledLabelInverse(unsigned a){
    unsigned value=1;
    for(unsigned i=0;i!=14;++i)value=scheduledLabelMultiply(value,a);
    return a?value:0;
}
// H_lambda=I+lambda*11^T+sum_k b_k*b_k^T, lambda=bc in AES GF256.
// Shared bit sums:26 XOR; correction ladder:15; output:16; one GFNI.
HDGROUP_INLINE void lambda16Column(__m512i* x){
    const auto p0=_mm512_xor_si512(x[0],x[1]),p1=_mm512_xor_si512(x[2],x[3]);
    const auto p2=_mm512_xor_si512(x[4],x[5]),p3=_mm512_xor_si512(x[6],x[7]);
    const auto p4=_mm512_xor_si512(x[8],x[9]),p5=_mm512_xor_si512(x[10],x[11]);
    const auto p6=_mm512_xor_si512(x[12],x[13]),p7=_mm512_xor_si512(x[14],x[15]);
    const auto a0=_mm512_xor_si512(
        _mm512_xor_si512(_mm512_xor_si512(x[1],x[3]),_mm512_xor_si512(x[5],x[7])),
        _mm512_xor_si512(_mm512_xor_si512(x[9],x[11]),_mm512_xor_si512(x[13],x[15])));
    const auto a1=_mm512_xor_si512(_mm512_xor_si512(p1,p3),_mm512_xor_si512(p5,p7));
    const auto q0=_mm512_xor_si512(p0,p1),q1=_mm512_xor_si512(p2,p3);
    const auto q2=_mm512_xor_si512(p4,p5),q3=_mm512_xor_si512(p6,p7);
    const auto a2=_mm512_xor_si512(q1,q3),a3=_mm512_xor_si512(q2,q3);
    const auto s=_mm512_xor_si512(_mm512_xor_si512(q0,q1),a3);
    const auto common=_mm512_gf2p8mul_epi8(s,_mm512_set1_epi8(static_cast<char>(0xbc)));
    const auto c1=_mm512_xor_si512(common,a0),c2=_mm512_xor_si512(common,a1);
    const auto c3=_mm512_xor_si512(c1,a1),c4=_mm512_xor_si512(common,a2);
    const auto c5=_mm512_xor_si512(c1,a2),c6=_mm512_xor_si512(c2,a2);
    const auto c7=_mm512_xor_si512(c3,a2),c8=_mm512_xor_si512(common,a3);
    const auto c9=_mm512_xor_si512(c1,a3),c10=_mm512_xor_si512(c2,a3);
    const auto c11=_mm512_xor_si512(c3,a3),c12=_mm512_xor_si512(c4,a3);
    const auto c13=_mm512_xor_si512(c5,a3),c14=_mm512_xor_si512(c6,a3);
    const auto c15=_mm512_xor_si512(c7,a3);
    x[0]=_mm512_xor_si512(x[0],common);x[1]=_mm512_xor_si512(x[1],c1);
    x[2]=_mm512_xor_si512(x[2],c2);x[3]=_mm512_xor_si512(x[3],c3);
    x[4]=_mm512_xor_si512(x[4],c4);x[5]=_mm512_xor_si512(x[5],c5);
    x[6]=_mm512_xor_si512(x[6],c6);x[7]=_mm512_xor_si512(x[7],c7);
    x[8]=_mm512_xor_si512(x[8],c8);x[9]=_mm512_xor_si512(x[9],c9);
    x[10]=_mm512_xor_si512(x[10],c10);x[11]=_mm512_xor_si512(x[11],c11);
    x[12]=_mm512_xor_si512(x[12],c12);x[13]=_mm512_xor_si512(x[13],c13);
    x[14]=_mm512_xor_si512(x[14],c14);x[15]=_mm512_xor_si512(x[15],c15);
}

template<unsigned Offset,unsigned Local,bool Stream>
HDGROUP_INLINE void scheduledLambdaStore(const __m512i* values,Block* planes,std::size_t stride,
                                        const Coefficient* before,const Coefficient* after){
    static_assert((Offset==0||Offset==8)&&Local<8);
    constexpr unsigned j=Offset+Local;
    alignas(64) __m512i x[16];
    // P_j(i)=inv(i xor j). The inverse lookup is i=inv(r) xor j.
    // D1 remains indexed by ORIGINAL i: its normalized pivot moves to P_j(0).
#define HDGROUP_LAMBDA_LOAD(R) do { \
        constexpr unsigned i=scheduledLabelInverse(R)^j; \
        x[R]=values[i*8+Local]; \
        if constexpr(HDGROUP_NORMALIZE){ \
            if constexpr(i!=0)x[R]=diagonalProduct(x[R],before+(j*15+i-1)*4); \
        }else x[R]=diagonalProduct(x[R],before+(j*16+i)*4); \
    } while(false)
    HDGROUP_LAMBDA_LOAD(0);HDGROUP_LAMBDA_LOAD(1);HDGROUP_LAMBDA_LOAD(2);HDGROUP_LAMBDA_LOAD(3);
    HDGROUP_LAMBDA_LOAD(4);HDGROUP_LAMBDA_LOAD(5);HDGROUP_LAMBDA_LOAD(6);HDGROUP_LAMBDA_LOAD(7);
    HDGROUP_LAMBDA_LOAD(8);HDGROUP_LAMBDA_LOAD(9);HDGROUP_LAMBDA_LOAD(10);HDGROUP_LAMBDA_LOAD(11);
    HDGROUP_LAMBDA_LOAD(12);HDGROUP_LAMBDA_LOAD(13);HDGROUP_LAMBDA_LOAD(14);HDGROUP_LAMBDA_LOAD(15);
#undef HDGROUP_LAMBDA_LOAD
    lambda16Column(x);
    // Q_j(r)=r xor j. Preserve physical store order and physical D2 keys.
#define HDGROUP_LAMBDA_STORE(R) do { \
        const auto y=diagonalProduct(x[(R)^j],after+(j*16+(R))*4); \
        auto* out=planes+(j*16+(R))*stride; \
        if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),y); \
        else _mm512_storeu_si512(out,y); \
    } while(false)
    HDGROUP_LAMBDA_STORE(0);HDGROUP_LAMBDA_STORE(1);HDGROUP_LAMBDA_STORE(2);HDGROUP_LAMBDA_STORE(3);
    HDGROUP_LAMBDA_STORE(4);HDGROUP_LAMBDA_STORE(5);HDGROUP_LAMBDA_STORE(6);HDGROUP_LAMBDA_STORE(7);
    HDGROUP_LAMBDA_STORE(8);HDGROUP_LAMBDA_STORE(9);HDGROUP_LAMBDA_STORE(10);HDGROUP_LAMBDA_STORE(11);
    HDGROUP_LAMBDA_STORE(12);HDGROUP_LAMBDA_STORE(13);HDGROUP_LAMBDA_STORE(14);HDGROUP_LAMBDA_STORE(15);
#undef HDGROUP_LAMBDA_STORE
}
#endif

template<unsigned Group,bool Stream>
HDGROUP_INLINE void mixStore(__m512i* column,unsigned j,Block* planes,std::size_t stride,
                            const Coefficient* before,const Coefficient* after){
    if constexpr(HDGROUP_XOR_MIXER){
        // New ensembles, not certificate-preserving replacements: modes
        // 1/2 use H=I+11^T; modes 3/4/5 use two-sum H. Modes 2/3/4/5
        // normalize the first D1 coefficient into the output diagonal.
        c8Unroll(std::make_index_sequence<Group>{},[&](auto r){
            if constexpr(HDGROUP_NORMALIZE){
                if constexpr(r!=0)column[r]=diagonalProduct(column[r],before+(j*(Group-1)+r-1)*4);
            }else column[r]=diagonalProduct(column[r],before+(j*Group+r)*4);
        });
#if HDGROUP_XOR_MIXER==5
        static_assert(Group==16,"The GF16-coset mixer has exactly sixteen rows");
        twoSum16Column(column);
#elif HDGROUP_XOR_MIXER>=3
        static_assert(Group==8,"The two-sum mixer has exactly eight rows");
        twoSumColumn(column);
#else
        auto sum=_mm512_setzero_si512();
        c8Unroll(std::make_index_sequence<Group>{},[&](auto r){sum=_mm512_xor_si512(sum,column[r]);});
        c8Unroll(std::make_index_sequence<Group>{},[&](auto r){column[r]=_mm512_xor_si512(column[r],sum);});
#endif
    }else if constexpr(HDGROUP_TOEPLITZ){
#if HDGROUP_TOEPLITZ_COMPACT7
        toeplitzColumn<Group>(column,before+j*toeplitzStoredCoefficients<Group>*4);
#else
        toeplitzColumn<Group>(column,before+j*toeplitzProducts<Group>*4);
#endif
    }else if constexpr(HDGROUP_FIELD32&&Group==2){
        // Scalar a+b*v in F[v]/(v^2+v+0x2000). No lane repacking.
        const auto* p=before+j*8;const auto* q=after+j*8;
        const auto x0=column[0],x1=column[1];
        if constexpr(HDGROUP_FIELD32==2){
            const auto a=diagonalProduct(x0,p),b=diagonalProduct(x1,p+4);
            const auto cross=diagonalProduct(_mm512_xor_si512(x0,x1),q+4);
            const auto swapped=medium32::selectBytes<1,0,3,2>(b);
            const auto deltaB=_mm512_gf2p8mul_epi8(_mm512_ternarylogic_epi32(swapped,b,_mm512_set1_epi16(0xff00),0x78),_mm512_set1_epi16(0x206c));
            column[0]=_mm512_xor_si512(a,deltaB);column[1]=_mm512_xor_si512(cross,a);
        }else{
            column[0]=_mm512_xor_si512(diagonalProduct(x0,p),diagonalProduct(x1,q));
            column[1]=_mm512_xor_si512(diagonalProduct(x0,p+4),diagonalProduct(x1,q+4));
        }
    }else{
        c8Unroll(std::make_index_sequence<Group>{},[&](auto r){
            if constexpr(HDGROUP_NORMALIZE&&r==0){} // exact gauge-normalized identity
            else if constexpr(HDGROUP_NORMALIZE)column[r]=diagonalProduct(column[r],before+(j*(Group-1)+r-1)*4);
            else column[r]=diagonalProduct(column[r],before+(j*Group+r)*4);
        });
        tinyFft<Group,0,true>(column);tinyFft<Group,Group>(column);
    }
#if HDGROUP_XOR_MIXER==5
    // Keep the sixteen diagonal/stores inside this always-inlined kernel.
    // GCC otherwise outlines c8Unroll's large closure, forcing the complete
    // mixed column through a 1 KiB stack temporary before reloading it.
    const auto* diagonals=after+j*Group*4;
    auto* destination=planes+j*Group*stride;
#define HDGROUP_TWO_SUM_STORE(R) do { \
        const auto y=diagonalProduct(column[R],diagonals+(R)*4); \
        auto* out=destination+(R)*stride; \
        if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),y); \
        else _mm512_storeu_si512(out,y); \
    } while(false)
    HDGROUP_TWO_SUM_STORE(0);HDGROUP_TWO_SUM_STORE(1);HDGROUP_TWO_SUM_STORE(2);HDGROUP_TWO_SUM_STORE(3);
    HDGROUP_TWO_SUM_STORE(4);HDGROUP_TWO_SUM_STORE(5);HDGROUP_TWO_SUM_STORE(6);HDGROUP_TWO_SUM_STORE(7);
    HDGROUP_TWO_SUM_STORE(8);HDGROUP_TWO_SUM_STORE(9);HDGROUP_TWO_SUM_STORE(10);HDGROUP_TWO_SUM_STORE(11);
    HDGROUP_TWO_SUM_STORE(12);HDGROUP_TWO_SUM_STORE(13);HDGROUP_TWO_SUM_STORE(14);HDGROUP_TWO_SUM_STORE(15);
#undef HDGROUP_TWO_SUM_STORE
#else
    c8Unroll(std::make_index_sequence<Group>{},[&](auto r){
        const auto y=[&](){if constexpr(HDGROUP_XOR_MIXER==1||HDGROUP_TOEPLITZ||(HDGROUP_FIELD32&&Group==2))return column[r];
            else return diagonalProduct(column[r],after+(j*Group+r)*4);}();
        auto* out=planes+(j*Group+r)*stride;
        if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),y);
        else _mm512_storeu_si512(out,y);
    });
#endif
}

#if HDGROUP_CROSS_LAYERS
// Every key pack contains four independent nonzero ratios, one per group.
// The second pack is 2*alpha, prepared once to avoid a dependent GFNI chain.
inline constexpr unsigned crossDelta1(unsigned column){
    if constexpr(HDGROUP_CROSS_VARYING&&HDGROUP_CROSS_LAYERS==1)return 8+(column>>1);
    else return 1;
}
inline constexpr unsigned crossDelta2(unsigned column){
    if constexpr(HDGROUP_CROSS_VARYING){
        constexpr unsigned deltas[4]={2,4,8,14};return deltas[column>>2];
    }else return 2;
}
inline constexpr unsigned crossPairIndex(unsigned layer,unsigned column,unsigned row){
    return layer==0?(column>>1)*16+row:(column>>2)*32+(column&1)*16+row;
}

template<unsigned Local,unsigned Absolute>
HDGROUP_INLINE void crossLoadMix(const __m512i* values,__m512i* column,const Coefficient* before){
    c8Unroll(std::make_index_sequence<16>{},[&](auto r){
        column[r]=values[r*8+Local];
        if constexpr(HDGROUP_NORMALIZE){
            if constexpr(r!=0)column[r]=diagonalProduct(column[r],before+(Absolute*15+r-1)*4);
        }else column[r]=diagonalProduct(column[r],before+(Absolute*16+r)*4);
    });
    twoSum16Column(column);
}
template<unsigned Column,unsigned Row,bool Stream>
HDGROUP_INLINE void crossFinalStore(__m512i value,Block* planes,std::size_t stride,const Coefficient* after){
    const auto y=diagonalProduct(value,after+(Column*16+Row)*4);
    auto* out=planes+(Column*16+Row)*stride;
    if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),y);
    else _mm512_storeu_si512(out,y);
}
template<unsigned Layer,unsigned Column,unsigned Row>
HDGROUP_INLINE void crossPair(__m512i a,__m512i b,const Coefficient* before,__m512i& low,__m512i& high){
    constexpr unsigned key=beforeBaseRegions<16>+2*(Layer*128+crossPairIndex(Layer,Column,Row));
    // Independent products, followed by the two MDS additions.
    low=_mm512_xor_si512(a,diagonalProduct(b,before+key*4));
    high=_mm512_xor_si512(a,diagonalProduct(b,before+(key+1)*4));
}

template<unsigned Offset,unsigned Local,bool Stream>
HDGROUP_INLINE void crossBlock2(__m512i* values,Block* planes,std::size_t stride,
                                const Coefficient* before,const Coefficient* after){
    constexpr unsigned column0=Offset+Local,delta=crossDelta1(column0);
    alignas(64) __m512i column[16];
    crossLoadMix<Local,column0>(values,column,before);
    c8Unroll(std::make_index_sequence<16>{},[&](auto r){values[r*8+Local]=column[r];});
    crossLoadMix<Local+1,column0+1>(values,column,before);
    // Keep this drain explicitly unrolled: no large closure may be outlined.
#define HDGROUP_CROSS_PAIR_STORE(R) do { \
        __m512i low,high; \
        crossPair<0,column0,R>(values[(R)*8+Local],column[(R)^delta],before,low,high); \
        crossFinalStore<column0,R,Stream>(low,planes,stride,after); \
        crossFinalStore<column0+1,(R)^delta,Stream>(high,planes,stride,after); \
    } while(false)
    HDGROUP_CROSS_PAIR_STORE(0);HDGROUP_CROSS_PAIR_STORE(1);HDGROUP_CROSS_PAIR_STORE(2);HDGROUP_CROSS_PAIR_STORE(3);
    HDGROUP_CROSS_PAIR_STORE(4);HDGROUP_CROSS_PAIR_STORE(5);HDGROUP_CROSS_PAIR_STORE(6);HDGROUP_CROSS_PAIR_STORE(7);
    HDGROUP_CROSS_PAIR_STORE(8);HDGROUP_CROSS_PAIR_STORE(9);HDGROUP_CROSS_PAIR_STORE(10);HDGROUP_CROSS_PAIR_STORE(11);
    HDGROUP_CROSS_PAIR_STORE(12);HDGROUP_CROSS_PAIR_STORE(13);HDGROUP_CROSS_PAIR_STORE(14);HDGROUP_CROSS_PAIR_STORE(15);
#undef HDGROUP_CROSS_PAIR_STORE
}

template<unsigned Offset,unsigned Local,bool Stream>
HDGROUP_INLINE void crossBlock4(__m512i* values,Block* planes,std::size_t stride,
                                const Coefficient* before,const Coefficient* after){
    constexpr unsigned column0=Offset+Local,delta1=crossDelta1(column0),delta2=crossDelta2(column0);
    static_assert(crossDelta1(column0+2)==delta1,"Two-layer first matching must be constant within each quartet");
    alignas(64) __m512i column[16];
    crossLoadMix<Local,column0>(values,column,before);
    c8Unroll(std::make_index_sequence<16>{},[&](auto r){values[r*8+Local]=column[r];});
    crossLoadMix<Local+1,column0+1>(values,column,before);
#define HDGROUP_CROSS_SAVE_PAIR(R) do { \
        __m512i low,high; \
        crossPair<0,column0,R>(values[(R)*8+Local],column[(R)^delta1],before,low,high); \
        values[(R)*8+Local]=low;values[((R)^delta1)*8+Local+1]=high; \
    } while(false)
    HDGROUP_CROSS_SAVE_PAIR(0);HDGROUP_CROSS_SAVE_PAIR(1);HDGROUP_CROSS_SAVE_PAIR(2);HDGROUP_CROSS_SAVE_PAIR(3);
    HDGROUP_CROSS_SAVE_PAIR(4);HDGROUP_CROSS_SAVE_PAIR(5);HDGROUP_CROSS_SAVE_PAIR(6);HDGROUP_CROSS_SAVE_PAIR(7);
    HDGROUP_CROSS_SAVE_PAIR(8);HDGROUP_CROSS_SAVE_PAIR(9);HDGROUP_CROSS_SAVE_PAIR(10);HDGROUP_CROSS_SAVE_PAIR(11);
    HDGROUP_CROSS_SAVE_PAIR(12);HDGROUP_CROSS_SAVE_PAIR(13);HDGROUP_CROSS_SAVE_PAIR(14);HDGROUP_CROSS_SAVE_PAIR(15);
#undef HDGROUP_CROSS_SAVE_PAIR
    crossLoadMix<Local+2,column0+2>(values,column,before);
    c8Unroll(std::make_index_sequence<16>{},[&](auto r){values[r*8+Local+2]=column[r];});
    crossLoadMix<Local+3,column0+3>(values,column,before);
#define HDGROUP_CROSS_QUARTET_STORE(R) do { \
        __m512i pairLow,pairHigh,out0,out2,out1,out3; \
        crossPair<0,column0+2,R>(values[(R)*8+Local+2],column[(R)^delta1],before,pairLow,pairHigh); \
        crossPair<1,column0,(R)^delta2>(values[((R)^delta2)*8+Local],pairLow,before,out0,out2); \
        crossFinalStore<column0,(R)^delta2,Stream>(out0,planes,stride,after); \
        crossFinalStore<column0+2,R,Stream>(out2,planes,stride,after); \
        crossPair<1,column0+1,(R)^delta1^delta2>(values[((R)^delta1^delta2)*8+Local+1],pairHigh,before,out1,out3); \
        crossFinalStore<column0+1,(R)^delta1^delta2,Stream>(out1,planes,stride,after); \
        crossFinalStore<column0+3,(R)^delta1,Stream>(out3,planes,stride,after); \
    } while(false)
    HDGROUP_CROSS_QUARTET_STORE(0);HDGROUP_CROSS_QUARTET_STORE(1);HDGROUP_CROSS_QUARTET_STORE(2);HDGROUP_CROSS_QUARTET_STORE(3);
    HDGROUP_CROSS_QUARTET_STORE(4);HDGROUP_CROSS_QUARTET_STORE(5);HDGROUP_CROSS_QUARTET_STORE(6);HDGROUP_CROSS_QUARTET_STORE(7);
    HDGROUP_CROSS_QUARTET_STORE(8);HDGROUP_CROSS_QUARTET_STORE(9);HDGROUP_CROSS_QUARTET_STORE(10);HDGROUP_CROSS_QUARTET_STORE(11);
    HDGROUP_CROSS_QUARTET_STORE(12);HDGROUP_CROSS_QUARTET_STORE(13);HDGROUP_CROSS_QUARTET_STORE(14);HDGROUP_CROSS_QUARTET_STORE(15);
#undef HDGROUP_CROSS_QUARTET_STORE
}
#endif

template<unsigned Group,unsigned Offset,bool Stream=true>
HD35_NOINLINE void coset(const Block* source,Block* planes,std::size_t stride,
                        const Coefficient* before,const Coefficient* after,std::size_t groups){
    constexpr unsigned dim=128/Group;
    alignas(64) __m512i values[128];
    // One ZMM carries four independent groups, not four uncoupled constituent
    // rows. Each column below mixes Group different ZMMs across those rows.
    for(unsigned r=0;r!=Group;++r){
#if HDGROUP_SOURCE_PREFIX
        if constexpr(Group==4){
            sourcePrefix32<Offset>(source+4*r*dim,values+r*dim);
        }else{
#endif
        c8Unroll(std::make_index_sequence<dim>{},[&](auto j){values[r*dim+j]=_mm512_loadu_si512(source+4*(r*dim+j));});
        if constexpr(HDGROUP_FUSE_LAST&&Group<=8)rsPrefix<dim,Offset>(values+r*dim);
        else rsFft<dim,Offset>(values+r*dim);
#if HDGROUP_SOURCE_PREFIX
        }
#endif
    }
#if HDGROUP_SCHEDULED_LAMBDA
    static_assert(Group==16);
    scheduledLambdaStore<Offset,0,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,1,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,2,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,3,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,4,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,5,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,6,Stream>(values,planes,stride,before,after);
    scheduledLambdaStore<Offset,7,Stream>(values,planes,stride,before,after);
#elif HDGROUP_CROSS_LAYERS
    static_assert(Group==16);
    if constexpr(HDGROUP_CROSS_LAYERS==1){
        crossBlock2<Offset,0,Stream>(values,planes,stride,before,after);
        crossBlock2<Offset,2,Stream>(values,planes,stride,before,after);
        crossBlock2<Offset,4,Stream>(values,planes,stride,before,after);
        crossBlock2<Offset,6,Stream>(values,planes,stride,before,after);
    }else{
        crossBlock4<Offset,0,Stream>(values,planes,stride,before,after);
        crossBlock4<Offset,4,Stream>(values,planes,stride,before,after);
    }
#else
    if constexpr(HDGROUP_FUSE_LAST&&Group<=8){
        for(unsigned j=0;j!=dim;j+=2){
            alignas(64) __m512i even[Group],odd[Group];
            const auto coefficient=_mm512_set1_epi8(static_cast<char>(Offset+j));
            c8Unroll(std::make_index_sequence<Group>{},[&](auto r){
                const auto a=values[r*dim+j],b=values[r*dim+j+1];
                even[r]=_mm512_xor_si512(a,_mm512_gf2p8mul_epi8(b,coefficient));
                odd[r]=_mm512_xor_si512(even[r],b);
            });
            mixStore<Group,Stream>(even,Offset+j,planes,stride,before,after);
            mixStore<Group,Stream>(odd,Offset+j+1,planes,stride,before,after);
        }
    }else{
        for(unsigned j=0;j!=dim;++j){
            alignas(64) __m512i column[Group];
            c8Unroll(std::make_index_sequence<Group>{},[&](auto r){column[r]=values[r*dim+j];});
            mixStore<Group,Stream>(column,Offset+j,planes,stride,before,after);
        }
    }
#endif
}

#if HDGROUP_COSET_PIPELINE && !HDGROUP_LONG_RS
// Produce one complete coset across all groups. Keep the full planes base:
// mixStore already uses absolute region indices for both coset offsets.
template<unsigned Group,unsigned Half,bool Stream=true,unsigned TileGroups=64>
HD35_NOINLINE void forwardCoset(std::size_t k,std::size_t stride,const Block* input,Block* planes,
                               const Coefficient* before,const Coefficient* after){
    static_assert(Group==2||Group==4||Group==8||Group==16);
    static_assert(Half<2);
    static_assert(TileGroups>0&&TileGroups%4==0);
    const auto groups=k/128;
    if(!k||k%512||stride%4||stride<groups)throw std::invalid_argument("grouped outer dimensions");
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(planes)&63U){forwardCoset<Group,Half,false,TileGroups>(k,stride,input,planes,before,after);return;}
    }
    for(std::size_t first=0;first<groups;first+=TileGroups){
        const auto last=std::min(first+TileGroups,groups);
        for(auto r=first;r<last;r+=4)
            coset<Group,Half*128/Group,Stream>(input+r*128,planes+r,stride,
                before+r*beforeRegions<Group>,(HDGROUP_TOEPLITZ||HDGROUP_XOR_MIXER==1)?after:after+r*256,groups);
    }
    for(unsigned c=128*Half;c!=128*(Half+1);++c)
        std::memset(planes+c*stride+groups,0,(stride-groups)*sizeof(Block));
    if constexpr(Stream)_mm_sfence();
}
#endif

template<unsigned Group,bool Stream=true,unsigned TileGroups=64>
HD35_NOINLINE void forward(std::size_t k,std::size_t stride,const Block* input,Block* planes,
                          const Coefficient* before,const Coefficient* after){
#if HDGROUP_LONG_RS
    static_assert(Group==1,"Long RS uses one RS[256,128] constituent per group");
#else
    static_assert(Group==2||Group==4||Group==8||Group==16);
#endif
    static_assert(TileGroups>0&&TileGroups%4==0);
    const auto groups=k/128;
    if(!k||k%512||stride%4||stride<groups)throw std::invalid_argument("grouped outer dimensions");
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(planes)&63U){forward<Group,false,TileGroups>(k,stride,input,planes,before,after);return;}
    }
    for(std::size_t first=0;first<groups;first+=TileGroups){
        const auto last=std::min(first+TileGroups,groups);
#if HDGROUP_LONG_RS
        for(auto r=first;r<last;r+=4)longCoset128<0,Stream>(input+r*128,planes+r,stride,after+r*256);
        for(auto r=first;r<last;r+=4)longCoset128<128,Stream>(input+r*128,planes+128*stride+r,stride,after+r*256+512);
#else
        for(auto r=first;r<last;r+=4)coset<Group,0,Stream>(input+r*128,planes+r,stride,before+r*beforeRegions<Group>,(HDGROUP_TOEPLITZ||HDGROUP_XOR_MIXER==1)?after:after+r*256,groups);
        for(auto r=first;r<last;r+=4)coset<Group,128/Group,Stream>(input+r*128,planes+r,stride,before+r*beforeRegions<Group>,(HDGROUP_TOEPLITZ||HDGROUP_XOR_MIXER==1)?after:after+r*256,groups);
#endif
    }
    for(unsigned c=0;c!=256;++c)std::memset(planes+c*stride+groups,0,(stride-groups)*sizeof(Block));
    if constexpr(Stream)_mm_sfence();
}
#undef HDGROUP_INLINE
} // namespace spin::research::char2::hdgroup
