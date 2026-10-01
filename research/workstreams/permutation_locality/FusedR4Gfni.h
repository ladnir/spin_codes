#pragma once
// Isolated exact-map GFNI experiment. Coordinate padding changes storage only.
#include "Inner.h"

namespace spin::research::gfni_r4 {
using namespace spin::detail::kernel;
template<unsigned S> struct Row {
    static_assert(S==16 || S==19);
    static constexpr unsigned groups=(S+7)/8;
    std::uint64_t reduce[groups]{},expand[groups]{};
};
template<unsigned S> struct alignas(64) Packed {
    // One byte per payload bit, eight state coordinates per byte. Each group
    // contains all 128 payload positions in two ZMMs; unused coordinates are 0.
    __m512i v[2*Row<S>::groups];
};
static_assert(sizeof(Row<19>)==48 && sizeof(Row<16>)==32);
static_assert(sizeof(Packed<19>)==384 && sizeof(Packed<16>)==256);

// Original forward map M=T4 T3 T2 T1=I+sum ai vi^T,
// ai=T4 ... T(i+1) ui. Transpose reduces by ai and expands by vi.
template<unsigned S> inline Row<S> makeRow(const unsigned* masks,bool transpose=true) {
    unsigned a[4],v[4];
    for(unsigned i=0;i<4;++i) {
        a[i]=masks[2*i];v[i]=masks[2*i+1];
        for(unsigned r=i+1;r<4;++r)
            if(std::popcount(a[i]&masks[2*r+1])&1)a[i]^=masks[2*r];
    }
    Row<S> row;
    for(unsigned g=0;g<Row<S>::groups;++g) {
        for(unsigned i=0;i<4;++i)
            row.reduce[g]|=std::uint64_t(((transpose?a[i]:v[i])>>(8*g))&255)<<(8*(7-i));
        for(unsigned j=0;j<8 && 8*g+j<S;++j) {
            unsigned bits=0;
            for(unsigned i=0;i<4;++i)bits|=(((transpose?v[i]:a[i])>>(8*g+j))&1)<<i;
            row.expand[g]|=std::uint64_t(bits)<<(8*(7-j));
        }
    }
    return row;
}

// Explicit lane-local byte transpose. Read every input before writing any
// reference, so calling this with the eight distinct arguments reversed is safe.
SPIN_FORCEINLINE void byteTranspose(__m128i& x0,__m128i& x1,__m128i& x2,__m128i& x3,
    __m128i& x4,__m128i& x5,__m128i& x6,__m128i& x7) {
    const auto t0=_mm_unpacklo_epi8(x0,x1),t1=_mm_unpackhi_epi8(x0,x1);
    const auto t2=_mm_unpacklo_epi8(x2,x3),t3=_mm_unpackhi_epi8(x2,x3);
    const auto t4=_mm_unpacklo_epi8(x4,x5),t5=_mm_unpackhi_epi8(x4,x5);
    const auto t6=_mm_unpacklo_epi8(x6,x7),t7=_mm_unpackhi_epi8(x6,x7);
    const auto u0=_mm_unpacklo_epi16(t0,t2),u1=_mm_unpackhi_epi16(t0,t2);
    const auto u2=_mm_unpacklo_epi16(t1,t3),u3=_mm_unpackhi_epi16(t1,t3);
    const auto u4=_mm_unpacklo_epi16(t4,t6),u5=_mm_unpackhi_epi16(t4,t6);
    const auto u6=_mm_unpacklo_epi16(t5,t7),u7=_mm_unpackhi_epi16(t5,t7);
    x0=_mm_unpacklo_epi32(u0,u4);x1=_mm_unpackhi_epi32(u0,u4);
    x2=_mm_unpacklo_epi32(u1,u5);x3=_mm_unpackhi_epi32(u1,u5);
    x4=_mm_unpacklo_epi32(u2,u6);x5=_mm_unpackhi_epi32(u2,u6);
    x6=_mm_unpacklo_epi32(u3,u7);x7=_mm_unpackhi_epi32(u3,u7);
}
SPIN_FORCEINLINE __m512i join(__m128i a,__m128i b,__m128i c,__m128i d) {
    auto x=_mm512_castsi128_si512(a);
    x=_mm512_inserti32x4(x,b,1);x=_mm512_inserti32x4(x,c,2);
    return _mm512_inserti32x4(x,d,3);
}
template<unsigned S,unsigned J> SPIN_FORCEINLINE __m128i coordinate(const __m128i* state) {
    if constexpr(J<S)return state[J];
    else return _mm_setzero_si128();
}
template<unsigned S,unsigned G=0> SPIN_FORCEINLINE void pack(const __m128i* state,Packed<S>& out) {
    auto x0=coordinate<S,8*G>(state),x1=coordinate<S,8*G+1>(state);
    auto x2=coordinate<S,8*G+2>(state),x3=coordinate<S,8*G+3>(state);
    auto x4=coordinate<S,8*G+4>(state),x5=coordinate<S,8*G+5>(state);
    auto x6=coordinate<S,8*G+6>(state),x7=coordinate<S,8*G+7>(state);
    byteTranspose(x7,x6,x5,x4,x3,x2,x1,x0);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    out.v[2*G]=_mm512_gf2p8affine_epi64_epi8(basis,join(x7,x6,x5,x4),0);
    out.v[2*G+1]=_mm512_gf2p8affine_epi64_epi8(basis,join(x3,x2,x1,x0),0);
    if constexpr(G+1<Row<S>::groups)pack<S,G+1>(state,out);
}
template<unsigned S,unsigned J> SPIN_FORCEINLINE void storeCoordinate(__m128i* out,__m128i value) {
    if constexpr(J<S)out[J]=value;
}
// pack/unpack require distinct source and destination storage. The final group
// loads/stores only its three real coordinates for S=19, never a twentieth word.
template<unsigned S,unsigned G=0> SPIN_FORCEINLINE void unpack(const Packed<S>& state,__m128i* out) {
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto lo=_mm512_gf2p8affine_epi64_epi8(basis,state.v[2*G],0);
    const auto hi=_mm512_gf2p8affine_epi64_epi8(basis,state.v[2*G+1],0);
    auto x0=_mm512_castsi512_si128(lo),x1=_mm512_extracti32x4_epi32(lo,1);
    auto x2=_mm512_extracti32x4_epi32(lo,2),x3=_mm512_extracti32x4_epi32(lo,3);
    auto x4=_mm512_castsi512_si128(hi),x5=_mm512_extracti32x4_epi32(hi,1);
    auto x6=_mm512_extracti32x4_epi32(hi,2),x7=_mm512_extracti32x4_epi32(hi,3);
    byteTranspose(x0,x1,x2,x3,x4,x5,x6,x7);
    storeCoordinate<S,8*G>(out,_mm_unpacklo_epi8(x0,x4));
    storeCoordinate<S,8*G+1>(out,_mm_unpackhi_epi8(x0,x4));
    storeCoordinate<S,8*G+2>(out,_mm_unpacklo_epi8(x1,x5));
    storeCoordinate<S,8*G+3>(out,_mm_unpackhi_epi8(x1,x5));
    storeCoordinate<S,8*G+4>(out,_mm_unpacklo_epi8(x2,x6));
    storeCoordinate<S,8*G+5>(out,_mm_unpackhi_epi8(x2,x6));
    storeCoordinate<S,8*G+6>(out,_mm_unpacklo_epi8(x3,x7));
    storeCoordinate<S,8*G+7>(out,_mm_unpackhi_epi8(x3,x7));
    if constexpr(G+1<Row<S>::groups)unpack<S,G+1>(state,out);
}
template<unsigned S,bool Add,unsigned G=0> SPIN_FORCEINLINE void expand(
    Packed<S>& state,const Row<S>& row,__m512i d0,__m512i d1,const Packed<S>* syndrome) {
    const auto matrix=_mm512_set1_epi64(row.expand[G]);
    const auto a=_mm512_gf2p8affine_epi64_epi8(d0,matrix,0);
    const auto b=_mm512_gf2p8affine_epi64_epi8(d1,matrix,0);
    if constexpr(Add) {
        state.v[2*G]=_mm512_ternarylogic_epi64(state.v[2*G],a,syndrome->v[2*G],0x96);
        state.v[2*G+1]=_mm512_ternarylogic_epi64(state.v[2*G+1],b,syndrome->v[2*G+1],0x96);
    } else {
        state.v[2*G]=_mm512_xor_si512(state.v[2*G],a);
        state.v[2*G+1]=_mm512_xor_si512(state.v[2*G+1],b);
    }
    if constexpr(G+1<Row<S>::groups)expand<S,Add,G+1>(state,row,d0,d1,syndrome);
}
template<unsigned S,bool Add=false> SPIN_FORCEINLINE void step(
    Packed<S>& state,const Row<S>& row,const Packed<S>* syndrome=nullptr) {
    const auto m0=_mm512_set1_epi64(row.reduce[0]);
    const auto m1=_mm512_set1_epi64(row.reduce[1]);
    auto d0=_mm512_gf2p8affine_epi64_epi8(state.v[0],m0,0);
    auto d1=_mm512_gf2p8affine_epi64_epi8(state.v[1],m0,0);
    const auto e0=_mm512_gf2p8affine_epi64_epi8(state.v[2],m1,0);
    const auto e1=_mm512_gf2p8affine_epi64_epi8(state.v[3],m1,0);
    if constexpr(S==19) {
        const auto m2=_mm512_set1_epi64(row.reduce[2]);
        d0=_mm512_ternarylogic_epi64(d0,e0,_mm512_gf2p8affine_epi64_epi8(state.v[4],m2,0),0x96);
        d1=_mm512_ternarylogic_epi64(d1,e1,_mm512_gf2p8affine_epi64_epi8(state.v[5],m2,0),0x96);
    } else {d0=_mm512_xor_si512(d0,e0);d1=_mm512_xor_si512(d1,e1);}
    // Every reduction reads the entering state. Mutation begins only here.
    expand<S,Add>(state,row,d0,d1,syndrome);
}

// Optional S=16 core only: four 8x8 blocks represent ANY fixed binary 16x16
// matrix. This is not an S=16 encoder or a distance/distribution claim.
struct Dense16Row {std::uint64_t matrix[4]{};};
inline Dense16Row makeDense16(const unsigned* rows) {
    Dense16Row result;
    for(unsigned out=0;out<2;++out)for(unsigned in=0;in<2;++in)
        for(unsigned j=0;j<8;++j)
            result.matrix[2*out+in]|=std::uint64_t((rows[8*out+j]>>(8*in))&255)<<(8*(7-j));
    return result;
}
template<bool Add=false> SPIN_FORCEINLINE void denseStep16(
    Packed<16>& state,const Dense16Row& row,const Packed<16>* syndrome=nullptr) {
    const auto m00=_mm512_set1_epi64(row.matrix[0]),m01=_mm512_set1_epi64(row.matrix[1]);
    const auto m10=_mm512_set1_epi64(row.matrix[2]),m11=_mm512_set1_epi64(row.matrix[3]);
    const auto a0=_mm512_gf2p8affine_epi64_epi8(state.v[0],m00,0);
    const auto b0=_mm512_gf2p8affine_epi64_epi8(state.v[2],m01,0);
    const auto a1=_mm512_gf2p8affine_epi64_epi8(state.v[1],m00,0);
    const auto b1=_mm512_gf2p8affine_epi64_epi8(state.v[3],m01,0);
    const auto a2=_mm512_gf2p8affine_epi64_epi8(state.v[0],m10,0);
    const auto b2=_mm512_gf2p8affine_epi64_epi8(state.v[2],m11,0);
    const auto a3=_mm512_gf2p8affine_epi64_epi8(state.v[1],m10,0);
    const auto b3=_mm512_gf2p8affine_epi64_epi8(state.v[3],m11,0);
    if constexpr(Add) {
        state.v[0]=_mm512_ternarylogic_epi64(a0,b0,syndrome->v[0],0x96);
        state.v[1]=_mm512_ternarylogic_epi64(a1,b1,syndrome->v[1],0x96);
        state.v[2]=_mm512_ternarylogic_epi64(a2,b2,syndrome->v[2],0x96);
        state.v[3]=_mm512_ternarylogic_epi64(a3,b3,syndrome->v[3],0x96);
    } else {
        state.v[0]=_mm512_xor_si512(a0,b0);state.v[1]=_mm512_xor_si512(a1,b1);
        state.v[2]=_mm512_xor_si512(a2,b2);state.v[3]=_mm512_xor_si512(a3,b3);
    }
}
}
