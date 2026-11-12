// Forward RS8 outer: the exact binary transpose of the tuned transposed outer
// in Fast.cpp. Each linear step is replaced by its adjoint: the message
// unpack inverts storePair, the parity network is transposed, each GF(2^32)
// column map applies the transposed flat9 tower circuit with adjoint byte
// multiplications (GFNI affine matrices prepared at setup), and the column
// unpack inverts packEight. Output is lanes contiguous scratch planes laid
// out exactly as forwardScalar's, so the s22 forward inner gathers routes.
#include "Rs8.h"
#include "Field.h"
#include "../packet/PacketInnerFast.h"
#include <utility>

namespace spin::detail::rs8 {
namespace {
using packet::fast::unpackIndex;

constexpr std::uint64_t transposeAffine(std::uint64_t m) {
    std::uint64_t result=0;
    for(unsigned i=0;i<8;++i) for(unsigned j=0;j<8;++j)
        result|=((m>>(8*(7-i)+j))&1ULL)<<(8*(7-j)+i);
    return result;
}
// Transposes of the parity network's adjoint nibble multiplications.
constexpr std::uint64_t m6t=transposeAffine(0x060c0b0560c0b050ULL);
constexpr std::uint64_t m8t=transposeAffine(0x0803060c803060c0ULL);

// Load every physical cache line once even for wide records; lane
// separation is register-local.
template<unsigned L> struct Quads { __m512i v[L]; };
template<unsigned L> static SPIN_FORCEINLINE Quads<L> loadQuads(const Block* input) {
    if constexpr(L==1) return {{_mm512_loadu_si512(input)}};
    else if constexpr(L==2) {
        const auto a=_mm512_loadu_si512(input),b=_mm512_loadu_si512(input+4);
        return {{_mm512_shuffle_i32x4(a,b,0x88),_mm512_shuffle_i32x4(a,b,0xdd)}};
    } else {
        const auto a=_mm512_loadu_si512(input),b=_mm512_loadu_si512(input+4);
        const auto c=_mm512_loadu_si512(input+8),d=_mm512_loadu_si512(input+12);
        const auto ab0=_mm512_shuffle_i32x4(a,b,0x44),ab1=_mm512_shuffle_i32x4(a,b,0xee);
        const auto cd0=_mm512_shuffle_i32x4(c,d,0x44),cd1=_mm512_shuffle_i32x4(c,d,0xee);
        return {{_mm512_shuffle_i32x4(ab0,cd0,0x88),_mm512_shuffle_i32x4(ab0,cd0,0xdd),
                 _mm512_shuffle_i32x4(ab1,cd1,0x88),_mm512_shuffle_i32x4(ab1,cd1,0xdd)}};
    }
}
template<unsigned L> struct alignas(64) Packed { __m512i v[L][64]; };

// Inverse of storePair<Symbol=S,Plane=P>: message blocks 16P+4S.. and
// 16(P+1)+4S.. become packed systematic column S at planes P and P+1.
template<unsigned S,unsigned P>
static SPIN_FORCEINLINE void packMessage(__m512i x,__m512i y,__m512i* packed) {
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    packed[8*S+P]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(x,_mm512_setr_epi64(0x0010203040506070ULL,0x0111213141516171ULL,0x0212223242526272ULL,0x0313233343536373ULL,0x0414243444546474ULL,0x0515253545556575ULL,0x0616263646566676ULL,0x0717273747576777ULL),y),0);
    packed[8*S+P+1]=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(x,_mm512_setr_epi64(0x0818283848586878ULL,0x0919293949596979ULL,0x0a1a2a3a4a5a6a7aULL,0x0b1b2b3b4b5b6b7bULL,0x0c1c2c3c4c5c6c7cULL,0x0d1d2d3d4d5d6d7dULL,0x0e1e2e3e4e5e6e7eULL,0x0f1f2f3f4f5f6f7fULL),y),0);
}
// Packs record lanes First+Lane of L-lane records into packed[Lane].
template<unsigned S,unsigned P,unsigned L,unsigned First,unsigned C,unsigned... Lane>
static SPIN_FORCEINLINE void packLanes(const Block* input,Packed<C>& packed,
    std::integer_sequence<unsigned,Lane...>) {
    const auto x=loadQuads<L>(input+L*(16*P+4*S));
    const auto y=loadQuads<L>(input+L*(16*(P+1)+4*S));
    (packMessage<S,P>(x.v[First+Lane],y.v[First+Lane],packed.v[Lane]),...);
}

// Transpose of parityValues<Plane>: systematic columns 0..3 at this plane
// produce parity columns 4..7.
template<unsigned Plane>
static SPIN_FORCEINLINE void parityForward(__m512i* packed) {
    const auto m6=_mm512_set1_epi64(std::int64_t(m6t));
    const auto m8=_mm512_set1_epi64(std::int64_t(m8t));
    auto y0=packed[0+Plane],y1=packed[8+Plane],y2=packed[16+Plane],y3=packed[24+Plane];
    // Reverse of: y0^=y1; y2^=y3; y0^=y2; y1^=y3.
    y3=_mm512_xor_si512(y3,y1);y2=_mm512_xor_si512(y2,y0);
    y3=_mm512_xor_si512(y3,y2);y1=_mm512_xor_si512(y1,y0);
    // Adjoint of y0=x0, y1=x1+A8x0, y2=x2+A6x0, y3=x3+x0+A6x0+A6x1+A8x2.
    auto x0=_mm512_ternarylogic_epi64(y0,y3,
        _mm512_gf2p8affine_epi64_epi8(y1,m8,0),0x96);
    x0=_mm512_ternarylogic_epi64(x0,
        _mm512_gf2p8affine_epi64_epi8(y2,m6,0),
        _mm512_gf2p8affine_epi64_epi8(y3,m6,0),0x96);
    auto x1=_mm512_xor_si512(y1,_mm512_gf2p8affine_epi64_epi8(y3,m6,0));
    auto x2=_mm512_xor_si512(y2,_mm512_gf2p8affine_epi64_epi8(y3,m8,0));
    auto x3=y3;
    // Reverse of: x0^=x1; x2^=x3; x0^=x2; x1^=x3.
    x3=_mm512_xor_si512(x3,x1);x2=_mm512_xor_si512(x2,x0);
    x3=_mm512_xor_si512(x3,x2);x1=_mm512_xor_si512(x1,x0);
    packed[32+Plane]=x0;packed[40+Plane]=x1;packed[48+Plane]=x2;packed[56+Plane]=x3;
}
static SPIN_FORCEINLINE void parityPlanes(__m512i* packed) {
    parityForward<0>(packed);parityForward<1>(packed);parityForward<2>(packed);parityForward<3>(packed);
    parityForward<4>(packed);parityForward<5>(packed);parityForward<6>(packed);parityForward<7>(packed);
}

// Transposed flat9: the adjoint of the column multiplication by c.
static SPIN_FORCEINLINE void flat9Adjoint(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,
    const __m512i* m) {
    const auto y0=x0,y1=x1,y2=x2,y3=x3;
    const auto y01=_mm512_xor_si512(y0,y1),y23=_mm512_xor_si512(y2,y3);
    const auto a0=_mm512_gf2p8affine_epi64_epi8(_mm512_xor_si512(y01,y23),m[0],0);
    const auto a1=_mm512_gf2p8affine_epi64_epi8(_mm512_xor_si512(y0,y2),m[1],0);
    const auto a2=_mm512_gf2p8affine_epi64_epi8(_mm512_xor_si512(y1,y3),m[2],0);
    const auto a3=_mm512_gf2p8affine_epi64_epi8(y01,m[3],0);
    const auto a4=_mm512_gf2p8affine_epi64_epi8(y0,m[4],0);
    const auto a5=_mm512_gf2p8affine_epi64_epi8(y1,m[5],0);
    const auto a6=_mm512_gf2p8affine_epi64_epi8(y23,m[6],0);
    const auto a7=_mm512_gf2p8affine_epi64_epi8(y2,m[7],0);
    const auto a8=_mm512_gf2p8affine_epi64_epi8(y3,m[8],0);
    const auto s68=_mm512_xor_si512(a6,a8),s78=_mm512_xor_si512(a7,a8);
    x0=_mm512_ternarylogic_epi64(a0,a2,s68,0x96);
    x1=_mm512_ternarylogic_epi64(a1,a2,s78,0x96);
    x2=_mm512_ternarylogic_epi64(a3,a5,s68,0x96);
    x3=_mm512_ternarylogic_epi64(a4,a5,s78,0x96);
}
// Exact inverse of packEight: eight 128-bit coordinates back to blocks.
static SPIN_FORCEINLINE void unpackPair(__m512i lo,__m512i hi,Block* output) {
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto a=_mm512_gf2p8affine_epi64_epi8(basis,lo,0);
    const auto b=_mm512_gf2p8affine_epi64_epi8(basis,hi,0);
    _mm512_store_si512(output,_mm512_permutex2var_epi8(a,unpackIndex<0>(),b));
    _mm512_store_si512(output+4,_mm512_permutex2var_epi8(a,unpackIndex<4>(),b));
}
static SPIN_FORCEINLINE void mixColumn(__m512i* v,Block* output,const __m512i* m) {
    flat9Adjoint(v[0],v[2],v[4],v[6],m);
    flat9Adjoint(v[1],v[3],v[5],v[7],m);
    unpackPair(v[0],v[1],output);unpackPair(v[2],v[3],output+8);
    unpackPair(v[4],v[5],output+16);unpackPair(v[6],v[7],output+24);
}
// Encodes record lanes First..First+C-1 of one group. Wider records may be
// split into passes so the packed state, input and output stay L1-resident.
template<unsigned L,unsigned First,unsigned C,unsigned... Lane>
static SPIN_NOINLINE void groupForward(const Block* __restrict input,Block* __restrict output,
    std::size_t plane,const std::array<std::uint64_t,9>* __restrict coeff,
    std::integer_sequence<unsigned,Lane...> lanes) {
    Packed<C> packed;
#define SPIN_RS8_PACK(S,P) packLanes<S,P,L,First,C>(input,packed,lanes)
    SPIN_RS8_PACK(0,0);SPIN_RS8_PACK(1,0);SPIN_RS8_PACK(2,0);SPIN_RS8_PACK(3,0);
    SPIN_RS8_PACK(0,2);SPIN_RS8_PACK(1,2);SPIN_RS8_PACK(2,2);SPIN_RS8_PACK(3,2);
    SPIN_RS8_PACK(0,4);SPIN_RS8_PACK(1,4);SPIN_RS8_PACK(2,4);SPIN_RS8_PACK(3,4);
    SPIN_RS8_PACK(0,6);SPIN_RS8_PACK(1,6);SPIN_RS8_PACK(2,6);SPIN_RS8_PACK(3,6);
#undef SPIN_RS8_PACK
    (parityPlanes(packed.v[Lane]),...);
    for(unsigned column=0;column<8;++column) {
        __m512i m[9];
        for(unsigned i=0;i<9;++i) m[i]=_mm512_set1_epi64(std::int64_t(coeff[column][i]));
        (mixColumn(packed.v[Lane]+8*column,output+(First+Lane)*plane+32*column,m),...);
    }
}
#ifndef SPIN_RS8_OUTER_SPLIT2
#define SPIN_RS8_OUTER_SPLIT2 0
#endif
#ifndef SPIN_RS8_OUTER_SPLIT
#define SPIN_RS8_OUTER_SPLIT 1
#endif
template<unsigned L>
static void outer(const Block* input,Block* scratch,const Plan& plan) {
    const auto plane=plan.scratchBlocks();
    for(std::size_t group=0;group<plan.groups;++group) {
        const auto* in=input+128*L*group;auto* out=scratch+groupStride*group;
        const auto* coeff=plan.adjoint.data()+8*group;
        if constexpr(L==2 && SPIN_RS8_OUTER_SPLIT2) {
            groupForward<2,0,1>(in,out,plane,coeff,std::make_integer_sequence<unsigned,1>{});
            groupForward<2,1,1>(in,out,plane,coeff,std::make_integer_sequence<unsigned,1>{});
        } else if constexpr(L==4 && SPIN_RS8_OUTER_SPLIT) {
            groupForward<4,0,2>(in,out,plane,coeff,std::make_integer_sequence<unsigned,2>{});
            groupForward<4,2,2>(in,out,plane,coeff,std::make_integer_sequence<unsigned,2>{});
        } else {
            groupForward<L,0,L>(in,out,plane,coeff,std::make_integer_sequence<unsigned,L>{});
        }
    }
}
}

void outerForwardFast(const Block* input,Block* scratch,const Plan& plan,unsigned lanes) {
    if(lanes==1) outer<1>(input,scratch,plan);
    else if(lanes==2) outer<2>(input,scratch,plan);
    else outer<4>(input,scratch,plan);
}
void forwardFast(const Block* input,Block* output,Block* scratch,const Plan& plan,
                 unsigned lanes,bool stream) {
    outerForwardFast(input,scratch,plan,lanes);
    certified22::forwardInnerFast(scratch,output,*plan.certified22Plan,lanes,stream,plan.scratchBlocks());
}
}
