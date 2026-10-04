#include "Rs8.h"
#include "Field.h"
#include "../../../src/packet/PacketLargeInner.h"


namespace spin::research::twoms::rs8 {
namespace {
using Coeff=std::array<std::uint8_t,9>;
struct PlaneValues { __m512i v[4]; };

template<unsigned Index>
static SPIN_FORCEINLINE __m512i product(__m512i value,const std::uint8_t* c) {
    return _mm512_gf2p8mul_epi8(value,_mm512_set1_epi8(char(c[Index])));
}

// The exact nine-product tower map, flattened to eliminate intermediate
// quadratic-result XORs. Both payload halves use the same sampled scalar.
static SPIN_FORCEINLINE void flat9(__m512i& x0,__m512i& x1,
    __m512i& x2,__m512i& x3,const std::uint8_t* c) {
    const auto x01=_mm512_xor_si512(x0,x1),x23=_mm512_xor_si512(x2,x3);
    const auto p0=product<0>(x0,c),p1=product<1>(x1,c),p2=product<2>(x01,c);
    const auto q0=product<3>(x2,c),q1=product<4>(x3,c),q2=product<5>(x23,c);
    const auto r0=product<6>(_mm512_xor_si512(x0,x2),c);
    const auto r1=product<7>(_mm512_xor_si512(x1,x3),c);
    const auto r2=product<8>(_mm512_xor_si512(x01,x23),c);
    const auto pq=_mm512_xor_si512(p0,q0),pr=_mm512_xor_si512(p0,r0);
    x0=_mm512_ternarylogic_epi64(pq,p1,q1,0x96);
    x1=_mm512_ternarylogic_epi64(pq,p2,q2,0x96);
    x2=_mm512_ternarylogic_epi64(pr,p1,r1,0x96);
    x3=_mm512_ternarylogic_epi64(pr,p2,r2,0x96);
}


template<unsigned Offset>
static SPIN_FORCEINLINE void packEight(const Block* input,__m512i& low,__m512i& high) {
    using detail::packet::large::packIndex;
    const auto a=_mm512_loadu_si512(input+Offset);
    const auto b=_mm512_loadu_si512(input+Offset+4);
    const auto basis=_mm512_set1_epi64(0x0102040810204080ULL);
    low=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<0>(),b),0);
    high=_mm512_gf2p8affine_epi64_epi8(basis,_mm512_permutex2var_epi8(a,packIndex<8>(),b),0);
}

static SPIN_FORCEINLINE void mix(const Block* input,__m512i* packed,const Coeff& coeff) {
    __m512i a0,a1,b0,b1,c0,c1,d0,d1;
    packEight<0>(input,a0,a1);packEight<8>(input,b0,b1);
    packEight<16>(input,c0,c1);packEight<24>(input,d0,d1);
    flat9(a0,b0,c0,d0,coeff.data());
    flat9(a1,b1,c1,d1,coeff.data());
    _mm512_store_si512(packed+0,a0);_mm512_store_si512(packed+1,a1);
    _mm512_store_si512(packed+2,b0);_mm512_store_si512(packed+3,b1);
    _mm512_store_si512(packed+4,c0);_mm512_store_si512(packed+5,c1);
    _mm512_store_si512(packed+6,d0);_mm512_store_si512(packed+7,d1);
}

static SPIN_FORCEINLINE __m512i xor3(__m512i a,__m512i b,__m512i c) {
    return _mm512_ternarylogic_epi64(a,b,c,0x96);
}

template<unsigned Plane>
static SPIN_FORCEINLINE PlaneValues parityValues(const __m512i* packed) {
    static_assert(Plane<8);
    // Adjoint nibble multiplication, independently in the two nibbles.
    const auto m6=_mm512_set1_epi64(0x060c0b0560c0b050ULL);
    const auto m8=_mm512_set1_epi64(0x0803060c803060c0ULL);
    auto x0=_mm512_load_si512(packed+32+Plane);
    auto x1=_mm512_load_si512(packed+40+Plane);
    auto x2=_mm512_load_si512(packed+48+Plane);
    auto x3=_mm512_load_si512(packed+56+Plane);
    x0=_mm512_xor_si512(x0,x1);x2=_mm512_xor_si512(x2,x3);
    x0=_mm512_xor_si512(x0,x2);x1=_mm512_xor_si512(x1,x3);
    const auto t6=_mm512_gf2p8affine_epi64_epi8(x0,m6,0);
    const auto t8=_mm512_gf2p8affine_epi64_epi8(x0,m8,0);
    auto y0=x0;
    auto y1=_mm512_xor_si512(x1,t8);
    auto y2=_mm512_xor_si512(x2,t6);
    auto y3=xor3(x3,t6,x0);
    y3=xor3(y3,_mm512_gf2p8affine_epi64_epi8(x1,m6,0),
                  _mm512_gf2p8affine_epi64_epi8(x2,m8,0));
    y0=_mm512_xor_si512(y0,y1);y2=_mm512_xor_si512(y2,y3);
    y0=_mm512_xor_si512(y0,y2);y1=_mm512_xor_si512(y1,y3);
    return {{
        _mm512_xor_si512(_mm512_load_si512(packed+0+Plane),y0),
        _mm512_xor_si512(_mm512_load_si512(packed+8+Plane),y1),
        _mm512_xor_si512(_mm512_load_si512(packed+16+Plane),y2),
        _mm512_xor_si512(_mm512_load_si512(packed+24+Plane),y3)}};
}

template<bool Stream,unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storePair(const PlaneValues& a,const PlaneValues& b,Block* out) {
    static_assert(Symbol<4 && (Plane==0 || Plane==2 || Plane==4 || Plane==6));
    const auto i0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
    const auto i1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x0=_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0);
    const auto x1=_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0);
    const auto lo=_mm512_permutex2var_epi8(x0,i0,x1);
    const auto hi=_mm512_permutex2var_epi8(x0,i1,x1);
    if constexpr(Stream) {
        _mm512_stream_si512(reinterpret_cast<__m512i*>(out+16*Plane+4*Symbol),lo);
        _mm512_stream_si512(reinterpret_cast<__m512i*>(out+16*(Plane+1)+4*Symbol),hi);
    } else {
        _mm512_storeu_si512(out+16*Plane+4*Symbol,lo);
        _mm512_storeu_si512(out+16*(Plane+1)+4*Symbol,hi);
    }
}

template<bool Stream,unsigned Plane>
static SPIN_FORCEINLINE void finishPair(const __m512i* packed,Block* output) {
    const auto a=parityValues<Plane>(packed);
    const auto b=parityValues<Plane+1>(packed);
    storePair<Stream,0,Plane>(a,b,output);storePair<Stream,1,Plane>(a,b,output);
    storePair<Stream,2,Plane>(a,b,output);storePair<Stream,3,Plane>(a,b,output);
}

template<bool Stream>
static SPIN_NOINLINE void group(const Block* __restrict in,Block* __restrict out,
    const Coeff* __restrict coeff) {
    alignas(64) __m512i packed[64];
    mix(in+32*0,packed+8*0,coeff[0]);mix(in+32*1,packed+8*1,coeff[1]);
    mix(in+32*2,packed+8*2,coeff[2]);mix(in+32*3,packed+8*3,coeff[3]);
    mix(in+32*4,packed+8*4,coeff[4]);mix(in+32*5,packed+8*5,coeff[5]);
    mix(in+32*6,packed+8*6,coeff[6]);mix(in+32*7,packed+8*7,coeff[7]);
    finishPair<Stream,0>(packed,out);finishPair<Stream,2>(packed,out);
    finishPair<Stream,4>(packed,out);finishPair<Stream,6>(packed,out);
}

template<bool Stream>
static void outer(const Block* in,Block* out,const Plan& plan,const Tables32& tables) {
    for(std::size_t g=0;g<plan.groups;++g)
        group<Stream>(in+detail::certified22::groupStride*g,out+128*g,tables.field.data()+8*g);
    if constexpr(Stream) _mm_sfence();
}
}

void outerFast(const Block* in,Block* out,const Plan& plan,const Tables32& tables) {
    if(plan.k>262144 && !(reinterpret_cast<std::uintptr_t>(out)&63U))
        outer<true>(in,out,plan,tables);
    else
        outer<false>(in,out,plan,tables);
}

void outerCached(const Block* in,Block* out,const Plan& plan,const Tables32& tables) {
    outer<false>(in,out,plan,tables);
}

bool selfTest() noexcept {
    constexpr std::array<std::uint32_t,8> scalars{
        1U,2U,0x100U,0x10000U,0x80000000U,0xffffffffU,0x738abcd1U,0x39870f12U};
    std::uint32_t random=0xc4397af1U;
    alignas(64) std::uint8_t input[4][64],output[4][64];
    for(const auto scalar:scalars) {
        const auto c=field::coefficients(scalar);
        for(unsigned batch=0;batch<17;++batch) {
            for(unsigned lane=0;lane<64;++lane) {
                random^=random<<13;random^=random>>17;random^=random<<5;
                const auto value=batch?random:std::uint32_t{1}<<(lane%32);
                for(unsigned byte=0;byte<4;++byte)
                    input[byte][lane]=std::uint8_t(value>>(8*byte));
            }
            auto a=_mm512_load_si512(input[0]),b=_mm512_load_si512(input[1]);
            auto c2=_mm512_load_si512(input[2]),d=_mm512_load_si512(input[3]);
            flat9(a,b,c2,d,c.data());
            _mm512_store_si512(output[0],a);_mm512_store_si512(output[1],b);
            _mm512_store_si512(output[2],c2);_mm512_store_si512(output[3],d);
            for(unsigned lane=0;lane<64;++lane) {
                std::uint32_t value=0,actual=0;
                for(unsigned byte=0;byte<4;++byte) {
                    value|=std::uint32_t(input[byte][lane])<<(8*byte);
                    actual|=std::uint32_t(output[byte][lane])<<(8*byte);
                }
                if(actual!=field::multiply32(scalar,value)) return false;
            }
        }
    }
    return true;
}
}

