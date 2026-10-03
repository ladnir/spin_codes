#include "StoreVariants.h"
// Reuse the frozen shared-parity circuit without changing or duplicating its
// arithmetic. Rename its experiment wrappers to avoid external-name clashes.
#define outervariants packet8wide24scaling_reference
#include "../k16_codesign_100us/outer_variants/OuterVariants.cpp"
#undef outervariants

namespace spin::research::packet8wide24scaling {
namespace {
namespace reference=k16codesign::packet8wide24scaling_reference;
using k16codesign::nativeouter::PlaneValues;

template<unsigned Index>
static SPIN_FORCEINLINE __m512i product(__m512i value,const std::uint8_t* coefficients) {
    return _mm512_gf2p8mul_epi8(value,_mm512_set1_epi8(char(coefficients[Index])));
}
// Exact flat9 MUL winner, including its shared ternary-logical schedule.
static SPIN_FORCEINLINE void flat9(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,
                                   const std::uint8_t* c) {
    const auto x01=_mm512_xor_si512(x0,x1),x23=_mm512_xor_si512(x2,x3);
    const auto p0=product<0>(x0,c);
    const auto p1=product<1>(x1,c);
    const auto p2=product<2>(x01,c);
    const auto q0=product<3>(x2,c);
    const auto q1=product<4>(x3,c);
    const auto q2=product<5>(x23,c);
    const auto r0=product<6>(_mm512_xor_si512(x0,x2),c);
    const auto r1=product<7>(_mm512_xor_si512(x1,x3),c);
    const auto r2=product<8>(_mm512_xor_si512(x01,x23),c);
    const auto pq=_mm512_xor_si512(p0,q0),pr=_mm512_xor_si512(p0,r0);
    x0=_mm512_ternarylogic_epi64(pq,p1,q1,0x96);
    x1=_mm512_ternarylogic_epi64(pq,p2,q2,0x96);
    x2=_mm512_ternarylogic_epi64(pr,p1,r1,0x96);
    x3=_mm512_ternarylogic_epi64(pr,p2,r2,0x96);
}
static SPIN_FORCEINLINE void mixSymbol(const Block* input,__m512i* low,__m512i* high,
                                       const std::uint8_t* coefficients) {
    auto v0=_mm512_load_si512(input),v1=_mm512_load_si512(input+4);
    auto v2=_mm512_load_si512(input+8),v3=_mm512_load_si512(input+12);
    auto v4=_mm512_load_si512(input+16),v5=_mm512_load_si512(input+20);
    auto v6=_mm512_load_si512(input+24),v7=_mm512_load_si512(input+28);
    flat9(v0,v2,v4,v6,coefficients);
    flat9(v1,v3,v5,v7,coefficients);
    low[0]=v0;low[1]=v1;low[2]=v2;low[3]=v3;
    high[0]=v4;high[1]=v5;high[2]=v6;high[3]=v7;
}

template<unsigned Symbol,unsigned Plane>
static SPIN_FORCEINLINE void storePairStreaming(const PlaneValues& a,const PlaneValues& b,Block* out) {
    static_assert(Symbol<8 && (Plane==0 || Plane==2));
    const auto i0=_mm512_setr_epi64(0x3830282018100800ULL,0x7870686058504840ULL,0x3931292119110901ULL,0x7971696159514941ULL,0x3a322a221a120a02ULL,0x7a726a625a524a42ULL,0x3b332b231b130b03ULL,0x7b736b635b534b43ULL);
    const auto i1=_mm512_setr_epi64(0x3c342c241c140c04ULL,0x7c746c645c544c44ULL,0x3d352d251d150d05ULL,0x7d756d655d554d45ULL,0x3e362e261e160e06ULL,0x7e766e665e564e46ULL,0x3f372f271f170f07ULL,0x7f776f675f574f47ULL);
    const auto basis=_mm512_set1_epi64(0x8040201008040201ULL);
    const auto x0=_mm512_gf2p8affine_epi64_epi8(basis,a.v[Symbol],0);
    const auto x1=_mm512_gf2p8affine_epi64_epi8(basis,b.v[Symbol],0);
    _mm512_stream_si512(reinterpret_cast<__m512i*>(out+32*Plane+4*Symbol),
        _mm512_permutex2var_epi8(x0,i0,x1));
    _mm512_stream_si512(reinterpret_cast<__m512i*>(out+32*(Plane+1)+4*Symbol),
        _mm512_permutex2var_epi8(x0,i1,x1));
}
template<unsigned Plane,bool Stream>
static SPIN_FORCEINLINE void finish(const __m512i* packed,Block* output) {
    if constexpr(!Stream)reference::finish<Plane,true>(packed,output);
    else {
        const auto a=reference::sharedParity<Plane>(packed);
        const auto b=reference::sharedParity<Plane+1>(packed);
        storePairStreaming<0,Plane>(a,b,output);storePairStreaming<1,Plane>(a,b,output);
        storePairStreaming<2,Plane>(a,b,output);storePairStreaming<3,Plane>(a,b,output);
        storePairStreaming<4,Plane>(a,b,output);storePairStreaming<5,Plane>(a,b,output);
        storePairStreaming<6,Plane>(a,b,output);storePairStreaming<7,Plane>(a,b,output);
    }
}
template<bool Stream>
static SPIN_NOINLINE void group(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,9>* __restrict coefficients) {
    alignas(64) __m512i low[64],high[64];
    for(unsigned symbol=0;symbol<16;++symbol)
        mixSymbol(input+32*symbol,low+4*symbol,high+4*symbol,coefficients[symbol].data());
    finish<0,Stream>(low,output);finish<2,Stream>(low,output);
    finish<0,Stream>(high,output+128);finish<2,Stream>(high,output+128);
}
template<bool Stream>
void run(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t g=0;g<plan.groups;++g)
        group<Stream>(scratch+packet8wide24::groupStride*g,output+256*g,
                      plan.outerCoefficients.data()+16*g);
    if constexpr(Stream)_mm_sfence();
}
}

void outer(const Block* scratch,Block* output,const Plan& plan,bool stream) {
    if(stream && (reinterpret_cast<std::uintptr_t>(output)&63U)==0)
        run<true>(scratch,output,plan);
    else run<false>(scratch,output,plan);
}
}
