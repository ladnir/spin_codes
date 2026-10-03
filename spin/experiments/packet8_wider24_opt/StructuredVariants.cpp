#include "StructuredVariants.h"
#include "../../src/kernels/SetupRandom.h"
#define outervariants packet8wide24structured_reference
#include "../k16_codesign_100us/outer_variants/OuterVariants.cpp"
#undef outervariants
#include <bit>
#include <stdexcept>

namespace spin::research::packet8wide24opt {
namespace {
constexpr std::uint8_t multiply4(std::uint8_t a,std::uint8_t b) {
    unsigned result=0;
    for(unsigned i=0;i<4;++i) {
        if(b&1)result^=a;
        a=std::uint8_t((a<<1)^((a&8)?0x13:0));b>>=1;
    }
    return std::uint8_t(result);
}
constexpr std::uint8_t multiplyTower8(std::uint8_t a,std::uint8_t b) {
    const auto a0=std::uint8_t(a&15),a1=std::uint8_t(a>>4);
    const auto b0=std::uint8_t(b&15),b1=std::uint8_t(b>>4);
    const auto p=multiply4(a0,b0),q=multiply4(a1,b1);
    const auto r=multiply4(std::uint8_t(a0^a1),std::uint8_t(b0^b1));
    return std::uint8_t((p^multiply4(8,q))|((r^p)<<4));
}
constexpr std::uint64_t ordinaryMatrix(std::uint8_t scalar) {
    std::uint64_t packed=0;
    for(unsigned row=0;row<8;++row) {
        unsigned mask=0;
        for(unsigned column=0;column<8;++column)
            mask|=((multiplyTower8(scalar,std::uint8_t(1U<<column))>>row)&1U)<<column;
        packed|=std::uint64_t(mask)<<(8*(7-row));
    }
    return packed;
}
constexpr auto xtimeMatrix=ordinaryMatrix(2);
std::uint8_t scalarAffine(std::uint64_t matrix,std::uint8_t input) {
    unsigned value=0;
    for(unsigned row=0;row<8;++row)
        value|=(std::popcount(unsigned(input&std::uint8_t(matrix>>(8*(7-row)))))&1U)<<row;
    return std::uint8_t(value);
}
std::uint32_t scalarMap(std::uint32_t input,const std::array<std::uint8_t,7>& coefficients) {
    std::uint8_t x[4]={std::uint8_t(input),std::uint8_t(input>>8),
        std::uint8_t(input>>16),std::uint8_t(input>>24)};
    for(unsigned i=1;i<4;++i)x[i]=multiplyTower8(coefficients[i-1],x[i]);
    constexpr std::uint8_t h[4][4]={{2,3,1,1},{1,2,3,1},{1,1,2,3},{3,1,1,2}};
    std::uint32_t result=0;
    for(unsigned row=0;row<4;++row) {
        std::uint8_t sum=0;
        for(unsigned column=0;column<4;++column)sum^=multiplyTower8(h[row][column],x[column]);
        result|=std::uint32_t(multiplyTower8(coefficients[3+row],sum))<<(8*row);
    }
    return result;
}

template<unsigned Index>
static SPIN_FORCEINLINE __m512i diagonal(__m512i value,const std::uint64_t* matrices) {
    return _mm512_gf2p8affine_epi64_epi8(value,_mm512_set1_epi64(matrices[Index]),0);
}
template<bool Bitwise>
static SPIN_FORCEINLINE __m512i xtime(__m512i value) {
    if constexpr(Bitwise) {
        const auto low=_mm512_and_si512(_mm512_slli_epi16(value,1),_mm512_set1_epi8(char(0xee)));
        const auto carry=_mm512_and_si512(_mm512_srli_epi16(value,3),_mm512_set1_epi8(0x11));
        return _mm512_ternarylogic_epi64(low,carry,_mm512_slli_epi16(carry,1),0x96);
    } else {
        return _mm512_gf2p8affine_epi64_epi8(value,_mm512_set1_epi64(xtimeMatrix),0);
    }
}
template<bool Bitwise>
static SPIN_FORCEINLINE void mix(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,
                                 const std::uint64_t* matrices) {
    x1=diagonal<0>(x1,matrices);x2=diagonal<1>(x2,matrices);x3=diagonal<2>(x3,matrices);
    const auto x01=_mm512_xor_si512(x0,x1),x23=_mm512_xor_si512(x2,x3);
    const auto all=_mm512_xor_si512(x01,x23);
    const auto d0=xtime<Bitwise>(x01);
    const auto d1=xtime<Bitwise>(_mm512_xor_si512(x1,x2));
    const auto d2=xtime<Bitwise>(x23);
    // The four cyclic differences sum to zero, and xtime is linear.
    const auto d3=_mm512_ternarylogic_epi64(d0,d1,d2,0x96);
    x0=diagonal<3>(_mm512_ternarylogic_epi64(x0,all,d0,0x96),matrices);
    x1=diagonal<4>(_mm512_ternarylogic_epi64(x1,all,d1,0x96),matrices);
    x2=diagonal<5>(_mm512_ternarylogic_epi64(x2,all,d2,0x96),matrices);
    x3=diagonal<6>(_mm512_ternarylogic_epi64(x3,all,d3,0x96),matrices);
}
template<bool Bitwise>
static SPIN_FORCEINLINE void mixSymbol(const Block* input,__m512i* low,__m512i* high,
                                       const std::uint64_t* matrices) {
    auto v0=_mm512_load_si512(input),v1=_mm512_load_si512(input+4);
    auto v2=_mm512_load_si512(input+8),v3=_mm512_load_si512(input+12);
    auto v4=_mm512_load_si512(input+16),v5=_mm512_load_si512(input+20);
    auto v6=_mm512_load_si512(input+24),v7=_mm512_load_si512(input+28);
    mix<Bitwise>(v0,v2,v4,v6,matrices);mix<Bitwise>(v1,v3,v5,v7,matrices);
    low[0]=v0;low[1]=v1;low[2]=v2;low[3]=v3;
    high[0]=v4;high[1]=v5;high[2]=v6;high[3]=v7;
}
template<bool Bitwise>
static SPIN_NOINLINE void group(const Block* __restrict input,Block* __restrict output,
                                const std::array<std::uint64_t,7>* __restrict matrices) {
    alignas(64) __m512i low[64],high[64];
    for(unsigned symbol=0;symbol<16;++symbol)
        mixSymbol<Bitwise>(input+32*symbol,low+4*symbol,high+4*symbol,matrices[symbol].data());
    k16codesign::packet8wide24structured_reference::finish<0,true>(low,output);
    k16codesign::packet8wide24structured_reference::finish<2,true>(low,output);
    k16codesign::packet8wide24structured_reference::finish<0,true>(high,output+128);
    k16codesign::packet8wide24structured_reference::finish<2,true>(high,output+128);
}
template<bool Bitwise>
void run(const Block* scratch,Block* output,const Plan& plan,const StructuredPlan& tables) {
    for(std::size_t g=0;g<plan.groups;++g)
        group<Bitwise>(scratch+packet8wide24::groupStride*g,output+256*g,tables.diagonal.data()+16*g);
}
}

StructuredPlan::StructuredPlan(const Plan& plan,std::uint64_t seed) {
    if(plan.outerRows.size()!=16*plan.groups)
        throw std::invalid_argument("structured randomizer geometry mismatch");
    coefficients.resize(plan.outerRows.size());diagonal.resize(plan.outerRows.size());
    reverseRows.resize(plan.outerRows.size());
    detail::kernel::setup::Words random(seed^0x43e091acf257d6b8ULL);
    for(std::size_t symbol=0;symbol<coefficients.size();++symbol) {
        for(unsigned i=0;i<7;++i) {
            do { coefficients[symbol][i]=std::uint8_t(random()); } while(!coefficients[symbol][i]);
            diagonal[symbol][i]=ordinaryMatrix(coefficients[symbol][i]);
        }
        for(unsigned column=0;column<32;++column) {
            const auto image=scalarMap(std::uint32_t(1U<<column),coefficients[symbol]);
            for(unsigned row=0;row<32;++row)reverseRows[symbol][row]|=((image>>row)&1U)<<column;
        }
    }
}

void applyStructuredRows(Plan& plan,const StructuredPlan& tables) {
    if(tables.reverseRows.size()!=16*plan.groups || tables.diagonal.size()!=tables.reverseRows.size())
        throw std::invalid_argument("structured oracle row geometry mismatch");
    plan.outerRows=tables.reverseRows;
}

bool structuredTablesMatch(const StructuredPlan& tables) {
    if(tables.coefficients.size()!=tables.diagonal.size() || tables.reverseRows.size()!=tables.diagonal.size())return false;
    for(unsigned value=0;value<256;++value) {
        const auto carry=(value>>3)&0x11;
        const auto bits=((value<<1)&0xee)^carry^(carry<<1);
        if(bits!=multiplyTower8(2,std::uint8_t(value)) || bits!=scalarAffine(xtimeMatrix,std::uint8_t(value)))return false;
    }
    for(std::size_t symbol=0;symbol<tables.coefficients.size();++symbol) {
        for(unsigned i=0;i<7;++i)for(unsigned bit=0;bit<8;++bit)
            if(!tables.coefficients[symbol][i] || scalarAffine(tables.diagonal[symbol][i],std::uint8_t(1U<<bit))
                !=multiplyTower8(tables.coefficients[symbol][i],std::uint8_t(1U<<bit)))return false;
        for(unsigned test=0;test<34;++test) {
            const auto input=test<32?std::uint32_t(1U<<test):
                (test==32?0x17aec983U:0xd97f1ba6U)^std::uint32_t(symbol*0x9e3779b9U);
            std::uint32_t value=0;
            for(unsigned row=0;row<32;++row)value|=(std::popcount(tables.reverseRows[symbol][row]&input)&1U)<<row;
            if(value!=scalarMap(input,tables.coefficients[symbol]))return false;
        }
    }
    return true;
}

void outerStructured(const Block* scratch,Block* output,const Plan& plan,
                     const StructuredPlan& tables,unsigned variant) {
    switch(variant) {
    case 0:run<false>(scratch,output,plan,tables);return;
    case 1:run<true>(scratch,output,plan,tables);return;
    default:throw std::invalid_argument("structured randomizer variant must be 0 or 1");
    }
}
}
