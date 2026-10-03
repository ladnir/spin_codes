#include "RandomizerVariants.h"
// Keep the certified shared-parity/output circuit byte-for-byte unchanged.
// Isolate its externally named experiment wrappers so this translation unit
// can coexist with the frozen baseline Outer.cpp in one benchmark binary.
#define outervariants packet8wide24randomizer_reference
#include "../k16_codesign_100us/outer_variants/OuterVariants.cpp"
#undef outervariants
#include <bit>
#include <stdexcept>

namespace spin::research::packet8wide24opt {
namespace {
std::uint64_t ordinaryMatrix(std::uint8_t scalar) {
    std::uint64_t packed=0;
    for(unsigned row=0;row<8;++row) {
        unsigned mask=0;
        for(unsigned column=0;column<8;++column)
            mask|=((packet8wide24::multiply8(scalar,std::uint8_t(1U<<column))>>row)&1U)<<column;
        packed|=std::uint64_t(mask)<<(8*(7-row));
    }
    return packed;
}
std::uint8_t scalarAffine(std::uint64_t matrix,std::uint8_t input) {
    unsigned value=0;
    for(unsigned row=0;row<8;++row)
        value|=(std::popcount(unsigned(input&std::uint8_t(matrix>>(8*(7-row)))))&1U)<<row;
    return std::uint8_t(value);
}
template<bool Affine>
std::uint32_t scalarFlat(std::uint32_t input,const std::array<std::uint8_t,9>& coefficients,
                         const std::array<std::uint64_t,9>& matrices) {
    const auto x0=std::uint8_t(input),x1=std::uint8_t(input>>8);
    const auto x2=std::uint8_t(input>>16),x3=std::uint8_t(input>>24);
    auto product=[&](unsigned i,std::uint8_t value) {
        if constexpr(Affine)return scalarAffine(matrices[i],value);
        else return packet8wide24::multiply8(coefficients[i],value);
    };
    const auto p0=product(0,x0),p1=product(1,x1),p2=product(2,std::uint8_t(x0^x1));
    const auto q0=product(3,x2),q1=product(4,x3),q2=product(5,std::uint8_t(x2^x3));
    const auto r0=product(6,std::uint8_t(x0^x2)),r1=product(7,std::uint8_t(x1^x3));
    const auto r2=product(8,std::uint8_t(x0^x1^x2^x3));
    return std::uint32_t(p0^q0^p1^q1)|(std::uint32_t(p0^q0^p2^q2)<<8)
        |(std::uint32_t(p0^r0^p1^r1)<<16)|(std::uint32_t(p0^r0^p2^r2)<<24);
}

template<unsigned Index,bool Affine>
static SPIN_FORCEINLINE __m512i product(__m512i value,const std::uint8_t* coefficients,
                                        const std::uint64_t* matrices) {
    if constexpr(Affine)
        return _mm512_gf2p8affine_epi64_epi8(value,_mm512_set1_epi64(matrices[Index]),0);
    else
        return _mm512_gf2p8mul_epi8(value,_mm512_set1_epi8(char(coefficients[Index])));
}
template<bool Affine>
static SPIN_FORCEINLINE void flat9(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,
                                   const std::uint8_t* c,const std::uint64_t* m) {
    const auto x01=_mm512_xor_si512(x0,x1),x23=_mm512_xor_si512(x2,x3);
    const auto p0=product<0,Affine>(x0,c,m);
    const auto p1=product<1,Affine>(x1,c,m);
    const auto p2=product<2,Affine>(x01,c,m);
    const auto q0=product<3,Affine>(x2,c,m);
    const auto q1=product<4,Affine>(x3,c,m);
    const auto q2=product<5,Affine>(x23,c,m);
    const auto r0=product<6,Affine>(_mm512_xor_si512(x0,x2),c,m);
    const auto r1=product<7,Affine>(_mm512_xor_si512(x1,x3),c,m);
    const auto r2=product<8,Affine>(_mm512_xor_si512(x01,x23),c,m);
    // Four products per output, sharing two partial sums. This is six
    // post-product logical instructions instead of the nested circuit's ten.
    const auto pq=_mm512_xor_si512(p0,q0),pr=_mm512_xor_si512(p0,r0);
    x0=_mm512_ternarylogic_epi64(pq,p1,q1,0x96);
    x1=_mm512_ternarylogic_epi64(pq,p2,q2,0x96);
    x2=_mm512_ternarylogic_epi64(pr,p1,r1,0x96);
    x3=_mm512_ternarylogic_epi64(pr,p2,r2,0x96);
}
template<unsigned Row>
static SPIN_FORCEINLINE __m512i directRow(__m512i x0,__m512i x1,__m512i x2,__m512i x3,
                                         const std::uint64_t* m) {
    const auto p0=_mm512_gf2p8affine_epi64_epi8(x0,_mm512_set1_epi64(m[4*Row]),0);
    const auto p1=_mm512_gf2p8affine_epi64_epi8(x1,_mm512_set1_epi64(m[4*Row+1]),0);
    const auto p2=_mm512_gf2p8affine_epi64_epi8(x2,_mm512_set1_epi64(m[4*Row+2]),0);
    const auto p3=_mm512_gf2p8affine_epi64_epi8(x3,_mm512_set1_epi64(m[4*Row+3]),0);
    return _mm512_xor_si512(_mm512_ternarylogic_epi64(p0,p1,p2,0x96),p3);
}
static SPIN_FORCEINLINE void direct16(__m512i& x0,__m512i& x1,__m512i& x2,__m512i& x3,
                                      const std::uint64_t* m) {
    const auto y0=directRow<0>(x0,x1,x2,x3,m),y1=directRow<1>(x0,x1,x2,x3,m);
    const auto y2=directRow<2>(x0,x1,x2,x3,m),y3=directRow<3>(x0,x1,x2,x3,m);
    x0=y0;x1=y1;x2=y2;x3=y3;
}

template<unsigned Variant>
static SPIN_FORCEINLINE void mixSymbol(const Block* input,__m512i* low,__m512i* high,
    const std::uint8_t* coefficients,const std::uint64_t* affine,const std::uint64_t* direct) {
    auto v0=_mm512_load_si512(input),v1=_mm512_load_si512(input+4);
    auto v2=_mm512_load_si512(input+8),v3=_mm512_load_si512(input+12);
    auto v4=_mm512_load_si512(input+16),v5=_mm512_load_si512(input+20);
    auto v6=_mm512_load_si512(input+24),v7=_mm512_load_si512(input+28);
    if constexpr(Variant==1) {
        flat9<false>(v0,v2,v4,v6,coefficients,affine);
        flat9<false>(v1,v3,v5,v7,coefficients,affine);
    } else if constexpr(Variant==2) {
        flat9<true>(v0,v2,v4,v6,coefficients,affine);
        flat9<true>(v1,v3,v5,v7,coefficients,affine);
    } else {
        static_assert(Variant==3);
        direct16(v0,v2,v4,v6,direct);direct16(v1,v3,v5,v7,direct);
    }
    low[0]=v0;low[1]=v1;low[2]=v2;low[3]=v3;
    high[0]=v4;high[1]=v5;high[2]=v6;high[3]=v7;
}
template<unsigned Variant>
static SPIN_NOINLINE void group(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,9>* __restrict coefficients,
    const std::array<std::uint64_t,9>* __restrict affine,
    const std::array<std::uint64_t,16>* __restrict direct) {
    alignas(64) __m512i low[64],high[64];
    if constexpr(Variant==4) {
        // Retain the precise byteplane layout, but omit all nine products.
        low[0]=_mm512_load_si512(input);low[1]=_mm512_load_si512(input+4);
        low[2]=_mm512_load_si512(input+8);low[3]=_mm512_load_si512(input+12);
        high[0]=_mm512_load_si512(input+16);high[1]=_mm512_load_si512(input+20);
        high[2]=_mm512_load_si512(input+24);high[3]=_mm512_load_si512(input+28);
    }
    for(unsigned symbol=Variant==4?1:0;symbol<16;++symbol)
        mixSymbol<Variant==4?1:Variant>(input+32*symbol,low+4*symbol,high+4*symbol,
            coefficients[symbol].data(),affine[symbol].data(),direct[symbol].data());
    k16codesign::packet8wide24randomizer_reference::finish<0,true>(low,output);
    k16codesign::packet8wide24randomizer_reference::finish<2,true>(low,output);
    k16codesign::packet8wide24randomizer_reference::finish<0,true>(high,output+128);
    k16codesign::packet8wide24randomizer_reference::finish<2,true>(high,output+128);
}
template<unsigned Variant>
void run(const Block* scratch,Block* output,const Plan& plan,const RandomizerPlan& tables) {
    for(std::size_t g=0;g<plan.groups;++g)
        group<Variant>(scratch+packet8wide24::groupStride*g,output+256*g,
            plan.outerCoefficients.data()+16*g,tables.affine9.data()+16*g,tables.direct16.data()+16*g);
}
}

RandomizerPlan::RandomizerPlan(const Plan& plan) {
    if(plan.outerRows.size()!=16*plan.groups || plan.outerCoefficients.size()!=plan.outerRows.size())
        throw std::invalid_argument("randomizer table geometry mismatch");
    affine9.resize(plan.outerRows.size());direct16.resize(plan.outerRows.size());
    for(std::size_t symbol=0;symbol<plan.outerRows.size();++symbol) {
        for(unsigned i=0;i<9;++i)affine9[symbol][i]=ordinaryMatrix(plan.outerCoefficients[symbol][i]);
        for(unsigned out=0;out<4;++out)for(unsigned in=0;in<4;++in) {
            std::uint64_t matrix=0;
            for(unsigned row=0;row<8;++row)
                matrix|=std::uint64_t((plan.outerRows[symbol][8*out+row]>>(8*in))&255U)<<(8*(7-row));
            direct16[symbol][4*out+in]=matrix;
        }
    }
}

bool randomizerTablesMatch(const Plan& plan,const RandomizerPlan& tables) {
    if(tables.affine9.size()!=plan.outerRows.size() || tables.direct16.size()!=plan.outerRows.size()
        || plan.outerCoefficients.size()!=plan.outerRows.size())return false;
    for(std::size_t symbol=0;symbol<plan.outerRows.size();++symbol)
        for(unsigned test=0;test<34;++test) {
            const auto input=test<32?std::uint32_t(1U<<test):
                (test==32?0xa573c96eU:0xf0192bd7U)^std::uint32_t(symbol*0x9e3779b9U);
            std::uint32_t expected=0,actual=0;
            for(unsigned row=0;row<32;++row)
                expected|=(std::popcount(plan.outerRows[symbol][row]&input)&1U)<<row;
            for(unsigned out=0;out<4;++out) {
                std::uint8_t value=0;
                for(unsigned in=0;in<4;++in)
                    value^=scalarAffine(tables.direct16[symbol][4*out+in],std::uint8_t(input>>(8*in)));
                actual|=std::uint32_t(value)<<(8*out);
            }
            if(actual!=expected || scalarFlat<false>(input,plan.outerCoefficients[symbol],tables.affine9[symbol])!=expected
                || scalarFlat<true>(input,plan.outerCoefficients[symbol],tables.affine9[symbol])!=expected)return false;
        }
    return true;
}

void applyOmittedRows(Plan& plan) {
    if(plan.outerRows.size()!=16*plan.groups)
        throw std::invalid_argument("omitted randomizer geometry mismatch");
    for(std::size_t group=0;group<plan.groups;++group)
        for(unsigned row=0;row<32;++row)plan.outerRows[16*group][row]=std::uint32_t(1U<<row);
}

void outerRandomizer(const Block* scratch,Block* output,const Plan& plan,
                     const RandomizerPlan& tables,unsigned variant) {
    switch(variant) {
    case 0:packet8wide24::outerFast(scratch,output,plan);return;
    case 1:run<1>(scratch,output,plan,tables);return;
    case 2:run<2>(scratch,output,plan,tables);return;
    case 3:run<3>(scratch,output,plan,tables);return;
    case 4:run<4>(scratch,output,plan,tables);return;
    default:throw std::invalid_argument("randomizer variant must be 0..4");
    }
}
}
