#pragma once
#include "MediumField32.h"

namespace spin::research::char2::wide32 {

using Record128=medium32::Record128;

// The same GF(2^32) tower and byte basis as MediumField32.h. Multiplication
// by one scalar is a4x4 GF256 matrix M. Store its four cyclic diagonals:
// bytes of diagonal[k] are M[row,(row+k)%4], row=0..3.
// One record shares this scalar across four32-bit field coordinates.
struct Coefficient { std::array<std::uint32_t,4> diagonal{}; };
static_assert(sizeof(Coefficient)==16);

// Independent existing field oracle supplies the four matrix columns.
// This function is setup/checking code, never part of the timed product.
inline Coefficient prepareCoefficientDirect(std::uint32_t scalar) noexcept {
    std::uint32_t columns[4]{};
    for(unsigned column=0;column!=4;++column)
        columns[column]=medium32::towerMultiplyOracle(scalar,1U<<(8*column));
    Coefficient result{};
    for(unsigned diagonal=0;diagonal!=4;++diagonal)for(unsigned row=0;row!=4;++row)
        result.diagonal[diagonal]|=columns[(row+diagonal)%4]&(0xffU<<(8*row));
    return result;
}

using ByteTables=std::array<std::array<Coefficient,256>,4>;
inline const ByteTables& byteTables(){
    // Multiplication is GF256-linear in the scalar. Four byte contributions
    // allow compact bulk setup after only4096 oracle column evaluations.
    static const ByteTables tables=[] {
        ByteTables result{};
        for(unsigned byte=0;byte!=4;++byte)for(unsigned value=0;value!=256;++value)
            result[byte][value]=prepareCoefficientDirect(std::uint32_t(value)<<(8*byte));
        return result;
    }();
    return tables;
}
inline Coefficient prepareFromTables(std::uint32_t scalar,const ByteTables& table) noexcept {
    const auto& a=table[0][scalar&255];const auto& b=table[1][(scalar>>8)&255];
    const auto& c=table[2][(scalar>>16)&255];const auto& d=table[3][scalar>>24];
    Coefficient result{};
    for(unsigned j=0;j!=4;++j)result.diagonal[j]=a.diagonal[j]^b.diagonal[j]^c.diagonal[j]^d.diagonal[j];
    return result;
}
inline Coefficient prepareCoefficient(std::uint32_t scalar){return prepareFromTables(scalar,byteTables());}
inline void prepare(const std::uint32_t* scalars,Coefficient* output,std::size_t count){
    const auto& table=byteTables();
    for(std::size_t i=0;i!=count;++i)output[i]=prepareFromTables(scalars[i],table);
}

template<unsigned Index>
inline __m512i expandDiagonal(__m512i packed) noexcept {
    static_assert(Index<4);
    constexpr unsigned first=4*Index;
    constexpr std::uint32_t indices=first|((first+1)<<8)|((first+2)<<16)|((first+3)<<24);
    return _mm512_shuffle_epi8(packed,_mm512_set1_epi32(static_cast<int>(indices)));
}
inline __m512i multiply512(__m512i value,const Coefficient* coefficients) noexcept {
    const auto packed=_mm512_loadu_si512(coefficients);
    const auto p0=_mm512_gf2p8mul_epi8(value,expandDiagonal<0>(packed));
    const auto p1=_mm512_gf2p8mul_epi8(_mm512_ror_epi32(value,8),expandDiagonal<1>(packed));
    const auto p2=_mm512_gf2p8mul_epi8(_mm512_ror_epi32(value,16),expandDiagonal<2>(packed));
    const auto p3=_mm512_gf2p8mul_epi8(_mm512_ror_epi32(value,24),expandDiagonal<3>(packed));
    return _mm512_xor_si512(_mm512_ternarylogic_epi32(p0,p1,p2,0x96),p3);
}
inline std::uint32_t recoverScalar(const Coefficient& coefficient) noexcept {
    // Column zero is the product of the scalar with the field identity.
    return (coefficient.diagonal[0]&0xffU)|(coefficient.diagonal[3]&0xff00U)|
        (coefficient.diagonal[2]&0xff0000U)|(coefficient.diagonal[1]&0xff000000U);
}

#if defined(_MSC_VER)
#define WIDE32_NOINLINE __declspec(noinline)
#else
#define WIDE32_NOINLINE __attribute__((noinline))
#endif

// Four GFNI, three dword rotations, four coefficient shuffles, two Boolean
// operations per4 records. Setup reads64 bytes per batch, not16 as in the
// existing tower kernel. Exact input/output alias and all legal record
// alignments are supported. Coefficients cannot overlap output.
WIDE32_NOINLINE inline void diagonal(const Record128* input,const Coefficient* coefficients,
    Record128* output,std::size_t count) noexcept {
    std::size_t i=0;
    for(;i+4<=count;i+=4)
        _mm512_storeu_si512(output+i,multiply512(_mm512_loadu_si512(input+i),coefficients+i));
    for(;i!=count;++i){
        const auto scalar=recoverScalar(coefficients[i]);
        for(unsigned j=0;j!=4;++j)output[i].lane[j]=medium32::towerMultiplyScalar(input[i].lane[j],scalar);
    }
}

#undef WIDE32_NOINLINE
} // namespace spin::research::char2::wide32
