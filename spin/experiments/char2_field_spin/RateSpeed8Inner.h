#pragma once
#include "RateClimb8Inner.h"

namespace spin::research::char2 {

// GF2P8AFFINEQB consumes matrix row i in byte7-i of each 64-bit lane.
// This matrix is exactly multiplication by c in the AES polynomial basis.
constexpr std::uint64_t rs8AffineMatrix(std::uint8_t c)noexcept{
    std::uint64_t matrix=0;
    for(unsigned row=0;row!=8;++row){
        unsigned bits=0;
        for(unsigned column=0;column!=8;++column)
            bits|=unsigned((sfMultiply(c,std::uint8_t(1U<<column))>>row)&1U)<<column;
        matrix|=std::uint64_t(bits)<<(8*(7-row));
    }
    return matrix;
}
static_assert(rs8AffineMatrix(1)==0x0102040810204080ULL);
constexpr std::uint8_t rs8AffineApply(std::uint8_t value,std::uint64_t matrix)noexcept{
    unsigned result=0;
    for(unsigned row=0;row!=8;++row){
        unsigned bits=value&unsigned((matrix>>(8*(7-row)))&255U);
        bits^=bits>>4;bits^=bits>>2;bits^=bits>>1;
        result|=(bits&1U)<<row;
    }
    return std::uint8_t(result);
}
// Optional untimed scalar check of all65,536 constant/input byte pairs.
inline bool rs8CheckAffineConstants()noexcept{
    for(unsigned c=0;c!=256;++c){
        const auto matrix=rs8AffineMatrix(std::uint8_t(c));
        for(unsigned value=0;value!=256;++value)
            if(rs8AffineApply(std::uint8_t(value),matrix)!=sfMultiply(std::uint8_t(c),std::uint8_t(value)))
                return false;
    }
    return true;
}

#if defined(_MSC_VER)
#define RS8_INLINE __forceinline
#else
#define RS8_INLINE inline __attribute__((always_inline))
#endif

template<bool Affine,std::uint8_t Constant>
RS8_INLINE __m512i rs8ConstantMultiply(__m512i value)noexcept{
    if constexpr(Affine)
        return _mm512_gf2p8affine_epi64_epi8(value,
            _mm512_set1_epi64(static_cast<long long>(rs8AffineMatrix(Constant))),0);
    else return _mm512_gf2p8mul_epi8(value,_mm512_set1_epi8(static_cast<char>(Constant)));
}

// Same nineteen-product subset circuit and register order as p8Circulant.
// Only the fixed scalar multiplication instruction changes; the field/map do not.
template<bool Emission,bool Affine>
RS8_INLINE void rs8Circulant(__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3,
                            __m512i& z4,__m512i& z5,__m512i& z6,__m512i& z7)noexcept{
    if constexpr(!Affine)p8Circulant<Emission>(z0,z1,z2,z3,z4,z5,z6,z7);
    else{
        static_assert(p8ZetaCoefficient<Emission,0>()==1);
        p8Zeta(z0,z1,z2,z3,z4,z5,z6,z7);
        z7=_mm512_ternarylogic_epi32(z7,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,1>()>(z6),
            rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,2>()>(z5),0x96);
        z7=_mm512_ternarylogic_epi32(z7,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,3>()>(z4),
            rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,4>()>(z3),0x96);
        z7=_mm512_ternarylogic_epi32(z7,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,5>()>(z2),
            rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,6>()>(z1),0x96);
        z7=_mm512_xor_si512(z7,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,7>()>(z0));
        z6=_mm512_ternarylogic_epi32(z6,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,2>()>(z4),
            rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,4>()>(z2),0x96);
        z6=_mm512_xor_si512(z6,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,6>()>(z0));
        z5=_mm512_ternarylogic_epi32(z5,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,1>()>(z4),
            rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,4>()>(z1),0x96);
        z5=_mm512_xor_si512(z5,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,5>()>(z0));
        z4=_mm512_xor_si512(z4,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,4>()>(z0));
        z3=_mm512_ternarylogic_epi32(z3,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,1>()>(z2),
            rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,2>()>(z1),0x96);
        z3=_mm512_xor_si512(z3,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,3>()>(z0));
        z2=_mm512_xor_si512(z2,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,2>()>(z0));
        z1=_mm512_xor_si512(z1,rs8ConstantMultiply<true,p8ZetaCoefficient<Emission,1>()>(z0));
        p8Zeta(z0,z1,z2,z3,z4,z5,z6,z7);
    }
}

template<bool AffineScale>
RS8_INLINE void rs8History(__m512i& moment0,__m512i& moment1,__m512i& h0,__m512i& h1,
                           const std::uint8_t* diagonal)noexcept{
    const auto old0=h0,old1=h1;
    h0=_mm512_gf2p8mul_epi8(_mm512_xor_si512(old0,moment0),sfBroadcastFour(diagonal));
    h1=_mm512_gf2p8mul_epi8(_mm512_xor_si512(old1,moment1),sfBroadcastFour(diagonal+4));
    moment0=_mm512_xor_si512(old0,rs8ConstantMultiply<AffineScale,242>(moment0));
    moment1=_mm512_xor_si512(old1,rs8ConstantMultiply<AffineScale,242>(moment1));
}

// Two A/B candidates: AffineMaps=true,AffineScale=false (maps only), or
// true,true (all constants). false,false is an exact scheduling control.
// All have54 GFNI sites before the external eight input-diagonal sites.
// Compact input layout and folded output are identical to RateClimb8Inner.
template<class Profile,bool AffineMaps=true,bool AffineScale=true,bool Folded=true>
RS8_INLINE void rs8InnerTransposeFour(const Block* input,Block* output,__m512i& h0,__m512i& h1,
                                      const std::uint8_t* diagonal)noexcept{
    static_assert(4%Profile::period==0&&Profile::tileRecords*Profile::period==4*Profile::count);
    static_assert(Profile::commonRetained>=5&&Profile::commonRetained<=7);
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_loadu_si512(input+12);
    auto z4=_mm512_loadu_si512(input+16);
    auto z5=_mm512_setzero_si512(),z6=_mm512_setzero_si512(),z7=_mm512_setzero_si512();
    if constexpr(Profile::commonRetained>=6)z5=_mm512_loadu_si512(input+20);
    if constexpr(Profile::commonRetained>=7)z6=_mm512_loadu_si512(input+24);
    if constexpr(Profile::tailLayout==R8TailLayout::Phase3)
        z6=_mm512_inserti32x4(z6,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+24)),3);
    else if constexpr(Profile::tailLayout==R8TailLayout::OddPhases){
        z5=_mm512_inserti32x4(z5,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+20)),1);
        z5=_mm512_inserti32x4(z5,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+21)),3);
    }
    rs8Circulant<true,AffineMaps>(z0,z1,z2,z3,z4,z5,z6,z7);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    rs8History<AffineScale>(z3,z7,h0,h1,diagonal+24);
    rs8History<AffineScale>(z2,z6,h0,h1,diagonal+16);
    rs8History<AffineScale>(z1,z5,h0,h1,diagonal+8);
    rs8History<AffineScale>(z0,z4,h0,h1,diagonal);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    rs8Circulant<false,AffineMaps>(z0,z1,z2,z3,z4,z5,z6,z7);
    if constexpr(Folded){
        _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
        _mm512_storeu_si512(output+8,z2);_mm512_storeu_si512(output+12,z3);
        _mm512_storeu_si512(output+16,z4);_mm512_storeu_si512(output+20,z5);
        _mm512_storeu_si512(output+24,z6);_mm512_storeu_si512(output+28,z7);
    }else{
        c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
        c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
        _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z4);
        _mm512_storeu_si512(output+8,z1);_mm512_storeu_si512(output+12,z5);
        _mm512_storeu_si512(output+16,z2);_mm512_storeu_si512(output+20,z6);
        _mm512_storeu_si512(output+24,z3);_mm512_storeu_si512(output+28,z7);
    }
}

#undef RS8_INLINE
} // namespace spin::research::char2
