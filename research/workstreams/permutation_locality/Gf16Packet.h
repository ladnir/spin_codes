#pragma once
#include <cstdint>
#include <immintrin.h>

namespace spin::research::gf16 {
// Polynomial basis modulo X^4+X+1. Setup/reference only.
constexpr unsigned multiply(unsigned a,unsigned b) {
    unsigned result=0;
    for(unsigned i=0;i<4;++i) {
        if(b&1)result^=a;
        b>>=1;a<<=1;if(a&16)a^=0x13;
    }
    return result;
}

// One byte per input lane: two mask bits select each 128-bit output lane.
// For the adjoint, output j takes the dot product with column j of M_a.
constexpr std::uint32_t transposeControl(unsigned a) {
    std::uint32_t result=0;
    for(unsigned in=0;in<4;++in)for(unsigned out=0;out<4;++out)
        if((multiply(a,1U<<out)>>in)&1)result|=std::uint32_t(3U<<(2*out))<<(8*in);
    return result;
}

// Four fixed shuffles with lane masks, no table gathers or data-dependent
// branches. The coefficient is public setup data; the elements remain opaque.
inline __m512i transpose(__m512i x,std::uint32_t control) {
    const auto a=_mm512_maskz_permutexvar_epi64(__mmask8(control),_mm512_setr_epi64(0,1,0,1,0,1,0,1),x);
    const auto b=_mm512_maskz_permutexvar_epi64(__mmask8(control>>8),_mm512_setr_epi64(2,3,2,3,2,3,2,3),x);
    const auto c=_mm512_maskz_permutexvar_epi64(__mmask8(control>>16),_mm512_setr_epi64(4,5,4,5,4,5,4,5),x);
    const auto d=_mm512_maskz_permutexvar_epi64(__mmask8(control>>24),_mm512_setr_epi64(6,7,6,7,6,7,6,7),x);
    return _mm512_xor_si512(_mm512_xor_si512(a,b),_mm512_xor_si512(c,d));
}
}
