#pragma once
#include "MediumField32.h"

namespace spin::research::char2::wide16 {

using Record128=medium32::Record128;

// GF(2^16)=AES-GF256[u]/(u^2+u+delta), delta=0x20. Each128-bit record
// contains eight coordinates in this field and shares one scalar across
// those coordinates. Coefficients are genuine16-bit field elements.
// Fixed embedded GF256 scalar maps remain bytewise GFNI in this basis.
struct Coefficient {
    std::uint8_t a,aPlusB,deltaB,b;
};
static_assert(sizeof(Coefficient)==4);
static_assert(medium32::byteTrace(medium32::delta)==1);

// Setup only. Zero is supported arithmetically; the benchmark ensemble
// samples from all65,535 nonzero scalars before entering the timed loop.
constexpr Coefficient prepareCoefficient(std::uint16_t value) noexcept {
    const auto a=static_cast<std::uint8_t>(value),b=static_cast<std::uint8_t>(value>>8);
    return {a,static_cast<std::uint8_t>(a^b),medium32::byteMultiply(medium32::delta,b),b};
}
inline void prepare(const std::uint16_t* values,Coefficient* output,std::size_t count) noexcept {
    for(std::size_t i=0;i!=count;++i)output[i]=prepareCoefficient(values[i]);
}

template<unsigned Offset>
inline __m512i expandPairs(__m512i packed) noexcept {
    static_assert(Offset==0||Offset==2);
    constexpr auto repeated=[](unsigned low){
        return static_cast<long long>((std::uint64_t(low)|std::uint64_t(low+1)<<8)*0x0001000100010001ULL);
    };
    const auto select=_mm512_setr_epi64(repeated(Offset),repeated(Offset),
        repeated(4+Offset),repeated(4+Offset),repeated(8+Offset),repeated(8+Offset),
        repeated(12+Offset),repeated(12+Offset));
    return _mm512_shuffle_epi8(packed,select);
}

// Four records, one exact16-byte setup load, two coefficient expansions,
// one byte-pair swap, two GFNI products, and one XOR. For x=x0+x1*u and
// c=a+b*u, result=(a*x0+delta*b*x1)+(b*x0+(a+b)*x1)*u.
inline __m512i multiply512(__m512i value,const Coefficient* coefficients) noexcept {
    const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(coefficients)));
    const auto direct=_mm512_gf2p8mul_epi8(value,expandPairs<0>(packed));
    const auto swapped=medium32::selectBytes<1,0,3,2>(value);
    return _mm512_xor_si512(direct,_mm512_gf2p8mul_epi8(swapped,expandPairs<2>(packed)));
}

#if defined(_MSC_VER)
#define WIDE16_NOINLINE __declspec(noinline)
#else
#define WIDE16_NOINLINE __attribute__((noinline))
#endif

// Exact input/output alias is supported; coefficient storage must be
// disjoint from output. Unaligned512-bit payload/setup loads are intentional.
// Tail reads touch only the declared records and coefficients.
WIDE16_NOINLINE inline void diagonal(const Record128* input,const Coefficient* coefficients,
    Record128* output,std::size_t count) noexcept {
    std::size_t i=0;
    for(;i+4<=count;i+=4)
        _mm512_storeu_si512(output+i,multiply512(_mm512_loadu_si512(input+i),coefficients+i));
    for(;i!=count;++i){
        const auto c=static_cast<std::uint16_t>(coefficients[i].a|(unsigned(coefficients[i].b)<<8));
        for(unsigned j=0;j!=4;++j){
            const auto x=input[i].lane[j];
            output[i].lane[j]=medium32::quadraticMultiply(static_cast<std::uint16_t>(x),c)|
                (std::uint32_t(medium32::quadraticMultiply(static_cast<std::uint16_t>(x>>16),c))<<16);
        }
    }
}

#undef WIDE16_NOINLINE
} // namespace spin::research::char2::wide16
