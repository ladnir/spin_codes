#pragma once
#include "Punctured8Inner.h"

namespace spin::research::char2 {

// count is the total retained output count over one chronological period.
// Every supported period divides a four-packet tile and the 64-row region.
enum class R8TailLayout : std::uint8_t { None, Phase3, OddPhases };
struct R8Fixed6 {
    static constexpr unsigned dimension=80,period=1,count=6,commonRetained=6,tileRecords=24;
    static constexpr R8TailLayout tailLayout=R8TailLayout::None;
    static constexpr const char* name="keep6";
    static constexpr unsigned retained(std::size_t)noexcept{return 6;}
};
struct R8Periodic84 {
    static constexpr unsigned dimension=84,period=4,count=25,commonRetained=6,tileRecords=25;
    static constexpr R8TailLayout tailLayout=R8TailLayout::Phase3;
    static constexpr const char* name="periodic84";
    static constexpr unsigned retained(std::size_t step)noexcept{return (step&3)==3?7:6;}
};
struct R8Periodic74 {
    static constexpr unsigned dimension=74,period=2,count=11,commonRetained=5,tileRecords=22;
    static constexpr R8TailLayout tailLayout=R8TailLayout::OddPhases;
    static constexpr const char* name="periodic74";
    static constexpr unsigned retained(std::size_t step)noexcept{return 5+unsigned(step&1);}
};
// Exact-map control: identical setup, source basis, and serialization to p8.
struct R8Fixed7 {
    static constexpr unsigned dimension=84,period=1,count=7,commonRetained=7,tileRecords=28;
    static constexpr R8TailLayout tailLayout=R8TailLayout::None;
    static constexpr const char* name="keep7-control";
    static constexpr unsigned retained(std::size_t)noexcept{return 7;}
};

#if defined(_MSC_VER)
#define R8_PACKET_INLINE __forceinline
#else
#define R8_PACKET_INLINE inline __attribute__((always_inline))
#endif

// Each full coordinate occupies four 128-bit records, in forward packet
// order. Remaining retained coordinates are packed by increasing phase:
// periodic84 has one record (coordinate6,phase3); periodic74 has two
// records (coordinate5,phases1,3). No load crosses the compact input tile.
// Baseline output is the same 32 packet-major raw records used by p8.
// Folded output stays coordinate-major: output[4*j+phase]. Its fixed
// permutation is absorbed into the input diagonal and route during setup.
template<class Profile,bool Folded=false>
R8_PACKET_INLINE void r8InnerTransposeFour(const Block* input,Block* output,
                            __m512i& h0,__m512i& h1,const std::uint8_t* diagonal)noexcept{
    static_assert(4%Profile::period==0&&Profile::tileRecords*Profile::period==4*Profile::count);
    static_assert(Profile::commonRetained>=5&&Profile::commonRetained<=7);
    auto z0=_mm512_loadu_si512(input),z1=_mm512_loadu_si512(input+4);
    auto z2=_mm512_loadu_si512(input+8),z3=_mm512_loadu_si512(input+12);
    auto z4=_mm512_loadu_si512(input+16);
    auto z5=_mm512_setzero_si512(),z6=_mm512_setzero_si512(),z7=_mm512_setzero_si512();
    if constexpr(Profile::commonRetained>=6)z5=_mm512_loadu_si512(input+20);
    if constexpr(Profile::commonRetained>=7)z6=_mm512_loadu_si512(input+24);
    if constexpr(Profile::tailLayout==R8TailLayout::Phase3){
        z6=_mm512_inserti32x4(z6,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+24)),3);
    }else if constexpr(Profile::tailLayout==R8TailLayout::OddPhases){
        z5=_mm512_inserti32x4(z5,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+20)),1);
        z5=_mm512_inserti32x4(z5,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+21)),3);
    }
    p8Circulant<true>(z0,z1,z2,z3,z4,z5,z6,z7);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    p8History(z3,z7,h0,h1,diagonal+24);
    p8History(z2,z6,h0,h1,diagonal+16);
    p8History(z1,z5,h0,h1,diagonal+8);
    p8History(z0,z4,h0,h1,diagonal);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    c8TransposeRecords(z4,z5,z6,z7,z4,z5,z6,z7);
    p8Circulant<false>(z0,z1,z2,z3,z4,z5,z6,z7);
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

#undef R8_PACKET_INLINE
} // namespace spin::research::char2
