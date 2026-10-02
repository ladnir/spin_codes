// Exact adjoint-multiplier forward family, ordinary GF(2^64) transpose.
#pragma once
#include "../rs16x8/Tower64Randomizer.h"
#include "../rs16x8/Tower32ByteRandomizer.h"

namespace spin::research::rs::tower64byte {
inline void coefficients32(std::uint32_t scalar,std::uint8_t* out) {
    const auto c0=std::uint16_t(scalar),c1=std::uint16_t(scalar>>16);
    tower32byte::coefficients16(c0,out);
    tower32byte::coefficients16(tower32::multiply16(tower32::Nu,c1),out+3);
    tower32byte::coefficients16(c0^c1,out+6);
}
inline std::array<std::uint8_t,27> coefficients(std::uint64_t scalar) {
    if(!scalar)throw std::invalid_argument("tower64byte requires a nonzero scalar");
    const auto c0=std::uint32_t(scalar),c1=std::uint32_t(scalar>>32);
    std::array<std::uint8_t,27> result{};
    coefficients32(c0,result.data());
    coefficients32(tower32::multiply32(tower64::Theta,c1),result.data()+9);
    coefficients32(c0^c1,result.data()+18);
    return result;
}
inline std::array<std::uint64_t,64> multiplyRows(std::uint64_t scalar) {
    if(!scalar)throw std::invalid_argument("tower64byte requires a nonzero scalar");
    std::array<std::uint64_t,64> rows{};
    for(unsigned column=0;column<64;++column) {
        const auto image=tower64::multiply64(scalar,std::uint64_t{1}<<column);
        for(unsigned row=0;row<64;++row)rows[row]|=((image>>row)&1ULL)<<column;
    }
    return rows;
}
#if defined(__AVX512F__) || defined(SPIN_PACKET_ENABLE_VBMI)
// One 64-bit payload half: 27 GFNI byte MULs and 57 vector XORs.
static SPIN_FORCEINLINE void applyMultiply(__m512i& y0,__m512i& y1,__m512i& y2,__m512i& y3,
    __m512i& y4,__m512i& y5,__m512i& y6,__m512i& y7,const std::uint8_t* coeff) {
    auto p0=y0,p1=y1,p2=y2,p3=y3;
    auto q0=y4,q1=y5,q2=y6,q3=y7;
    auto r0=_mm512_xor_si512(y0,y4),r1=_mm512_xor_si512(y1,y5);
    auto r2=_mm512_xor_si512(y2,y6),r3=_mm512_xor_si512(y3,y7);
    tower32byte::applyMultiply(p0,p1,p2,p3,coeff);
    tower32byte::applyMultiply(q0,q1,q2,q3,coeff+9);
    tower32byte::applyMultiply(r0,r1,r2,r3,coeff+18);
    y0=_mm512_xor_si512(p0,q0);y1=_mm512_xor_si512(p1,q1);
    y2=_mm512_xor_si512(p2,q2);y3=_mm512_xor_si512(p3,q3);
    y4=_mm512_xor_si512(p0,r0);y5=_mm512_xor_si512(p1,r1);
    y6=_mm512_xor_si512(p2,r2);y7=_mm512_xor_si512(p3,r3);
}
#endif
}
