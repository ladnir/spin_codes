#pragma once
#include <array>
#include <cstdint>

namespace spin::research::twoms::rs8::field {
// Fixed tower basis: bytes (1,u,v,uv), least significant byte first.
// GF(256) uses z^8+z^4+z^3+z+1 (0x11b); u^2+u+0x20=0;
// v^2+v+0x2000=0 over GF(2^16). The two constants have absolute trace one.
// This basis is separate from the retained Plan's GF(2^16) outer.
inline std::uint8_t multiply8(std::uint8_t a, std::uint8_t b) noexcept {
    std::uint8_t result=0;
    for(unsigned i=0;i<8;++i) {
        if(b&1U) result^=a;
        a=std::uint8_t((unsigned(a)<<1)^((a&128U)?0x11bU:0U));
        b>>=1;
    }
    return result;
}
inline std::uint16_t multiply16(std::uint16_t a,std::uint16_t b) noexcept {
    const auto a0=std::uint8_t(a),a1=std::uint8_t(a>>8);
    const auto b0=std::uint8_t(b),b1=std::uint8_t(b>>8);
    const auto lo=multiply8(a0,b0)^multiply8(0x20,multiply8(a1,b1));
    const auto hi=multiply8(a0,b1)^multiply8(a1,b0)^multiply8(a1,b1);
    return std::uint16_t(lo|(unsigned(hi)<<8));
}
inline std::uint32_t multiply32(std::uint32_t a,std::uint32_t b) noexcept {
    const auto a0=std::uint16_t(a),a1=std::uint16_t(a>>16);
    const auto b0=std::uint16_t(b),b1=std::uint16_t(b>>16);
    const auto lo=multiply16(a0,b0)^multiply16(0x2000,multiply16(a1,b1));
    const auto hi=multiply16(a0,b1)^multiply16(a1,b0)^multiply16(a1,b1);
    return std::uint32_t(lo)|(std::uint32_t(hi)<<16);
}
inline void coefficients16(std::uint16_t c,std::uint8_t* result) noexcept {
    const auto c0=std::uint8_t(c),c1=std::uint8_t(c>>8);
    result[0]=c0;
    result[1]=multiply8(0x20,c1);
    result[2]=c0^c1;
}
inline std::array<std::uint8_t,9> coefficients(std::uint32_t c) noexcept {
    const auto c0=std::uint16_t(c),c1=std::uint16_t(c>>16);
    std::array<std::uint8_t,9> result{};
    coefficients16(c0,result.data());
    coefficients16(multiply16(0x2000,c1),result.data()+3);
    coefficients16(c0^c1,result.data()+6);
    return result;
}
}

