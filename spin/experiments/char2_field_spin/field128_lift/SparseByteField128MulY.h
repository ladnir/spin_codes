#pragma once
#include "SparseByteField128.h"

namespace spin::research::char2::field128_lift::sparse_byte {

// Isolated fixed-scalar experiment. Y has degree16 over AES GF256 in
// this representation. No production encoder or sampled key is changed.
inline constexpr Element towerY{0x100,0};

constexpr std::uint8_t aesXtime(std::uint8_t value) noexcept {
    const unsigned x=value;
    return static_cast<std::uint8_t>((x<<1)^((0U-(x>>7))&0x1bU));
}

// 0x32 = X^5+X^4+X in AES GF256. Fixed straight-line arithmetic;
// no lookup table, general field product, or data-dependent branch.
constexpr std::uint8_t multiplyByteBy32(std::uint8_t value) noexcept {
    const auto x2=aesXtime(value);
    const auto x4=aesXtime(x2);
    const auto x8=aesXtime(x4);
    const auto x16=aesXtime(x8);
    const auto x32=aesXtime(x16);
    return static_cast<std::uint8_t>(x2^x16^x32);
}

// For a=sum a_i Y^i, shift by one coefficient and fold a_15 Y^16
// using Y^16=Y^5+Y^3+0x32. The words store little-endian coefficients.
constexpr Element multiplyY(Element value) noexcept {
    const auto top=static_cast<std::uint8_t>(value.hi>>56);
    return {(value.lo<<8)^multiplyByteBy32(top)^(std::uint64_t(top)<<24)^(std::uint64_t(top)<<40),
            (value.hi<<8)|(value.lo>>56)};
}

} // namespace spin::research::char2::field128_lift::sparse_byte
