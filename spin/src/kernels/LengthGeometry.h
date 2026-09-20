#pragma once
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace spin::detail::kernel::length_geometry {
// Routes and encoded-coordinate counts use uint32_t: a representation limit,
// not a performance or certificate policy.
inline constexpr std::size_t maxMessageBlocks =
    std::numeric_limits<std::uint32_t>::max() / std::size_t{2};
constexpr bool valid(std::size_t k, unsigned step) noexcept {
    return step && k && k <= maxMessageBlocks && k % (std::size_t{128} * step) == 0;
}
inline void check(std::size_t k, unsigned step) {
    if (!valid(k, step))
        throw std::invalid_argument("K must be a positive multiple of 128*t within the 32-bit routing range");
}
inline std::size_t fromExponent(unsigned exponent) {
    if (exponent >= std::numeric_limits<std::size_t>::digits)
        throw std::invalid_argument("message exponent overflows size_t");
    return std::size_t{1} << exponent;
}
}
