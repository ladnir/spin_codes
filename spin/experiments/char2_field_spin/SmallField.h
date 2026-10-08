#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>

namespace spin::research::char2 {

// Native polynomial bases. Coefficient c is the polynomial whose coefficient
// bits are c, not an element of an embedded GF(256) when the field is larger.
constexpr std::uint8_t reduction(unsigned bits) noexcept {
    return bits == 8 ? 0x1b : bits == 32 ? 0x8d : bits == 64 ? 0x1b : 0x87;
}
constexpr bool validBits(unsigned bits) noexcept {
    return bits == 8 || bits == 32 || bits == 64 || bits == 128;
}

constexpr std::uint16_t polynomialProduct8(std::uint8_t x,
                                           std::uint8_t c) noexcept {
    std::uint16_t result = 0;
    for (unsigned bit = 0; bit != 8; ++bit)
        if ((c >> bit) & 1U) result ^= std::uint16_t(x) << bit;
    return result;
}

// GFNI's byte 7 is the row for output bit 0; the matrix acts over GF(2).
constexpr std::uint64_t affineMatrix(std::uint8_t c, bool high) noexcept {
    std::uint64_t result = 0;
    for (unsigned output = 0; output != 8; ++output) {
        std::uint8_t row = 0;
        for (unsigned input = 0; input != 8; ++input) {
            const int coefficientBit = int(output + (high ? 8 : 0)) - int(input);
            if (coefficientBit >= 0 && coefficientBit < 8 &&
                ((c >> coefficientBit) & 1U)) row |= std::uint8_t(1U << input);
        }
        result |= std::uint64_t(row) << (8 * (7 - output));
    }
    return result;
}

struct Coeff128 {
    std::uint64_t low, high;
};

struct alignas(64) AffineTable {
    std::array<std::uint64_t, 256> low{}, high{};
    constexpr AffineTable() noexcept {
        for (unsigned c = 0; c != 256; ++c) {
            low[c] = affineMatrix(std::uint8_t(c), false);
            high[c] = affineMatrix(std::uint8_t(c), true);
        }
    }
};
inline constexpr AffineTable affineTable{};

inline Coeff128 coefficient128(std::uint8_t c) noexcept {
    return {affineTable.low[c], affineTable.high[c]};
}

struct alignas(64) Coeff4 {
    std::array<std::uint64_t, 8> low{}, high{};
};
inline Coeff4 prepare4(const std::uint8_t* coefficients) noexcept {
    Coeff4 result;
    for (unsigned i = 0; i != 4; ++i) {
        const auto c = coefficient128(coefficients[i]);
        result.low[2 * i] = result.low[2 * i + 1] = c.low;
        result.high[2 * i] = result.high[2 * i + 1] = c.high;
    }
    return result;
}

// Independent schoolbook carryless product and long division oracle. Input
// and output may alias. Field elements are little-endian polynomial bits.
void multiplyOracle(const std::uint8_t* input, std::uint8_t coefficient,
                    std::uint8_t* output, unsigned bits);
void multiplyPortable(const std::uint8_t* input, const std::uint8_t* coefficients,
                      std::uint8_t* output, std::size_t count, unsigned bits);

bool gfni128Available() noexcept;
bool gfni512Available() noexcept;

// Array entry points do not check CPU features. Call the predicates first.
void multiplyGfni128(const std::uint8_t* input, const std::uint8_t* coefficients,
                     std::uint8_t* output, std::size_t count);
void multiplyGfni512(const std::uint8_t* input, const std::uint8_t* coefficients,
                     std::uint8_t* output, std::size_t count, unsigned bits);
void multiplyGfniPrepared(const std::uint8_t* input, const Coeff4* coefficients,
                          std::uint8_t* output, std::size_t count);

// Latency-chain modes have no input-array traffic. They are benchmark controls,
// not complete encoders. setup and coefficient preparation are caller-owned.
void chainGfni128(std::uint8_t* state, const std::uint8_t* coefficients,
                  std::size_t count);
void chainGfni128Prepared(std::uint8_t* state, const Coeff128* coefficients,
                          std::size_t count);

} // namespace spin::research::char2
