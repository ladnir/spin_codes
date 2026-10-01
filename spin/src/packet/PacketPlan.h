#pragma once
#include "../kernels/Block.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::detail::packet {
inline constexpr std::size_t tileStride = 1028;
bool validMessageSize(std::size_t k) noexcept;

// One immutable setup, shared by the scalar and fused implementations.
// All counts and route offsets are in 128-bit records, not bytes.
struct Plan {
    explicit Plan(std::size_t k, std::uint64_t seed);
    Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed);
    Plan(const Plan&)=delete;
    Plan& operator=(const Plan&)=delete;
    Plan(Plan&&)=delete;
    Plan& operator=(Plan&&)=delete;
    std::size_t n = 0;
    std::vector<std::uint32_t> route;
    std::vector<std::array<std::uint16_t, 16>> reverseMatrices;
    std::vector<std::array<std::uint32_t, 32>> outerMatrices;
    std::vector<std::uint64_t> compactStorage;
    std::vector<std::uint64_t> composedUpdates;

    const std::uint64_t* compactCoefficients() const noexcept;
    std::size_t scratchBlocks() const noexcept { return (n / 1024) * tileStride; }
    std::size_t setupBytes() const noexcept;
};

// Allocation-free. The whole input is consumed before any output is written,
// so output may equal input. Scratch is disjoint and holds scratchBlocks().
void transposeScalar(const storage::block* input, storage::block* output,
                     storage::block* scratch, const Plan&);
void transposeFast(const storage::block* input, storage::block* output,
                   storage::block* scratch, std::size_t n,
                   const std::uint32_t* routeBases,
                   const std::uint64_t* compactCoefficients,
                   const std::uint64_t* composedUpdateWords);
}
