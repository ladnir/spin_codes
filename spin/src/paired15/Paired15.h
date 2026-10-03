#pragma once
#include "../kernels/Block.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::detail::paired15 {
using Block = storage::block;
inline constexpr std::size_t groupStride = 260;

// Frozen K16 mode-52 ensemble: four GF(16) RS[16,8] rows per
// 128-to-256-coordinate group, native GF(2^16) symbol maps, four-bit
// routing through 64 regions, and the paired t64/s15 inner with GL15.
// Setup resets the same three independent streams as the research profile:
// routeSeed, routeSeed ^ 0x75a1dc09, innerSeed ^ 0x3f625a92.
struct Plan {
    Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed);
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    std::size_t k = 0, n = 0, groups = 0;
    std::vector<std::uint32_t> route, routeBytes;
    std::vector<std::array<std::uint16_t,16>> reverseMatrices, outerMatrices16;
    std::vector<std::uint64_t> updates, forwardUpdates;
    std::vector<std::array<std::uint8_t,3>> outerField;
    std::vector<std::array<std::uint64_t,3>> outerAdjoint;
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Lanes is 1, 2, or 4, representing interleaved 128/256/512-bit records.
// Scratch has scratchBlocks()*lanes blocks, in contiguous lane planes.
// It is 64-byte aligned, disjoint from message/output. No entry point
// allocates. The entire message is consumed before output is written;
// exact input/output aliasing is supported, with 16-byte minimum alignment.
void forwardScalar(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
void forwardFast(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes, bool stream = true);
// Internal diagnostic/control variant; identical arithmetic with cached output.
void forwardFastCached(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
void transposeScalar(const Block*, Block*, Block* scratch, const Plan&);
void transposeFast(const Block*, Block*, Block* scratch, const Plan&);

// Internal stage boundary shared by the fast translation units.
void forwardOuterFast(const Block*, Block* scratch, const Plan&, unsigned lanes);
}
