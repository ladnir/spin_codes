#pragma once
#include "../kernels/Block.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::detail::certified22 {
using Block = storage::block;
inline constexpr std::size_t groupStride = 260;

bool validMessageSize(std::size_t k) noexcept;

// Four systematic GF(16) RS[16,8] rows per 128-to-256-coordinate group,
// native GF(2^16) symbol maps, four-bit routing through 64 regions, and
// the exact t64/s22 inner. State runs continuously across region boundaries.
struct Plan {
    Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed,
         bool withOuter = true);
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    std::size_t k = 0, n = 0, groups = 0;
    std::vector<std::uint32_t> route, routeBytes;
    std::vector<std::array<std::uint32_t,22>> reverseMatrices;
    std::vector<std::uint64_t> updates;
    std::vector<std::array<std::uint16_t,16>> outerMatrices16;
    std::vector<std::array<std::uint8_t,3>> outerField;
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Scratch is 64-byte aligned and disjoint from input/output. No entry point
// allocates. Exact input/output aliasing is supported at 16-byte alignment.
// Forward lanes=1,2,4 means interleaved 128/256/512-bit records; scratch
// consists of lanes contiguous planes of scratchBlocks() 128-bit blocks.
void transposeScalar(const Block*, Block*, Block* scratch, const Plan&);
void forwardScalar(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
void transposeFast(const Block*, Block*, Block* scratch, const Plan&);

// Internal stage entry points; also used for differential kernel checks.
// withOuter=false retains identical route/inner streams without RS16 tables.
void reverseRouteScalar(const Block*, Block* scratch, const Plan&);
void forwardInnerScalar(const Block* scratch, Block*, const Plan&, unsigned lanes);
void reverseRouteFast(const Block*, Block* scratch, const Plan&);
void outerFastCached(const Block* scratch, Block*, const Plan&);
void outerFastCachedUnrolled(const Block* scratch, Block*, const Plan&);
void outerFastStreaming(const Block* scratch, Block*, const Plan&);
}
