#pragma once
#include "../paired15/Paired15.h"
#include "../certified22/Plan.h"
#include <optional>

namespace spin::detail::rs8 {
using Block = storage::block;
inline constexpr std::size_t groupStride = 260;
static_assert(groupStride == paired15::groupStride);
static_assert(groupStride == certified22::groupStride);

bool validMessageSize(std::size_t k) noexcept;

// Eight systematic GF(16) RS[8,4] rows per 128-to-256-coordinate group,
// independent nonzero GF(2^32) reverse symbol maps, and four-bit routing
// through 64 regions. The paired t64/s15 inner is used through K16; larger
// messages use the exact t64/s22 inner. State crosses region boundaries.
struct Plan {
    Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed);
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    std::size_t k = 0, n = 0, groups = 0;
    std::optional<paired15::Plan> paired15Plan;
    std::optional<certified22::Plan> certified22Plan;
    std::vector<std::array<std::uint8_t,9>> field;
    // Column j of the literal reverse symbol map M_c, for scalar encoding.
    std::vector<std::array<std::uint32_t,32>> reverseColumns;
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Scratch is 64-byte aligned and disjoint from input/output. No entry point
// allocates. Exact input/output aliasing is supported at 16-byte alignment.
// Forward lanes=1,2,4 selects interleaved 128/256/512-bit records; scratch
// consists of lanes contiguous planes of scratchBlocks() 128-bit blocks.
void transposeScalar(const Block*, Block*, Block* scratch, const Plan&);
void forwardScalar(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
void transposeFast(const Block*, Block*, Block* scratch, const Plan&);

// Internal stage entry points and independent literal outer oracles. Fast
// outer stores stream above K18 when output is 64-byte aligned.
void outerFast(const Block* scratch, Block*, const Plan&);
void outerCached(const Block* scratch, Block*, const Plan&);
void outerScalar(const Block* scratch, Block*, const Plan&);
void outerForwardScalar(const Block*, Block* scratch, const Plan&, unsigned lanes);
bool selfTest() noexcept;
}
