#pragma once
#include "../../../src/certified22/Plan.h"

namespace spin::research::twoms::rs8 {
using Block = detail::certified22::Block;
using Plan = detail::certified22::Plan;

// Eight GF(16) RS[8,4] rows per 128-to-256-coordinate group. Each column
// contains 32 contiguous coordinates, ordered 4*row+bit, hence eight packets.
// Reverse symbol maps are multiplication by independent nonzero GF(2^32)
// scalars; the forward maps are their binary adjoints.
struct Tables32 {
    std::vector<std::array<std::uint8_t,9>> field;
    // Column j of the literal reverse map M_c; used only by the scalar oracle.
    std::vector<std::array<std::uint32_t,32>> reverseColumns;
};

// Leaves the Plan's route, inner, and all other fields unchanged. Only its
// 128-to-256 geometry is consumed; the retained certified22 inner stays s=22.
void setupField32(const Plan&, std::uint64_t seed, Tables32&);

// All stage entry points require scratch disjoint from output. outerFast
// uses NT stores above K18 when output is 64-byte aligned; otherwise it uses
// cached unaligned stores. No entry point allocates. Scratch is 64-byte aligned.
void outerFast(const Block* scratch, Block* output, const Plan&, const Tables32&);
void outerCached(const Block* scratch, Block* output, const Plan&, const Tables32&);

// Independent literal Lagrange and binary-matrix oracles.
void outerScalar(const Block* scratch, Block* output, const Plan&, const Tables32&);
void outerForwardScalar(const Block* message, Block* scratch, const Plan&, const Tables32&);

// Requires the same ISA as outerFast. Checks flattened arithmetic against the
// literal scalar field multiplication for basis and deterministic random data.
bool selfTest() noexcept;
}

