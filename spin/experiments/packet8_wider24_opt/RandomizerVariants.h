#pragma once
#include "../packet8_wider24/Packet8Wide24.h"
#include <array>
#include <cstdint>
#include <vector>

namespace spin::research::packet8wide24opt {
using Block=packet8wide24::Block;
using Plan=packet8wide24::Plan;

struct RandomizerPlan {
    explicit RandomizerPlan(const Plan&);
    std::vector<std::array<std::uint64_t,9>> affine9;
    std::vector<std::array<std::uint64_t,16>> direct16;
};

// Variants 0..3 apply exactly the baseline Plan's outer map.
// 0: frozen baseline; 1: flat9 MUL with ternary XORs;
// 2: flat9 AFFINE; 3: direct16 AFFINE blocks from literal outerRows.
// 4: flat9 MUL, but with symbol 0's randomizer omitted in every group.
// Dispatch is outside the hot group loop; no allocations occur in encoding.
void outerRandomizer(const Block* scratch,Block* output,const Plan&,
                     const RandomizerPlan&,unsigned variant);

// Setup-time algebra check: all32 basis inputs and two dense words per
// configured symbol. This is not called in the hot encoder.
bool randomizerTablesMatch(const Plan&,const RandomizerPlan&);

// Variant 4 uses baseline tables for symbols 1..15. Apply this to an oracle
// Plan copy only; its outerCoefficients intentionally remain unchanged.
void applyOmittedRows(Plan&);
}
