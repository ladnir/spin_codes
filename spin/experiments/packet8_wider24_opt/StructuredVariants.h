#pragma once
#include "RandomizerVariants.h"

namespace spin::research::packet8wide24opt {
// Byte field F16[u]/(u^2+u+8), with F16=F2[x]/(x^4+x+1).
// R = D_out H D_in, where H is the four-byte MixColumns matrix,
// D_in[0]=1 and the other seven diagonal entries are uniform nonzero.
struct StructuredPlan {
    StructuredPlan(const Plan&,std::uint64_t seed);
    // First three entries are input scales 1..3, final four output scales.
    std::vector<std::array<std::uint8_t,7>> coefficients;
    std::vector<std::array<std::uint64_t,7>> diagonal;
    std::vector<std::array<std::uint32_t,32>> reverseRows;
};

// Use a copy of the baseline Plan for literal forward/transpose checks.
// Its baseline outerCoefficients are deliberately not used by this kernel.
void applyStructuredRows(Plan&,const StructuredPlan&);
bool structuredTablesMatch(const StructuredPlan&);

// 0: seven random affine products + three fixed affine products for H.
// 1: seven random affine products + nibble xtime shifts/XORs for H.
// Allocation and variant dispatch are outside the hot group loop.
void outerStructured(const Block* scratch,Block* output,const Plan&,
                     const StructuredPlan&,unsigned variant);
}
