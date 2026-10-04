#pragma once
#include "../k18_certified22/K18Certified22.h"

namespace spin::research::k20small {
void reverseRoutePackedCached(const rs::Block* input,rs::Block* scratch,
                             const rs::Plan& plan,const k18cert22::Tables& tables);
// Exact certified s22 reverse/route with non-temporal stores. Input holds
// plan.n records; scratch holds plan.scratchBlocks() records, is 64-byte
// aligned, and is disjoint from input. Tables come from k18cert22::prepare.
// The required store fence completes before this function returns.
void reverseRouteNt(const rs::Block* input,rs::Block* scratch,
    const rs::Plan& plan,const k18cert22::Tables& tables);
// Same map and storage contract; retain the quadratic groups' interleaved
// update output between steps to reduce cross-vector state repacking.
void reverseRoutePackedNt(const rs::Block* input,rs::Block* scratch,
    const rs::Plan& plan,const k18cert22::Tables& tables);
// Untimed AVX-512/GFNI check of the new packing conversion against the
// original unpack-then-pack path on basis bits, random values, and padding.
bool packedCarrySelfTest() noexcept;
}
