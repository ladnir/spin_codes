#pragma once
#include "K18Certified22.h"
#include "K20Route.h"
#include "K20Outer.h"
#include "outer_variants/OuterVariants.h"

namespace spin::research::k20small {
// Compile-time experiment selectors. Route: retained, streaming, packed streaming.
// Outer: retained, streaming, split finish, unrolled mixing. No packet-level dispatch.
template<unsigned Kind>
inline void route(const rs::Block* input, rs::Block* scratch,
                  const rs::Plan& plan, const k18cert22::Tables& tables) {
    static_assert(Kind <= 3);
    if constexpr (Kind == 3) reverseRoutePackedCached(input, scratch, plan, tables);
    else if constexpr (Kind == 2) reverseRoutePackedNt(input, scratch, plan, tables);
    else if constexpr (Kind == 1) reverseRouteNt(input, scratch, plan, tables);
    else k18cert22::reverseRouteFast(input, scratch, plan, tables);
}
template<unsigned Kind>
inline void finish(const rs::Block* scratch, rs::Block* output,
                   const rs::Plan& plan, const k16codesign::NativeOuterTables& tables) {
    static_assert(Kind <= 4);
    if constexpr (Kind == 4) k16codesign::outervariants::fieldSharedParity(scratch, output, plan, tables);
    else if constexpr (Kind == 3) fieldUnrolledSharedParityNt(scratch, output, plan, tables);
    else if constexpr (Kind == 2) fieldSplitSharedParityNt(scratch, output, plan, tables);
    else if constexpr (Kind == 1) fieldLoopSharedParityNt(scratch, output, plan, tables);
    else k16codesign::outervariants::fieldLoopSharedParity(scratch, output, plan, tables);
}
}
