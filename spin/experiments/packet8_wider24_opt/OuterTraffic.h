#pragma once
#include "../packet8_wider24/Packet8Wide24.h"

namespace spin::research::packet8wide24opt {
using Block=packet8wide24::Block;
using Plan=packet8wide24::Plan;

// Same map and buffer contract as packet8wide24::outerFast.
// 0: frozen baseline, 8 KiB mixed-symbol workspace.
// 1: one payload half at a time, 4 KiB workspace; masked output stores.
// 2: parity first, then fused systematic mix/output, 4 KiB workspace.
// 3: parity first within each payload half, 2 KiB workspace; masked stores.
// 4: parity first with the exact flat9 MUL circuit, 4 KiB workspace.
// 5: as 4, but symbol 0 is identity in every group (changed sampled code).
// For variant 5, compare against an oracle Plan copy with identity outerRows
// at symbol 0; do not change the stored coefficients for symbols 1..15.
// Dispatch occurs once per call, never inside a group or symbol loop.
void outerTraffic(const Block* scratch,Block* output,const Plan&,unsigned variant);
}
