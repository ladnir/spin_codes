#pragma once
#include "../k16_codesign_100us/kernel/K16NativeOuter.h"

namespace spin::research::k20small {
using rs::Block;
using k16codesign::NativeOuterTables;

// Exact fieldLoopSharedParity map and schedule. A 64-byte-aligned output uses
// non-temporal stores and one completion SFENCE. Otherwise this calls the
// retained cached kernel. Scratch/output non-overlap requirements are unchanged.
void fieldLoopSharedParityNt(const Block* scratch,Block* output,
    const rs::Plan&,const NativeOuterTables&);
// Same map/store contract, with separate noinline finish calls per plane pair.
void fieldSplitSharedParityNt(const Block* scratch,Block* output,
    const rs::Plan&,const NativeOuterTables&);
// Same map/store contract, with the sixteen symbol-mix calls explicitly unrolled.
void fieldUnrolledSharedParityNt(const Block* scratch,Block* output,
    const rs::Plan&,const NativeOuterTables&);
}
