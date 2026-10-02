#pragma once
#include "../kernel/K16NativeOuter.h"

namespace spin::research::k16codesign::outervariants {
// Four allocation-free, compile-time-specialized alternatives for the exact
// outerNativeField16 map. All consume unchanged NativeOuterTables::field.
// Scratch/output alignment and non-overlap requirements are unchanged.
void fieldInline(const Block*, Block*, const rs::Plan&, const NativeOuterTables&);
void fieldUnrolled(const Block*, Block*, const rs::Plan&, const NativeOuterTables&);
void fieldSharedParity(const Block*, Block*, const rs::Plan&, const NativeOuterTables&);
void fieldBatch2(const Block*, Block*, const rs::Plan&, const NativeOuterTables&);
// Isolate the shared-parity circuit without fully unrolling symbol mixing.
void fieldLoopSharedParity(const Block*, Block*, const rs::Plan&, const NativeOuterTables&);
}
