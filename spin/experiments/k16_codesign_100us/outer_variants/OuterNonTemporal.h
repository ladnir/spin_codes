#pragma once
#include "OuterVariants.h"
namespace spin::research::k16codesign::outervariants {
// Exact shared-loop map. Uses NT stores only for64-byte-aligned output;
// otherwise falls back to ordinary stores. The required SFENCE is included.
void fieldLoopSharedNt(const Block*,Block*,const rs::Plan&,const NativeOuterTables&);
}

