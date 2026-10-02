#pragma once
#include "K16Paired15Shuffle.h"
namespace spin::research::k16codesign {
// Same GL15 raw-shuffle path/tables, explicit disjoint-buffer contract.
// 0: restrict only. 1: also peel initial/final steps out of the main loop.
// Reuse customizePaired15Shuffle and paired15 scalar oracles.
void reverseRoutePaired15NoAlias(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
}
