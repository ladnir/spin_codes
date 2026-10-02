#pragma once
#include "K16PairedShuffle.h"
namespace spin::research::k16codesign {
// Same GL16 raw-shuffle path/tables; makes the routing buffer's non-overlap
// contract explicit to the compiler. Input, scratch, and tables must be disjoint.
// 0: original loop with restrict pointers. 1: additionally peel both boundaries.
void reverseRoutePairedNoAlias(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
}
