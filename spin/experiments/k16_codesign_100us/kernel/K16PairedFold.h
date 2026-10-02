#pragma once
#include "K16PairedNoAlias.h"
#include "K16Paired15NoAlias.h"
namespace spin::research::k16codesign {
// Exact paired peeled/raw kernels with a lower-shuffle four-moment reduction.
// Reuse shuffle-basis setup and each dimension's existing scalar oracle.
void reverseRoutePairedFold(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&);
void reverseRoutePaired15Fold(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&);
}
