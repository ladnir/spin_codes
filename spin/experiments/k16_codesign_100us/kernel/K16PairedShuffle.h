#pragma once
#include "K16PairedBasis.h"
namespace spin::research::k16codesign {
// Same certified paired code and routes, with a different internal state order.
// Uses immediate 128-bit-lane feedback shuffles, rather than variable permutes.
void preparePairedShuffle(const rs::Plan&,PairedTables&,PairedOptimizedTables&);
void customizePairedField16Shuffle(rs::Plan&,std::uint64_t seed,PairedTables&,PairedOptimizedTables&);
// 0/1: emitted-feedback dense/field. 2/3: single-load raw-feedback dense/field.
// Existing paired scalar oracles consume the unchanged literal coordinates.
void reverseRoutePairedShuffle(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
}
