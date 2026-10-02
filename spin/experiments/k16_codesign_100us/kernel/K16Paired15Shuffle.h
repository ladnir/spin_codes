#pragma once
#include "K16Paired15.h"
namespace spin::research::k16codesign {
// Uniform GL15 in compact literal coordinates, with zero fast coordinate 9.
// Keeps the successful shuffle-only state order, except for the omitted row.
void customizePaired15Shuffle(rs::Plan&,std::uint64_t seed,PairedTables&,PairedOptimizedTables&);
// 0: emitted feedback; 1: single-load raw feedback. Reuse paired15 scalar oracles.
void reverseRoutePaired15Shuffle(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
}
