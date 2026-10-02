#pragma once
#include "K16PairedOptimized.h"
namespace spin::research::k16codesign {
// Swap literal state coordinates 3 and 7. This places (1,x0,x1,x0*x1)
// together, giving one four-lane butterfly for packet0's ANF polynomial.
// No encoded bit, route address, or sampled literal matrix changes.
void preparePairedBasis(const rs::Plan&,PairedTables& basisTables,PairedOptimizedTables&);
// Uniform field family conjugated into the literal paired state basis.
// Updates Plan literal matrices for the existing scalar oracle.
void customizePairedField16Basis(rs::Plan&,std::uint64_t seed,PairedTables& basisTables,PairedOptimizedTables&);
// 0: dense; 1: six-GFNI field. Tables must match the preparation above.
void reverseRoutePairedBasis(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables& basisTables,const PairedOptimizedTables&,unsigned variant);
}
