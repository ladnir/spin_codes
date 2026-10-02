#pragma once
#include "K16PairedBasis.h"
namespace spin::research::k16codesign {
// Same construction, tables, scalar oracle and coordinate basis as PairedBasis.
// Variant 0: dense GL16. Variant 1: transitive six-GFNI field multiplier.
// Reduce raw-input feedback while loading once for expansion and emission.
void reverseRoutePairedBasisRaw(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables& basisTables,const PairedOptimizedTables&,unsigned variant);
}
