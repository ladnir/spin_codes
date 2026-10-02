#pragma once
#include "K16PairedOptimized.h"
namespace spin::research::k16codesign {
// Exact paired feedback remains in four ordinary-word ZMMs throughout.
// 0: dense update; 1: field update; 2: raw-feedback dense; 3: raw-feedback field.
// Field variants require customizePairedField16; byte routes always prepared.
void reverseRoutePairedWide(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
}
