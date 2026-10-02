#pragma once
#include "K16Paired.h"
namespace spin::research::k16codesign {
struct PairedOptimizedTables {
    std::vector<std::uint32_t> routeBytes;
    // Three ascending-row GF8 multiplication matrices per physical step.
    std::vector<std::uint64_t> fieldUpdates;
};
void preparePairedOptimized(const rs::Plan&,PairedOptimizedTables&);
// Changes only the independent GL16 family to uniform nonzero GF(2^16)
// multiplication. Updates literal rows and the ordinary paired table so both
// existing scalar and dense kernels remain independent implementation oracles.
void customizePairedField16(rs::Plan&,std::uint64_t seed,PairedTables&,PairedOptimizedTables&);
// 0: byte route; 1: raw feedback + early next-state; 2: explicit XOR3 circuit;
// 3: three-map field update; 4: raw feedback + three-map field update.
// 3/4 require customizePairedField16; 0/1/2 accept either prepared family.
void reverseRoutePairedOptimized(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
}
