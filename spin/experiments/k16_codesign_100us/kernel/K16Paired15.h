#pragma once
#include "K16PairedBasis.h"
namespace spin::research::k16codesign {
// Uniform GL15 in compact literal row order, omitting paired row 10.
// Installs diag(GL15,1) in Plan.reverseMatrices; coordinate 15 stays zero.
// Fast tables use the same 3<->7 permutation as PairedBasis.
void customizePaired15(rs::Plan&,std::uint64_t seed,PairedTables&,PairedOptimizedTables&);
// 0: feedback from emitted packets; 1: raw-feedback single-load variant.
void reverseRoutePaired15(const Block*,Block* scratch,const rs::Plan&,
    const PairedTables&,const PairedOptimizedTables&,unsigned variant);
void reverseRoutePaired15Scalar(const Block*,Block* scratch,const rs::Plan&);
void forwardInnerPaired15Scalar(const Block* routed,Block* encoded,const rs::Plan&);
}
