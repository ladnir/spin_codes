#pragma once
#include "K16CodeSign.h"
namespace spin::research::k16codesign {
struct T128Tables { std::vector<detail::packet::large::Dense16Row> updates; };
// Uses the first n/128 independent GL16 matrices from the retained Plan,
// expressed in the literal prefix16 basis of the pinned t128/s20 expansion.
// Plan route, outer, and old dense tables remain unchanged.
void prepareT128(const rs::Plan&,T128Tables&);
void reverseRouteT128(const Block*,Block* scratch,const rs::Plan&,const T128Tables&);
void reverseRouteT128Scalar(const Block*,Block* scratch,const rs::Plan&);
// Validation-only literal forward inner; input is group-order routed scratch.
void forwardInnerT128Scalar(const Block* routed,Block* encoded,const rs::Plan&);
}
