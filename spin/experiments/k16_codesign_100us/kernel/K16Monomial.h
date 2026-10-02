#pragma once
#include "K16CodeSign.h"
namespace spin::research::k16codesign {
struct MonomialTables { std::vector<std::uint64_t> updates; };
// State order: 0,1,2,4,8,16,32,3,17,33,10,34,12,20,36,24.
// Reuses Plan's literal uniform GL16 rows with identity state basis.
void prepareMonomial(const rs::Plan&,MonomialTables&);
void reverseRouteMonomial(const Block*,Block* scratch,const rs::Plan&,const MonomialTables&);
void reverseRouteMonomialScalar(const Block*,Block* scratch,const rs::Plan&);
void forwardInnerMonomialScalar(const Block* routed,Block* encoded,const rs::Plan&);
}
