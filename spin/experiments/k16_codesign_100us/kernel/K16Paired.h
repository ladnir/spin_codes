#pragma once
#include "K16CodeSign.h"
namespace spin::research::k16codesign {
struct PairedTables { std::vector<std::uint64_t> updates; };
// State: constant + six linear terms + the nine disjoint groups in the
// proof proposal disjoint-pairs-proposal.json. Only six feedback XORs.
// Reuses Plan's literal uniform GL16 rows with identity state basis.
void preparePaired(const rs::Plan&,PairedTables&);
void reverseRoutePaired(const Block*,Block* scratch,const rs::Plan&,const PairedTables&);
void reverseRoutePairedScalar(const Block*,Block* scratch,const rs::Plan&);
void forwardInnerPairedScalar(const Block* routed,Block* encoded,const rs::Plan&);
}


