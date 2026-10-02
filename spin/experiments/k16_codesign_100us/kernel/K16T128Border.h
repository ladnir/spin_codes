#pragma once
#include "K16CodeSign.h"
namespace spin::research::k16codesign {
struct T128BorderTables {
    unsigned stateBits=20;
    // Uniform GL19/GL20 in literal expansion coordinates; upper unused bits zero.
    std::vector<std::array<std::uint32_t,20>> reverseMatrices;
    std::vector<std::array<std::uint64_t,9>> packedUpdates;
};
void prepareT128Border(const rs::Plan&,std::uint64_t seed,T128BorderTables&,unsigned stateBits=20);
void reverseRouteT128Border(const Block*,Block* scratch,const rs::Plan&,const T128BorderTables&);
void reverseRouteT128BorderScalar(const Block*,Block* scratch,const rs::Plan&,const T128BorderTables&);
void forwardInnerT128BorderScalar(const Block* routed,Block* encoded,const rs::Plan&,const T128BorderTables&);
}
