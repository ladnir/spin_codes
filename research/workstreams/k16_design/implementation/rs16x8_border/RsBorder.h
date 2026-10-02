#pragma once
#include "../rs16x8/RsWide.h"
namespace spin::research::rsborder {
using Block=rswide::Block;
struct BorderRow {
    detail::packet::large::Dense16Row base;
    std::uint64_t upper[2]{},lower[2]{},corner=0;
};
struct Plan {
    Plan(std::size_t k,std::uint64_t routeSeed,std::uint64_t innerSeed,unsigned stateBits,
         rswide::Randomizer randomizer=rswide::Randomizer::TowerByte32);
    rswide::Plan outer;
    unsigned stateBits;
    std::vector<std::array<std::uint32_t,20>> reverseMatrices;
    std::vector<BorderRow> updates;
};
void reverseRoute(const Block*,Block*,const Plan&);
void transposeFast(const Block*,Block*,Block*,const Plan&);
void transposeScalar(const Block*,Block*,Block*,const Plan&);
void forwardScalar(const Block*,Block*,const Plan&);
}
