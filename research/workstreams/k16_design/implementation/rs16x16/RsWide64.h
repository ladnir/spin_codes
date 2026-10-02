#pragma once
#include "../RsPrototype.h"

namespace spin::research::rswide64 {
using Block=rs::Block;
inline constexpr std::size_t groupStride=1028;
struct Plan {
    Plan(std::size_t messageBits,std::uint64_t routeSeed,std::uint64_t innerSeed);
    std::size_t k=0,n=0,groups=0;
    std::vector<std::uint32_t> route;
    std::vector<std::array<std::uint16_t,16>> reverseMatrices;
    std::vector<detail::packet::large::Dense16Row> denseUpdates;
    std::vector<std::array<std::uint64_t,64>> outerMatrices;
    std::vector<std::uint8_t> coefficientBytes;
    std::size_t scratchBlocks() const noexcept {return groups*groupStride;}
};
// Sixteen parallel GF16 RS[16,8] rows, independently sampled adjoints of
// uniform nonzero GF(2^64) multipliers on aligned 64-bit symbols. Each group
// shuffles its 256 four-bit packets; each region independently shuffles all
// groups. Uses the retained t64/s16 inner. No new numerical proof claim.
void reverseRoute(const Block*,Block*,const Plan&);
void outerFast(const Block*,Block*,const Plan&);
void transposeFast(const Block*,Block*,Block*,const Plan&);
void transposeScalar(const Block*,Block*,Block*,const Plan&);
void forwardScalar(const Block*,Block*,const Plan&);
void outerScalar(const Block*,Block*,const Plan&);
}
