#pragma once
#include "../../../research/workstreams/k16_design/implementation/RsPrototype.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::research::k18cert22 {
using Block = rs::Block;
struct Tables {
    // Literal binary matrices in the certified (not monomial) state basis.
    std::vector<std::array<std::uint32_t,22>> reverseMatrices;
    // P U P^-1, in output-byte/input-byte order, low row first per uint64.
    // Coordinates 22 and 23 are zero padding, never sampled state bits.
    std::vector<std::uint64_t> updates;
    std::vector<std::uint32_t> routeBytes;
};
void prepare(rs::Plan&, std::uint64_t seed, Tables&);
void reverseRouteScalar(const Block* input, Block* scratch, const rs::Plan&, const Tables&);
void forwardInnerScalar(const Block* scratch, Block* output, const rs::Plan&, const Tables&);
void reverseRouteFast(const Block* input, Block* scratch, const rs::Plan&, const Tables&);
}
