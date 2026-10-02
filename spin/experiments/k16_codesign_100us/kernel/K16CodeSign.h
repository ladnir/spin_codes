#pragma once
#include "../../../../research/workstreams/k16_design/implementation/RsPrototype.h"

// Isolated K65536 construction/implementation experiment. No production dispatch.
namespace spin::research::k16codesign {
using rs::Block;
struct Tables {
    std::vector<std::array<std::uint8_t,3>> outerField;
    std::vector<std::array<std::uint8_t,3>> innerField;
};

// Require Rs16Gf16; the inner customizer also requires RetainedStreaming.
// These also update the literal matrices
// and dense coefficients, allowing the retained scalar and GL kernels to serve
// as independent oracles for the new field-multiplication implementation.
void customizeOuterField16(rs::Plan&, std::uint64_t seed, Tables&);
void customizeInnerField16(rs::Plan&, std::uint64_t seed, Tables&);
void reverseRouteCachedGl16(const Block*, Block* scratch, const rs::Plan&);
void reverseRouteCachedField16(const Block*, Block* scratch, const rs::Plan&, const Tables&);
void outerFusedGl16(const Block* scratch, Block*, const rs::Plan&);
void outerFusedField16(const Block* scratch, Block*, const rs::Plan&, const Tables&);
}
