#pragma once
#include "../../src/paired15/Paired15.h"
#include "../../src/certified22/Plan.h"

namespace spin::research::rs8sizes {
using Block = detail::certified22::Block;
// Only needed outside K16, whose production paired15 Plan is used unchanged.
void resizePaired15(detail::paired15::Plan&, const detail::certified22::Plan&,
                    std::uint64_t innerSeed);
void reverseRoute15(const Block*, Block*, const detail::paired15::Plan&);
void outer15(const Block*, Block*, const detail::paired15::Plan&);
}
