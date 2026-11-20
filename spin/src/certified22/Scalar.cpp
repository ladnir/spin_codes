#include "Plan.h"
#include "productionMap.h"
#include <bit>
#include <cstring>

namespace spin::detail::certified22 {
namespace {
template<bool Forward>
void inner(const Block* input, Block* output, const Plan& p, unsigned lanes, unsigned lane) {
    Block state[22]{}, next[22], feedback[22];
    for(std::size_t step = 0; step < p.n/64; ++step) {
        const auto epoch = Forward ? step : p.n/64 - 1 - step;
        for(auto& value : feedback) value = Block::zero();
        for(unsigned j = 0; j < 64; ++j) {
            const auto i = 64*epoch+j, routed = p.route[i/4]+(i&3);
            const auto raw = input[Forward ? routed : i];
            auto value = raw;
            for(auto mask = map::oldColumns[j]; mask; mask &= mask-1) {
                const auto c = std::countr_zero(mask);
                value = (value^ state[c]);
                feedback[c] = (feedback[c]^ raw);
            }
            output[Forward ? i*lanes+lane : routed] = Block(value);
        }
        if(step+1 == p.n/64) break;
        for(unsigned j = 0; j < 22; ++j) {
            auto value = feedback[j];
            if constexpr(Forward) {
                for(unsigned c = 0; c < 22; ++c)
                    if((p.reverseMatrices[epoch][c] >> j) & 1U)
                        value = (value^ state[c]);
            } else {
                for(unsigned mask = p.reverseMatrices[epoch][j]; mask; mask &= mask-1)
                    value = (value^ state[std::countr_zero(mask)]);
            }
            next[j] = value;
        }
        std::memcpy(state, next, sizeof(state));
    }
}
}
void reverseRouteScalar(const Block* input, Block* scratch, const Plan& p) {
    inner<false>(input, scratch, p, 1, 0);
}
void forwardInnerScalar(const Block* scratch, Block* output, const Plan& p, unsigned lanes, std::size_t plane) {
    if(!plane) plane = p.scratchBlocks();
    for(unsigned lane = 0; lane < lanes; ++lane)
        inner<true>(scratch+lane*plane, output, p, lanes, lane);
}
}
