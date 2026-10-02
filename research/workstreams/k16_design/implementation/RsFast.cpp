#include "RsPrototype.h"
#include "../../../../spin/src/packet/PacketInnerFast.h"
#include "RsOuterFast.h"
#include "Rs16OuterFast.h"

namespace spin::research::rs {
namespace {
struct CachedRoute {
    Block* scratch;
    const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) const {
        _mm512_store_si512(scratch + bases[packet], value);
    }
};
}

void reverseRouteFast(const Block* input, Block* scratch, const Plan& plan) {
    CachedRoute emit{scratch, plan.route.data()};
    detail::packet::fast::reverse(input, plan.n, plan.composedUpdates.data(), emit);
}

void outerFast(const Block* scratch, Block* output, const Plan& plan) {
    const auto* coefficients = plan.compactCoefficients();
    if(plan.variant == Variant::Rs8Gf256) {
        for(std::size_t group = 0; group < plan.groups; ++group)
            fast::outerGroup(scratch + groupStride * group, output + 128 * group,
                             coefficients + 128 * group);
    } else {
        for(std::size_t group = 0; group < plan.groups; ++group)
            fast16::outerGroup(scratch + groupStride * group, output + 128 * group,
                               coefficients + 64 * group);
    }
}

void transposeFast(const Block* input, Block* output, Block* scratch, const Plan& plan) {
    reverseRouteFast(input, scratch, plan);
    outerFast(scratch, output, plan);
}
}
