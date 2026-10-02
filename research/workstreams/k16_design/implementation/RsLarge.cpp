// Import the complete retained inner schedule, without rewriting its hot loops.
// Only route addresses and the subsequent outer belong to the RS construction.
#include "RsPrototype.h"
#include "../../../../spin/src/packet/PacketLargeInner.h"

namespace spin::research::rs {
namespace {
struct StreamingRoute {
    Block* scratch;
    const std::uint32_t* bases;
    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) {
        _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch + bases[packet]), value);
    }
};
}

void reverseRouteLarge(const Block* input, Block* scratch, const Plan& plan) {
    StreamingRoute emit{scratch, plan.route.data()};
    detail::packet::large::reverse(input, plan.n, plan.denseUpdates.data(), emit);
    _mm_sfence();
}

void transposeLarge(const Block* input, Block* output, Block* scratch, const Plan& plan) {
    reverseRouteLarge(input, scratch, plan);
    outerFast(scratch, output, plan);
}
}
