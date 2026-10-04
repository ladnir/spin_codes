#include "Kernels.h"
#include "../../src/paired15/Paired15FastCommon.h"
#include "../../src/paired15/Paired15TransposeOuter.h"

namespace spin::research::rs8sizes {
namespace {
using namespace detail::paired15::fast;
struct ByteRoute {
    std::byte* scratch;
    const std::uint32_t* offsets;
    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) {
        _mm512_store_si512(reinterpret_cast<__m512i*>(scratch+offsets[packet]), value);
    }
};
}

// Exact Paired15Fast.cpp reverseRoute body and unchanged production helpers.
// Only the entry point name/visibility differs; all size adaptation is setup.
SPIN_NOINLINE void reverseRoute15(const Block* __restrict input,
                                 Block* __restrict scratch,
                                 const detail::paired15::Plan& p) {
    const auto epochs = p.n/64;
    ByteRoute emit{reinterpret_cast<std::byte*>(scratch), p.routeBytes.data()};
    alignas(64) __m512i packets[16], high[16];
    PackedState before;
    WordPackets16 state, feedback;
    const auto first = epochs-1;
    loadPackets(input+64*first, packets, std::make_index_sequence<16>{});
    for(unsigned h=16; h-->0;) emit(16*first+h, packets[h]);
    highFromPackets(packets, high);
    finishWide(high, state);
    for(std::size_t epoch=first; --epoch;) {
        packWordGroup<0>(state, before);
        packWordGroup<1>(state, before);
        streamHigh<false,true>(state, input+64*epoch, nullptr, 16*epoch, emit, high);
        finishWide(high, feedback);
        updateWide<false>(before, p.updates.data()+4*epoch, feedback, state);
    }
    streamHigh<false,true>(state, input, nullptr, 0, emit, high);
}

void outer15(const Block* scratch, Block* output, const detail::paired15::Plan& p) {
    detail::paired15::outerfast::run(scratch, output, p);
}
}
