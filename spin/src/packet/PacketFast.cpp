// All AVX-512/VBMI/GFNI implementation details stay in this translation unit.
// Public and setup translation units see only PacketPlan.h's plain-data API.
#include "PacketPlan.h"
#include "PacketInnerFast.h"
#include "PacketOuterFast.h"

namespace spin::detail::packet {
namespace {
struct CachedRoute {
    storage::block* scratch;
    const std::uint32_t* bases;

    SPIN_FORCEINLINE void operator()(std::size_t packet, __m512i value) const {
        _mm512_store_si512(scratch + bases[packet], value);
    }
};
}

void transposeFast(const storage::block* input, storage::block* output,
                   storage::block* scratch, std::size_t n,
                   const std::uint32_t* routeBases,
                   const std::uint64_t* compactCoefficients,
                   const std::uint64_t* composedUpdateWords) {
    // Validation, CPU selection, preparation, and ownership belong to the caller.
    // Scratch is 64-byte aligned and every routed packet fills one cache line.
    // Compact coefficient groups are 32-byte aligned. Input/output may be only
    // block-aligned; the aligned workspace path avoids split cache-line access.
    // Complete the route before output writes, including when output == input.
    CachedRoute emit{scratch, routeBases};
    fast::reverse(input, n, composedUpdateWords, emit);
    for(std::size_t tile = 0; tile < n / 1024; ++tile)
        fast::outerCompact(scratch + 1028 * tile, output + 512 * tile,
                           compactCoefficients + 512 * tile);
}
} // namespace spin::detail::packet
