#pragma once
#include "HdBlockwiseSharedInner.h"

namespace spin::research::char2::hdbwrandom {

// Experimental randomized outer = D * GRS[256,128], followed by the existing
// regional permutation and beta64 shared inner. D has one independent nonzero
// GF65536 tower-field scalar per outer output RECORD, shared by all eight
// field lanes of that record. The caller samples/prepares D once at setup.
// D[c*rows+row] is indexed BEFORE the regional route; after routing the same
// map uses the conjugated scalar D[c*rows+route[j]]. This is not a distance
// certificate, a sampler for cryptographic keys, or a GF128 multiplier.
//
// Fuse D into the already-required sequential region staging. The fixed-width
// four-vector schedule matches the old copy control; the identity template
// removes coefficient accesses and arithmetic entirely. count is a multiple
// of 16 and payload/coefficient loads may be unaligned. Exact payload alias
// is supported; coefficient storage must be disjoint from output.
template<bool Diagonal>
inline void stageRegion(const Block* source, const wide16::Coefficient* diagonal,
                        Block* output, std::size_t count) noexcept {
    for (std::size_t h = 0; h != count; h += 16) {
        auto a = _mm512_loadu_si512(source + h);
        auto b = _mm512_loadu_si512(source + h + 4);
        auto c = _mm512_loadu_si512(source + h + 8);
        auto d = _mm512_loadu_si512(source + h + 12);
        if constexpr (Diagonal) {
            a = wide16::multiply512(a, diagonal + h);
            b = wide16::multiply512(b, diagonal + h + 4);
            c = wide16::multiply512(c, diagonal + h + 8);
            d = wide16::multiply512(d, diagonal + h + 12);
        }
        _mm512_storeu_si512(output + h, a);
        _mm512_storeu_si512(output + h + 4, b);
        _mm512_storeu_si512(output + h + 8, c);
        _mm512_storeu_si512(output + h + 12, d);
    }
}

} // namespace spin::research::char2::hdbwrandom
