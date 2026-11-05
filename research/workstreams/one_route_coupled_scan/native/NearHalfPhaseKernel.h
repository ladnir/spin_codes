#pragma once

#include "ToeplitzWindow1024Nine.h"

// Measurement-only split of the frozen nine-leaf coordinate/raw-scaling
// feedback kernel. This reproduces its operation order; it does not replace
// the encoder or introduce a new construction. The profiling harness must
// compare complete outputs against the frozen encoder and compare timed and
// untimed split runs before interpreting phase measurements.
static_assert(SPIN_OR_WINDOW_COORDINATES == 1,
              "Phase diagnostic requires the coordinate variant");
static_assert(SPIN_OR_WINDOW_RAW_SCALE == 1,
              "Phase diagnostic requires deferred raw scaling");

namespace spin::research::one_route::phase {

using PreparedBatch = ToeplitzWindow1024Nine::PreparedBatch;
using Scratch = ToeplitzWindow1024Nine::Scratch;

#if defined(_MSC_VER)
#define SPIN_OR_PHASE_INLINE __forceinline
#else
#define SPIN_OR_PHASE_INLINE inline __attribute__((always_inline))
#endif

namespace detail {

// Reproduce the two private helpers used by the frozen feedback p-loop.
// Field arithmetic, the eight-point transform, and raw scaling are shared
// with that implementation. All half/node choices remain compile-time.
template<std::size_t... Half>
SPIN_OR_PHASE_INLINE void loadMixed(Element* mixed,
    const Element* const* halves, unsigned p,
    std::index_sequence<Half...>) noexcept {
    (_mm512_store_si512(mixed + 4 * Half,
        _mm512_load_si512(halves[Half] + 4 * p)), ...);
}

template<unsigned Node>
SPIN_OR_PHASE_INLINE void accumulateNodeDeferred(const PreparedBatch& keys,
    const Element* mixed, unsigned p,
    char2::field128_lift::sparse_byte::deferred::Product& all,
    char2::field128_lift::sparse_byte::deferred::Product& bit0,
    char2::field128_lift::sparse_byte::deferred::Product& bit1,
    char2::field128_lift::sparse_byte::deferred::Product& bit2) noexcept {
    namespace deferred = char2::field128_lift::sparse_byte::deferred;
    const auto value = deferred::multiplyRaw4(
        _mm512_load_si512(mixed + 4 * Node),
        _mm512_load_si512(keys.values[p][Node]));
    all = deferred::sum(all, value);
    if constexpr (Node & 1) bit0 = deferred::sum(bit0, value);
    if constexpr (Node & 2) bit1 = deferred::sum(bit1, value);
    if constexpr (Node & 4) bit2 = deferred::sum(bit2, value);
}

} // namespace detail

// state, keys, and scratch must be mutually disjoint. This phase includes
// ring-address selection, all key reads, the eight-point block mixing, nine
// full-field products per leaf frequency, raw scaling, and both reductions.
// It writes both complete spectra without modifying the cached state.
SPIN_OR_PHASE_INLINE void spectra(const ToeplitzWindow1024Nine::State& state,
    const PreparedBatch& keys, Scratch& scratch) noexcept {
    const auto newest = (state.oldest + 3) & 3U;
    const auto preceding = (state.oldest + 2) & 3U;
    const auto earlier = (state.oldest + 1) & 3U;
    const Element* halves[8] = {
        state.transformed[newest] + 256, state.transformed[newest],
        state.transformed[preceding] + 256, state.transformed[preceding],
        state.transformed[earlier] + 256, state.transformed[earlier],
        state.transformed[state.oldest] + 256, state.transformed[state.oldest]
    };
    constexpr auto beta1 = Transform256Cantor::evaluationPoint(2);
    constexpr auto beta2 = Transform256Cantor::evaluationPoint(4);
    for (unsigned p = 0; p != 64; ++p) {
        detail::loadMixed(scratch.mixed, halves, p,
                          std::make_index_sequence<8>{});
        TransformCantor<8>::evaluateFour(scratch.mixed);
        const auto zero = _mm512_setzero_si512();
        namespace deferred = char2::field128_lift::sparse_byte::deferred;
        deferred::Product all{zero, zero, zero, zero};
        deferred::Product bit0{zero, zero, zero, zero};
        deferred::Product bit1{zero, zero, zero, zero};
        deferred::Product bit2{zero, zero, zero, zero};
        detail::accumulateNodeDeferred<0>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<1>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<2>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<3>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<4>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<5>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<6>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        detail::accumulateNodeDeferred<7>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
        bit0 = deferred::sum(bit0, deferred::multiplyRaw4(
            _mm512_load_si512(halves[7] + 4 * p),
            _mm512_load_si512(keys.values[p][8])));
        _mm512_store_si512(scratch.spectra[0] + 4 * p, deferred::reduce(all));
        bit0 = deferred::sum(bit0, window_nine_detail::scaleRaw<beta1>(bit1));
        bit0 = deferred::sum(bit0, window_nine_detail::scaleRaw<beta2>(bit2));
        _mm512_store_si512(scratch.spectra[1] + 4 * p, deferred::reduce(bit0));
    }
}

// Consume fresh spectra once. output addresses256 Elements and must not
// overlap scratch, state, or keys. As in the frozen feedback method, these
// two inverse transforms are followed by extraction of coefficients128..255.
// Routed-input addition and the two append transforms belong to the caller.
SPIN_OR_PHASE_INLINE void finish(Scratch& scratch, Element* output) noexcept {
    Transform256Cantor::interpolateNovel(scratch.spectra[0]);
    Transform256Cantor::interpolateNovel(scratch.spectra[1]);
    for (unsigned j = 0; j != 128; j += 4) {
        _mm512_storeu_si512(output + j,
            _mm512_loadu_si512(scratch.spectra[0] + 128 + j));
        _mm512_storeu_si512(output + 128 + j,
            _mm512_loadu_si512(scratch.spectra[1] + 128 + j));
    }
}

#undef SPIN_OR_PHASE_INLINE
} // namespace spin::research::one_route::phase
