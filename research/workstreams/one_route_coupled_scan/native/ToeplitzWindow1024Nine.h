#pragma once

#include "ToeplitzWindow1024.h"

#ifndef SPIN_OR_WINDOW_RAW_SCALE
#define SPIN_OR_WINDOW_RAW_SCALE 0
#endif

namespace spin::research::one_route {

#if defined(_MSC_VER)
#define SPIN_OR_NINE_INLINE __forceinline
#else
#define SPIN_OR_NINE_INLINE inline __attribute__((always_inline))
#endif

namespace window_nine_detail {

// The first eight Cantor points form V3, whose vanishing polynomial is
// Z8(X)=X^8+X^4+X^2+X. Its derivative is1. Consequently the Lagrange
// polynomial for point a is Z8(X)/(X+a), with the coefficients below.
// These are the ADJOINT interpolation weights, not evaluations of the key.
constexpr auto lagrangeWeights() noexcept {
    std::array<std::array<std::uint8_t, 8>, 8> result{};
    for (unsigned node = 0; node != 8; ++node) {
        const auto a = Transform256Cantor::evaluationPoint(node);
        result[node][7] = 1;
        for (unsigned degree = 7; degree != 0; --degree) {
            const auto coefficient = static_cast<std::uint8_t>(
                degree == 1 || degree == 2 || degree == 4);
            result[node][degree - 1] = coefficient ^
                transform_detail::mul(a, result[node][degree]);
        }
    }
    return result;
}
inline constexpr auto weights = lagrangeWeights();

template<unsigned Coefficient>
SPIN_OR_NINE_INLINE __m512i scale(__m512i value) noexcept {
    if constexpr (Coefficient == 0) return _mm512_setzero_si512();
    else if constexpr (Coefficient == 1) return value;
    else return _mm512_gf2p8mul_epi8(value,
        _mm512_set1_epi8(static_cast<char>(Coefficient)));
}

template<unsigned Coefficient>
SPIN_OR_NINE_INLINE char2::field128_lift::sparse_byte::deferred::Product scaleRaw(
    char2::field128_lift::sparse_byte::deferred::Product value) noexcept {
    // A16-bit slot reduces to lo+0x1b*hi in AES. Scaling BOTH bytes by
    // Coefficient therefore commutes with AES and sparse-extension reduction.
    // The high byte may now use its top bit; reduceAndPack also handles that
    // degree15 representative. Only XOR and reduction follow this operation.
    return {scale<Coefficient>(value.p0), scale<Coefficient>(value.p1),
            scale<Coefficient>(value.p2), scale<Coefficient>(value.p3)};
}

template<bool Tail>
SPIN_OR_NINE_INLINE __m512i loadKey(const Element* raw, unsigned offset) noexcept {
    // Tail is coefficients252..254 and a zero fourth coefficient. The final
    // source leaf ends at raw[1278], so even an unused load of1279 is invalid.
    if constexpr (Tail) return _mm512_maskz_loadu_epi64(0x3f, raw + offset);
    else return _mm512_loadu_si512(raw + offset);
}

template<unsigned Node, unsigned Leaf, bool Tail>
SPIN_OR_NINE_INLINE void addKey(const Element* raw, unsigned offset, __m512i& sum) noexcept {
    constexpr auto coefficient = weights[Node][Leaf];
    if constexpr (coefficient != 0)
        sum = _mm512_xor_si512(sum,
            scale<coefficient>(loadKey<Tail>(raw, 128 * Leaf + offset)));
}

template<unsigned Node, bool Tail, std::size_t... Leaf>
SPIN_OR_NINE_INLINE __m512i nodeKey(const Element* raw, unsigned offset,
                                   std::index_sequence<Leaf...>) noexcept {
    auto result = _mm512_setzero_si512();
    (addKey<Node, Leaf, Tail>(raw, offset, result), ...);
    return result;
}

template<bool Tail>
SPIN_OR_NINE_INLINE __m512i infinityKey(const Element* raw, unsigned offset) noexcept {
    const auto a = _mm512_xor_si512(loadKey<Tail>(raw, 8 * 128 + offset),
                                   loadKey<Tail>(raw, 4 * 128 + offset));
    const auto b = _mm512_xor_si512(loadKey<Tail>(raw, 2 * 128 + offset),
                                   loadKey<Tail>(raw, 128 + offset));
    return _mm512_xor_si512(a, b);
}

} // namespace window_nine_detail

// Exact rank-nine block factorization of the SAME256x1024 Toeplitz matrix.
// Define128-square leaves A_l[t]=raw[128*l+t], l=0..8,t=0..254. Reverse
// the eight128-symbol state halves, z_l=x_(7-l). The two output halves are
// y_i=sum_l A_(i+l) z_l, i=0,1: the adjoint of multiplying
// S(X)=sum_l z_l X^l by u0+u1 X. Evaluate that degree-eight product at the
// eight Cantor points and retain its leading coefficient separately.
//
// Let W_a=sum_l<8 [X^l](Z8/(X+a))*A_l and W_inf=A8+A4+A2+A1. Then
// y0=sum_a W_a*S(a), y1=sum_a a*W_a*S(a)+W_inf*z7.
// All sums/scalings here are over GF256. They commute with the exact128-square
// leaf kernels. No newly sampled spectral coefficients or distribution
// substitution occurs: every W is computed from the original1279 diagonals.
// With SPIN_OR_WINDOW_COORDINATES=1, the explicitly labeled alternate map is
// diag(U,U)*T*diag(V,...,V), as in the12-leaf control. That is a different
// matrix distribution with the same inverse-uniformity proof interface, not
// an implementation of the original Toeplitz map on the same coordinates.
struct ToeplitzWindow1024Nine {
    using Field = ToeplitzWindow1024::Field;
    static constexpr unsigned memory = 1024;
    static constexpr unsigned batch = 256;
    static constexpr unsigned rawCoefficients = 1279;
    static constexpr unsigned preparedCoefficients = 9 * 256;
    static constexpr unsigned productsPerBatch = preparedCoefficients;

    struct alignas(64) PreparedBatch {
        Field values[64][9][4];
    };

    struct alignas(64) Scratch {
        Field spectra[2][256];
        // Four independent degree-seven polynomials in point-major layout.
        // Reuse this512-byte tile instead of keeping eight vectors live during
        // the full-field products or materializing eight whole mixed spectra.
        Field mixed[8 * 4];
    };

private:
    template<unsigned Node>
    static void prepareNode(const Field* raw, PreparedBatch& destination,
                            Field* temporary) noexcept {
        for (unsigned k = 0; k != 252; k += 4)
            _mm512_store_si512(temporary + k,
                window_nine_detail::nodeKey<Node, false>(raw, k, std::make_index_sequence<8>{}));
        _mm512_store_si512(temporary + 252,
            window_nine_detail::nodeKey<Node, true>(raw, 252, std::make_index_sequence<8>{}));
        Transform256Cantor::evaluate(temporary);
        ToeplitzWindow1024::prepareSpectrum(temporary);
        for (unsigned p = 0; p != 64; ++p)
            _mm512_store_si512(destination.values[p][Node], _mm512_load_si512(temporary + 4 * p));
    }

    template<std::size_t... Node>
    static void prepareNodes(const Field* raw, PreparedBatch& destination,
                             Field* temporary, std::index_sequence<Node...>) noexcept {
        (prepareNode<Node>(raw, destination, temporary), ...);
    }

    template<bool Gfni>
    static SPIN_OR_NINE_INLINE __m512i product(__m512i value, __m512i key) noexcept {
        if constexpr (Gfni) return char2::field128_lift::sparse_byte::multiply4Gfni(value, key);
        else return char2::field128_lift::sparse_byte::multiply4(value, key);
    }

    template<std::size_t... Half>
    static SPIN_OR_NINE_INLINE void loadMixed(Field* mixed, const Field* const* halves,
                                             unsigned p, std::index_sequence<Half...>) noexcept {
        (_mm512_store_si512(mixed + 4 * Half, _mm512_load_si512(halves[Half] + 4 * p)), ...);
    }

    template<unsigned Node, bool Gfni>
    static SPIN_OR_NINE_INLINE void accumulateNode(const PreparedBatch& keys,
        const Field* mixed, unsigned p, __m512i& all, __m512i& bit0,
        __m512i& bit1, __m512i& bit2) noexcept {
        const auto value = product<Gfni>(_mm512_load_si512(mixed + 4 * Node),
            _mm512_load_si512(keys.values[p][Node]));
        all = _mm512_xor_si512(all, value);
        if constexpr (Node & 1) bit0 = _mm512_xor_si512(bit0, value);
        if constexpr (Node & 2) bit1 = _mm512_xor_si512(bit1, value);
        if constexpr (Node & 4) bit2 = _mm512_xor_si512(bit2, value);
    }

    template<unsigned Node>
    static SPIN_OR_NINE_INLINE void accumulateNodeDeferred(const PreparedBatch& keys,
        const Field* mixed, unsigned p,
        char2::field128_lift::sparse_byte::deferred::Product& all,
        char2::field128_lift::sparse_byte::deferred::Product& bit0,
        char2::field128_lift::sparse_byte::deferred::Product& bit1,
        char2::field128_lift::sparse_byte::deferred::Product& bit2) noexcept {
        namespace deferred = char2::field128_lift::sparse_byte::deferred;
        const auto value = deferred::multiplyRaw4(_mm512_load_si512(mixed + 4 * Node),
            _mm512_load_si512(keys.values[p][Node]));
        all = deferred::sum(all, value);
        if constexpr (Node & 1) bit0 = deferred::sum(bit0, value);
        if constexpr (Node & 2) bit1 = deferred::sum(bit1, value);
        if constexpr (Node & 4) bit2 = deferred::sum(bit2, value);
    }

public:
    // Input has exactly1279 field values and must not overlap destination.
    // Preparation explicitly zero-pads each leaf coefficient255, including
    // the tail masked loads. The256-point product then has no alias at127.
    static void prepare(const Field* raw, PreparedBatch& destination) noexcept {
        alignas(64) Field temporary[256];
        prepareNodes(raw, destination, temporary, std::make_index_sequence<8>{});
        for (unsigned k = 0; k != 252; k += 4)
            _mm512_store_si512(temporary + k, window_nine_detail::infinityKey<false>(raw, k));
        _mm512_store_si512(temporary + 252, window_nine_detail::infinityKey<true>(raw, 252));
        Transform256Cantor::evaluate(temporary);
        ToeplitzWindow1024::prepareSpectrum(temporary);
        for (unsigned p = 0; p != 64; ++p)
            _mm512_store_si512(destination.values[p][8], _mm512_load_si512(temporary + 4 * p));
    }

    // Reuse only the frozen cache/ring implementation: no virtual dispatch.
    // reset(), load(raw1024), append(out256), and the point ordering are
    // identical to the12-leaf control. Two forward transforms per append.
    struct alignas(64) State : ToeplitzWindow1024::State {
        template<bool Gfni = false, bool Deferred = false>
        void feedback(const PreparedBatch& keys, Field* output, Scratch& scratch) const noexcept {
            static_assert(!(Gfni && Deferred), "Deferred is the CLMUL arithmetic alternative");
            // output, keys, this State, and scratch must be mutually disjoint.
            const auto newest = (oldest + 3) & 3U;
            const auto preceding = (oldest + 2) & 3U;
            const auto earlier = (oldest + 1) & 3U;
            const Field* halves[8] = {
                transformed[newest] + 256, transformed[newest],
                transformed[preceding] + 256, transformed[preceding],
                transformed[earlier] + 256, transformed[earlier],
                transformed[oldest] + 256, transformed[oldest]
            };
            constexpr auto beta1 = Transform256Cantor::evaluationPoint(2);
            constexpr auto beta2 = Transform256Cantor::evaluationPoint(4);
            for (unsigned p = 0; p != 64; ++p) {
                loadMixed(scratch.mixed, halves, p, std::make_index_sequence<8>{});
                // Fixed8-point monomial transform; not a dense Vandermonde
                // loop. Cantor conversion is XOR-only, and this FFT has five
                // nonzero GF256 scaling operations per four leaf frequencies.
                TransformCantor<8>::evaluateFour(scratch.mixed);
                const auto zero = _mm512_setzero_si512();
                if constexpr (Deferred) {
                    namespace deferred = char2::field128_lift::sparse_byte::deferred;
                    deferred::Product all{zero, zero, zero, zero};
                    deferred::Product bit0{zero, zero, zero, zero};
                    deferred::Product bit1{zero, zero, zero, zero};
                    deferred::Product bit2{zero, zero, zero, zero};
                    accumulateNodeDeferred<0>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<1>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<2>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<3>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<4>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<5>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<6>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNodeDeferred<7>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    // Infinity contributes only to y1, with scalar weight1.
                    bit0 = deferred::sum(bit0, deferred::multiplyRaw4(
                        _mm512_load_si512(halves[7] + 4 * p),
                        _mm512_load_si512(keys.values[p][8])));
                    _mm512_store_si512(scratch.spectra[0] + 4 * p, deferred::reduce(all));
                    if constexpr (SPIN_OR_WINDOW_RAW_SCALE) {
                        // This alternate is meaningful only for Deferred.
                        // Eight bytewise scales plus two final reductions use
                        //20 GFNI multiplies, versus26 in the retained path.
                        bit0 = deferred::sum(bit0, window_nine_detail::scaleRaw<beta1>(bit1));
                        bit0 = deferred::sum(bit0, window_nine_detail::scaleRaw<beta2>(bit2));
                        _mm512_store_si512(scratch.spectra[1] + 4 * p, deferred::reduce(bit0));
                    } else {
                        // Retained control: reduce each subset before scaling.
                        auto high = deferred::reduce(bit0);
                        high = _mm512_xor_si512(high,
                            window_nine_detail::scale<beta1>(deferred::reduce(bit1)));
                        high = _mm512_xor_si512(high,
                            window_nine_detail::scale<beta2>(deferred::reduce(bit2)));
                        _mm512_store_si512(scratch.spectra[1] + 4 * p, high);
                    }
                } else {
                    auto all = zero, bit0 = zero, bit1 = zero, bit2 = zero;
                    accumulateNode<0, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<1, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<2, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<3, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<4, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<5, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<6, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    accumulateNode<7, Gfni>(keys, scratch.mixed, p, all, bit0, bit1, bit2);
                    bit0 = _mm512_xor_si512(bit0, product<Gfni>(
                        _mm512_load_si512(halves[7] + 4 * p),
                        _mm512_load_si512(keys.values[p][8])));
                    auto high = _mm512_xor_si512(bit0, window_nine_detail::scale<beta1>(bit1));
                    high = _mm512_xor_si512(high, window_nine_detail::scale<beta2>(bit2));
                    _mm512_store_si512(scratch.spectra[0] + 4 * p, all);
                    _mm512_store_si512(scratch.spectra[1] + 4 * p, high);
                }
            }
            if constexpr (SPIN_OR_WINDOW_COORDINATES) {
                Transform256Cantor::interpolateNovel(scratch.spectra[0]);
                Transform256Cantor::interpolateNovel(scratch.spectra[1]);
            } else {
                Transform256Cantor::interpolate(scratch.spectra[0]);
                Transform256Cantor::interpolate(scratch.spectra[1]);
            }
            constexpr unsigned extraction = SPIN_OR_WINDOW_COORDINATES ? 128 : 127;
            for (unsigned j = 0; j != 128; j += 4) {
                _mm512_storeu_si512(output + j,
                    _mm512_loadu_si512(scratch.spectra[0] + extraction + j));
                _mm512_storeu_si512(output + 128 + j,
                    _mm512_loadu_si512(scratch.spectra[1] + extraction + j));
            }
        }
    };
};

static_assert(sizeof(ToeplitzWindow1024Nine::PreparedBatch) == 9 * 256 * sizeof(Element));
static_assert(alignof(ToeplitzWindow1024Nine::PreparedBatch) == 64);
static_assert(alignof(ToeplitzWindow1024Nine::Scratch) == 64);
static_assert(alignof(ToeplitzWindow1024Nine::State) == 64);

#undef SPIN_OR_NINE_INLINE
} // namespace spin::research::one_route
