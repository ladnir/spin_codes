#pragma once

#include "Transform256.h"
#include "WideDiagonal32.h"
#include <cstddef>
#include <utility>

namespace spin::research::one_route {

#if defined(_MSC_VER)
#define SPIN_OR_PERSISTENT_INLINE __forceinline
#else
#define SPIN_OR_PERSISTENT_INLINE inline __attribute__((always_inline))
#endif

// Research kernel for PERSISTENT_MDS.md. Element is only a 16-byte carrier
// here: each record contains four coordinates of WideDiagonal32's GF32
// tower. This is not multiplication in SparseByteField128's field basis.
// Each prepared scalar is shared across the four payload lanes.
//
// Setup supplies 2*Batch independent uniform GF32 feedback scalars (zero
// included), and Batch fixed scalars 1,alpha,...,alpha^(Batch-1). This
// header neither samples setup nor selects or validates primitive alpha.
template<unsigned Batch>
struct Persistent32 {
    static_assert(Batch == 16 || Batch == 32 || Batch == 64 || Batch == 128);
    static_assert(sizeof(Element) == 16);
    using Coefficient = char2::wide32::Coefficient;
    using Transform = TransformCantor<Batch>;
    static constexpr unsigned batchSize = Batch;
    static constexpr unsigned levels = Transform::levels;
    static constexpr unsigned keyCount = 2 * Batch;

private:
    using Blocks = std::make_index_sequence<Batch / 4>;

    template<std::size_t I>
    static SPIN_OR_PERSISTENT_INLINE void copyBlock(
        const Element* input, Element* output) noexcept {
        _mm512_storeu_si512(output + 4 * I, _mm512_loadu_si512(input + 4 * I));
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void copy(
        const Element* input, Element* output, std::index_sequence<I...>) noexcept {
        (copyBlock<I>(input, output), ...);
    }

    template<std::size_t I>
    static SPIN_OR_PERSISTENT_INLINE void duplicateBlock(
        const Element* state, Element* low, Element* high) noexcept {
        const auto value = _mm512_loadu_si512(state + 4 * I);
        _mm512_storeu_si512(low + 4 * I, value);
        _mm512_storeu_si512(high + 4 * I, value);
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void duplicate(
        const Element* state, Element* low, Element* high,
        std::index_sequence<I...>) noexcept {
        (duplicateBlock<I>(state, low, high), ...);
    }

    template<std::size_t I>
    static SPIN_OR_PERSISTENT_INLINE void multiplyBlock(
        Element* low, Element* high, const Coefficient* keys) noexcept {
        const auto a = char2::wide32::multiply512(
            _mm512_loadu_si512(low + 4 * I), keys + 4 * I);
        const auto b = char2::wide32::multiply512(
            _mm512_loadu_si512(high + 4 * I), keys + Batch + 4 * I);
        _mm512_storeu_si512(low + 4 * I, a);
        _mm512_storeu_si512(high + 4 * I, b);
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void multiply(
        Element* low, Element* high, const Coefficient* keys,
        std::index_sequence<I...>) noexcept {
        (multiplyBlock<I>(low, high, keys), ...);
    }

    template<std::size_t I>
    static SPIN_OR_PERSISTENT_INLINE void combineBlock(
        const Element* low, const Element* high, Element* output) noexcept {
        _mm512_storeu_si512(output + 4 * I,
            _mm512_xor_si512(_mm512_loadu_si512(low + 4 * I),
                             _mm512_loadu_si512(high + 4 * I)));
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void combine(
        const Element* low, const Element* high, Element* output,
        std::index_sequence<I...>) noexcept {
        (combineBlock<I>(low, high, output), ...);
    }

    template<std::size_t I>
    static SPIN_OR_PERSISTENT_INLINE void reverseBlockPair(Element* values) noexcept {
        constexpr unsigned other = Batch / 4 - 1 - I;
        const auto low = _mm512_loadu_si512(values + 4 * I);
        const auto high = _mm512_loadu_si512(values + 4 * other);
        _mm512_storeu_si512(values + 4 * I, _mm512_shuffle_i32x4(high, high, 0x1b));
        _mm512_storeu_si512(values + 4 * other, _mm512_shuffle_i32x4(low, low, 0x1b));
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void reverse(
        Element* values, std::index_sequence<I...>) noexcept {
        (reverseBlockPair<I>(values), ...);
    }

    template<std::size_t Block, std::size_t J>
    static constexpr bool wholeMomentTerm = J >= 2 && (J & ~levels) == 0
        && 4 * Block >= Batch - (1U << J);

    template<std::size_t Block, std::size_t J>
    static SPIN_OR_PERSISTENT_INLINE __m512i addWholeMomentTerm(
        __m512i value, const Element* values) noexcept {
        if constexpr (wholeMomentTerm<Block, J>) {
            constexpr unsigned source = 4 * Block - Batch + (1U << J);
            return _mm512_xor_si512(value, _mm512_loadu_si512(values + source));
        } else {
            return value;
        }
    }

    template<std::size_t Block, std::size_t... J>
    static SPIN_OR_PERSISTENT_INLINE void correctWholeMomentBlock(
        Element* values, std::index_sequence<J...>) noexcept {
        if constexpr ((wholeMomentTerm<Block, J> || ...)) {
            auto value = _mm512_loadu_si512(values + 4 * Block);
            ((value = addWholeMomentTerm<Block, J>(value, values)), ...);
            _mm512_storeu_si512(values + 4 * Block, value);
        }
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void correctMoments(
        Element* values, std::index_sequence<I...>) noexcept {
        // On the nonzero Cantor coset indexed by Batch, the vanishing
        // polynomial is Delta^levels(X)+1, with derivative one. After
        // reversal, coefficient c_i equals M_i plus M_(i-Batch+2^j)
        // for every proper j subset of levels with i>=Batch-2^j.
        // All correction sources lie in the unchanged first half.
        (correctWholeMomentBlock<I>(values, std::make_index_sequence<levels>{}), ...);
        if constexpr ((1U & ~levels) == 0) {
            auto* destination = reinterpret_cast<__m256i*>(values + Batch - 2);
            const auto* source = reinterpret_cast<const __m256i*>(values);
            _mm256_storeu_si256(destination,
                _mm256_xor_si256(_mm256_loadu_si256(destination), _mm256_loadu_si256(source)));
        }
        // j=0 occurs in every Delta power and contributes only to M_(Batch-1).
        auto* destination = reinterpret_cast<__m128i*>(values + Batch - 1);
        const auto* source = reinterpret_cast<const __m128i*>(values);
        _mm_storeu_si128(destination,
            _mm_xor_si128(_mm_loadu_si128(destination), _mm_loadu_si128(source)));
    }

    static SPIN_OR_PERSISTENT_INLINE void momentsInPlace(Element* values) noexcept {
        Transform::template interpolateNovelCoset<Batch>(values);
        transform_cantor_detail::convert<levels, true, false>(values);
        reverse(values, std::make_index_sequence<Batch / 8>{});
        correctMoments(values, Blocks{});
    }

    template<std::size_t I>
    static SPIN_OR_PERSISTENT_INLINE void updateBlock(
        Element* state, const Element* moments, const Coefficient* fixedP) noexcept {
        const auto shifted = char2::wide32::multiply512(
            _mm512_loadu_si512(state + 4 * I), fixedP + 4 * I);
        _mm512_storeu_si512(state + 4 * I,
            _mm512_xor_si512(shifted, _mm512_loadu_si512(moments + 4 * I)));
    }

    template<std::size_t... I>
    static SPIN_OR_PERSISTENT_INLINE void updateBlocks(
        Element* state, const Element* moments, const Coefficient* fixedP,
        std::index_sequence<I...>) noexcept {
        (updateBlock<I>(state, moments, fixedP), ...);
    }

public:
    // Compute T*state only; the caller adds X to obtain Y. Keys are ordered
    // by Cantor point index 0,...,2*Batch-1. Exact state/result alias is safe.
    // Coefficient storage must not overlap writable Element storage.
    static SPIN_OR_PERSISTENT_INLINE void feedback(
        const Element* state, const Coefficient* keys, Element* result) noexcept {
        alignas(64) Element low[Batch], high[Batch];
        duplicate(state, low, high, Blocks{});
        Transform::template evaluateNovelCoset<0>(low);
        Transform::template evaluateNovelCoset<Batch>(high);
        multiply(low, high, keys, Blocks{});
        Transform::template interpolateNovelCoset<0>(low);
        Transform::template interpolateNovelCoset<Batch>(high);
        // The top-level inverse butterfly has zero offset and normalization
        // one. Its highest Batch novel coefficients are low XOR high.
        combine(low, high, result, Blocks{});
    }

    // result[i] = sum_j beta_j^i * output[j], beta_j=point(Batch+j).
    // Input and result may be disjoint or exactly aliased, not partly overlapped.
    static SPIN_OR_PERSISTENT_INLINE void moments(
        const Element* output, Element* result) noexcept {
        copy(output, result, Blocks{});
        momentsInPlace(result);
    }

    // Replace state by P*state+C*output. Output remains unchanged unless it
    // aliases state. No heap allocation; the only scratch is Batch records.
    static SPIN_OR_PERSISTENT_INLINE void update(
        Element* state, const Element* output, const Coefficient* fixedP) noexcept {
        alignas(64) Element projected[Batch];
        copy(output, projected, Blocks{});
        momentsInPlace(projected);
        updateBlocks(state, projected, fixedP, Blocks{});
    }
};

#undef SPIN_OR_PERSISTENT_INLINE

} // namespace spin::research::one_route
