#pragma once

#include "Transform256.h"
#include "WideDiagonal32.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

#ifndef SPIN_OR_FULL_CANTOR
#define SPIN_OR_FULL_CANTOR 0
#endif
static_assert(SPIN_OR_FULL_CANTOR == 0 || SPIN_OR_FULL_CANTOR == 1);

namespace spin::research::one_route {

#if defined(_MSC_VER)
#define SPIN_OR_ADDITIVE_INLINE __forceinline
#else
#define SPIN_OR_ADDITIVE_INLINE inline __attribute__((always_inline))
#endif

// Research-only translation between additive evaluation cosets. Each Element
// carries four GF32 tower coordinates; no SparseByteField128 products occur.
// D_delta=FFT_Batch(delta)*IFFT_Batch(0), in the fixed Cantor novel basis.
// The input and output points are point(j) and point(j)+delta, respectively.
template<unsigned Batch>
struct AdditivePersistent32 {
    static_assert(Batch == 32 || Batch == 64 || Batch == 128);
    static_assert(sizeof(Element) == 16);
    using Coefficient = char2::wide32::Coefficient;
    using Transform = TransformCantor<Batch>;
    static constexpr unsigned batchSize = Batch;
    static constexpr unsigned levels = Transform::levels;
    static constexpr unsigned preparedCount = 3 * Batch / 4 - 1;

    struct alignas(64) PackedTwiddle {
        std::array<Coefficient, 4> lanes;
    };
    static_assert(sizeof(PackedTwiddle) == 64);

    struct Shift {
        unsigned fieldBit = 0;
        std::uint32_t delta = 0;
        // Nodes of levels>=2 have heap order. The final Batch/4 entries
        // combine the two length-two butterflies within each SIMD vector.
        std::array<PackedTwiddle, preparedCount> twiddles;
    };

private:
    using Blocks = std::make_index_sequence<Batch / 4>;

#if SPIN_OR_FULL_CANTOR
    static const std::array<std::uint32_t, 32>& fullCantorBasis() {
        // Setup only. Delta has kernel {0,1}; fixing input bit0 to zero
        // leaves 31 independent columns. Gaussian elimination therefore
        // selects one deterministic root of Delta(x)=target, when it exists.
        static const auto result = [] {
            std::array<std::uint64_t, 32> equation{};
            for (unsigned column = 1; column != 32; ++column) {
                const auto value = 1U << column;
                const auto image = char2::medium32::towerMultiplyScalar(value, value) ^ value;
                for (unsigned row = 0; row != 32; ++row)
                    if ((image >> row) & 1U) equation[row] |= std::uint64_t{1} << column;
            }
            const auto root = [&equation](std::uint32_t target) {
                auto rows = equation;
                for (unsigned row = 0; row != 32; ++row)
                    if ((target >> row) & 1U) rows[row] |= std::uint64_t{1} << 32;
                std::array<unsigned, 32> pivotRow{};
                unsigned used = 0;
                for (unsigned column = 1; column != 32; ++column) {
                    unsigned pivot = used;
                    while (pivot != 32 && ((rows[pivot] >> column) & 1U) == 0) ++pivot;
                    if (pivot == 32)
                        throw std::logic_error("GF32 Delta must have rank31 with bit0 fixed");
                    std::swap(rows[used], rows[pivot]);
                    for (unsigned row = 0; row != 32; ++row)
                        if (row != used && ((rows[row] >> column) & 1U)) rows[row] ^= rows[used];
                    pivotRow[column] = used++;
                }
                for (const auto row : rows)
                    if ((row & 0xffffffffULL) == 0 && (row >> 32))
                        throw std::logic_error("GF32 Cantor chain has no next root");
                std::uint32_t value = 0;
                for (unsigned column = 1; column != 32; ++column)
                    if ((rows[pivotRow[column]] >> 32) & 1U) value |= 1U << column;
                return value;
            };

            std::array<std::uint32_t, 32> beta{};
            for (unsigned i = 0; i != 8; ++i) beta[i] = transform_cantor_detail::basis.beta[i];
            for (unsigned i = 8; i != 32; ++i) beta[i] = root(beta[i - 1]);
            if (beta[0] != 1) throw std::logic_error("GF32 Cantor chain must start at one");
            for (unsigned i = 1; i != 32; ++i)
                if ((char2::medium32::towerMultiplyScalar(beta[i], beta[i]) ^ beta[i]) != beta[i - 1])
                    throw std::logic_error("GF32 Cantor chain relation failed");

            // Verify binary independence separately from the chain identities.
            std::array<std::uint32_t, 32> pivots{};
            for (auto value : beta) {
                bool inserted = false;
                for (unsigned bit = 32; bit-- != 0;) {
                    if (((value >> bit) & 1U) == 0) continue;
                    if (pivots[bit]) value ^= pivots[bit];
                    else {
                        pivots[bit] = value;
                        inserted = true;
                        break;
                    }
                }
                if (!inserted) throw std::logic_error("GF32 Cantor directions are dependent");
            }
            return beta;
        }();
        return result;
    }
#endif

    static constexpr unsigned node(unsigned level, unsigned offset) noexcept {
        return (1U << (levels - level)) - 1 + (offset >> level);
    }

    static constexpr unsigned leafNode(unsigned offset) noexcept {
        return Batch / 2 - 1 + offset / 4;
    }

    static void prepareFour(PackedTwiddle& destination,
                            std::uint32_t a, std::uint32_t b,
                            std::uint32_t c, std::uint32_t d) noexcept {
        destination.lanes[0] = char2::wide32::prepareCoefficientDirect(a);
        destination.lanes[1] = char2::wide32::prepareCoefficientDirect(b);
        destination.lanes[2] = char2::wide32::prepareCoefficientDirect(c);
        destination.lanes[3] = char2::wide32::prepareCoefficientDirect(d);
    }

    template<std::size_t I>
    static SPIN_OR_ADDITIVE_INLINE void copyBlock(
        const Element* input, Element* output) noexcept {
        _mm512_storeu_si512(output + 4 * I, _mm512_loadu_si512(input + 4 * I));
    }

    template<std::size_t... I>
    static SPIN_OR_ADDITIVE_INLINE void copy(
        const Element* input, Element* output, std::index_sequence<I...>) noexcept {
        (copyBlock<I>(input, output), ...);
    }

    template<unsigned Level, unsigned ExtensionDepth>
    static SPIN_OR_ADDITIVE_INLINE __m512i multiplyTwiddle(
        __m512i value, const PackedTwiddle& twiddle) noexcept {
        if constexpr (Level > ExtensionDepth) {
            // Full Cantor beta_k enters GF256 after k-7 Delta applications.
            // The raw-basis control has the same property at k8 only.
            // Byte twiddles have only their zeroth cyclic diagonal nonzero.
            const auto packed = _mm512_loadu_si512(twiddle.lanes.data());
            const auto coefficient = char2::wide32::expandDiagonal<0>(packed);
            return _mm512_gf2p8mul_epi8(value, coefficient);
        } else {
            return char2::wide32::multiply512(value, twiddle.lanes.data());
        }
    }

    template<unsigned Level, unsigned Offset, unsigned ExtensionDepth, std::size_t I>
    static SPIN_OR_ADDITIVE_INLINE void stageBlock(
        const Shift& shift, Element* values) noexcept {
        constexpr unsigned half = 1U << (Level - 1);
        auto low = _mm512_loadu_si512(values + 4 * I);
        auto high = _mm512_loadu_si512(values + half + 4 * I);
        low = _mm512_xor_si512(low, multiplyTwiddle<Level, ExtensionDepth>(
            high, shift.twiddles[node(Level, Offset)]));
        high = _mm512_xor_si512(high, low);
        _mm512_storeu_si512(values + 4 * I, low);
        _mm512_storeu_si512(values + half + 4 * I, high);
    }

    template<unsigned Level, unsigned Offset, unsigned ExtensionDepth, std::size_t... I>
    static SPIN_OR_ADDITIVE_INLINE void stage(
        const Shift& shift, Element* values, std::index_sequence<I...>) noexcept {
        (stageBlock<Level, Offset, ExtensionDepth, I>(shift, values), ...);
    }

    template<unsigned Offset, unsigned ExtensionDepth>
    static SPIN_OR_ADDITIVE_INLINE void leaf4(const Shift& shift, Element* values) noexcept {
        auto value = _mm512_loadu_si512(values);
        auto low = _mm512_shuffle_i32x4(value, value, 0x44);
        auto high = _mm512_shuffle_i32x4(value, value, 0xee);
        value = _mm512_xor_si512(low, multiplyTwiddle<2, ExtensionDepth>(
            high, shift.twiddles[node(2, Offset)]));
        low = _mm512_shuffle_i32x4(value, value, 0xa0);
        high = _mm512_shuffle_i32x4(value, value, 0xf5);
        value = _mm512_xor_si512(low, multiplyTwiddle<1, ExtensionDepth>(
            high, shift.twiddles[leafNode(Offset)]));
        _mm512_storeu_si512(values, value);
    }

    template<unsigned Level, unsigned Offset, unsigned ExtensionDepth>
    static SPIN_OR_ADDITIVE_INLINE void forward(const Shift& shift, Element* values) noexcept {
        if constexpr (Level == 2) {
            leaf4<Offset, ExtensionDepth>(shift, values);
        } else {
            constexpr unsigned half = 1U << (Level - 1);
            stage<Level, Offset, ExtensionDepth>(
                shift, values, std::make_index_sequence<half / 4>{});
            forward<Level - 1, Offset, ExtensionDepth>(shift, values);
            forward<Level - 1, Offset ^ half, ExtensionDepth>(shift, values + half);
        }
    }

public:
    // Both variants preserve the existing first eight GF256 directions.
    // The optional full chain changes the later deterministic cosets, not
    // the independent random feedback diagonals or the route distribution.
    static std::uint32_t shiftDelta(unsigned fieldBit) {
        if (fieldBit >= 32) throw std::invalid_argument("additive shift bit must be below32");
#if SPIN_OR_FULL_CANTOR
        return fullCantorBasis()[fieldBit];
#else
        return fieldBit < 8 ? Transform::evaluationPoint(1U << fieldBit) : (1U << fieldBit);
#endif
    }

    // This setup routine runs once per used direction, never per batch.
    static Shift prepare(unsigned fieldBit) {
        if (fieldBit < levels || fieldBit >= 32)
            throw std::invalid_argument("additive shift bit must lie in [log2(batch),32)");
        Shift result;
        result.fieldBit = fieldBit;
        result.delta = shiftDelta(fieldBit);
        if (fieldBit < 8) return result;

        std::array<std::uint32_t, levels> deltaIterate;
        deltaIterate[0] = result.delta;
        for (unsigned i = 1; i != levels; ++i) {
            const auto previous = deltaIterate[i - 1];
            deltaIterate[i] = char2::medium32::towerMultiplyScalar(previous, previous) ^ previous;
        }
        for (unsigned level = 2; level <= levels; ++level) {
            for (unsigned offset = 0; offset != Batch; offset += 1U << level) {
                const auto coefficient = deltaIterate[level - 1]
                    ^ transform_cantor_detail::twiddle(level - 1, offset);
                auto& destination = result.twiddles[node(level, offset)];
                if (level == 2)
                    prepareFour(destination, coefficient, coefficient, coefficient ^ 1, coefficient ^ 1);
                else
                    prepareFour(destination, coefficient, coefficient, coefficient, coefficient);
            }
        }
        for (unsigned offset = 0; offset != Batch; offset += 4) {
            const auto a = result.delta ^ Transform::evaluationPoint(offset);
            const auto b = result.delta ^ Transform::evaluationPoint(offset ^ 2);
            prepareFour(result.twiddles[leafNode(offset)], a, a ^ 1, b, b ^ 1);
        }
        return result;
    }

    // Shift must come from prepare(). Exact input/output alias is supported;
    // otherwise ranges must be disjoint. No allocations, scratch buffers,
    // coefficient expansion, or field-oracle calls occur in this hot path.
    static SPIN_OR_ADDITIVE_INLINE void apply(
        const Shift& shift, const Element* input, Element* output) noexcept {
        copy(input, output, Blocks{});
        Transform::template interpolateNovelCoset<0>(output);
        switch (shift.fieldBit) {
        case 5:
            if constexpr (Batch == 32) {
                Transform::template evaluateNovelCoset<32>(output);
                return;
            }
            break; // Unreachable for a Shift returned by prepare() at larger batches.
        case 6:
            if constexpr (Batch <= 64) {
                Transform::template evaluateNovelCoset<64>(output);
                return;
            }
            break; // Unreachable for a Shift returned by prepare() at Batch128.
        case 7:
            Transform::template evaluateNovelCoset<128>(output);
            return;
        case 8:
            forward<levels, 0, 1>(shift, output);
            return;
#if SPIN_OR_FULL_CANTOR
        case 9:
            forward<levels, 0, 2>(shift, output);
            return;
        case 10:
            forward<levels, 0, 3>(shift, output);
            return;
        case 11:
            forward<levels, 0, 4>(shift, output);
            return;
        case 12:
            if constexpr (levels > 5) {
                forward<levels, 0, 5>(shift, output);
                return;
            }
            break;
        case 13:
            if constexpr (levels > 6) {
                forward<levels, 0, 6>(shift, output);
                return;
            }
            break;
#endif
        default:
            break;
        }
        // All higher directions require full GF32 products at every level;
        // they share one instantiation rather than duplicating code per bit.
        forward<levels, 0, levels>(shift, output);
    }
};

static_assert((char2::medium32::towerMultiplyScalar(1U << 8, 1U << 8) ^ (1U << 8))
              == char2::medium32::delta,
              "fieldBit8 fast path requires the declared first tower generator");

#undef SPIN_OR_ADDITIVE_INLINE

} // namespace spin::research::one_route
