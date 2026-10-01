#pragma once

// Reusable exact t64/s16 inner kernel. This header does not select a code,
// change the setup distribution, allocate a route, or certify a new length.
// Its three primary recipes reproduce the measured packet5 / fusion17 /19
// inner implementations; StreamLegacy retains fusion6 for size tuning.
#include "InnerPacketMaps.h"
#include "InnerPacketFeedback.h"
#include "InnerPacketStream.h"
#if defined(__AVX512VBMI__) || (defined(SPIN_PACKET_ENABLE_VBMI) && SPIN_PACKET_ENABLE_VBMI)
#include "InnerVbmiState.h"
#define SPIN_PACKET_KERNEL_VBMI_COMPILED 1
#else
#define SPIN_PACKET_KERNEL_VBMI_COMPILED 0
#endif
#include <array>
#include <bit>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace spin::research::packet_inner {

enum class InnerRecipe {
    Packet,        // Materialized packets; retained state conversion.
    PacketVbmi,    // Materialized packets; VBMI pack and unpack.
    StreamVbmi,    // Grouped four-packet stream; VBMI pack and unpack.
    StreamLegacy,  // Incremental four-packet stream; retained conversion.
    StreamVbmiEarly, // Compute the independent next-state product before emission.
    StreamVbmiPruned, // Omit algebraically unused intermediate XORs.
    StreamVbmiPrunedEarly
};

inline constexpr bool packetInnerVbmiCompiled = SPIN_PACKET_KERNEL_VBMI_COMPILED != 0;

class PreparedPacketUpdates16 {
public:
    // Rows encode a matrix's output coordinates as masks of its inputs.
    // The constructor takes ORIGINAL-COORDINATE REVERSE matrices: M^T,
    // already transposed by the caller. It must not transpose them again.
    using OriginalReverseMatrix = std::array<unsigned,16>;
    using PackedRow = spin::research::gfni_r4::Dense16Row;

    // Only a prepared owner can construct this typed view. It is cheap to
    // copy, owns no storage, and must not outlive or survive modification of
    // the owner. No mutable pointer escapes the preparation interface.
    class ConjugatedRows {
        friend class PreparedPacketUpdates16;
        const PackedRow* rows_;
        std::size_t count_;
        ConjugatedRows(const PackedRow* rows, std::size_t count) noexcept
            : rows_(rows), count_(count) {}
    public:
        [[nodiscard]] const PackedRow* data() const noexcept { return rows_; }
        [[nodiscard]] std::size_t size() const noexcept { return count_; }
    };

    // Preparation is outside the encode path. Mask width is always checked;
    // full-rank checking is the default but can be skipped when the caller
    // has already validated its sampled GL16 matrices. Neither setting
    // samples randomness or establishes a setup-failure certificate.
    explicit PreparedPacketUpdates16(std::span<const OriginalReverseMatrix> originalReverse,
                                     bool checkFullRank = true)
        : rows_(originalReverse.size()) {
        const auto& p = Maps<true>::basisRows;
        const auto& inverse = Maps<true>::inverseRows;
        for(std::size_t epoch = 0; epoch < originalReverse.size(); ++epoch) {
            const auto& original = originalReverse[epoch];
            for(auto row : original)
                if(row & ~0xffffU)
                    throw std::invalid_argument("packet update matrix exceeds sixteen coordinates");
            if(checkFullRank && !fullRank(original))
                throw std::invalid_argument("packet update matrix is not invertible");
            OriginalReverseMatrix transformed{};
            for(unsigned output = 0; output < 16; ++output) {
                unsigned originalMask = 0;
                for(unsigned mask = p[output]; mask; mask &= mask - 1)
                    originalMask ^= original[std::countr_zero(mask)];
                for(unsigned mask = originalMask; mask; mask &= mask - 1)
                    transformed[output] ^= inverse[std::countr_zero(mask)];
            }
            // This is P M^T P^-1, in the same basis as Maps<true>::finish
            // (P A^T) and Maps<true>::emission (A P^-1).
            rows_[epoch] = spin::research::gfni_r4::makeDense16(transformed.data());
        }
    }

    [[nodiscard]] ConjugatedRows conjugatedRows() const & noexcept {
        return ConjugatedRows(rows_.data(), rows_.size());
    }
    ConjugatedRows conjugatedRows() const && = delete;
    [[nodiscard]] std::size_t epochs() const noexcept { return rows_.size(); }

private:
    std::vector<PackedRow> rows_;

    static bool fullRank(OriginalReverseMatrix rows) noexcept {
        unsigned rank = 0;
        for(unsigned column = 0; column < 16; ++column) {
            unsigned pivot = rank;
            while(pivot < 16 && !(rows[pivot] & (1U << column))) ++pivot;
            if(pivot == 16) return false;
            std::swap(rows[rank], rows[pivot]);
            for(unsigned row = rank + 1; row < 16; ++row)
                if(rows[row] & (1U << column)) rows[row] ^= rows[rank];
            ++rank;
        }
        return true;
    }
};

namespace packet_kernel_detail {
using Block = spin::detail::kernel::block;
namespace gf = spin::research::gfni_r4;

template<InnerRecipe Recipe>
struct RecipeTraits {
    static constexpr bool early = Recipe == InnerRecipe::StreamVbmiEarly || Recipe == InnerRecipe::StreamVbmiPrunedEarly;
    static constexpr bool pruned = Recipe == InnerRecipe::StreamVbmiPruned || Recipe == InnerRecipe::StreamVbmiPrunedEarly;
    static constexpr bool vbmi = Recipe == InnerRecipe::PacketVbmi || Recipe == InnerRecipe::StreamVbmi || early || pruned;
    static constexpr bool stream = Recipe == InnerRecipe::StreamVbmi || Recipe == InnerRecipe::StreamLegacy || early || pruned;
    static constexpr bool incremental = Recipe == InnerRecipe::StreamLegacy;
    static_assert(Recipe == InnerRecipe::Packet || Recipe == InnerRecipe::PacketVbmi ||
                  Recipe == InnerRecipe::StreamVbmi || Recipe == InnerRecipe::StreamLegacy || early || pruned,
                  "unsupported packet inner recipe");
    static_assert(!vbmi || packetInnerVbmiCompiled,
                  "VBMI packet recipe requires a VBMI-enabled translation unit");
};

template<InnerRecipe Recipe>
static SPIN_FORCEINLINE void unpack(const gf::Packed<16>& state, __m128i* words) {
    if constexpr(RecipeTraits<Recipe>::vbmi) {
#if SPIN_PACKET_KERNEL_VBMI_COMPILED
        wideUnpackVbmi(state, words);
#endif
    } else gf::unpack(state, words);
}

template<InnerRecipe Recipe>
static SPIN_FORCEINLINE void pack(const __m128i* words, gf::Packed<16>& state) {
    if constexpr(RecipeTraits<Recipe>::vbmi) {
#if SPIN_PACKET_KERNEL_VBMI_COMPILED
        widePackVbmi(words, state);
#endif
    } else gf::pack<16>(words, state);
}

template<std::size_t... I>
static SPIN_FORCEINLINE void loadPackets(const Block* raw, __m512i* packets,
                                         std::index_sequence<I...>) {
    ((packets[I] = _mm512_loadu_si512(raw + 4 * I)), ...);
}

template<std::size_t H, class Emit>
static SPIN_FORCEINLINE void emitOutput(const Block* raw, __m512i* packets,
                                        std::size_t packetBase, Emit& emit) {
    packets[H] = _mm512_xor_si512(packets[H], _mm512_loadu_si512(raw + 4 * H));
    emit(packetBase + H, packets[H]);
}

template<class Emit, std::size_t... I>
static SPIN_FORCEINLINE void emitOutputs(const Block* raw, __m512i* packets,
                                         std::size_t packetBase, Emit& emit,
                                         std::index_sequence<I...>) {
    (emitOutput<15 - I>(raw, packets, packetBase, emit), ...);
}

template<InnerRecipe Recipe, class Emit>
static SPIN_FORCEINLINE void reverseUnchecked(const Block* input, std::size_t n,
                                               const gf::Dense16Row* rows, Emit& emit) {
    using Traits = RecipeTraits<Recipe>;
    // Ordinary-function SIMD scratch. Do not move these arrays into a
    // coroutine frame. No setup allocation or runtime dispatch occurs here.
    alignas(64) __m128i words[16], moments[64], syndrome[16];
    alignas(64) __m512i packets[16];
    gf::Packed<16> state{}, feedback;
    for(std::size_t epoch = n / 64; epoch-- > 0;) {
        const auto* raw = input + 64 * epoch;
        const bool first = epoch + 1 == n / 64;
        if(first) {
            loadPackets(raw, packets, std::make_index_sequence<16>{});
            for(unsigned h = 16; h-- > 0;) emit(16 * epoch + h, packets[h]);
        } else {
            unpack<Recipe>(state, words);
            if constexpr(Traits::early) if(epoch) gf::denseStep16<false>(state, rows[epoch]);
            if constexpr(Traits::stream) {
                streamStep<Traits::pruned,Traits::incremental>(words, raw, 16 * epoch, emit, moments);
            } else {
                Maps<true>::template emission<false>(words, packets);
                emitOutputs(raw, packets, 16 * epoch, emit, std::make_index_sequence<16>{});
            }
        }
        if(!epoch) break; // No flush or unused update after the final output.
        if constexpr(Traits::stream) {
            if(first) packetMoments(packets, moments);
        } else packetMoments(packets, moments);
        Maps<true>::finish(moments, syndrome);
        pack<Recipe>(syndrome, feedback);
        if(first) state = feedback; // Entering state was zero; M^T 0 is zero.
        else if constexpr(Traits::early) {
            state.v[0] = _mm512_xor_si512(state.v[0], feedback.v[0]);
            state.v[1] = _mm512_xor_si512(state.v[1], feedback.v[1]);
            state.v[2] = _mm512_xor_si512(state.v[2], feedback.v[2]);
            state.v[3] = _mm512_xor_si512(state.v[3], feedback.v[3]);
        }
        else gf::denseStep16<true>(state, rows[epoch], &feedback);
    }
}

} // namespace packet_kernel_detail

// Encode the transpose inner using original physical packet indices.
// input has n block elements; n is a multiple of64. The prepared view has
// one update per physical step, including the unused boundary entries.
// emit(packetIndex, __m512i) receives four consecutive output elements in
// one vector, in descending packet order. The emitter must not overwrite
// unread input or the prepared update table. It controls cached/NT stores;
// if it uses non-temporal stores, the caller supplies the required fence
// before consuming the routed scratch. No fence is hidden in this kernel.
// Select ISA/size policy outside the hot loop; a VBMI recipe additionally
// requires a supporting CPU, not just a VBMI-enabled build.
template<InnerRecipe Recipe, class Emit>
static SPIN_FORCEINLINE void transposeInner(const spin::detail::kernel::block* input,
    std::size_t n, PreparedPacketUpdates16::ConjugatedRows updates, Emit&& emit) {
    using Traits = packet_kernel_detail::RecipeTraits<Recipe>;
    static_assert(!Traits::vbmi || packetInnerVbmiCompiled);
    if(n % 64 || updates.size() != n / 64)
        throw std::invalid_argument("packet inner length and prepared update count disagree");
    if(!n) return;
    if(!input) throw std::invalid_argument("packet inner input is null");
    packet_kernel_detail::reverseUnchecked<Recipe>(input, n, updates.data(), emit);
}

} // namespace spin::research::packet_inner

#undef SPIN_PACKET_KERNEL_VBMI_COMPILED
