#pragma once

// Experimental exact-map composition: store persistent state as ordinary
// words, but do NOT apply its random matrix with word lookup tables.
// A packed old state is retained during emission. Swapping the GFNI operand
// roles then fuses the random matrix product with the inverse bit transpose.
#include "PacketInnerKernel.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace spin::research::packet_inner {

enum class ComposedRecipe { Packet, Stream };

struct ComposedWordRow16 {
    // Four matrix blocks. Unlike Dense16Row, bytes contain output rows in
    // ASCENDING order: row0 occupies the low byte. These qwords are the DATA
    // operands of GFNI; the changing packed state is its MATRIX operand.
    std::uint64_t data[4]{};
};
static_assert(sizeof(ComposedWordRow16) == sizeof(spin::research::gfni_r4::Dense16Row));

static constexpr std::uint64_t composedByteReverse(std::uint64_t x) noexcept {
    x = ((x & 0x00ff00ff00ff00ffULL) << 8) | ((x >> 8) & 0x00ff00ff00ff00ffULL);
    x = ((x & 0x0000ffff0000ffffULL) << 16) | ((x >> 16) & 0x0000ffff0000ffffULL);
    return (x << 32) | (x >> 32);
}

static inline ComposedWordRow16 makeComposedWordRow16(
    const spin::research::gfni_r4::Dense16Row& row) noexcept {
    return {{composedByteReverse(row.matrix[0]), composedByteReverse(row.matrix[1]),
             composedByteReverse(row.matrix[2]), composedByteReverse(row.matrix[3])}};
}

// Owns exactly32 bytes per physical epoch, like PreparedPacketUpdates16.
// Construction changes only the byte orientation of its already conjugated
// matrices. It neither samples matrices nor transposes M a second time.
class PreparedComposedUpdates16 {
public:
    class Rows {
        friend class PreparedComposedUpdates16;
        const ComposedWordRow16* data_;
        std::size_t size_;
        Rows(const ComposedWordRow16* data, std::size_t size) noexcept : data_(data), size_(size) {}
    public:
        [[nodiscard]] const ComposedWordRow16* data() const noexcept { return data_; }
        [[nodiscard]] std::size_t size() const noexcept { return size_; }
    };

    explicit PreparedComposedUpdates16(PreparedPacketUpdates16::ConjugatedRows original)
        : rows_(original.size()) {
        for(std::size_t i = 0; i < rows_.size(); ++i)
            rows_[i] = makeComposedWordRow16(original.data()[i]);
    }

    [[nodiscard]] Rows rows() const & noexcept { return Rows(rows_.data(), rows_.size()); }
    Rows rows() const && = delete;
private:
    std::vector<ComposedWordRow16> rows_;
};

#if defined(__AVX512VBMI__) || (defined(SPIN_PACKET_ENABLE_VBMI) && SPIN_PACKET_ENABLE_VBMI)
namespace composed_detail {
namespace gf = spin::research::gfni_r4;
using Block = spin::detail::kernel::block;

template<unsigned First>
static SPIN_FORCEINLINE void storeFourAdd(__m128i* out, __m512i value, const __m128i* feedback) {
    out[First] = _mm_xor_si128(_mm512_castsi512_si128(value), feedback[First]);
    out[First + 1] = _mm_xor_si128(_mm512_extracti32x4_epi32(value, 1), feedback[First + 1]);
    out[First + 2] = _mm_xor_si128(_mm512_extracti32x4_epi32(value, 2), feedback[First + 2]);
    out[First + 3] = _mm_xor_si128(_mm512_extracti32x4_epi32(value, 3), feedback[First + 3]);
}

// For one8x8 state byte matrix Q, the retained route computes Q*M^T and
// then transposes it for word output. Instead GFNI(M_rows, Q) directly
// computes that transposed result. Ascending M row bytes account for the
// instruction's reversed matrix-row convention. This identity applies to
// arbitrary binary M, not just invertible matrices.
// Eight GFNI + four byte permutes produce ordinary next-state words;
// no additional inverse-transpose GFNI instructions are required.
static SPIN_FORCEINLINE void applyPackedToWords(const gf::Packed<16>& before,
    const ComposedWordRow16& row, const __m128i* feedback, __m128i* words) {
    const auto m00 = _mm512_set1_epi64(row.data[0]);
    const auto m01 = _mm512_set1_epi64(row.data[1]);
    const auto m10 = _mm512_set1_epi64(row.data[2]);
    const auto m11 = _mm512_set1_epi64(row.data[3]);
    const auto lo0 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m00, before.v[0], 0),
        _mm512_gf2p8affine_epi64_epi8(m01, before.v[2], 0));
    const auto hi0 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m00, before.v[1], 0),
        _mm512_gf2p8affine_epi64_epi8(m01, before.v[3], 0));
    const auto lo1 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m10, before.v[0], 0),
        _mm512_gf2p8affine_epi64_epi8(m11, before.v[2], 0));
    const auto hi1 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m10, before.v[1], 0),
        _mm512_gf2p8affine_epi64_epi8(m11, before.v[3], 0));
    storeFourAdd<0>(words, _mm512_permutex2var_epi8(lo0, vbmi_state_detail::unpackIndex<0>(), hi0), feedback);
    storeFourAdd<4>(words, _mm512_permutex2var_epi8(lo0, vbmi_state_detail::unpackIndex<4>(), hi0), feedback);
    storeFourAdd<8>(words, _mm512_permutex2var_epi8(lo1, vbmi_state_detail::unpackIndex<0>(), hi1), feedback);
    storeFourAdd<12>(words, _mm512_permutex2var_epi8(lo1, vbmi_state_detail::unpackIndex<4>(), hi1), feedback);
}

template<std::size_t... I>
static SPIN_FORCEINLINE void copyWords(const __m128i* in, __m128i* out, std::index_sequence<I...>) {
    ((out[I] = in[I]), ...);
}

struct alignas(64) WordPackets16 { __m512i v[4]; };

template<unsigned Group>
static SPIN_FORCEINLINE void packWordGroup(const WordPackets16& words, gf::Packed<16>& packed) {
    const auto lo = _mm512_permutex2var_epi8(words.v[2 * Group],
        vbmi_state_detail::packIndex<0>(), words.v[2 * Group + 1]);
    const auto hi = _mm512_permutex2var_epi8(words.v[2 * Group],
        vbmi_state_detail::packIndex<8>(), words.v[2 * Group + 1]);
    const auto basis = _mm512_set1_epi64(0x0102040810204080ULL);
    packed.v[2 * Group] = _mm512_gf2p8affine_epi64_epi8(basis, lo, 0);
    packed.v[2 * Group + 1] = _mm512_gf2p8affine_epi64_epi8(basis, hi, 0);
}

static SPIN_FORCEINLINE void wordPackets(const __m128i* words, WordPackets16& packets) {
    packets.v[0] = gf::join(words[0], words[1], words[2], words[3]);
    packets.v[1] = gf::join(words[4], words[5], words[6], words[7]);
    packets.v[2] = gf::join(words[8], words[9], words[10], words[11]);
    packets.v[3] = gf::join(words[12], words[13], words[14], words[15]);
}

static SPIN_FORCEINLINE void coordinateWords(const WordPackets16& packets, __m128i* words) {
    vbmi_state_detail::storeFour<0>(words, packets.v[0]);
    vbmi_state_detail::storeFour<4>(words, packets.v[1]);
    vbmi_state_detail::storeFour<8>(words, packets.v[2]);
    vbmi_state_detail::storeFour<12>(words, packets.v[3]);
}

// Keeping ordinary-coordinate state in four ZMMs moves, rather than adds,
// the twelve boundary lane inserts: the input pack needs none; feedback
// joins need twelve. Four ZMM XORs then replace sixteen XMM feedback XORs.
static SPIN_FORCEINLINE void applyPackedToPackets(const gf::Packed<16>& before,
    const ComposedWordRow16& row, const __m128i* feedback, WordPackets16& words) {
    const auto m00 = _mm512_set1_epi64(row.data[0]);
    const auto m01 = _mm512_set1_epi64(row.data[1]);
    const auto m10 = _mm512_set1_epi64(row.data[2]);
    const auto m11 = _mm512_set1_epi64(row.data[3]);
    const auto lo0 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m00, before.v[0], 0),
        _mm512_gf2p8affine_epi64_epi8(m01, before.v[2], 0));
    const auto hi0 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m00, before.v[1], 0),
        _mm512_gf2p8affine_epi64_epi8(m01, before.v[3], 0));
    const auto lo1 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m10, before.v[0], 0),
        _mm512_gf2p8affine_epi64_epi8(m11, before.v[2], 0));
    const auto hi1 = _mm512_xor_si512(
        _mm512_gf2p8affine_epi64_epi8(m10, before.v[1], 0),
        _mm512_gf2p8affine_epi64_epi8(m11, before.v[3], 0));
    words.v[0] = _mm512_xor_si512(_mm512_permutex2var_epi8(lo0, vbmi_state_detail::unpackIndex<0>(), hi0),
        gf::join(feedback[0], feedback[1], feedback[2], feedback[3]));
    words.v[1] = _mm512_xor_si512(_mm512_permutex2var_epi8(lo0, vbmi_state_detail::unpackIndex<4>(), hi0),
        gf::join(feedback[4], feedback[5], feedback[6], feedback[7]));
    words.v[2] = _mm512_xor_si512(_mm512_permutex2var_epi8(lo1, vbmi_state_detail::unpackIndex<0>(), hi1),
        gf::join(feedback[8], feedback[9], feedback[10], feedback[11]));
    words.v[3] = _mm512_xor_si512(_mm512_permutex2var_epi8(lo1, vbmi_state_detail::unpackIndex<4>(), hi1),
        gf::join(feedback[12], feedback[13], feedback[14], feedback[15]));
}

template<ComposedRecipe Recipe, class Emit>
static SPIN_FORCEINLINE void reverse(const Block* input, std::size_t n,
    const ComposedWordRow16* rows, Emit& emit) {
    alignas(64) __m128i words[16], moments[64], syndrome[16];
    alignas(64) __m512i packets[16];
    gf::Packed<16> before;
    for(std::size_t epoch = n / 64; epoch-- > 0;) {
        const auto* raw = input + 64 * epoch;
        const bool first = epoch + 1 == n / 64;
        if(first) {
            packet_kernel_detail::loadPackets(raw, packets, std::make_index_sequence<16>{});
            for(unsigned h = 16; h-- > 0;) emit(16 * epoch + h, packets[h]);
        } else {
            // Preserve old state in only four packed ZMMs before expansion;
            // sixteen coordinate words need not remain live through feedback.
            // The final physical epoch needs no state matrix update or pack.
            if(epoch) widePackVbmi(words, before);
            if constexpr(Recipe == ComposedRecipe::Stream)
                streamStep<false,false>(words, raw, 16 * epoch, emit, moments);
            else {
                Maps<true>::template emission<false>(words, packets);
                packet_kernel_detail::emitOutputs(raw, packets, 16 * epoch, emit, std::make_index_sequence<16>{});
            }
        }
        if(!epoch) break;
        if constexpr(Recipe == ComposedRecipe::Stream) {
            if(first) packetMoments(packets, moments);
        } else packetMoments(packets, moments);
        Maps<true>::finish(moments, syndrome);
        if(first) copyWords(syndrome, words, std::make_index_sequence<16>{});
        else applyPackedToWords(before, rows[epoch], syndrome, words);
    }
}

template<ComposedRecipe Recipe, class Emit>
static SPIN_FORCEINLINE void reverseWide(const Block* input, std::size_t n,
    const ComposedWordRow16* rows, Emit& emit) {
    alignas(64) __m128i words[16], moments[64], syndrome[16];
    alignas(64) __m512i packets[16];
    gf::Packed<16> before;
    WordPackets16 state;
    for(std::size_t epoch = n / 64; epoch-- > 0;) {
        const auto* raw = input + 64 * epoch;
        const bool first = epoch + 1 == n / 64;
        if(first) {
            packet_kernel_detail::loadPackets(raw, packets, std::make_index_sequence<16>{});
            for(unsigned h = 16; h-- > 0;) emit(16 * epoch + h, packets[h]);
        } else {
            if(epoch) { packWordGroup<0>(state, before); packWordGroup<1>(state, before); }
            coordinateWords(state, words);
            if constexpr(Recipe == ComposedRecipe::Stream)
                streamStep<false,false>(words, raw, 16 * epoch, emit, moments);
            else {
                Maps<true>::template emission<false>(words, packets);
                packet_kernel_detail::emitOutputs(raw, packets, 16 * epoch, emit, std::make_index_sequence<16>{});
            }
        }
        if(!epoch) break;
        if constexpr(Recipe == ComposedRecipe::Stream) {
            if(first) packetMoments(packets, moments);
        } else packetMoments(packets, moments);
        Maps<true>::finish(moments, syndrome);
        if(first) wordPackets(syndrome, state);
        else applyPackedToPackets(before, rows[epoch], syndrome, state);
    }
}

} // namespace composed_detail
#endif

// Same physical inputs, output indices, transformed basis, and typed-view
// lifecycle as transposeInner. The only setup change is byte-reversing each
// update block once. Routing-store policy and any required NT fence remain
// the caller's responsibility. Requires VBMI-capable compilation and CPU.
template<ComposedRecipe Recipe, class Emit>
static SPIN_FORCEINLINE void transposeInnerComposed(const spin::detail::kernel::block* input,
    std::size_t n, PreparedComposedUpdates16::Rows updates, Emit&& emit) {
    static_assert(Recipe == ComposedRecipe::Packet || Recipe == ComposedRecipe::Stream);
    static_assert(packetInnerVbmiCompiled || Recipe != Recipe,
                  "composed packet inner requires a VBMI-enabled translation unit");
    if(n % 64 || updates.size() != n / 64)
        throw std::invalid_argument("composed packet inner length and update count disagree");
    if(!n) return;
    if(!input) throw std::invalid_argument("composed packet inner input is null");
#if defined(__AVX512VBMI__) || (defined(SPIN_PACKET_ENABLE_VBMI) && SPIN_PACKET_ENABLE_VBMI)
    composed_detail::reverse<Recipe>(input, n, updates.data(), emit);
#endif
}

template<ComposedRecipe Recipe, class Emit>
static SPIN_FORCEINLINE void transposeInnerComposedWide(const spin::detail::kernel::block* input,
    std::size_t n, PreparedComposedUpdates16::Rows updates, Emit&& emit) {
    static_assert(Recipe == ComposedRecipe::Packet || Recipe == ComposedRecipe::Stream);
    static_assert(packetInnerVbmiCompiled || Recipe != Recipe,
                  "composed packet inner requires a VBMI-enabled translation unit");
    if(n % 64 || updates.size() != n / 64)
        throw std::invalid_argument("composed packet inner length and update count disagree");
    if(!n) return;
    if(!input) throw std::invalid_argument("composed packet inner input is null");
#if defined(__AVX512VBMI__) || (defined(SPIN_PACKET_ENABLE_VBMI) && SPIN_PACKET_ENABLE_VBMI)
    composed_detail::reverseWide<Recipe>(input, n, updates.data(), emit);
#endif
}

// Setup-only independent scalar check. All2,048 state input bits are tested
// under eight explicit binary matrices, including nonsymmetric and singular
// matrices. A separate2,048-bit feedback basis checks the added word map.
// The oracle uses ordinary uint64 XOR, never a packed update or conversion.
[[nodiscard]] inline bool composedStateSelfCheck() {
#if defined(__AVX512VBMI__) || (defined(SPIN_PACKET_ENABLE_VBMI) && SPIN_PACKET_ENABLE_VBMI)
    namespace gf = spin::research::gfni_r4;
    using Matrix = std::array<unsigned,16>;
    using ScalarWords = std::array<std::array<std::uint64_t,2>,16>;
    alignas(64) __m128i words[16], feedback[16], actual[16];
    gf::Packed<16> packed, packedWide;
    composed_detail::WordPackets16 wordPacketState, actualPacketState;
    std::array<Matrix,8> matrices{};
    for(unsigned j = 0; j < 16; ++j) {
        matrices[0][j] = 1U << j;
        matrices[1][j] = 1U << ((j + 5) & 15);
        matrices[2][j] = (1U << j) ^ (j ? (1U << (j - 1)) : 0U);
        matrices[3][j] = 0xffffU;
        matrices[4][j] = (0xa6d3U << (j & 3)) & 0xffffU;
        matrices[5][j] = 1U << (15 - j);
        matrices[6][j] = (0x8b29U ^ (j * 0x2479U)) & 0xffffU;
        // Matrix7 remains zero.
    }
    for(const auto& matrix : matrices) {
        const auto row = makeComposedWordRow16(gf::makeDense16(matrix.data()));
        for(unsigned bit = 0; bit < 2048; ++bit) {
            ScalarWords scalar{}, expected{};
            scalar[bit / 128][(bit % 128) / 64] = std::uint64_t(1) << (bit % 64);
            std::memcpy(words, scalar.data(), sizeof(words));
            std::memset(feedback, 0, sizeof(feedback));
            for(unsigned j = 0; j < 16; ++j)
                for(unsigned k = 0; k < 16; ++k)
                    if((matrix[j] >> k) & 1U) {
                        expected[j][0] ^= scalar[k][0]; expected[j][1] ^= scalar[k][1];
                    }
            widePackVbmi(words, packed);
            composed_detail::wordPackets(words, wordPacketState);
            composed_detail::packWordGroup<0>(wordPacketState, packedWide);
            composed_detail::packWordGroup<1>(wordPacketState, packedWide);
            if(std::memcmp(&packed, &packedWide, sizeof(packed))) return false;
            composed_detail::applyPackedToWords(packed, row, feedback, actual);
            if(std::memcmp(actual, expected.data(), sizeof(actual))) return false;
            composed_detail::applyPackedToPackets(packedWide, row, feedback, actualPacketState);
            composed_detail::coordinateWords(actualPacketState, actual);
            if(std::memcmp(actual, expected.data(), sizeof(actual))) return false;
        }
    }
    const auto row = makeComposedWordRow16(gf::makeDense16(matrices[6].data()));
    std::memset(words, 0, sizeof(words));
    widePackVbmi(words, packed);
    for(unsigned bit = 0; bit < 2048; ++bit) {
        ScalarWords expected{};
        expected[bit / 128][(bit % 128) / 64] = std::uint64_t(1) << (bit % 64);
        std::memcpy(feedback, expected.data(), sizeof(feedback));
        composed_detail::applyPackedToWords(packed, row, feedback, actual);
        if(std::memcmp(actual, expected.data(), sizeof(actual))) return false;
        composed_detail::applyPackedToPackets(packed, row, feedback, actualPacketState);
        composed_detail::coordinateWords(actualPacketState, actual);
        if(std::memcmp(actual, expected.data(), sizeof(actual))) return false;
    }
    return true;
#else
    return false;
#endif
}

} // namespace spin::research::packet_inner
