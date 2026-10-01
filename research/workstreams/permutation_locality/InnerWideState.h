#pragma once

#include "FusedR4Gfni.h"
#include <cstdint>
#include <cstring>

namespace spin::research::packet_inner {
namespace wide_state_detail {

// Use the broadly available epi32 constructor for two adjacent word
// indices; no _mm512_setr_epi16 convenience intrinsic is required.
static constexpr int indexPair(unsigned lo, unsigned hi) {
    return static_cast<int>(lo | (hi << 16));
}

// An eight-coordinate group contains 8 rows of 16 bytes. Gather one
// adjacent payload-byte pair from all eight rows into each 128-bit lane,
// then separate its even/odd bytes into the two GFNI matrices. Coordinates
// are reversed here to match gfni_r4::pack's matrix-row convention.
template<unsigned Group>
static SPIN_FORCEINLINE void packGroup(const __m128i* state,
    spin::research::gfni_r4::Packed<16>& out) {
    static_assert(Group < 2);
    const auto a = spin::research::gfni_r4::join(
        state[8 * Group], state[8 * Group + 1],
        state[8 * Group + 2], state[8 * Group + 3]);
    const auto b = spin::research::gfni_r4::join(
        state[8 * Group + 4], state[8 * Group + 5],
        state[8 * Group + 6], state[8 * Group + 7]);
    const auto loIndex = _mm512_setr_epi32(
        indexPair(56,48),indexPair(40,32),indexPair(24,16),indexPair(8,0),
        indexPair(57,49),indexPair(41,33),indexPair(25,17),indexPair(9,1),
        indexPair(58,50),indexPair(42,34),indexPair(26,18),indexPair(10,2),
        indexPair(59,51),indexPair(43,35),indexPair(27,19),indexPair(11,3));
    const auto hiIndex = _mm512_setr_epi32(
        indexPair(60,52),indexPair(44,36),indexPair(28,20),indexPair(12,4),
        indexPair(61,53),indexPair(45,37),indexPair(29,21),indexPair(13,5),
        indexPair(62,54),indexPair(46,38),indexPair(30,22),indexPair(14,6),
        indexPair(63,55),indexPair(47,39),indexPair(31,23),indexPair(15,7));
    const auto separate = _mm512_broadcast_i32x4(_mm_setr_epi8(
        0,2,4,6,8,10,12,14, 1,3,5,7,9,11,13,15));
    const auto lo = _mm512_shuffle_epi8(
        _mm512_permutex2var_epi16(a, loIndex, b), separate);
    const auto hi = _mm512_shuffle_epi8(
        _mm512_permutex2var_epi16(a, hiIndex, b), separate);
    const auto basis = _mm512_set1_epi64(0x0102040810204080ULL);
    out.v[2 * Group] = _mm512_gf2p8affine_epi64_epi8(basis, lo, 0);
    out.v[2 * Group + 1] = _mm512_gf2p8affine_epi64_epi8(basis, hi, 0);
}

// Keep the API's coordinate words scalarized at the caller boundary. A
// wide store followed by sixteen narrow reads can otherwise create a
// stack bridge even though all words were already register resident.
template<unsigned First>
static SPIN_FORCEINLINE void storeFour(__m128i* out, __m512i x) {
    out[First] = _mm512_castsi512_si128(x);
    out[First + 1] = _mm512_extracti32x4_epi32(x, 1);
    out[First + 2] = _mm512_extracti32x4_epi32(x, 2);
    out[First + 3] = _mm512_extracti32x4_epi32(x, 3);
}

template<unsigned Group>
static SPIN_FORCEINLINE void unpackGroup(
    const spin::research::gfni_r4::Packed<16>& state, __m128i* out) {
    static_assert(Group < 2);
    const auto basis = _mm512_set1_epi64(0x8040201008040201ULL);
    // The retained inverse GFNI produces coordinates in ascending order,
    // unlike the reversed matrix rows consumed by packGroup.
    const auto lo = _mm512_gf2p8affine_epi64_epi8(basis, state.v[2 * Group], 0);
    const auto hi = _mm512_gf2p8affine_epi64_epi8(basis, state.v[2 * Group + 1], 0);
    const auto interleave = _mm512_broadcast_i32x4(_mm_setr_epi8(
        0,8,1,9,2,10,3,11, 4,12,5,13,6,14,7,15));
    const auto a = _mm512_shuffle_epi8(lo, interleave);
    const auto b = _mm512_shuffle_epi8(hi, interleave);
    const auto firstIndex = _mm512_setr_epi32(
        indexPair(0,8),indexPair(16,24),indexPair(32,40),indexPair(48,56),
        indexPair(1,9),indexPair(17,25),indexPair(33,41),indexPair(49,57),
        indexPair(2,10),indexPair(18,26),indexPair(34,42),indexPair(50,58),
        indexPair(3,11),indexPair(19,27),indexPair(35,43),indexPair(51,59));
    const auto secondIndex = _mm512_setr_epi32(
        indexPair(4,12),indexPair(20,28),indexPair(36,44),indexPair(52,60),
        indexPair(5,13),indexPair(21,29),indexPair(37,45),indexPair(53,61),
        indexPair(6,14),indexPair(22,30),indexPair(38,46),indexPair(54,62),
        indexPair(7,15),indexPair(23,31),indexPair(39,47),indexPair(55,63));
    storeFour<8 * Group>(out, _mm512_permutex2var_epi16(a, firstIndex, b));
    storeFour<8 * Group + 4>(out, _mm512_permutex2var_epi16(a, secondIndex, b));
}

} // namespace wide_state_detail

// Exact byte-compatible replacements for gfni_r4::pack<16>/unpack<16>.
// Requires AVX-512F/BW and GFNI, but not VBMI. Source and destination must
// be disjoint, as in the retained helpers. Ordinary __m128i alignment is
// sufficient for the word array; Packed<16> retains its 64-byte alignment.
// Per direction: four word permutes, four byte shuffles, four GFNI, and
// twelve boundary lane inserts (pack) or extracts (unpack), before compiler
// simplification. No runtime loops, byteTranspose calls, or scratch arrays.
static SPIN_FORCEINLINE void widePack(const __m128i* state,
    spin::research::gfni_r4::Packed<16>& out) {
    wide_state_detail::packGroup<0>(state, out);
    wide_state_detail::packGroup<1>(state, out);
}

static SPIN_FORCEINLINE void wideUnpack(
    const spin::research::gfni_r4::Packed<16>& state, __m128i* out) {
    wide_state_detail::unpackGroup<0>(state, out);
    wide_state_detail::unpackGroup<1>(state, out);
}

// Setup-only exhaustive validation, active under NDEBUG. Compare pack on
// every coordinate-word input bit against the retained pack, and compare
// unpack on every independently chosen packed input bit against retained
// unpack. These are two separate 2048-bit bases, not only a round trip.
// Invoke from an ordinary function, not a coroutine with SIMD scratch.
[[nodiscard]] inline bool wideStateSelfCheck() {
    namespace gf = spin::research::gfni_r4;
    alignas(64) __m128i words[16], expectedWords[16], actualWords[16];
    gf::Packed<16> expectedPacked, actualPacked, inputPacked;
    alignas(64) std::uint64_t bitWords[8];
    for (unsigned inputBit = 0; inputBit < 2048; ++inputBit) {
        for (auto& word : words) word = _mm_setzero_si128();
        for (auto& word : bitWords) word = 0;
        bitWords[(inputBit & 127) >> 6] = std::uint64_t(1) << (inputBit & 63);
        words[inputBit >> 7] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(bitWords));
        gf::pack<16>(words, expectedPacked);
        widePack(words, actualPacked);
        if (std::memcmp(&actualPacked, &expectedPacked, sizeof(expectedPacked))) return false;

        for (auto& vector : inputPacked.v) vector = _mm512_setzero_si512();
        for (auto& word : bitWords) word = 0;
        bitWords[(inputBit & 511) >> 6] = std::uint64_t(1) << (inputBit & 63);
        inputPacked.v[inputBit >> 9] = _mm512_load_si512(bitWords);
        gf::unpack(inputPacked, expectedWords);
        wideUnpack(inputPacked, actualWords);
        if (std::memcmp(actualWords, expectedWords, sizeof(expectedWords))) return false;
    }
    return true;
}

} // namespace spin::research::packet_inner
