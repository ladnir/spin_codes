#pragma once

#include "FusedR4Gfni.h"
#include <array>
#include <cstdint>
#include <cstring>

namespace spin::research::packet_inner {
namespace vbmi_state_detail {

// Pack eight coordinate words, each sixteen bytes. GFNI's matrix rows
// consume coordinates in reverse order. Each output qword contains one
// payload byte from all eight coordinates; a direct VBMI byte permutation
// performs both the wide gather and the previous lane-local byte shuffle.
template<unsigned FirstPayloadByte>
static SPIN_FORCEINLINE __m512i packIndex() {
    static_assert(FirstPayloadByte == 0 || FirstPayloadByte == 8);
    alignas(64) static constexpr auto indices = [] {
        std::array<std::uint8_t,64> result{};
        for(unsigned byte = 0; byte < 8; ++byte)
            for(unsigned coordinate = 0; coordinate < 8; ++coordinate)
                result[8 * byte + coordinate] = static_cast<std::uint8_t>(
                    16 * (7 - coordinate) + FirstPayloadByte + byte);
        return result;
    }();
    return _mm512_load_si512(indices.data());
}

// After inverse GFNI, concatenating lo/hi gives sixteen eight-byte
// payload columns with ascending coordinate order. Gather one complete
// coordinate word into each output 128-bit lane.
template<unsigned FirstCoordinate>
static SPIN_FORCEINLINE __m512i unpackIndex() {
    static_assert(FirstCoordinate == 0 || FirstCoordinate == 4);
    alignas(64) static constexpr auto indices = [] {
        std::array<std::uint8_t,64> result{};
        for(unsigned coordinate = 0; coordinate < 4; ++coordinate)
            for(unsigned byte = 0; byte < 16; ++byte)
                result[16 * coordinate + byte] = static_cast<std::uint8_t>(
                    8 * byte + FirstCoordinate + coordinate);
        return result;
    }();
    return _mm512_load_si512(indices.data());
}

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
    const auto lo = _mm512_permutex2var_epi8(a, packIndex<0>(), b);
    const auto hi = _mm512_permutex2var_epi8(a, packIndex<8>(), b);
    const auto basis = _mm512_set1_epi64(0x0102040810204080ULL);
    out.v[2 * Group] = _mm512_gf2p8affine_epi64_epi8(basis, lo, 0);
    out.v[2 * Group + 1] = _mm512_gf2p8affine_epi64_epi8(basis, hi, 0);
}

template<unsigned First>
static SPIN_FORCEINLINE void storeFour(__m128i* out, __m512i value) {
    out[First] = _mm512_castsi512_si128(value);
    out[First + 1] = _mm512_extracti32x4_epi32(value, 1);
    out[First + 2] = _mm512_extracti32x4_epi32(value, 2);
    out[First + 3] = _mm512_extracti32x4_epi32(value, 3);
}

template<unsigned Group>
static SPIN_FORCEINLINE void unpackGroup(
    const spin::research::gfni_r4::Packed<16>& state, __m128i* out) {
    static_assert(Group < 2);
    const auto basis = _mm512_set1_epi64(0x8040201008040201ULL);
    const auto lo = _mm512_gf2p8affine_epi64_epi8(basis, state.v[2 * Group], 0);
    const auto hi = _mm512_gf2p8affine_epi64_epi8(basis, state.v[2 * Group + 1], 0);
    storeFour<8 * Group>(out, _mm512_permutex2var_epi8(lo, unpackIndex<0>(), hi));
    storeFour<8 * Group + 4>(out, _mm512_permutex2var_epi8(lo, unpackIndex<4>(), hi));
}

} // namespace vbmi_state_detail

// Exact storage-compatible conversions. Requires AVX-512F/BW/VBMI and
// GFNI; callers select this only on a supporting CPU/build. Source and
// destination storage must be disjoint, as in the retained conversions.
// Each direction: four VPERMI2B, four GFNI, no VPSHUFB, and twelve
// boundary 128-bit lane inserts/extracts before compiler simplification.
// All index-building loops above are constexpr, never part of the hot path.
static SPIN_FORCEINLINE void widePackVbmi(const __m128i* state,
    spin::research::gfni_r4::Packed<16>& out) {
    vbmi_state_detail::packGroup<0>(state, out);
    vbmi_state_detail::packGroup<1>(state, out);
}

static SPIN_FORCEINLINE void wideUnpackVbmi(
    const spin::research::gfni_r4::Packed<16>& state, __m128i* out) {
    vbmi_state_detail::unpackGroup<0>(state, out);
    vbmi_state_detail::unpackGroup<1>(state, out);
}

// Setup-only validation, not a round-trip-only test. Exhaust all 2048
// coordinate-word bits for pack and all 2048 independent packed bits for
// unpack. Both directions must equal the retained implementation AND the
// explicit scalar storage convention. SIMD scratch belongs in an ordinary
// function, not a coroutine. This check remains active under NDEBUG.
[[nodiscard]] inline bool vbmiStateSelfCheck() {
    namespace gf = spin::research::gfni_r4;
    alignas(64) __m128i inputWords[16], actualWords[16], expectedWords[16];
    alignas(64) std::uint8_t scalarWords[256], scalarPacked[256];
    gf::Packed<16> inputPacked, actualPacked, expectedPacked;
    for(unsigned inputBit = 0; inputBit < 2048; ++inputBit) {
        std::memset(inputWords, 0, sizeof(inputWords));
        std::memset(scalarPacked, 0, sizeof(scalarPacked));
        const unsigned coordinate = inputBit / 128;
        const unsigned payloadBit = inputBit % 128;
        reinterpret_cast<unsigned char*>(inputWords)[inputBit / 8] =
            static_cast<unsigned char>(1U << (inputBit % 8));
        // One packed byte contains eight coordinate bits at one payload
        // position. Within each payload byte, GFNI stores bit positions
        // in descending order, hence payloadBit XOR 7.
        scalarPacked[128 * (coordinate / 8) + (payloadBit ^ 7)] =
            static_cast<std::uint8_t>(1U << (coordinate % 8));
        gf::pack<16>(inputWords, expectedPacked);
        widePackVbmi(inputWords, actualPacked);
        if(std::memcmp(&actualPacked, &expectedPacked, sizeof(actualPacked)) ||
           std::memcmp(&actualPacked, scalarPacked, sizeof(actualPacked))) return false;

        std::memset(&inputPacked, 0, sizeof(inputPacked));
        std::memset(scalarWords, 0, sizeof(scalarWords));
        reinterpret_cast<unsigned char*>(&inputPacked)[inputBit / 8] =
            static_cast<unsigned char>(1U << (inputBit % 8));
        const unsigned packedByte = inputBit / 8;
        const unsigned unpackedCoordinate = 8 * (packedByte / 128) + (inputBit % 8);
        const unsigned unpackedPayloadBit = (packedByte % 128) ^ 7;
        scalarWords[16 * unpackedCoordinate + unpackedPayloadBit / 8] =
            static_cast<std::uint8_t>(1U << (unpackedPayloadBit % 8));
        gf::unpack(inputPacked, expectedWords);
        wideUnpackVbmi(inputPacked, actualWords);
        if(std::memcmp(actualWords, expectedWords, sizeof(actualWords)) ||
           std::memcmp(actualWords, scalarWords, sizeof(actualWords))) return false;
    }
    return true;
}

} // namespace spin::research::packet_inner
