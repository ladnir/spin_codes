#pragma once

#include <bit>
#include <cstdint>
#include <immintrin.h>

#if defined(_MSC_VER)
#define SPIN_PACKET_FEEDBACK_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
#define SPIN_PACKET_FEEDBACK_INLINE inline __attribute__((always_inline))
#else
#define SPIN_PACKET_FEEDBACK_INLINE inline
#endif

namespace spin::research::packet_inner {
namespace feedback_detail {

struct FourMoments {
    __m512i z0, z1, z2, z3;
};

// Superset zeta on two packet-index bits. XOR3 combines the all-input
// moment without the fourth binary XOR of a two-stage butterfly.
template<unsigned First>
SPIN_PACKET_FEEDBACK_INLINE FourMoments firstStage(const __m512i* packets) {
    const auto x0 = _mm512_loadu_si512(packets + First);
    const auto x1 = _mm512_loadu_si512(packets + First + 1);
    const auto x2 = _mm512_loadu_si512(packets + First + 2);
    const auto x3 = _mm512_loadu_si512(packets + First + 3);
    const auto z1 = _mm512_xor_si512(x1, x3);
    const auto z2 = _mm512_xor_si512(x2, x3);
    const auto z0 = _mm512_ternarylogic_epi64(x0, x2, z1, 0x96);
    return {z0, z1, z2, x3};
}

// Fold four 128-bit lanes. Degree is the remaining low-two-bit degree;
// all lane extractions and stores are fixed at compile time.
template<unsigned Mask, unsigned Degree>
SPIN_PACKET_FEEDBACK_INLINE void storeLowMoments(__m512i x, __m128i* moments) {
    static_assert((Mask & 3) == 0 && Mask < 64 && Degree <= 2);
    const auto upper = _mm512_extracti64x4_epi64(x, 1);
    const auto folded = _mm256_xor_si256(_mm512_castsi512_si256(x), upper);
    const auto odd = _mm256_extracti128_si256(folded, 1);
    _mm_storeu_si128(moments + Mask,
        _mm_xor_si128(_mm256_castsi256_si128(folded), odd));
    if constexpr (Degree >= 1) {
        const auto last = _mm256_extracti128_si256(upper, 1);
        _mm_storeu_si128(moments + Mask + 1, odd);
        _mm_storeu_si128(moments + Mask + 2,
            _mm_xor_si128(_mm256_castsi256_si128(upper), last));
        if constexpr (Degree == 2)
            _mm_storeu_si128(moments + Mask + 3, last);
    }
}

} // namespace feedback_detail

// Input word p is 128-bit lane (p & 3) of packets[p >> 2]. For every
// mask m with popcount(m) <= 2, write
//   moments[m] = XOR { input[p] : (p & m) == m }.
// The other 42 entries are untouched; the t64/s16 finisher reads only the
// 22 written entries. No input unpacking or temporary moment array is used.
// Requires AVX2, AVX-512F. Input loads and output stores may be unaligned.
SPIN_PACKET_FEEDBACK_INLINE void packetMoments(
    const __m512i* packets, __m128i* moments) {
    using feedback_detail::firstStage;
    using feedback_detail::storeLowMoments;

    const auto a = firstStage<0>(packets);
    const auto b = firstStage<4>(packets);
    const auto c = firstStage<8>(packets);
    const auto d = firstStage<12>(packets);

    // The remaining high-index zeta is pruned to degree two. Together with
    // firstStage, this needs 15 binary XORs and 8 XOR3 operations.
    // Low lane folding then needs 27 logical operations and 27 extracts.
    const auto z4 = _mm512_xor_si512(b.z0, d.z0);
    const auto z8 = _mm512_xor_si512(c.z0, d.z0);
    const auto z0 = _mm512_ternarylogic_epi64(a.z0, c.z0, z4, 0x96);
    storeLowMoments<0, 2>(z0, moments);
    storeLowMoments<16, 1>(z4, moments);
    storeLowMoments<32, 1>(z8, moments);
    storeLowMoments<48, 0>(d.z0, moments);

    const auto z5 = _mm512_xor_si512(b.z1, d.z1);
    const auto z9 = _mm512_xor_si512(c.z1, d.z1);
    const auto z1 = _mm512_ternarylogic_epi64(a.z1, c.z1, z5, 0x96);
    storeLowMoments<4, 1>(z1, moments);
    storeLowMoments<20, 0>(z5, moments);
    storeLowMoments<36, 0>(z9, moments);

    const auto z6 = _mm512_xor_si512(b.z2, d.z2);
    const auto z10 = _mm512_xor_si512(c.z2, d.z2);
    const auto z2 = _mm512_ternarylogic_epi64(a.z2, c.z2, z6, 0x96);
    storeLowMoments<8, 1>(z2, moments);
    storeLowMoments<24, 0>(z6, moments);
    storeLowMoments<40, 0>(z10, moments);

    const auto z3 = _mm512_ternarylogic_epi64(a.z3, b.z3,
        _mm512_xor_si512(c.z3, d.z3), 0x96);
    storeLowMoments<12, 0>(z3, moments);
}

// Optional setup-only check, with no executable or timing harness attached.
// Exhaust all 64*128 input-bit bases against the scalar superset formula.
// Since the implementation is linear, these bases determine its entire map.
// This has ordinary stack scratch and must not be moved into a coroutine.
[[nodiscard]] inline bool packetMomentsSelfCheck() {
    alignas(64) std::uint64_t words[128]{};
    alignas(64) __m512i packets[16];
    alignas(64) __m128i moments[64];
    std::uint64_t actual[2];
    for (unsigned p = 0; p < 64; ++p) {
        for (unsigned bit = 0; bit < 128; ++bit) {
            const auto value = std::uint64_t(1) << (bit & 63);
            words[2 * p + (bit >> 6)] = value;
            for (unsigned packet = 0; packet < 16; ++packet)
                packets[packet] = _mm512_load_si512(words + 8 * packet);
            packetMoments(packets, moments);
            for (unsigned mask = 0; mask < 64; ++mask) {
                if (std::popcount(mask) > 2) continue;
                _mm_storeu_si128(reinterpret_cast<__m128i*>(actual), moments[mask]);
                const auto expected = (p & mask) == mask ? value : 0;
                if (actual[0] != (bit < 64 ? expected : 0) ||
                    actual[1] != (bit >= 64 ? expected : 0)) return false;
            }
            words[2 * p + (bit >> 6)] = 0;
        }
    }
    return true;
}

} // namespace spin::research::packet_inner

#undef SPIN_PACKET_FEEDBACK_INLINE
