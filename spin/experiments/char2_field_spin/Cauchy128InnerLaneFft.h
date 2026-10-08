#pragma once
#include "Cauchy256InnerFft.h"

namespace spin::research::char2 {

// Exact representation conjugation of Cauchy256InnerFft, NOT a new state basis.
// Ordinary ZMM layout is four 128-bit records. P transposes their 4x4 matrix
// of 32-bit chunks, so each 128-bit lane holds one chunk from all four records.
// P is self-inverse. All field operations still act independently on bytes.
//
// Hook contracts:
//   input4/output4: ordinary record layout; exact alias is supported.
//   state2: persistent P layout across EVERY packet and region.
//   diagonal2: P layout, produced by c128LaneDiagonalFour (NOT sfBroadcastFour).
// Zero initial state needs no conversion. For a nonzero scalar-reference state,
// apply c128LanePermutation to its two ordinary packed vectors, and invert P
// again only when observing the state. State/diagonal storage must not overlap
// input/output. External input-coordinate diagonals remain outside the hooks.
//
// Per packet: unchanged 26 fixed-map GFNI, 2 eta GFNI, 2 state GFNI;
// plus 4 external input GFNI = 34 full-packet GFNI. The maps use 46 lane-local
// VPSHUFD shuffles instead of 46 cross-128-bit shuffles. Eight boundary VPERMD
// operations convert the four inputs and outputs. No heap or large scratch.

inline __m512i c128LanePermutation(__m512i x) noexcept {
    const auto index = _mm512_setr_epi32(0, 4, 8, 12, 1, 5, 9, 13,
                                       2, 6, 10, 14, 3, 7, 11, 15);
    return _mm512_permutexvar_epi32(index, x);
}

inline __m512i c128LaneDiagonalFour(const std::uint8_t* coefficients) noexcept {
    std::uint32_t packed;
    std::memcpy(&packed, coefficients, 4);
    const auto repeated = _mm512_broadcast_i32x4(_mm_cvtsi32_si128(static_cast<int>(packed)));
    const auto selectors = _mm512_setr_epi32(
        0, 0x01010101, 0x02020202, 0x03030303,
        0, 0x01010101, 0x02020202, 0x03030303,
        0, 0x01010101, 0x02020202, 0x03030303,
        0, 0x01010101, 0x02020202, 0x03030303);
    return _mm512_shuffle_epi8(repeated, selectors);
}

namespace c128_lane_fft_detail {

using Kind = c256_inner_fft_detail::Kind;

template<std::uint8_t A, std::uint8_t B, std::uint8_t C, std::uint8_t D>
inline __m512i coefficients() noexcept {
    constexpr auto repeat = [](std::uint8_t x) {
        return static_cast<int>(0x01010101U*x);
    };
    return _mm512_setr_epi32(
        repeat(A), repeat(B), repeat(C), repeat(D),
        repeat(A), repeat(B), repeat(C), repeat(D),
        repeat(A), repeat(B), repeat(C), repeat(D),
        repeat(A), repeat(B), repeat(C), repeat(D));
}

template<unsigned Half, std::uint8_t Offset, unsigned First, Kind K>
inline __m512i smallStage(__m512i x) noexcept {
    static_assert(Half == 1 || Half == 2);
    static_assert(First == 0 || First == 4);
    constexpr auto c0 = c256_inner_fft_detail::normalized(
        Half == 1 ? 0 : 1, std::uint8_t(Offset ^ c256_inner_fft_detail::w[First]));
    constexpr auto c1 = c256_inner_fft_detail::normalized(Half == 1 ? 0 : 1,
        std::uint8_t(Offset ^ c256_inner_fft_detail::w[First + (Half == 1 ? 2 : 0)]));
    constexpr auto lowShuffle = static_cast<_MM_PERM_ENUM>(Half == 1 ? 0xa0 : 0x44);
    constexpr auto highShuffle = static_cast<_MM_PERM_ENUM>(Half == 1 ? 0xf5 : 0xee);
    constexpr __mmask16 lowMask = Half == 1 ? 0x5555 : 0x3333;
    constexpr __mmask16 highMask = static_cast<__mmask16>(~lowMask);
    if constexpr (c0 == 0 && c1 == 0) {
        if constexpr (K == Kind::Forward || K == Kind::Inverse) {
            const auto a = _mm512_shuffle_epi32(x, lowShuffle);
            return _mm512_mask_xor_epi32(x, highMask, x, a);
        } else {
            const auto b = _mm512_shuffle_epi32(x, highShuffle);
            return _mm512_mask_xor_epi32(x, lowMask, x, b);
        }
    } else {
        const auto a = _mm512_shuffle_epi32(x, lowShuffle);
        const auto b = _mm512_shuffle_epi32(x, highShuffle);
        if constexpr (K == Kind::Forward || K == Kind::InverseTranspose) {
            constexpr unsigned flip = K == Kind::Forward ? 0 : 1;
            const auto c = [&]() {
                if constexpr (Half == 1)
                    return coefficients<c0^flip, c0^flip^1, c1^flip, c1^flip^1>();
                else return coefficients<c0^flip, c0^flip, c0^flip^1, c0^flip^1>();
            }();
            if constexpr (K == Kind::Forward)
                return _mm512_xor_si512(a, _mm512_gf2p8mul_epi8(b, c));
            else return _mm512_xor_si512(b, _mm512_gf2p8mul_epi8(a, c));
        } else {
            const auto c = [&]() {
                if constexpr (Half == 1) return coefficients<c0, c0, c1, c1>();
                else return coefficients<c0, c0, c0, c0>();
            }();
            const auto sum = _mm512_xor_si512(a, b);
            const auto product = _mm512_gf2p8mul_epi8(sum, c);
            if constexpr (K == Kind::Inverse)
                return _mm512_mask_xor_epi32(sum, lowMask, a, product);
            else return _mm512_mask_xor_epi32(sum, highMask, b, product);
        }
    }
}

template<std::uint8_t Offset, Kind K>
inline void transform(__m512i& a, __m512i& b) noexcept {
    // The top stage pairs corresponding bytes across two ZMMs. Its scalar
    // coefficients commute with P, so the original helper is unchanged.
    if constexpr (K == Kind::Forward || K == Kind::InverseTranspose) {
        c256_inner_fft_detail::topStage<Offset, K>(a, b);
        a = smallStage<2, Offset, 0, K>(a);
        b = smallStage<2, Offset, 4, K>(b);
        a = smallStage<1, Offset, 0, K>(a);
        b = smallStage<1, Offset, 4, K>(b);
    } else {
        a = smallStage<1, Offset, 0, K>(a);
        b = smallStage<1, Offset, 4, K>(b);
        a = smallStage<2, Offset, 0, K>(a);
        b = smallStage<2, Offset, 4, K>(b);
        c256_inner_fft_detail::topStage<Offset, K>(a, b);
    }
}

template<unsigned Phase>
inline void emit(__m512i s0, __m512i s1, __m512i* output4) noexcept {
    static_assert(Phase < 4);
    transform<0, Kind::Inverse>(s0, s1);
    auto t0 = s0, t1 = s1;
    transform<c256Labels[Phase][0]^16, Kind::Forward>(s0, s1);
    transform<c256Labels[Phase][8]^16, Kind::Forward>(t0, t1);
    output4[0] = s0; output4[1] = s1;
    output4[2] = t0; output4[3] = t1;
}

template<unsigned Phase>
inline void inject(const __m512i* input4, __m512i& s0, __m512i& s1) noexcept {
    static_assert(Phase < 4);
    auto a0 = input4[0], a1 = input4[1], b0 = input4[2], b1 = input4[3];
    transform<c256Labels[Phase][0]^16, Kind::ForwardTranspose>(a0, a1);
    transform<c256Labels[Phase][8]^16, Kind::ForwardTranspose>(b0, b1);
    s0 = _mm512_xor_si512(a0, b0);
    s1 = _mm512_xor_si512(a1, b1);
    transform<0, Kind::InverseTranspose>(s0, s1);
}

} // namespace c128_lane_fft_detail

template<unsigned Phase>
inline void c128InnerTransposeLaneFftPacked(const __m512i* input4, __m512i* output4,
                                           __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {
        c128LanePermutation(input4[0]), c128LanePermutation(input4[1]),
        c128LanePermutation(input4[2]), c128LanePermutation(input4[3])};
    const auto h0 = state2[0], h1 = state2[1];
    __m512i m0, m1, emitted[4];
    c128_lane_fft_detail::inject<Phase>(input, m0, m1);
    m0 = sfConstant4<2>(m0); m1 = sfConstant4<2>(m1);
    c128_lane_fft_detail::emit<Phase>(h0, h1, emitted);
    output4[0] = c128LanePermutation(_mm512_xor_si512(input[0], emitted[0]));
    output4[1] = c128LanePermutation(_mm512_xor_si512(input[1], emitted[1]));
    output4[2] = c128LanePermutation(_mm512_xor_si512(input[2], emitted[2]));
    output4[3] = c128LanePermutation(_mm512_xor_si512(input[3], emitted[3]));
    state2[0] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h0, m0), diagonal2[0]);
    state2[1] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h1, m1), diagonal2[1]);
}

template<unsigned Phase>
inline void c128InnerForwardLaneFftPacked(const __m512i* input4, __m512i* output4,
                                         __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {
        c128LanePermutation(input4[0]), c128LanePermutation(input4[1]),
        c128LanePermutation(input4[2]), c128LanePermutation(input4[3])};
    const auto r0 = _mm512_gf2p8mul_epi8(state2[0], diagonal2[0]);
    const auto r1 = _mm512_gf2p8mul_epi8(state2[1], diagonal2[1]);
    __m512i m0, m1, emitted[4];
    c128_lane_fft_detail::inject<Phase>(input, m0, m1);
    c128_lane_fft_detail::emit<Phase>(sfConstant4<2>(r0), sfConstant4<2>(r1), emitted);
    output4[0] = c128LanePermutation(_mm512_xor_si512(input[0], emitted[0]));
    output4[1] = c128LanePermutation(_mm512_xor_si512(input[1], emitted[1]));
    output4[2] = c128LanePermutation(_mm512_xor_si512(input[2], emitted[2]));
    output4[3] = c128LanePermutation(_mm512_xor_si512(input[3], emitted[3]));
    state2[0] = _mm512_xor_si512(r0, m0);
    state2[1] = _mm512_xor_si512(r1, m1);
}

} // namespace spin::research::char2
