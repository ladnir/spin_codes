#pragma once
#include "Cauchy256Inner.h"

namespace spin::research::char2 {

// Exact-map alternative in the ORIGINAL certified residue-state coordinates.
// A = (FFT8(label0^16); FFT8(label8^16)) IFFT8(0), C = A^T.
// State refresh is unchanged; no coefficient-basis diagonal is substituted.
// Full sixteen-record packet, including the external input diagonal:
// 34 GFNI512 + 53 XOR/masked-XOR + 46 128-bit-lane shuffles.
// These hooks exclude input D: 30 GFNI512. The two fixed maps alone use 26.
// Counts include a separate two-vector eta scale, not an unproved free fold.
// All butterflies are compile-time specialized, with two live packed vectors.
// Same alias/diagonal contracts as Cauchy256Inner.h; no heap or large scratch.
namespace c256_inner_fft_detail {

inline constexpr std::uint8_t w[8] = {0, 1, 92, 93, 224, 225, 188, 189};

constexpr std::uint8_t power(std::uint8_t a, unsigned exponent) noexcept {
    std::uint8_t result = 1;
    while (exponent) {
        if (exponent & 1U) result = sfMultiply(result, a);
        a = sfMultiply(a, a);
        exponent >>= 1;
    }
    return result;
}

constexpr std::uint8_t normalized(unsigned level, std::uint8_t offset) noexcept {
    std::uint8_t value = 1, scale = 1;
    for (unsigned i = 0; i != (1U << level); ++i) {
        value = sfMultiply(value, std::uint8_t(offset ^ w[i]));
        scale = sfMultiply(scale, std::uint8_t(w[1U << level] ^ w[i]));
    }
    return sfMultiply(value, power(scale, 254));
}

template<std::uint8_t A, std::uint8_t B, std::uint8_t C, std::uint8_t D>
inline __m512i coefficients() noexcept {
    constexpr auto repeat = [](std::uint8_t x) {
        return static_cast<long long>(0x0101010101010101ULL*x);
    };
    return _mm512_setr_epi64(repeat(A), repeat(A), repeat(B), repeat(B),
                            repeat(C), repeat(C), repeat(D), repeat(D));
}

enum class Kind { Forward, Inverse, ForwardTranspose, InverseTranspose };

template<unsigned Half, std::uint8_t Offset, unsigned First, Kind K>
inline __m512i smallStage(__m512i x) noexcept {
    static_assert(Half == 1 || Half == 2);
    static_assert(First == 0 || First == 4);
    constexpr auto c0 = normalized(Half == 1 ? 0 : 1, std::uint8_t(Offset ^ w[First]));
    constexpr auto c1 = normalized(Half == 1 ? 0 : 1,
        std::uint8_t(Offset ^ w[First + (Half == 1 ? 2 : 0)]));
    constexpr int lowShuffle = Half == 1 ? 0xa0 : 0x44;
    constexpr int highShuffle = Half == 1 ? 0xf5 : 0xee;
    constexpr __mmask8 lowMask = Half == 1 ? 0x33 : 0x0f;
    constexpr __mmask8 highMask = static_cast<__mmask8>(~lowMask);
    if constexpr (c0 == 0 && c1 == 0) {
        if constexpr (K == Kind::Forward || K == Kind::Inverse) {
            const auto a = _mm512_shuffle_i32x4(x, x, lowShuffle);
            return _mm512_mask_xor_epi64(x, highMask, x, a);
        } else {
            const auto b = _mm512_shuffle_i32x4(x, x, highShuffle);
            return _mm512_mask_xor_epi64(x, lowMask, x, b);
        }
    } else {
        const auto a = _mm512_shuffle_i32x4(x, x, lowShuffle);
        const auto b = _mm512_shuffle_i32x4(x, x, highShuffle);
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
                return _mm512_mask_xor_epi64(sum, lowMask, a, product);
            else return _mm512_mask_xor_epi64(sum, highMask, b, product);
        }
    }
}

template<std::uint8_t Offset, Kind K>
inline void topStage(__m512i& a, __m512i& b) noexcept {
    constexpr auto c = normalized(2, Offset);
    if constexpr (c == 0) {
        if constexpr (K == Kind::Forward || K == Kind::Inverse)
            b = _mm512_xor_si512(a, b);
        else a = _mm512_xor_si512(a, b);
    } else if constexpr (K == Kind::Forward) {
        a = _mm512_xor_si512(a, sfConstant4<c>(b));
        b = _mm512_xor_si512(a, b);
    } else if constexpr (K == Kind::Inverse) {
        b = _mm512_xor_si512(a, b);
        a = _mm512_xor_si512(a, sfConstant4<c>(b));
    } else if constexpr (K == Kind::ForwardTranspose) {
        a = _mm512_xor_si512(a, b);
        b = _mm512_xor_si512(b, sfConstant4<c>(a));
    } else {
        b = _mm512_xor_si512(b, sfConstant4<c>(a));
        a = _mm512_xor_si512(a, b);
    }
}

template<std::uint8_t Offset, Kind K>
inline void transform(__m512i& a, __m512i& b) noexcept {
    if constexpr (K == Kind::Forward || K == Kind::InverseTranspose) {
        topStage<Offset, K>(a, b);
        a = smallStage<2, Offset, 0, K>(a);
        b = smallStage<2, Offset, 4, K>(b);
        a = smallStage<1, Offset, 0, K>(a);
        b = smallStage<1, Offset, 4, K>(b);
    } else {
        a = smallStage<1, Offset, 0, K>(a);
        b = smallStage<1, Offset, 4, K>(b);
        a = smallStage<2, Offset, 0, K>(a);
        b = smallStage<2, Offset, 4, K>(b);
        topStage<Offset, K>(a, b);
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

} // namespace c256_inner_fft_detail

template<unsigned Phase>
inline void c256InnerTransposeFftPacked(const __m512i* input4, __m512i* output4,
                                       __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {input4[0], input4[1], input4[2], input4[3]};
    const auto h0 = state2[0], h1 = state2[1];
    __m512i m0, m1, emitted[4];
    c256_inner_fft_detail::inject<Phase>(input, m0, m1);
    m0 = sfConstant4<2>(m0); m1 = sfConstant4<2>(m1);
    c256_inner_fft_detail::emit<Phase>(h0, h1, emitted);
    output4[0] = _mm512_xor_si512(input[0], emitted[0]);
    output4[1] = _mm512_xor_si512(input[1], emitted[1]);
    output4[2] = _mm512_xor_si512(input[2], emitted[2]);
    output4[3] = _mm512_xor_si512(input[3], emitted[3]);
    state2[0] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h0, m0), diagonal2[0]);
    state2[1] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h1, m1), diagonal2[1]);
}

template<unsigned Phase>
inline void c256InnerForwardFftPacked(const __m512i* input4, __m512i* output4,
                                     __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {input4[0], input4[1], input4[2], input4[3]};
    const auto r0 = _mm512_gf2p8mul_epi8(state2[0], diagonal2[0]);
    const auto r1 = _mm512_gf2p8mul_epi8(state2[1], diagonal2[1]);
    __m512i m0, m1, emitted[4];
    c256_inner_fft_detail::inject<Phase>(input, m0, m1);
    c256_inner_fft_detail::emit<Phase>(sfConstant4<2>(r0), sfConstant4<2>(r1), emitted);
    output4[0] = _mm512_xor_si512(input[0], emitted[0]);
    output4[1] = _mm512_xor_si512(input[1], emitted[1]);
    output4[2] = _mm512_xor_si512(input[2], emitted[2]);
    output4[3] = _mm512_xor_si512(input[3], emitted[3]);
    state2[0] = _mm512_xor_si512(r0, m0);
    state2[1] = _mm512_xor_si512(r1, m1);
}

} // namespace spin::research::char2
