#pragma once

// Keep the original AES-polynomial-basis implementation available unchanged
// arithmetically as an A/B control. Selection below is entirely compile-time.
#include "Transform256Polynomial.h"
#include <type_traits>

namespace spin::research::one_route {

#if defined(_MSC_VER)
#define SPIN_OR_CANTOR_INLINE __forceinline
#else
#define SPIN_OR_CANTOR_INLINE inline __attribute__((always_inline))
#endif

namespace transform_cantor_detail {

using transform_detail::mul;
using transform_detail::scale;
using transform_detail::laneConstants;
using transform_detail::pair;

// Delta(x)=x^2+x. In GF(2^8), Delta^8=0 and its kernel is GF2.
// Consequently it has one nilpotent Jordan chain of length eight. Starting
// with beta0=1, the root choices below extend through beta7; choosing the
// smaller root changes only the basis, not this existence argument.
// Construction is constexpr, bounded by 7*256 AES-byte candidates.
struct CantorBasis {
    std::array<std::uint8_t, 8> beta{};
    bool valid = false;

    constexpr CantorBasis() noexcept {
        beta[0] = 1;
        for (unsigned i = 1; i != 8; ++i) {
            bool found = false;
            for (unsigned candidate = 0; candidate != 256; ++candidate) {
                const auto value = static_cast<std::uint8_t>(candidate);
                if (std::uint8_t(mul(value, value) ^ value) == beta[i - 1]) {
                    beta[i] = value;
                    found = true;
                    break;
                }
            }
            if (!found) return;
        }
        valid = true;
    }
};
inline constexpr CantorBasis basis{};
static_assert(basis.valid, "AES GF256 Cantor chain must have eight entries");

constexpr std::uint8_t point(unsigned index) noexcept {
    std::uint8_t result = 0;
    for (unsigned bit = 0; bit != 8; ++bit)
        if (index & (1U << bit)) result ^= basis.beta[bit];
    return result;
}
constexpr std::uint8_t vanishing(unsigned level, std::uint8_t value) noexcept {
    for (unsigned i = 0; i != level; ++i)
        value = std::uint8_t(mul(value, value) ^ value);
    return value;
}
constexpr bool validateBasis() noexcept {
    if (basis.beta[0] != 1) return false;
    for (unsigned j = 0; j != 8; ++j) {
        for (unsigned i = 0; i != 8; ++i) {
            const auto expected = j < i ? 0 : basis.beta[j - i];
            if (vanishing(i, basis.beta[j]) != expected) return false;
        }
    }
    std::array<bool, 256> seen{};
    for (unsigned index = 0; index != 256; ++index) {
        const auto value = point(index);
        if (seen[value] || vanishing(8, value) != 0) return false;
        seen[value] = true;
    }
    return true;
}
static_assert(validateBasis(), "Cantor points and subspace normalizations must agree");

// Let V_i=span(beta0,...,beta_(i-1)). Delta^i kills V_i, has degree2^i,
// and is monic; it is therefore exactly s_i(X)=product_(v in V_i)(X+v).
// Moreover s_i(beta_i)=beta0=1. No normalization factors are needed.
// In characteristic two, Delta^i(X)=sum_(j subset i) X^(2^j), where
// j subset i means (j & ~i)==0. Thus every basis-conversion coefficient
// is 0 or1. This is polynomial-basis conversion, not a field-basis change.
constexpr std::uint8_t twiddle(unsigned level, unsigned offsetIndex) noexcept {
    return vanishing(level, point(offsetIndex));
}

template<unsigned Level, unsigned J>
SPIN_OR_CANTOR_INLINE void divisionTerm(Element* p, unsigned source, __m512i value) noexcept {
    if constexpr ((J & ~(Level - 1)) == 0) {
        constexpr unsigned half = 1U << (Level - 1);
        auto* target = p + source - half + (1U << J);
        _mm512_storeu_si512(target, _mm512_xor_si512(_mm512_loadu_si512(target), value));
    }
}
template<unsigned Level, unsigned J>
SPIN_OR_CANTOR_INLINE void divisionTerm(Element* p, unsigned source, __m128i value) noexcept {
    if constexpr ((J & ~(Level - 1)) == 0) {
        constexpr unsigned half = 1U << (Level - 1);
        auto* target = reinterpret_cast<__m128i*>(p + source - half + (1U << J));
        _mm_storeu_si128(target, _mm_xor_si128(_mm_loadu_si128(target), value));
    }
}
template<unsigned Level, unsigned J>
SPIN_OR_CANTOR_INLINE void divisionTermFour(Element* p, unsigned source, __m512i value) noexcept {
    if constexpr ((J & ~(Level - 1)) == 0) {
        constexpr unsigned half = 1U << (Level - 1);
        auto* target = p + 4 * (source - half + (1U << J));
        _mm512_storeu_si512(target, _mm512_xor_si512(_mm512_loadu_si512(target), value));
    }
}
template<unsigned Level, bool Four, class Vector, std::size_t... J>
SPIN_OR_CANTOR_INLINE void divisionTerms(Element* p, unsigned source, Vector value,
                                       std::index_sequence<J...>) noexcept {
    if constexpr (Four) (divisionTermFour<Level, J>(p, source, value), ...);
    else (divisionTerm<Level, J>(p, source, value), ...);
}
template<unsigned Level, bool Reverse, bool Four>
SPIN_OR_CANTOR_INLINE void split(Element* p) noexcept {
    if constexpr (Level > 1) {
        constexpr unsigned half = 1U << (Level - 1);
        constexpr unsigned step = Four ? 1 : (half >= 8 ? 4 : 1);
        // Every lower term has exponent<=half/2. Four-source batches do
        // not modify their own source values when half>=8.
        for (unsigned j = 0; j != half; j += step) {
            const unsigned source = Reverse ? half + j : 2 * half - step - j;
            if constexpr (Four) {
                const auto value = _mm512_loadu_si512(p + 4 * source);
                divisionTerms<Level, true>(p, source, value, std::make_index_sequence<Level - 1>{});
            } else if constexpr (step == 4) {
                const auto value = _mm512_loadu_si512(p + source);
                divisionTerms<Level, false>(p, source, value, std::make_index_sequence<Level - 1>{});
            } else {
                const auto value = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p + source));
                divisionTerms<Level, false>(p, source, value, std::make_index_sequence<Level - 1>{});
            }
        }
    }
}
template<unsigned Level, bool Reverse, bool Four>
SPIN_OR_CANTOR_INLINE void convert(Element* p) noexcept {
    if constexpr (Level > 0) {
        constexpr unsigned half = 1U << (Level - 1);
        // Descending long division gives (remainder,quotient). Reverse
        // traversal reconstructs f=remainder+s_(Level-1)*quotient.
        if constexpr (!Reverse) split<Level, false, Four>(p);
        convert<Level - 1, Reverse, Four>(p);
        convert<Level - 1, Reverse, Four>(p + (Four ? 4 : 1) * half);
        if constexpr (Reverse) split<Level, true, Four>(p);
    }
}

template<unsigned Half, unsigned Offset, bool Inverse>
SPIN_OR_CANTOR_INLINE __m512i smallStage(__m512i x) noexcept {
    static_assert(Half == 1 || Half == 2);
    constexpr unsigned c0 = twiddle(Half == 1 ? 0 : 1, Offset);
    constexpr unsigned c1 = twiddle(Half == 1 ? 0 : 1, Offset ^ (Half == 1 ? 2 : 0));
    constexpr int lowShuffle = Half == 1 ? 0xa0 : 0x44;
    constexpr int highShuffle = Half == 1 ? 0xf5 : 0xee;
    constexpr __mmask8 lowMask = Half == 1 ? 0x33 : 0x0f;
    constexpr __mmask8 highMask = static_cast<__mmask8>(~lowMask);
    const auto low = _mm512_shuffle_i32x4(x, x, lowShuffle);
    if constexpr (c0 == 0 && c1 == 0) {
        return _mm512_mask_xor_epi64(x, highMask, x, low);
    } else {
        const auto high = _mm512_shuffle_i32x4(x, x, highShuffle);
        if constexpr (!Inverse) {
            const auto coefficients = [] {
                if constexpr (Half == 1) return laneConstants<c0, c0 ^ 1, c1, c1 ^ 1>();
                else return laneConstants<c0, c0, c0 ^ 1, c0 ^ 1>();
            }();
            return _mm512_xor_si512(low, _mm512_gf2p8mul_epi8(high, coefficients));
        } else {
            const auto coefficients = [] {
                if constexpr (Half == 1) return laneConstants<c0, c0, c1, c1>();
                else return laneConstants<c0, c0, c0, c0>();
            }();
            const auto sum = _mm512_xor_si512(low, high);
            return _mm512_mask_xor_epi64(sum, lowMask, low,
                _mm512_gf2p8mul_epi8(sum, coefficients));
        }
    }
}
template<unsigned Offset, bool Inverse>
SPIN_OR_CANTOR_INLINE void leaf16(Element* p) noexcept {
    auto a = _mm512_loadu_si512(p);
    auto b = _mm512_loadu_si512(p + 4);
    auto c = _mm512_loadu_si512(p + 8);
    auto d = _mm512_loadu_si512(p + 12);
    if constexpr (!Inverse) {
        pair<twiddle(3, Offset), false>(a, c);
        pair<twiddle(3, Offset), false>(b, d);
        pair<twiddle(2, Offset), false>(a, b);
        pair<twiddle(2, Offset ^ 8), false>(c, d);
    }
#define SPIN_OR_CANTOR_SMALL(H) \
    a = smallStage<H, Offset, Inverse>(a); \
    b = smallStage<H, Offset ^ 4, Inverse>(b); \
    c = smallStage<H, Offset ^ 8, Inverse>(c); \
    d = smallStage<H, Offset ^ 12, Inverse>(d)
    if constexpr (Inverse) {
        SPIN_OR_CANTOR_SMALL(1);
        SPIN_OR_CANTOR_SMALL(2);
        pair<twiddle(2, Offset), true>(a, b);
        pair<twiddle(2, Offset ^ 8), true>(c, d);
        pair<twiddle(3, Offset), true>(a, c);
        pair<twiddle(3, Offset), true>(b, d);
    } else {
        SPIN_OR_CANTOR_SMALL(2);
        SPIN_OR_CANTOR_SMALL(1);
    }
#undef SPIN_OR_CANTOR_SMALL
    _mm512_storeu_si512(p, a);
    _mm512_storeu_si512(p + 4, b);
    _mm512_storeu_si512(p + 8, c);
    _mm512_storeu_si512(p + 12, d);
}
template<unsigned Level, unsigned Offset, bool Inverse, bool Four>
SPIN_OR_CANTOR_INLINE void stage(Element* p) noexcept {
    constexpr unsigned half = 1U << (Level - 1);
    constexpr unsigned c = twiddle(Level - 1, Offset);
    constexpr unsigned step = Four ? 1 : (half >= 4 ? 4 : 1);
    for (unsigned j = 0; j != half; j += step) {
        if constexpr (Four || step == 4) {
            auto* low = p + (Four ? 4 : 1) * j;
            auto* high = p + (Four ? 4 : 1) * (j + half);
            auto a = _mm512_loadu_si512(low), b = _mm512_loadu_si512(high);
            pair<c, Inverse>(a, b);
            _mm512_storeu_si512(low, a);
            _mm512_storeu_si512(high, b);
        } else {
            auto* low = reinterpret_cast<__m128i*>(p + j);
            auto* high = reinterpret_cast<__m128i*>(p + j + half);
            auto a = _mm_loadu_si128(low), b = _mm_loadu_si128(high);
            if constexpr (Inverse) b = _mm_xor_si128(b, a);
            if constexpr (c != 0) a = _mm_xor_si128(a, scale<c>(b));
            if constexpr (!Inverse) b = _mm_xor_si128(b, a);
            _mm_storeu_si128(low, a);
            _mm_storeu_si128(high, b);
        }
    }
}
template<unsigned Level, unsigned Offset, bool Inverse, bool Four>
SPIN_OR_CANTOR_INLINE void fft(Element* p) noexcept {
    if constexpr (Level == 4 && !Four) {
        leaf16<Offset, Inverse>(p);
    } else if constexpr (Level > 0) {
        constexpr unsigned half = 1U << (Level - 1);
        if constexpr (!Inverse) stage<Level, Offset, false, Four>(p);
        fft<Level - 1, Offset, Inverse, Four>(p);
        fft<Level - 1, Offset ^ half, Inverse, Four>(p + (Four ? 4 : 1) * half);
        if constexpr (Inverse) stage<Level, Offset, true, Four>(p);
    }
}

} // namespace transform_cantor_detail

// Same monomial coefficient interface and memory contract as the baseline.
// Output i is evaluation at evaluationPoint(i), not at numeric AES byte i.
// At N256 this is a permutation of all GF256 points. For smaller N it is
// a different additive subspace; consumers must use the declared points.
template<unsigned N>
struct TransformCantor {
    static_assert(N >= 1 && N <= 256 && (N & (N - 1)) == 0);
    static constexpr unsigned levels = transform_detail::logSize(N);
    // Precondition: index<256. Ordinary transforms use indices<N; coset
    // transforms use OffsetIndex+i. No ordering permutation is hidden.
    static constexpr std::uint8_t evaluationPoint(unsigned index) noexcept {
        return transform_cantor_detail::point(index);
    }
    // P_degree evaluated at evaluationPoint(pointIndex). Both indices<256.
    static constexpr std::uint8_t novelBasisValue(unsigned degree, unsigned pointIndex) noexcept {
        std::uint8_t result = 1;
        for (unsigned bit = 0; bit != 8; ++bit)
            if (degree & (1U << bit))
                result = transform_detail::mul(result,
                    transform_cantor_detail::twiddle(bit, pointIndex));
        return result;
    }
    static void evaluate(Element* values) noexcept {
        transform_cantor_detail::convert<levels, false, false>(values);
        transform_cantor_detail::fft<levels, 0, false, false>(values);
    }
    static void interpolate(Element* values) noexcept {
        transform_cantor_detail::fft<levels, 0, true, false>(values);
        transform_cantor_detail::convert<levels, true, false>(values);
    }
    static void evaluateFour(Element* values) noexcept {
        transform_cantor_detail::convert<levels, false, true>(values);
        transform_cantor_detail::fft<levels, 0, false, true>(values);
    }
    static void interpolateFour(Element* values) noexcept {
        transform_cantor_detail::fft<levels, 0, true, true>(values);
        transform_cantor_detail::convert<levels, true, true>(values);
    }

    // Explicit alternate interface: P_j=product_i s_i(X)^(bit_i(j)).
    // This basis is already normalized because s_i(beta_i)=1. These
    // calls skip monomial conversion; they do not silently reinterpret
    // monomial coefficients or change the Element's field coordinates.
    static void evaluateNovel(Element* values) noexcept {
        transform_cantor_detail::fft<levels, 0, false, false>(values);
    }
    static void interpolateNovel(Element* values) noexcept {
        transform_cantor_detail::fft<levels, 0, true, false>(values);
    }
    static void evaluateNovelFour(Element* values) noexcept {
        transform_cantor_detail::fft<levels, 0, false, true>(values);
    }
    static void interpolateNovelFour(Element* values) noexcept {
        transform_cantor_detail::fft<levels, 0, true, true>(values);
    }
    // The coefficient basis stays P_0,...,P_(N-1); only the evaluation
    // points change to evaluationPoint(OffsetIndex+i), 0<=i<N.
    template<unsigned OffsetIndex>
    static void evaluateNovelCoset(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_cantor_detail::fft<levels, OffsetIndex, false, false>(values);
    }
    template<unsigned OffsetIndex>
    static void interpolateNovelCoset(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_cantor_detail::fft<levels, OffsetIndex, true, false>(values);
    }
    template<unsigned OffsetIndex>
    static void evaluateNovelCosetFour(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_cantor_detail::fft<levels, OffsetIndex, false, true>(values);
    }
    template<unsigned OffsetIndex>
    static void interpolateNovelCosetFour(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_cantor_detail::fft<levels, OffsetIndex, true, true>(values);
    }
};

enum class TransformBasis { Polynomial, Cantor };
#ifndef SPIN_OR_TRANSFORM_CANTOR
#define SPIN_OR_TRANSFORM_CANTOR 0
#endif
static_assert(SPIN_OR_TRANSFORM_CANTOR == 0 || SPIN_OR_TRANSFORM_CANTOR == 1);
inline constexpr TransformBasis defaultTransformBasis = SPIN_OR_TRANSFORM_CANTOR
    ? TransformBasis::Cantor : TransformBasis::Polynomial;
template<unsigned N, TransformBasis Basis = defaultTransformBasis>
using Transform = std::conditional_t<Basis == TransformBasis::Cantor,
                                     TransformCantor<N>, TransformPolynomial<N>>;
using Transform256 = Transform<256>;
using Transform256Cantor = Transform<256, TransformBasis::Cantor>;

#undef SPIN_OR_CANTOR_INLINE
} // namespace spin::research::one_route
