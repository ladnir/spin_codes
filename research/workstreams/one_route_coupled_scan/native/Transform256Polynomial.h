#pragma once

#include "SparseByteField128Simd.h"
#include <array>
#include <cstdint>
#include <utility>

namespace spin::research::one_route {

using Element = char2::field128_lift::sparse_byte::Element;

#if defined(_MSC_VER)
#define SPIN_OR_FFT_INLINE __forceinline
#else
#define SPIN_OR_FFT_INLINE inline __attribute__((always_inline))
#endif

namespace transform_detail {

constexpr std::uint8_t mul(std::uint8_t a, std::uint8_t b) noexcept {
    return char2::field128_lift::sparse_byte::byteMultiply(a, b);
}
constexpr std::uint8_t inverse(std::uint8_t a) noexcept {
    std::uint8_t result = 1;
    for (unsigned exponent = 254; exponent; exponent >>= 1) {
        if (exponent & 1) result = mul(result, a);
        a = mul(a, a);
    }
    return result;
}

// s_i(X)=product_(v=0)^(2^i-1)(X+v), represented by coefficients of
// X,X^2,...,X^(2^i). Every constant is in the AES GF256 subfield.
struct Subspaces {
    std::array<std::array<std::uint8_t, 9>, 9> polynomial{};
    std::array<std::uint8_t, 8> scale{};

    constexpr Subspaces() noexcept {
        polynomial[0][0] = 1;
        for (unsigned i = 0; i != 8; ++i) {
            std::uint8_t power = static_cast<std::uint8_t>(1U << i);
            for (unsigned j = 0; j <= i; ++j) {
                scale[i] ^= mul(polynomial[i][j], power);
                power = mul(power, power);
            }
            for (unsigned j = 0; j <= i; ++j) {
                polynomial[i + 1][j] ^= mul(scale[i], polynomial[i][j]);
                polynomial[i + 1][j + 1] ^= mul(polynomial[i][j], polynomial[i][j]);
            }
        }
    }
};
inline constexpr Subspaces subspaces{};
static_assert(subspaces.polynomial[8][0] == 1 && subspaces.polynomial[8][8] == 1);

constexpr std::uint8_t normalized(unsigned level, std::uint8_t offset) noexcept {
    std::uint8_t value = 0, power = offset;
    for (unsigned j = 0; j <= level; ++j) {
        value ^= mul(subspaces.polynomial[level][j], power);
        power = mul(power, power);
    }
    return mul(value, inverse(subspaces.scale[level]));
}
constexpr std::uint8_t basisScale(unsigned index) noexcept {
    std::uint8_t result = 1;
    for (unsigned bit = 0; bit != 8; ++bit)
        if (index & (1U << bit)) result = mul(result, subspaces.scale[bit]);
    return result;
}

template<unsigned C>
SPIN_OR_FFT_INLINE __m512i scale(__m512i x) noexcept {
    if constexpr (C == 0) return _mm512_setzero_si512();
    else if constexpr (C == 1) return x;
    else return _mm512_gf2p8mul_epi8(x, _mm512_set1_epi8(static_cast<char>(C)));
}
template<unsigned C>
SPIN_OR_FFT_INLINE __m128i scale(__m128i x) noexcept {
    if constexpr (C == 0) return _mm_setzero_si128();
    else if constexpr (C == 1) return x;
    else return _mm_gf2p8mul_epi8(x, _mm_set1_epi8(static_cast<char>(C)));
}
template<unsigned A, unsigned B, unsigned C, unsigned D>
SPIN_OR_FFT_INLINE __m512i laneConstants() noexcept {
    constexpr auto repeat = [](unsigned x) {
        return static_cast<long long>(0x0101010101010101ULL * x);
    };
    return _mm512_setr_epi64(repeat(A), repeat(A), repeat(B), repeat(B),
                            repeat(C), repeat(C), repeat(D), repeat(D));
}

template<unsigned Level, unsigned J>
SPIN_OR_FFT_INLINE void divisionTerm(Element* p, unsigned source, __m512i value) noexcept {
    constexpr unsigned c = subspaces.polynomial[Level - 1][J];
    if constexpr (c != 0) {
        constexpr unsigned half = 1U << (Level - 1);
        auto* target = p + source - half + (1U << J);
        _mm512_storeu_si512(target, _mm512_xor_si512(_mm512_loadu_si512(target), scale<c>(value)));
    }
}
template<unsigned Level, unsigned J>
SPIN_OR_FFT_INLINE void divisionTerm(Element* p, unsigned source, __m128i value) noexcept {
    constexpr unsigned c = subspaces.polynomial[Level - 1][J];
    if constexpr (c != 0) {
        constexpr unsigned half = 1U << (Level - 1);
        auto* target = reinterpret_cast<__m128i*>(p + source - half + (1U << J));
        _mm_storeu_si128(target, _mm_xor_si128(_mm_loadu_si128(target), scale<c>(value)));
    }
}
template<unsigned Level, unsigned J>
SPIN_OR_FFT_INLINE void divisionTermFour(Element* p, unsigned source, __m512i value) noexcept {
    constexpr unsigned c = subspaces.polynomial[Level - 1][J];
    if constexpr (c != 0) {
        constexpr unsigned half = 1U << (Level - 1);
        auto* target = p + 4 * (source - half + (1U << J));
        _mm512_storeu_si512(target, _mm512_xor_si512(_mm512_loadu_si512(target), scale<c>(value)));
    }
}
template<unsigned Level, bool Four, class Vector, std::size_t... J>
SPIN_OR_FFT_INLINE void divisionTerms(Element* p, unsigned source, Vector value,
                                    std::index_sequence<J...>) noexcept {
    if constexpr (Four) (divisionTermFour<Level, J>(p, source, value), ...);
    else (divisionTerm<Level, J>(p, source, value), ...);
}

// Divide by the monic sparse s_(Level-1). Descending source order forms
// (remainder,quotient); ascending order reverses the same XOR operations.
// For half>=8, the nearest dependency is half/2>=4 records away, so a
// four-record load is safe before all its lower-coordinate updates.
template<unsigned Level, bool Reverse, bool Four>
SPIN_OR_FFT_INLINE void split(Element* p) noexcept {
    if constexpr (Level > 1) {
        constexpr unsigned half = 1U << (Level - 1);
        constexpr unsigned step = Four ? 1 : (half >= 8 ? 4 : 1);
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
SPIN_OR_FFT_INLINE void convert(Element* p) noexcept {
    if constexpr (Level > 0) {
        constexpr unsigned half = 1U << (Level - 1);
        if constexpr (!Reverse) split<Level, false, Four>(p);
        convert<Level - 1, Reverse, Four>(p);
        convert<Level - 1, Reverse, Four>(p + (Four ? 4 : 1) * half);
        if constexpr (Reverse) split<Level, true, Four>(p);
    }
}

template<unsigned Index, bool Reverse>
constexpr unsigned normalization() noexcept {
    if constexpr (Reverse) return inverse(basisScale(Index));
    else return basisScale(Index);
}
template<unsigned First, bool Reverse>
SPIN_OR_FFT_INLINE void normalizeFourRecords(Element* p) noexcept {
    constexpr unsigned a = normalization<First, Reverse>();
    constexpr unsigned b = normalization<First + 1, Reverse>();
    constexpr unsigned c = normalization<First + 2, Reverse>();
    constexpr unsigned d = normalization<First + 3, Reverse>();
    if constexpr (a != 1 || b != 1 || c != 1 || d != 1) {
        const auto value = _mm512_loadu_si512(p + First);
        _mm512_storeu_si512(p + First,
            _mm512_gf2p8mul_epi8(value, laneConstants<a, b, c, d>()));
    }
}
template<unsigned Index, bool Reverse, bool Four>
SPIN_OR_FFT_INLINE void normalizeOne(Element* p) noexcept {
    constexpr unsigned c = normalization<Index, Reverse>();
    if constexpr (c != 1) {
        if constexpr (Four) {
            _mm512_storeu_si512(p + 4 * Index, scale<c>(_mm512_loadu_si512(p + 4 * Index)));
        } else {
            auto* target = reinterpret_cast<__m128i*>(p + Index);
            _mm_storeu_si128(target, scale<c>(_mm_loadu_si128(target)));
        }
    }
}
template<unsigned N, bool Reverse, bool Four, std::size_t... J>
SPIN_OR_FFT_INLINE void normalize(Element* p, std::index_sequence<J...>) noexcept {
    if constexpr (!Four && N >= 4) (normalizeFourRecords<4 * J, Reverse>(p), ...);
    else (normalizeOne<J, Reverse, Four>(p), ...);
}

template<unsigned C, bool Inverse>
SPIN_OR_FFT_INLINE void pair(__m512i& low, __m512i& high) noexcept {
    if constexpr (Inverse) high = _mm512_xor_si512(high, low);
    if constexpr (C != 0) low = _mm512_xor_si512(low, scale<C>(high));
    if constexpr (!Inverse) high = _mm512_xor_si512(high, low);
}

// In-register h=2/1 butterflies for four consecutive Element records.
template<unsigned Half, unsigned Offset, bool Inverse>
SPIN_OR_FFT_INLINE __m512i smallStage(__m512i x) noexcept {
    static_assert(Half == 1 || Half == 2);
    constexpr unsigned c0 = normalized(Half == 1 ? 0 : 1, Offset);
    constexpr unsigned c1 = normalized(Half == 1 ? 0 : 1, Offset ^ (Half == 1 ? 2 : 0));
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
SPIN_OR_FFT_INLINE void leaf16(Element* p) noexcept {
    auto a = _mm512_loadu_si512(p);
    auto b = _mm512_loadu_si512(p + 4);
    auto c = _mm512_loadu_si512(p + 8);
    auto d = _mm512_loadu_si512(p + 12);
    if constexpr (!Inverse) {
        pair<normalized(3, Offset), false>(a, c);
        pair<normalized(3, Offset), false>(b, d);
        pair<normalized(2, Offset), false>(a, b);
        pair<normalized(2, Offset ^ 8), false>(c, d);
    }
#define SPIN_OR_SMALL_STAGE(H) \
    a = smallStage<H, Offset, Inverse>(a); \
    b = smallStage<H, Offset ^ 4, Inverse>(b); \
    c = smallStage<H, Offset ^ 8, Inverse>(c); \
    d = smallStage<H, Offset ^ 12, Inverse>(d)
    if constexpr (Inverse) {
        SPIN_OR_SMALL_STAGE(1);
        SPIN_OR_SMALL_STAGE(2);
        pair<normalized(2, Offset), true>(a, b);
        pair<normalized(2, Offset ^ 8), true>(c, d);
        pair<normalized(3, Offset), true>(a, c);
        pair<normalized(3, Offset), true>(b, d);
    } else {
        SPIN_OR_SMALL_STAGE(2);
        SPIN_OR_SMALL_STAGE(1);
    }
#undef SPIN_OR_SMALL_STAGE
    _mm512_storeu_si512(p, a);
    _mm512_storeu_si512(p + 4, b);
    _mm512_storeu_si512(p + 8, c);
    _mm512_storeu_si512(p + 12, d);
}

template<unsigned Level, unsigned Offset, bool Inverse, bool Four>
SPIN_OR_FFT_INLINE void stage(Element* p) noexcept {
    constexpr unsigned half = 1U << (Level - 1);
    constexpr unsigned c = normalized(Level - 1, Offset);
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
SPIN_OR_FFT_INLINE void fft(Element* p) noexcept {
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

constexpr unsigned logSize(unsigned n) noexcept {
    unsigned result = 0;
    while (n > 1) { n >>= 1; ++result; }
    return result;
}

} // namespace transform_detail

// In-place, unaligned-safe evaluation/interpolation on the AES subfield
// points 0,...,N-1. Both public coefficient interfaces are MONOMIAL.
// Each Element is one GF(2^128) value in sparse_byte coordinates. Subfield
// factors act on its sixteen AES bytes; no field-coordinate conversion is
// needed. There is no allocation or scratch array. Basis conversion uses
// O(N log^2 N) sparse operations, not a free or dense hidden transform.
template<unsigned N>
struct TransformPolynomial {
    static_assert(N >= 1 && N <= 256 && (N & (N - 1)) == 0);
    static constexpr unsigned levels = transform_detail::logSize(N);
    // Precondition: index<256. Ordinary transforms use indices<N; coset
    // transforms use OffsetIndex+i. Numeric AES-byte order is the baseline.
    static constexpr std::uint8_t evaluationPoint(unsigned index) noexcept {
        return static_cast<std::uint8_t>(index);
    }
    // P_degree evaluated at evaluationPoint(pointIndex). Both indices<256.
    static constexpr std::uint8_t novelBasisValue(unsigned degree, unsigned pointIndex) noexcept {
        std::uint8_t result = 1;
        for (unsigned bit = 0; bit != 8; ++bit)
            if (degree & (1U << bit))
                result = transform_detail::mul(result,
                    transform_detail::normalized(bit, evaluationPoint(pointIndex)));
        return result;
    }

    static void evaluate(Element* values) noexcept {
        transform_detail::convert<levels, false, false>(values);
        transform_detail::normalize<N, false, false>(values, std::make_index_sequence<(N >= 4 ? N / 4 : N)>{});
        transform_detail::fft<levels, 0, false, false>(values);
    }
    static void interpolate(Element* values) noexcept {
        transform_detail::fft<levels, 0, true, false>(values);
        transform_detail::normalize<N, true, false>(values, std::make_index_sequence<(N >= 4 ? N / 4 : N)>{});
        transform_detail::convert<levels, true, false>(values);
    }

    // Four independent polynomials, point-major: values[4*degree+lane].
    // Each SIMD register is one coefficient from all four polynomials.
    static void evaluateFour(Element* values) noexcept {
        transform_detail::convert<levels, false, true>(values);
        transform_detail::normalize<N, false, true>(values, std::make_index_sequence<N>{});
        transform_detail::fft<levels, 0, false, true>(values);
    }
    static void interpolateFour(Element* values) noexcept {
        transform_detail::fft<levels, 0, true, true>(values);
        transform_detail::normalize<N, true, true>(values, std::make_index_sequence<N>{});
        transform_detail::convert<levels, true, true>(values);
    }

    // Explicit alternate interface: coefficients are in the NORMALIZED
    // novel basis P_j=product_i (s_i(X)/s_i(2^i))^(bit_i(j)). These calls
    // intentionally perform neither monomial conversion nor normalization.
    static void evaluateNovel(Element* values) noexcept {
        transform_detail::fft<levels, 0, false, false>(values);
    }
    static void interpolateNovel(Element* values) noexcept {
        transform_detail::fft<levels, 0, true, false>(values);
    }
    static void evaluateNovelFour(Element* values) noexcept {
        transform_detail::fft<levels, 0, false, true>(values);
    }
    static void interpolateNovelFour(Element* values) noexcept {
        transform_detail::fft<levels, 0, true, true>(values);
    }
    // The coefficient basis stays P_0,...,P_(N-1); only the evaluation
    // points change to evaluationPoint(OffsetIndex+i), 0<=i<N.
    template<unsigned OffsetIndex>
    static void evaluateNovelCoset(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_detail::fft<levels, OffsetIndex, false, false>(values);
    }
    template<unsigned OffsetIndex>
    static void interpolateNovelCoset(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_detail::fft<levels, OffsetIndex, true, false>(values);
    }
    template<unsigned OffsetIndex>
    static void evaluateNovelCosetFour(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_detail::fft<levels, OffsetIndex, false, true>(values);
    }
    template<unsigned OffsetIndex>
    static void interpolateNovelCosetFour(Element* values) noexcept {
        static_assert(OffsetIndex <= 256 - N && OffsetIndex % N == 0);
        transform_detail::fft<levels, OffsetIndex, true, true>(values);
    }
};

using Transform256Polynomial = TransformPolynomial<256>;

#undef SPIN_OR_FFT_INLINE
} // namespace spin::research::one_route
