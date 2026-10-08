#pragma once
#include "Full.h"

namespace spin::research::char2 {

enum class SubfieldInnerKernel : std::uint8_t { Baseline, Prefetch128 };

// All encoder coefficients lie in AES GF(256). A record of 1,4,8,16 bytes is
// interpreted in a fixed GF(256) basis of GF(2^8), GF(2^32), GF(2^64), GF(2^128).
// Multiplication by a subfield scalar acts independently on these components.
// Thus the map is linear over each full extension field, for any chosen basis:
// it commutes with every GF(256)-linear transformation of the record components,
// including multiplication by an arbitrary extension-field element. This is
// NOT the native polynomial-bit representation used by Full.h's other family.
// Scalar extension preserves symbol distance exactly; larger records do not
// acquire a larger distance guarantee merely by having more components.
// For the full-extension tests byte i denotes Z^i, using the independently
// checked irreducible moduli:
//   4 bytes: Z^4  + d4 Z^3  + d1 Z + 3b;
//   8 bytes: Z^8  + 0f Z^7  + 48 Z + 57;
//  16 bytes: Z^16 + 41 Z^15 + a5 Z + 63. (All constants are hexadecimal.)

constexpr std::uint8_t sfMultiply(std::uint8_t a, std::uint8_t b) noexcept {
    std::uint8_t result = 0;
    for (unsigned i = 0; i != 8; ++i) {
        if (b & 1U) result ^= a;
        a = std::uint8_t((unsigned(a) << 1) ^ ((a & 128U) ? 0x11bU : 0U));
        b >>= 1;
    }
    return result;
}

inline std::uint8_t sfInverse(std::uint8_t a) {
    if (!a) throw std::invalid_argument("zero subfield inverse");
    std::uint8_t result = 1;
    for (unsigned exponent = 254; exponent; exponent >>= 1) {
        if (exponent & 1U) result = sfMultiply(result, a);
        a = sfMultiply(a, a);
    }
    return result;
}

constexpr std::uint8_t sfSubspace(unsigned j, unsigned x) noexcept {
    std::uint8_t result = 1;
    for (unsigned v = 0; v != (1U << j); ++v)
        result = sfMultiply(result, std::uint8_t(x ^ v));
    return result;
}

constexpr std::uint8_t sfBasis(unsigned j, unsigned x) noexcept {
    std::uint8_t result = 1;
    for (unsigned degree = 0; j; ++degree, j >>= 1)
        if (j & 1U) result = sfMultiply(result, sfSubspace(degree, x));
    return result;
}

struct SubfieldMatrix {
    std::uint8_t a, b, c, d; // row-major [[a,b],[c,d]], determinant ad+bc=1.
};

// For a nonzero first column (a,c), exactly 256 second columns have determinant
// one. This map is a bijection from (nonzero first column, free byte) to SL2.
inline SubfieldMatrix sfMatrixFromColumn(std::uint8_t a, std::uint8_t c,
                                         std::uint8_t free) {
    if (a) return {a, free, c, sfMultiply(std::uint8_t(1 ^ sfMultiply(free, c)), sfInverse(a))};
    if (c) return {0, sfInverse(c), c, free};
    throw std::invalid_argument("SL2 first column must be nonzero");
}

template<class Generator>
inline SubfieldMatrix sfSampleMatrix(Generator& generator) {
    std::uint16_t column;
    do { column = std::uint16_t(generator()); } while (!column);
    return sfMatrixFromColumn(std::uint8_t(column), std::uint8_t(column >> 8),
                              std::uint8_t(generator()));
}

template<class Generator>
inline std::uint8_t sfSampleNonzero(Generator& generator) {
    std::uint8_t value;
    do { value = std::uint8_t(generator()); } while (!value);
    return value;
}

struct alignas(64) SubfieldIndexCoefficients {
    std::array<std::uint64_t, 8> bytes;
};

// Exact-law research candidate: GRS outer, a fresh nonzero diagonal per outer
// output, fixed-coordinate regions with independent uniform group permutations,
// and independently uniform SL2 refresh BEFORE each emission. No finite
// distance certificate is asserted by this implementation. mt19937_64 supplies
// reproducible test fixtures; it is not a production setup random generator.
struct SubfieldPlan {
    std::size_t k, n, groups;
    unsigned width, outer;
    std::uint64_t seed;
    std::vector<std::uint32_t> route;
    std::vector<std::uint8_t> routedDiagonal;
    std::vector<SubfieldMatrix> refresh;
    std::array<SubfieldIndexCoefficients, 4> indexCoefficients;

    SubfieldPlan(std::size_t message, unsigned t, std::uint64_t setupSeed, unsigned b = 16)
        : k(message), n(2 * message), groups(0),
          width(t), outer(b), seed(setupSeed) {
        if ((b != 8 && b != 16) || (t != 8 && t != 16) || !k || k % (b / 2) ||
            k > std::numeric_limits<std::uint32_t>::max() / 2 || n % t)
            throw std::invalid_argument("subfield candidate requires natural K, t8/16, outer8/16");
        groups = k / (outer / 2);
        std::mt19937_64 routeRng(seed ^ 0xdf4a1b3285ec7096ULL);
        std::mt19937_64 diagonalRng(seed ^ 0x2e7d18c4f0b3a965ULL);
        std::mt19937_64 matrixRng(seed ^ 0x918b745cf62d03eaULL);
        route.resize(n);
        for (unsigned region = 0; region != outer; ++region) {
            auto begin = route.begin() + region * groups;
            std::iota(begin, begin + groups, std::uint32_t(region * groups));
            for (std::size_t i = groups - 1; i; --i) {
                std::uniform_int_distribution<std::size_t> draw(0, i);
                std::swap(begin[i], begin[draw(routeRng)]);
            }
        }
        std::vector<std::uint8_t> originalDiagonal(n);
        for (auto& d : originalDiagonal) d = sfSampleNonzero(diagonalRng);
        routedDiagonal.resize(n);
        for (std::size_t i = 0; i != n; ++i) routedDiagonal[i] = originalDiagonal[route[i]];
        refresh.resize(n / width);
        for (auto& matrix : refresh) matrix = sfSampleMatrix(matrixRng);
        for (unsigned batch = 0; batch != 4; ++batch)
            for (unsigned lane = 0; lane != 4; ++lane) {
                const auto repeated = 0x0101010101010101ULL * (4 * batch + lane);
                indexCoefficients[batch].bytes[2 * lane] = repeated;
                indexCoefficients[batch].bytes[2 * lane + 1] = repeated;
            }
    }

    std::size_t setupBytes() const noexcept {
        return route.size() * sizeof(route[0]) + routedDiagonal.size() +
               refresh.size() * sizeof(refresh[0]) + sizeof(indexCoefficients);
    }
};

inline bool sfValidWidth(unsigned bytes) noexcept {
    return bytes == 1 || bytes == 4 || bytes == 8 || bytes == 16;
}

inline std::vector<std::uint8_t> subfieldForwardOracle(
    const SubfieldPlan& plan, const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.k * bytes)
        throw std::invalid_argument("subfield forward oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.n * bytes);
    const unsigned dimension = plan.outer / 2;
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned coordinate = 0; coordinate != plan.outer; ++coordinate)
            for (unsigned j = 0; j != dimension; ++j)
                for (unsigned component = 0; component != bytes; ++component)
                    planes[(coordinate * plan.groups + group) * bytes + component] ^=
                        sfMultiply(input[(group * dimension + j) * bytes + component], sfBasis(j, coordinate));
    std::array<std::uint8_t, 16> q0{}, q1{};
    for (std::size_t step = 0; step != plan.refresh.size(); ++step) {
        const auto m = plan.refresh[step];
        for (unsigned component = 0; component != bytes; ++component) {
            const auto r0 = sfMultiply(m.a, q0[component]) ^ sfMultiply(m.b, q1[component]);
            const auto r1 = sfMultiply(m.c, q0[component]) ^ sfMultiply(m.d, q1[component]);
            q0[component] = std::uint8_t(r0);
            q1[component] = std::uint8_t(r1);
        }
        std::array<std::uint8_t, 16> sum0{}, sum1{};
        for (unsigned j = 0; j != plan.width; ++j) {
            const auto i = step * plan.width + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto x = sfMultiply(planes[std::size_t(plan.route[i]) * bytes + component], plan.routedDiagonal[i]);
                output[i * bytes + component] = x ^ q0[component] ^ sfMultiply(std::uint8_t(j), q1[component]);
                sum0[component] ^= x;
                sum1[component] ^= sfMultiply(std::uint8_t(j), x);
            }
        }
        for (unsigned component = 0; component != bytes; ++component) {
            q0[component] ^= sum0[component];
            q1[component] ^= sum1[component];
        }
    }
    return output;
}

inline std::vector<std::uint8_t> subfieldTransposeOracle(
    const SubfieldPlan& plan, const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.n * bytes)
        throw std::invalid_argument("subfield transpose oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.k * bytes);
    std::array<std::uint8_t, 16> r0{}, r1{};
    for (std::size_t step = plan.refresh.size(); step--;) {
        std::array<std::uint8_t, 16> sum0{}, sum1{};
        for (unsigned j = 0; j != plan.width; ++j) {
            const auto i = step * plan.width + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto v = input[i * bytes + component];
                const auto u = std::uint8_t(v ^ r0[component] ^ sfMultiply(std::uint8_t(j), r1[component]));
                planes[std::size_t(plan.route[i]) * bytes + component] = sfMultiply(u, plan.routedDiagonal[i]);
                sum0[component] ^= v;
                sum1[component] ^= sfMultiply(std::uint8_t(j), v);
            }
        }
        const auto m = plan.refresh[step];
        for (unsigned component = 0; component != bytes; ++component) {
            const auto s0 = std::uint8_t(r0[component] ^ sum0[component]);
            const auto s1 = std::uint8_t(r1[component] ^ sum1[component]);
            r0[component] = sfMultiply(m.a, s0) ^ sfMultiply(m.c, s1);
            r1[component] = sfMultiply(m.b, s0) ^ sfMultiply(m.d, s1);
        }
    }
    const unsigned dimension = plan.outer / 2;
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned j = 0; j != dimension; ++j)
            for (unsigned coordinate = 0; coordinate != plan.outer; ++coordinate)
                for (unsigned component = 0; component != bytes; ++component)
                    output[(group * dimension + j) * bytes + component] ^=
                        sfMultiply(planes[(coordinate * plan.groups + group) * bytes + component], sfBasis(j, coordinate));
    return output;
}

inline __m512i sfBroadcastFour(const std::uint8_t* coefficients) noexcept {
    std::uint32_t packed;
    std::memcpy(&packed, coefficients, 4);
    const auto repeated = _mm512_broadcast_i32x4(_mm_cvtsi32_si128(static_cast<int>(packed)));
    const auto selectors = _mm512_setr_epi64(0, 0, 0x0101010101010101LL, 0x0101010101010101LL,
        0x0202020202020202LL, 0x0202020202020202LL, 0x0303030303030303LL, 0x0303030303030303LL);
    return _mm512_shuffle_epi8(repeated, selectors);
}

inline __m128i sfSumRecords4(__m512i x) noexcept {
    return _mm_xor_si128(
        _mm_xor_si128(_mm512_castsi512_si128(x), _mm512_extracti32x4_epi32(x, 1)),
        _mm_xor_si128(_mm512_extracti32x4_epi32(x, 2), _mm512_extracti32x4_epi32(x, 3)));
}

template<unsigned Batch, unsigned Prefetch = 0>
inline void sfInnerBatch(const SubfieldPlan& plan, const Block* input, Block* planes,
                          std::size_t begin, __m512i state0, __m512i state1,
                          __m512i& sum0, __m512i& sum1) {
    const auto i = begin + 4 * Batch;
    if constexpr (Prefetch != 0) {
        if (i >= Prefetch) {
            _mm_prefetch(reinterpret_cast<const char*>(planes + plan.route[i - Prefetch]), _MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(planes + plan.route[i - Prefetch + 1]), _MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(planes + plan.route[i - Prefetch + 2]), _MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(planes + plan.route[i - Prefetch + 3]), _MM_HINT_T0);
        }
    }
    const auto index = _mm512_load_si512(plan.indexCoefficients[Batch].bytes.data());
    const auto v = _mm512_loadu_si512(input + i);
    sum0 = _mm512_xor_si512(sum0, v);
    sum1 = _mm512_xor_si512(sum1, _mm512_gf2p8mul_epi8(v, index));
    auto u = _mm512_xor_si512(v, _mm512_xor_si512(state0, _mm512_gf2p8mul_epi8(state1, index)));
    u = _mm512_gf2p8mul_epi8(u, sfBroadcastFour(plan.routedDiagonal.data() + i));
    planes[plan.route[i]] = block(_mm512_castsi512_si128(u));
    planes[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
    planes[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
    planes[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
}

template<unsigned Width, unsigned Prefetch = 0>
inline void subfieldInnerTranspose(const SubfieldPlan& plan, const Block* input, Block* planes) {
    static_assert(Width == 8 || Width == 16);
    auto r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128();
    for (std::size_t step = plan.n / Width; step--;) {
        const auto begin = step * Width;
        const auto state0 = _mm512_broadcast_i32x4(r0), state1 = _mm512_broadcast_i32x4(r1);
        auto sum0 = _mm512_setzero_si512(), sum1 = _mm512_setzero_si512();
        sfInnerBatch<0, Prefetch>(plan, input, planes, begin, state0, state1, sum0, sum1);
        sfInnerBatch<1, Prefetch>(plan, input, planes, begin, state0, state1, sum0, sum1);
        if constexpr (Width == 16) {
            sfInnerBatch<2, Prefetch>(plan, input, planes, begin, state0, state1, sum0, sum1);
            sfInnerBatch<3, Prefetch>(plan, input, planes, begin, state0, state1, sum0, sum1);
        }
        const auto s0 = _mm_xor_si128(r0, sfSumRecords4(sum0));
        const auto s1 = _mm_xor_si128(r1, sfSumRecords4(sum1));
        const auto m = plan.refresh[step];
        // Four independent products; unlike three serial shears this has a
        // single multiplication level on the state dependency path.
        r0 = _mm_xor_si128(_mm_gf2p8mul_epi8(s0, _mm_set1_epi8(static_cast<char>(m.a))),
                           _mm_gf2p8mul_epi8(s1, _mm_set1_epi8(static_cast<char>(m.c))));
        r1 = _mm_xor_si128(_mm_gf2p8mul_epi8(s0, _mm_set1_epi8(static_cast<char>(m.b))),
                           _mm_gf2p8mul_epi8(s1, _mm_set1_epi8(static_cast<char>(m.d))));
    }
}

template<std::uint8_t C> inline __m512i sfConstant4(__m512i x) noexcept {
    if constexpr (C == 0) return _mm512_setzero_si512();
    else if constexpr (C == 1) return x;
    else return _mm512_gf2p8mul_epi8(x, _mm512_set1_epi8(static_cast<char>(C)));
}

template<unsigned Count, unsigned Offset>
inline void sfEvalTranspose4(__m512i* values) {
    if constexpr (Count > 1) {
        constexpr unsigned half = Count / 2, degree = std::countr_zero(half);
        sfEvalTranspose4<half, Offset>(values);
        sfEvalTranspose4<half, Offset ^ half>(values + half);
        for (unsigned j = 0; j != half; ++j) {
            const auto sum = _mm512_xor_si512(values[j], values[j + half]);
            values[j + half] = _mm512_xor_si512(sfConstant4<sfSubspace(degree, Offset)>(sum),
                sfConstant4<sfSubspace(degree, half)>(values[j + half]));
            values[j] = sum;
        }
    }
}

template<unsigned Outer>
inline void subfieldOuterTranspose(const SubfieldPlan& plan, const Block* planes, Block* output) {
    static_assert(Outer == 8 || Outer == 16);
    constexpr unsigned dimension = Outer / 2;
    std::size_t group = 0;
    for (; group + 4 <= plan.groups; group += 4) {
        __m512i values[Outer];
        for (unsigned coordinate = 0; coordinate != Outer; ++coordinate)
            values[coordinate] = _mm512_loadu_si512(planes + coordinate * plan.groups + group);
        sfEvalTranspose4<dimension, 0>(values);
        sfEvalTranspose4<dimension, dimension>(values + dimension);
        for (unsigned j = 0; j != dimension; ++j)
            values[j] = _mm512_xor_si512(values[j], values[j + dimension]);
        for (unsigned j = 0; j != dimension; j += 4) {
            const auto a01 = _mm512_shuffle_i32x4(values[j], values[j + 1], 0x44);
            const auto b01 = _mm512_shuffle_i32x4(values[j + 2], values[j + 3], 0x44);
            const auto a23 = _mm512_shuffle_i32x4(values[j], values[j + 1], 0xee);
            const auto b23 = _mm512_shuffle_i32x4(values[j + 2], values[j + 3], 0xee);
            _mm512_storeu_si512(output + dimension * group + j,
                _mm512_shuffle_i32x4(a01, b01, 0x88));
            _mm512_storeu_si512(output + dimension * (group + 1) + j,
                _mm512_shuffle_i32x4(a01, b01, 0xdd));
            _mm512_storeu_si512(output + dimension * (group + 2) + j,
                _mm512_shuffle_i32x4(a23, b23, 0x88));
            _mm512_storeu_si512(output + dimension * (group + 3) + j,
                _mm512_shuffle_i32x4(a23, b23, 0xdd));
        }
    }
    // Natural-size remainder; one field-element GFNI multiply per component
    // vector. This is a correctness fallback, not used by power-of-two sizes.
    for (; group != plan.groups; ++group)
        for (unsigned j = 0; j != dimension; ++j) {
            auto result = _mm_setzero_si128();
            for (unsigned coordinate = 0; coordinate != Outer; ++coordinate)
                result = _mm_xor_si128(result, _mm_gf2p8mul_epi8(reg(planes[coordinate * plan.groups + group]),
                    _mm_set1_epi8(static_cast<char>(sfBasis(j, coordinate)))));
            output[dimension * group + j] = block(result);
        }
}

inline void subfieldInnerTransposeDispatch(const SubfieldPlan& plan, const Block* input, Block* scratch,
                                           SubfieldInnerKernel kernel = SubfieldInnerKernel::Baseline) {
    if (kernel == SubfieldInnerKernel::Prefetch128) {
        if (plan.width == 8) subfieldInnerTranspose<8, 128>(plan, input, scratch);
        else subfieldInnerTranspose<16, 128>(plan, input, scratch);
    } else {
        if (plan.width == 8) subfieldInnerTranspose<8, 0>(plan, input, scratch);
        else subfieldInnerTranspose<16, 0>(plan, input, scratch);
    }
}

inline void subfieldOuterTransposeDispatch(const SubfieldPlan& plan, const Block* scratch, Block* output) {
    if (plan.outer == 8) subfieldOuterTranspose<8>(plan, scratch, output);
    else subfieldOuterTranspose<16>(plan, scratch, output);
}

inline void subfieldTransposeFast128(const SubfieldPlan& plan, const Block* input,
                                     Block* output, Block* scratch,
                                     SubfieldInnerKernel kernel = SubfieldInnerKernel::Baseline) {
    subfieldInnerTransposeDispatch(plan, input, scratch, kernel);
    subfieldOuterTransposeDispatch(plan, scratch, output);
}

} // namespace spin::research::char2
