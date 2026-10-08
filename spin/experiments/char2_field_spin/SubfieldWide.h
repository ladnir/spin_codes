#pragma once
#include "Subfield16.h"

#if defined(_MSC_VER)
#define SPIN_WIDE_INLINE __forceinline
#else
#define SPIN_WIDE_INLINE inline __attribute__((always_inline))
#endif

namespace spin::research::char2 {

// Experimental six-state profile. The outer, embedded GF16 field, diagonal,
// and single regional route match Subfield16. Only the inner changes to
// A_j=(1,j,...,j^5), C=A^T, and independent uniform GL(6,16) pre-emission
// refreshes. This implementation does not assert a completed certificate.
struct SubfieldWideMatrix {
    // Six embedded coefficients per row plus two zero padding bytes. Padding
    // permits one bounded eight-byte load for the packed refresh operation.
    std::array<std::uint8_t, 48> entries{};
};

inline bool subfieldWideInvertible(std::array<std::uint8_t, 36> matrix) noexcept {
    for (unsigned column = 0; column != 6; ++column) {
        unsigned pivot = column;
        while (pivot != 6 && matrix[6 * pivot + column] == 0) ++pivot;
        if (pivot == 6) return false;
        if (pivot != column)
            for (unsigned j = column; j != 6; ++j)
                std::swap(matrix[6 * pivot + j], matrix[6 * column + j]);
        const auto inverse = sf16Inverse(matrix[7 * column]);
        for (unsigned row = column + 1; row != 6; ++row) {
            const auto factor = sf16Multiply(matrix[6 * row + column], inverse);
            for (unsigned j = column; j != 6; ++j)
                matrix[6 * row + j] ^= sf16Multiply(factor, matrix[6 * column + j]);
        }
    }
    return true;
}

template<class Generator>
inline SubfieldWideMatrix subfieldWideSampleMatrix(Generator& generator) {
    std::array<std::uint8_t, 36> abstract;
    do {
        std::uint64_t packed = 0;
        for (unsigned j = 0; j != 36; ++j) {
            if (j % 16 == 0) packed = std::uint64_t(generator());
            abstract[j] = std::uint8_t(packed & 15U);
            packed >>= 4;
        }
    } while (!subfieldWideInvertible(abstract));
    SubfieldWideMatrix matrix;
    for (unsigned row = 0; row != 6; ++row)
        for (unsigned column = 0; column != 6; ++column)
            matrix.entries[8 * row + column] = sf16Embed(abstract[6 * row + column]);
    return matrix;
}

struct SubfieldWidePlan {
    static constexpr unsigned width = 16, outer = 16, dimension = 8, state = 6;
    std::size_t k, n, groups;
    std::uint64_t seed;
    std::vector<std::uint32_t> route, inverseRoute;
    std::vector<std::uint8_t> routedDiagonal;
    std::vector<SubfieldWideMatrix> refresh;
    std::array<std::array<SubfieldIndexCoefficients, 5>, 4> coefficients;

    SubfieldWidePlan(std::size_t message, std::uint64_t setupSeed)
        : k(message), n(2 * message), groups(message / dimension), seed(setupSeed) {
        if (!k || k % dimension || k > std::numeric_limits<std::uint32_t>::max() / 2)
            throw std::invalid_argument("GF16 s6 needs nonzero K divisible by 8 within uint32 route limits");
        // Construct only this profile's setup: no throwaway Subfield16Plan or
        // GL4 sampling. Seeds are research fixtures, not production randomness.
        std::mt19937_64 routeRng(seed ^ 0xdf4a1b3285ec7096ULL);
        std::mt19937_64 diagonalRng(seed ^ 0x782d043af9165bc1ULL);
        std::mt19937_64 matrixRng(seed ^ 0xd489051f36a7cb2eULL);
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
        for (auto& entry : originalDiagonal) {
            std::uint8_t coefficient;
            do { coefficient = std::uint8_t(diagonalRng() & 15U); } while (!coefficient);
            entry = sf16Embed(coefficient);
        }
        inverseRoute.resize(n);
        routedDiagonal.resize(n);
        for (std::size_t i = 0; i != n; ++i) {
            inverseRoute[route[i]] = static_cast<std::uint32_t>(i);
            routedDiagonal[i] = originalDiagonal[route[i]];
        }
        refresh.resize(n / width);
        for (auto& matrix : refresh) matrix = subfieldWideSampleMatrix(matrixRng);
        for (unsigned batch = 0; batch != 4; ++batch)
            for (unsigned power = 1; power != 6; ++power)
                for (unsigned lane = 0; lane != 4; ++lane) {
                    const auto repeated = 0x0101010101010101ULL * sf16Embed(sf16Monomial(power, 4 * batch + lane));
                    coefficients[batch][power - 1].bytes[2 * lane] = repeated;
                    coefficients[batch][power - 1].bytes[2 * lane + 1] = repeated;
                }
    }

    bool regionAligned() const noexcept { return groups % width == 0; }
    std::size_t setupBytes() const noexcept {
        return (route.size() + inverseRoute.size()) * sizeof(std::uint32_t) + routedDiagonal.size() +
               refresh.size() * sizeof(SubfieldWideMatrix) + sizeof(coefficients);
    }
};

inline std::vector<std::uint8_t> subfieldWideForwardOracle(const SubfieldWidePlan& plan,
    const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.k * bytes)
        throw std::invalid_argument("GF16 s6 forward oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.n * bytes);
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned coordinate = 0; coordinate != 16; ++coordinate)
            for (unsigned j = 0; j != 8; ++j)
                for (unsigned component = 0; component != bytes; ++component)
                    planes[(coordinate * plan.groups + group) * bytes + component] ^=
                        sfMultiply(input[(group * 8 + j) * bytes + component], sf16Embed(sf16Basis(j, coordinate)));
    std::array<std::array<std::uint8_t, 16>, 6> q{};
    for (std::size_t step = 0; step != plan.refresh.size(); ++step) {
        std::array<std::array<std::uint8_t, 16>, 6> r{}, sum{};
        const auto& m = plan.refresh[step].entries;
        for (unsigned a = 0; a != 6; ++a)
            for (unsigned b = 0; b != 6; ++b)
                for (unsigned component = 0; component != bytes; ++component)
                    r[a][component] ^= sfMultiply(m[8 * a + b], q[b][component]);
        for (unsigned j = 0; j != 16; ++j) {
            const auto i = step * 16 + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto x = sfMultiply(planes[std::size_t(plan.route[i]) * bytes + component], plan.routedDiagonal[i]);
                auto y = x;
                for (unsigned a = 0; a != 6; ++a) {
                    const auto coefficient = sf16Embed(sf16Monomial(a, j));
                    y ^= sfMultiply(coefficient, r[a][component]);
                    sum[a][component] ^= sfMultiply(coefficient, x);
                }
                output[i * bytes + component] = y;
            }
        }
        for (unsigned a = 0; a != 6; ++a)
            for (unsigned component = 0; component != bytes; ++component)
                q[a][component] = r[a][component] ^ sum[a][component];
    }
    return output;
}

inline std::vector<std::uint8_t> subfieldWideTransposeOracle(const SubfieldWidePlan& plan,
    const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.n * bytes)
        throw std::invalid_argument("GF16 s6 transpose oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.k * bytes);
    std::array<std::array<std::uint8_t, 16>, 6> r{};
    for (std::size_t step = plan.refresh.size(); step--;) {
        std::array<std::array<std::uint8_t, 16>, 6> sum{}, previous{};
        for (unsigned j = 0; j != 16; ++j) {
            const auto i = step * 16 + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto v = input[i * bytes + component];
                auto u = v;
                for (unsigned a = 0; a != 6; ++a) {
                    const auto coefficient = sf16Embed(sf16Monomial(a, j));
                    u ^= sfMultiply(coefficient, r[a][component]);
                    sum[a][component] ^= sfMultiply(coefficient, v);
                }
                planes[std::size_t(plan.route[i]) * bytes + component] = sfMultiply(u, plan.routedDiagonal[i]);
            }
        }
        const auto& m = plan.refresh[step].entries;
        for (unsigned a = 0; a != 6; ++a)
            for (unsigned b = 0; b != 6; ++b)
                for (unsigned component = 0; component != bytes; ++component)
                    previous[a][component] ^= sfMultiply(m[8 * b + a], r[b][component] ^ sum[b][component]);
        r = previous;
    }
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned j = 0; j != 8; ++j)
            for (unsigned coordinate = 0; coordinate != 16; ++coordinate)
                for (unsigned component = 0; component != bytes; ++component)
                    output[(group * 8 + j) * bytes + component] ^=
                        sfMultiply(planes[(coordinate * plan.groups + group) * bytes + component], sf16Embed(sf16Basis(j, coordinate)));
    return output;
}

struct SubfieldWideState {
    __m128i r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128(), r2 = _mm_setzero_si128();
    __m128i r3 = _mm_setzero_si128(), r4 = _mm_setzero_si128(), r5 = _mm_setzero_si128();
};

enum class SubfieldWideMapKernel : std::uint8_t { Direct, AdditiveFft };

inline constexpr auto subfieldWideFftSubspaces = sf16Pack4({sf16Subspace(2, 0), sf16Subspace(2, 4),
                                                          sf16Subspace(2, 8), sf16Subspace(2, 12)});

template<unsigned Batch, bool Scatter>
SPIN_WIDE_INLINE void subfieldWideFftOutput(const SubfieldWidePlan& plan, Block* destination,
    __m512i value, __m512i diagonal, std::size_t begin, std::size_t outputBegin) {
    const auto i = begin + 4 * Batch;
    const auto u = _mm512_gf2p8mul_epi8(value, sf16ExpandRow<Batch>(diagonal));
    if constexpr (Scatter) {
        destination[plan.route[i]] = block(_mm512_castsi512_si128(u));
        destination[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
        destination[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
        destination[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
    } else _mm512_storeu_si512(destination + outputBegin + 4 * Batch, u);
}

// Exact factorization of the same monomial A and A^T. In abstract GF16,
// S2(X)=X^4+7X^2+6X. Use additive basis
// (1,X,X^2+X,X^3+X^2,S2(X),X*S2(X)). Its last two terms reduce to constants
// times (1,X) on each coset a+{0,1,2,3}, so only two extra products per map
// are needed beyond the four-state FFT. All constants below are embedded.
// Internal outputs 'moments' and input 'state' must be distinct objects.
template<bool Scatter>
SPIN_WIDE_INLINE void subfieldWideInnerFft(const SubfieldWidePlan& plan, const Block* input,
    Block* destination, std::size_t begin, std::size_t outputBegin,
    const SubfieldWideState& state, SubfieldWideState& moments) {
    auto v0 = _mm512_loadu_si512(input + begin), v1 = _mm512_loadu_si512(input + begin + 4);
    auto v2 = _mm512_loadu_si512(input + begin + 8), v3 = _mm512_loadu_si512(input + begin + 12);
    sf16TransposeRecords4(v0, v1, v2, v3);
    const auto a = _mm512_load_si512(sf16FftOffsets.bytes.data());
    const auto a2 = _mm512_load_si512(sf16FftSecondOffsets.bytes.data());
    const auto s1a = _mm512_load_si512(sf16FftSubspaces.bytes.data());
    const auto s2a = _mm512_load_si512(subfieldWideFftSubspaces.bytes.data());
    const auto six = _mm512_set1_epi8(static_cast<char>(sf16Embed(6)));
    const auto seven = _mm_set1_epi8(static_cast<char>(sf16Embed(7)));

    const auto v01 = _mm512_xor_si512(v0, v1), v23 = _mm512_xor_si512(v2, v3);
    const auto w1 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(a, v01), v1);
    const auto w3 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(a2, v23), v3);
    const auto b0 = _mm512_xor_si512(v01, v23), b1 = _mm512_xor_si512(w1, w3);
    const auto b2 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(s1a, b0), _mm512_gf2p8mul_epi8(six, v23));
    const auto b3 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(s1a, b1), _mm512_gf2p8mul_epi8(six, w3));
    const auto b4 = _mm512_gf2p8mul_epi8(s2a, b0), b5 = _mm512_gf2p8mul_epi8(s2a, b1);
    const auto m1 = sfSumRecords4(b1), m2 = sfSumRecords4(b2), m3 = sfSumRecords4(b3);
    moments.r0 = sfSumRecords4(b0);
    moments.r1 = m1;
    moments.r2 = _mm_xor_si128(m1, m2);
    moments.r3 = _mm_xor_si128(moments.r2, m3);
    moments.r4 = _mm_xor_si128(m1, _mm_xor_si128(_mm_gf2p8mul_epi8(m2, seven), sfSumRecords4(b4)));
    moments.r5 = _mm_xor_si128(moments.r2, _mm_xor_si128(_mm_gf2p8mul_epi8(m3, seven), sfSumRecords4(b5)));

    // Monomial-to-additive conversion:
    // b1=r1+r2+r3+r4+r5, b2=r2+r3+7*r4+r5, b3=r3+7*r5.
    const auto r23 = _mm_xor_si128(state.r2, state.r3);
    const auto basis1 = _mm_xor_si128(_mm_xor_si128(state.r1, r23), _mm_xor_si128(state.r4, state.r5));
    const auto basis2 = _mm_xor_si128(_mm_xor_si128(r23, state.r5), _mm_gf2p8mul_epi8(state.r4, seven));
    const auto basis3 = _mm_xor_si128(state.r3, _mm_gf2p8mul_epi8(state.r5, seven));
    const auto c0 = _mm512_xor_si512(_mm512_broadcast_i32x4(state.r0),
        _mm512_gf2p8mul_epi8(s2a, _mm512_broadcast_i32x4(state.r4)));
    const auto c1 = _mm512_xor_si512(_mm512_broadcast_i32x4(basis1),
        _mm512_gf2p8mul_epi8(s2a, _mm512_broadcast_i32x4(state.r5)));
    const auto c2 = _mm512_broadcast_i32x4(basis2), c3 = _mm512_broadcast_i32x4(basis3);
    const auto u0 = _mm512_xor_si512(c0, _mm512_gf2p8mul_epi8(s1a, c2));
    const auto u1 = _mm512_xor_si512(c1, _mm512_gf2p8mul_epi8(s1a, c3));
    const auto u2 = _mm512_xor_si512(u0, _mm512_gf2p8mul_epi8(six, c2));
    const auto u3 = _mm512_xor_si512(u1, _mm512_gf2p8mul_epi8(six, c3));
    const auto e0 = _mm512_xor_si512(u0, _mm512_gf2p8mul_epi8(a, u1));
    const auto e2 = _mm512_xor_si512(u2, _mm512_gf2p8mul_epi8(a2, u3));
    v0 = _mm512_xor_si512(v0, e0);
    v1 = _mm512_xor_si512(v1, _mm512_xor_si512(e0, u1));
    v2 = _mm512_xor_si512(v2, e2);
    v3 = _mm512_xor_si512(v3, _mm512_xor_si512(e2, u3));
    sf16TransposeRecords4(v0, v1, v2, v3);
    const auto diagonal = _mm512_broadcast_i32x4(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(plan.routedDiagonal.data() + begin)));
    subfieldWideFftOutput<0, Scatter>(plan, destination, v0, diagonal, begin, outputBegin);
    subfieldWideFftOutput<1, Scatter>(plan, destination, v1, diagonal, begin, outputBegin);
    subfieldWideFftOutput<2, Scatter>(plan, destination, v2, diagonal, begin, outputBegin);
    subfieldWideFftOutput<3, Scatter>(plan, destination, v3, diagonal, begin, outputBegin);
}

template<unsigned Batch, bool Scatter>
SPIN_WIDE_INLINE void subfieldWideInnerBatch(const SubfieldWidePlan& plan, const Block* input, Block* destination,
    std::size_t begin, std::size_t outputBegin, __m512i diagonal,
    __m512i r0, __m512i r1, __m512i r2, __m512i r3, __m512i r4, __m512i r5,
    __m512i& s0, __m512i& s1, __m512i& s2, __m512i& s3, __m512i& s4, __m512i& s5) {
    const auto i = begin + 4 * Batch;
    const auto a1 = _mm512_load_si512(plan.coefficients[Batch][0].bytes.data());
    const auto a2 = _mm512_load_si512(plan.coefficients[Batch][1].bytes.data());
    const auto a3 = _mm512_load_si512(plan.coefficients[Batch][2].bytes.data());
    const auto a4 = _mm512_load_si512(plan.coefficients[Batch][3].bytes.data());
    const auto a5 = _mm512_load_si512(plan.coefficients[Batch][4].bytes.data());
    const auto v = _mm512_loadu_si512(input + i);
    s0 = _mm512_xor_si512(s0, v);
    s1 = _mm512_xor_si512(s1, _mm512_gf2p8mul_epi8(v, a1));
    s2 = _mm512_xor_si512(s2, _mm512_gf2p8mul_epi8(v, a2));
    s3 = _mm512_xor_si512(s3, _mm512_gf2p8mul_epi8(v, a3));
    s4 = _mm512_xor_si512(s4, _mm512_gf2p8mul_epi8(v, a4));
    s5 = _mm512_xor_si512(s5, _mm512_gf2p8mul_epi8(v, a5));
    const auto p01 = _mm512_xor_si512(_mm512_xor_si512(v, r0), _mm512_gf2p8mul_epi8(r1, a1));
    const auto p23 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(r2, a2), _mm512_gf2p8mul_epi8(r3, a3));
    const auto p45 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(r4, a4), _mm512_gf2p8mul_epi8(r5, a5));
    const auto u = _mm512_gf2p8mul_epi8(_mm512_xor_si512(p01, _mm512_xor_si512(p23, p45)), sf16ExpandRow<Batch>(diagonal));
    if constexpr (Scatter) {
        destination[plan.route[i]] = block(_mm512_castsi512_si128(u));
        destination[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
        destination[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
        destination[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
    } else _mm512_storeu_si512(destination + outputBegin + 4 * Batch, u);
}

template<unsigned Row>
SPIN_WIDE_INLINE void subfieldWideRefreshAdd(const SubfieldWideMatrix& matrix, __m128i source,
                                            __m512i& low, __m256i& high) {
    const auto packed = _mm512_broadcast_i32x4(_mm_loadl_epi64(
        reinterpret_cast<const __m128i*>(matrix.entries.data() + 8 * Row)));
    const auto x = _mm512_broadcast_i32x4(source);
    low = _mm512_xor_si512(low, _mm512_gf2p8mul_epi8(x, sf16ExpandRow<0>(packed)));
    const auto tail = _mm256_shuffle_epi8(_mm512_castsi512_si256(packed), _mm256_setr_epi64x(
        0x0404040404040404LL, 0x0404040404040404LL, 0x0505050505050505LL, 0x0505050505050505LL));
    high = _mm256_xor_si256(high, _mm256_gf2p8mul_epi8(_mm512_castsi512_si256(x), tail));
}

template<bool Scatter, bool Fft = false>
inline void subfieldWideInnerRange(const SubfieldWidePlan& plan, const Block* input, Block* destination,
                                   std::size_t first, std::size_t past, std::size_t outputOffset,
                                   SubfieldWideState& state) {
    for (std::size_t step = past; step-- > first;) {
        const auto begin = step * 16, outputBegin = begin - outputOffset;
        __m128i u0, u1, u2, u3, u4, u5;
        if constexpr (Fft) {
            SubfieldWideState moments;
            subfieldWideInnerFft<Scatter>(plan, input, destination, begin, outputBegin, state, moments);
            u0 = _mm_xor_si128(state.r0, moments.r0);
            u1 = _mm_xor_si128(state.r1, moments.r1);
            u2 = _mm_xor_si128(state.r2, moments.r2);
            u3 = _mm_xor_si128(state.r3, moments.r3);
            u4 = _mm_xor_si128(state.r4, moments.r4);
            u5 = _mm_xor_si128(state.r5, moments.r5);
        } else {
            const auto r0 = _mm512_broadcast_i32x4(state.r0), r1 = _mm512_broadcast_i32x4(state.r1);
            const auto r2 = _mm512_broadcast_i32x4(state.r2), r3 = _mm512_broadcast_i32x4(state.r3);
            const auto r4 = _mm512_broadcast_i32x4(state.r4), r5 = _mm512_broadcast_i32x4(state.r5);
            auto s0 = _mm512_setzero_si512(), s1 = _mm512_setzero_si512(), s2 = _mm512_setzero_si512();
            auto s3 = _mm512_setzero_si512(), s4 = _mm512_setzero_si512(), s5 = _mm512_setzero_si512();
            const auto diagonal = _mm512_broadcast_i32x4(_mm_loadu_si128(
                reinterpret_cast<const __m128i*>(plan.routedDiagonal.data() + begin)));
            subfieldWideInnerBatch<0, Scatter>(plan, input, destination, begin, outputBegin, diagonal, r0,r1,r2,r3,r4,r5,s0,s1,s2,s3,s4,s5);
            subfieldWideInnerBatch<1, Scatter>(plan, input, destination, begin, outputBegin, diagonal, r0,r1,r2,r3,r4,r5,s0,s1,s2,s3,s4,s5);
            subfieldWideInnerBatch<2, Scatter>(plan, input, destination, begin, outputBegin, diagonal, r0,r1,r2,r3,r4,r5,s0,s1,s2,s3,s4,s5);
            subfieldWideInnerBatch<3, Scatter>(plan, input, destination, begin, outputBegin, diagonal, r0,r1,r2,r3,r4,r5,s0,s1,s2,s3,s4,s5);
            u0 = _mm_xor_si128(state.r0, sfSumRecords4(s0));
            u1 = _mm_xor_si128(state.r1, sfSumRecords4(s1));
            u2 = _mm_xor_si128(state.r2, sfSumRecords4(s2));
            u3 = _mm_xor_si128(state.r3, sfSumRecords4(s3));
            u4 = _mm_xor_si128(state.r4, sfSumRecords4(s4));
            u5 = _mm_xor_si128(state.r5, sfSumRecords4(s5));
        }
        auto low = _mm512_setzero_si512();
        auto high = _mm256_setzero_si256();
        const auto& matrix = plan.refresh[step];
        subfieldWideRefreshAdd<0>(matrix, u0, low, high);
        subfieldWideRefreshAdd<1>(matrix, u1, low, high);
        subfieldWideRefreshAdd<2>(matrix, u2, low, high);
        subfieldWideRefreshAdd<3>(matrix, u3, low, high);
        subfieldWideRefreshAdd<4>(matrix, u4, low, high);
        subfieldWideRefreshAdd<5>(matrix, u5, low, high);
        state.r0 = _mm512_castsi512_si128(low);
        state.r1 = _mm512_extracti32x4_epi32(low, 1);
        state.r2 = _mm512_extracti32x4_epi32(low, 2);
        state.r3 = _mm512_extracti32x4_epi32(low, 3);
        state.r4 = _mm256_castsi256_si128(high);
        state.r5 = _mm256_extracti128_si256(high, 1);
    }
}

template<bool Stream, bool Fft = false>
inline void subfieldWideRegionInner(const SubfieldWidePlan& plan, const Block* input,
                                    Block* planes, Block* temporary) {
    SubfieldWideState state;
    for (unsigned region = 16; region--;) {
        const auto offset = region * plan.groups;
        subfieldWideInnerRange<false, Fft>(plan, input, temporary, offset / 16,
            (offset + plan.groups) / 16, offset, state);
        for (std::size_t j = 0; j != plan.groups; j += 4) {
            const auto value = sf16GatherRegion4(temporary, plan.inverseRoute.data() + offset + j,
                                                static_cast<std::uint32_t>(offset));
            if constexpr (Stream) _mm512_stream_si512(reinterpret_cast<__m512i*>(planes + offset + j), value);
            else _mm512_storeu_si512(planes + offset + j, value);
        }
    }
    if constexpr (Stream) _mm_sfence();
}

// Same four-group outer circuit and layout as Subfield16. Its helpers are
// reused directly; this wrapper accepts the independent six-state setup type.
template<bool Stream>
inline void subfieldWideOuterTranspose(const SubfieldWidePlan& plan, const Block* planes, Block* output) {
    std::size_t group = 0;
    for (; group + 4 <= plan.groups; group += 4) {
        __m512i values[16];
        for (unsigned coordinate = 0; coordinate != 16; ++coordinate)
            values[coordinate] = _mm512_loadu_si512(planes + coordinate * plan.groups + group);
        sf16EvalTranspose4<8, 0>(values);
        sf16EvalTranspose4<8, 8>(values + 8);
        for (unsigned j = 0; j != 8; ++j) values[j] = _mm512_xor_si512(values[j], values[j + 8]);
        for (unsigned j = 0; j != 8; j += 4) {
            const auto a01 = _mm512_shuffle_i32x4(values[j], values[j + 1], 0x44);
            const auto b01 = _mm512_shuffle_i32x4(values[j + 2], values[j + 3], 0x44);
            const auto a23 = _mm512_shuffle_i32x4(values[j], values[j + 1], 0xee);
            const auto b23 = _mm512_shuffle_i32x4(values[j + 2], values[j + 3], 0xee);
            sf16Store<Stream>(output + 8 * group + j, _mm512_shuffle_i32x4(a01, b01, 0x88));
            sf16Store<Stream>(output + 8 * (group + 1) + j, _mm512_shuffle_i32x4(a01, b01, 0xdd));
            sf16Store<Stream>(output + 8 * (group + 2) + j, _mm512_shuffle_i32x4(a23, b23, 0x88));
            sf16Store<Stream>(output + 8 * (group + 3) + j, _mm512_shuffle_i32x4(a23, b23, 0xdd));
        }
    }
    for (; group != plan.groups; ++group)
        for (unsigned j = 0; j != 8; ++j) {
            auto result = _mm_setzero_si128();
            for (unsigned coordinate = 0; coordinate != 16; ++coordinate)
                result = _mm_xor_si128(result, _mm_gf2p8mul_epi8(reg(planes[coordinate * plan.groups + group]),
                    _mm_set1_epi8(static_cast<char>(sf16Embed(sf16Basis(j, coordinate))))));
            output[8 * group + j] = block(result);
        }
    if constexpr (Stream) _mm_sfence();
}

inline std::size_t subfieldWideTemporaryRecords(const SubfieldWidePlan& plan, Subfield16Kernel kernel) noexcept {
    return kernel != Subfield16Kernel::Scatter && plan.regionAligned() ? plan.groups : 0;
}

template<bool Fft>
inline void subfieldWideInnerMapDispatch(const SubfieldWidePlan& plan, const Block* input, Block* planes,
                                         Block* temporary, Subfield16Kernel kernel) {
    if (kernel == Subfield16Kernel::Scatter || !plan.regionAligned()) {
        SubfieldWideState state;
        subfieldWideInnerRange<true, Fft>(plan, input, planes, 0, plan.n / 16, 0, state);
    } else {
        if (!temporary) throw std::invalid_argument("GF16 s6 region schedule needs L temporary records");
        if ((reinterpret_cast<std::uintptr_t>(planes) & 63U) == 0)
            subfieldWideRegionInner<true, Fft>(plan, input, planes, temporary);
        else subfieldWideRegionInner<false, Fft>(plan, input, planes, temporary);
    }
}

inline void subfieldWideInnerDispatch(const SubfieldWidePlan& plan, const Block* input, Block* planes,
                                      Block* temporary, Subfield16Kernel kernel,
                                      SubfieldWideMapKernel map = SubfieldWideMapKernel::Direct) {
    if (map == SubfieldWideMapKernel::AdditiveFft)
        subfieldWideInnerMapDispatch<true>(plan, input, planes, temporary, kernel);
    else subfieldWideInnerMapDispatch<false>(plan, input, planes, temporary, kernel);
}

inline void subfieldWideOuterDispatch(const SubfieldWidePlan& plan, const Block* planes, Block* output,
                                      Subfield16Kernel kernel) {
    if (kernel == Subfield16Kernel::RegionAllNonTemporal &&
        (reinterpret_cast<std::uintptr_t>(output) & 63U) == 0)
        subfieldWideOuterTranspose<true>(plan, planes, output);
    else subfieldWideOuterTranspose<false>(plan, planes, output);
}

// In-place input==output changes only [0,k). Both scratch areas are caller
// owned, disjoint from input/output and one another; no allocation in encode.
inline void subfieldWideTransposeFast128(const SubfieldWidePlan& plan, const Block* input, Block* output,
                                        Block* planes, Block* temporary, Subfield16Kernel kernel,
                                        SubfieldWideMapKernel map = SubfieldWideMapKernel::Direct) {
    subfieldWideInnerDispatch(plan, input, planes, temporary, kernel, map);
    subfieldWideOuterDispatch(plan, planes, output, kernel);
}

} // namespace spin::research::char2

#undef SPIN_WIDE_INLINE
