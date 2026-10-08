#pragma once
#include "Subfield.h"

namespace spin::research::char2 {

// Separate research family: GF16 = GF2[X]/(X^4+X+1), embedded into AES GF256.
// Every stored scalar below is the image of an abstract nibble, not that nibble
// interpreted directly as an AES-field element. Records use the GF256 basis
// documented by Subfield.h; these subfield scalars therefore act linearly over
// the full GF(2^m) extension, m in {8,32,64,128}. The ideal setup ensemble has
// a K=2^16, 10% symbol-distance certificate with >42 bits of failure margin:
// research/workstreams/char2_field_spin/theory/gf16_certificate.py.
// The deterministic benchmark PRNG is not itself the ideal setup sampler.
constexpr std::uint8_t sf16Multiply(std::uint8_t a, std::uint8_t b) noexcept {
    std::uint8_t result = 0;
    for (unsigned i = 0; i != 4; ++i) {
        if (b & 1U) result ^= a;
        a = std::uint8_t((unsigned(a) << 1) ^ ((a & 8U) ? 0x13U : 0U));
        b >>= 1;
    }
    return result;
}

constexpr std::uint8_t sf16EmbeddingRoot() noexcept {
    for (unsigned value = 2; value != 256; ++value) {
        const auto a = std::uint8_t(value), a2 = sfMultiply(a, a);
        if ((sfMultiply(a2, a2) ^ a ^ 1U) == 0) return a;
    }
    return 0;
}

inline constexpr auto sf16Beta = sf16EmbeddingRoot();
static_assert(sf16Beta != 0);

constexpr std::uint8_t sf16Embed(std::uint8_t nibble) noexcept {
    std::uint8_t result = 0, power = 1;
    for (unsigned i = 0; i != 4; ++i) {
        if (nibble & (1U << i)) result ^= power;
        power = sfMultiply(power, sf16Beta);
    }
    return result;
}

constexpr std::uint8_t sf16Inverse(std::uint8_t a) noexcept {
    const auto a2 = sf16Multiply(a, a), a4 = sf16Multiply(a2, a2);
    return sf16Multiply(sf16Multiply(a4, a4), sf16Multiply(a4, a2));
}

constexpr std::uint8_t sf16Subspace(unsigned degree, unsigned point) noexcept {
    std::uint8_t result = 1;
    for (unsigned v = 0; v != (1U << degree); ++v)
        result = sf16Multiply(result, std::uint8_t(point ^ v));
    return result;
}

constexpr std::uint8_t sf16Basis(unsigned index, unsigned point) noexcept {
    std::uint8_t result = 1;
    for (unsigned degree = 0; index; ++degree, index >>= 1)
        if (index & 1U) result = sf16Multiply(result, sf16Subspace(degree, point));
    return result;
}

constexpr std::uint8_t sf16Monomial(unsigned degree, unsigned point) noexcept {
    std::uint8_t result = 1;
    for (unsigned i = 0; i != degree; ++i) result = sf16Multiply(result, std::uint8_t(point));
    return result;
}

struct Subfield16Matrix {
    std::array<std::uint8_t, 16> entries; // Row-major, embedded AES-field bytes.
};

inline bool sf16Invertible(std::array<std::uint8_t, 16> matrix) noexcept {
    for (unsigned column = 0; column != 4; ++column) {
        unsigned pivot = column;
        while (pivot != 4 && matrix[4 * pivot + column] == 0) ++pivot;
        if (pivot == 4) return false;
        if (pivot != column)
            for (unsigned j = column; j != 4; ++j)
                std::swap(matrix[4 * pivot + j], matrix[4 * column + j]);
        const auto inverse = sf16Inverse(matrix[5 * column]);
        for (unsigned row = column + 1; row != 4; ++row) {
            const auto factor = sf16Multiply(matrix[4 * row + column], inverse);
            for (unsigned j = column; j != 4; ++j)
                matrix[4 * row + j] ^= sf16Multiply(factor, matrix[4 * column + j]);
        }
    }
    return true;
}

template<class Generator>
inline Subfield16Matrix sf16SampleMatrix(Generator& generator) {
    std::array<std::uint8_t, 16> abstract;
    do {
        auto packed = std::uint64_t(generator());
        for (auto& entry : abstract) { entry = std::uint8_t(packed & 15U); packed >>= 4; }
    } while (!sf16Invertible(abstract));
    Subfield16Matrix matrix;
    for (unsigned i = 0; i != 16; ++i) matrix.entries[i] = sf16Embed(abstract[i]);
    return matrix;
}

enum class Subfield16Kernel : std::uint8_t { Scatter, RegionNonTemporal, RegionAllNonTemporal };
enum class Subfield16MapKernel : std::uint8_t { Direct, AdditiveFft, PackedDirect };

// GRS[16,8] outer; independent nonzero GF16 diagonal; fixed outer-coordinate
// regions independently permuted over individual groups; width-16 inner with
// A_j=(1,j,j^2,j^3), C=A^T, and a fresh uniform GL4 refresh before each emission.
// Setup uses unbiased rejection/Fisher-Yates given uniform generator words.
// mt19937_64 is a reproducible fixture, not a production setup randomness source.
struct Subfield16Plan {
    static constexpr unsigned width = 16, outer = 16, dimension = 8, state = 4;
    std::size_t k, n, groups;
    std::uint64_t seed;
    std::vector<std::uint32_t> route, inverseRoute;
    std::vector<std::uint8_t> routedDiagonal;
    std::vector<Subfield16Matrix> refresh;
    std::array<std::array<SubfieldIndexCoefficients, 3>, 4> coefficients;

    Subfield16Plan(std::size_t message, std::uint64_t setupSeed)
        : k(message), n(2 * message), groups(message / dimension), seed(setupSeed) {
        if (!k || k % dimension || k > std::numeric_limits<std::uint32_t>::max() / 2)
            throw std::invalid_argument("GF16 candidate needs nonzero K divisible by 8 within uint32 route limits");
        std::mt19937_64 routeRng(seed ^ 0xdf4a1b3285ec7096ULL);
        std::mt19937_64 diagonalRng(seed ^ 0x782d043af9165bc1ULL);
        std::mt19937_64 matrixRng(seed ^ 0x621bdfea9037c845ULL);
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
        for (auto& matrix : refresh) matrix = sf16SampleMatrix(matrixRng);
        for (unsigned batch = 0; batch != 4; ++batch)
            for (unsigned power = 1; power != 4; ++power)
                for (unsigned lane = 0; lane != 4; ++lane) {
                    const auto repeated = 0x0101010101010101ULL * sf16Embed(sf16Monomial(power, 4 * batch + lane));
                    coefficients[batch][power - 1].bytes[2 * lane] = repeated;
                    coefficients[batch][power - 1].bytes[2 * lane + 1] = repeated;
                }
    }

    bool regionAligned() const noexcept { return groups % width == 0; }
    std::size_t setupBytes() const noexcept {
        return (route.size() + inverseRoute.size()) * sizeof(std::uint32_t) + routedDiagonal.size() +
               refresh.size() * sizeof(Subfield16Matrix) + sizeof(coefficients);
    }
};

inline std::vector<std::uint8_t> subfield16ForwardOracle(const Subfield16Plan& plan,
    const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.k * bytes)
        throw std::invalid_argument("GF16 forward oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.n * bytes);
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned coordinate = 0; coordinate != 16; ++coordinate)
            for (unsigned j = 0; j != 8; ++j)
                for (unsigned component = 0; component != bytes; ++component)
                    planes[(coordinate * plan.groups + group) * bytes + component] ^=
                        sfMultiply(input[(group * 8 + j) * bytes + component], sf16Embed(sf16Basis(j, coordinate)));
    std::array<std::array<std::uint8_t, 16>, 4> q{};
    for (std::size_t step = 0; step != plan.refresh.size(); ++step) {
        std::array<std::array<std::uint8_t, 16>, 4> r{}, sum{};
        const auto& m = plan.refresh[step].entries;
        for (unsigned a = 0; a != 4; ++a)
            for (unsigned b = 0; b != 4; ++b)
                for (unsigned component = 0; component != bytes; ++component)
                    r[a][component] ^= sfMultiply(m[4 * a + b], q[b][component]);
        for (unsigned j = 0; j != 16; ++j) {
            const auto i = step * 16 + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto x = sfMultiply(planes[std::size_t(plan.route[i]) * bytes + component], plan.routedDiagonal[i]);
                auto y = x;
                for (unsigned a = 0; a != 4; ++a) {
                    const auto coefficient = sf16Embed(sf16Monomial(a, j));
                    y ^= sfMultiply(coefficient, r[a][component]);
                    sum[a][component] ^= sfMultiply(coefficient, x);
                }
                output[i * bytes + component] = y;
            }
        }
        for (unsigned a = 0; a != 4; ++a)
            for (unsigned component = 0; component != bytes; ++component)
                q[a][component] = r[a][component] ^ sum[a][component];
    }
    return output;
}

inline std::vector<std::uint8_t> subfield16TransposeOracle(const Subfield16Plan& plan,
    const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.n * bytes)
        throw std::invalid_argument("GF16 transpose oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.k * bytes);
    std::array<std::array<std::uint8_t, 16>, 4> r{};
    for (std::size_t step = plan.refresh.size(); step--;) {
        std::array<std::array<std::uint8_t, 16>, 4> sum{}, previous{};
        for (unsigned j = 0; j != 16; ++j) {
            const auto i = step * 16 + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto v = input[i * bytes + component];
                auto u = v;
                for (unsigned a = 0; a != 4; ++a) {
                    const auto coefficient = sf16Embed(sf16Monomial(a, j));
                    u ^= sfMultiply(coefficient, r[a][component]);
                    sum[a][component] ^= sfMultiply(coefficient, v);
                }
                planes[std::size_t(plan.route[i]) * bytes + component] = sfMultiply(u, plan.routedDiagonal[i]);
            }
        }
        const auto& m = plan.refresh[step].entries;
        for (unsigned a = 0; a != 4; ++a)
            for (unsigned b = 0; b != 4; ++b)
                for (unsigned component = 0; component != bytes; ++component)
                    previous[a][component] ^= sfMultiply(m[4 * b + a], r[b][component] ^ sum[b][component]);
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

// One compact 16-byte load supplies all four coefficient rows. Each 128-bit
// lane receives its own selected byte, with no expanded setup table.
template<unsigned Row>
inline __m512i sf16ExpandRow(__m512i packed) noexcept {
    constexpr auto repeated = [](unsigned x) { return std::int64_t(x * 0x0101010101010101ULL); };
    const auto select = _mm512_setr_epi64(repeated(4 * Row), repeated(4 * Row),
        repeated(4 * Row + 1), repeated(4 * Row + 1),
        repeated(4 * Row + 2), repeated(4 * Row + 2),
        repeated(4 * Row + 3), repeated(4 * Row + 3));
    return _mm512_shuffle_epi8(packed, select);
}

template<unsigned Batch, bool Scatter, bool Packed = false>
inline void sf16InnerBatch(const Subfield16Plan& plan, const Block* input, Block* destination,
                           std::size_t begin, std::size_t outputBegin,
                           __m512i r0, __m512i r1, __m512i r2, __m512i r3,
                           __m512i& sum0, __m512i& sum1, __m512i& sum2, __m512i& sum3,
                           __m512i packedDiagonal = _mm512_setzero_si512()) {
    const auto i = begin + 4 * Batch;
    const auto a1 = _mm512_load_si512(plan.coefficients[Batch][0].bytes.data());
    const auto a2 = _mm512_load_si512(plan.coefficients[Batch][1].bytes.data());
    const auto a3 = _mm512_load_si512(plan.coefficients[Batch][2].bytes.data());
    const auto v = _mm512_loadu_si512(input + i);
    sum0 = _mm512_xor_si512(sum0, v);
    sum1 = _mm512_xor_si512(sum1, _mm512_gf2p8mul_epi8(v, a1));
    sum2 = _mm512_xor_si512(sum2, _mm512_gf2p8mul_epi8(v, a2));
    sum3 = _mm512_xor_si512(sum3, _mm512_gf2p8mul_epi8(v, a3));
    const auto low = _mm512_xor_si512(v, _mm512_xor_si512(r0, _mm512_gf2p8mul_epi8(r1, a1)));
    const auto high = _mm512_xor_si512(_mm512_gf2p8mul_epi8(r2, a2), _mm512_gf2p8mul_epi8(r3, a3));
    __m512i diagonal;
    if constexpr (Packed) diagonal = sf16ExpandRow<Batch>(packedDiagonal);
    else diagonal = sfBroadcastFour(plan.routedDiagonal.data() + i);
    const auto u = _mm512_gf2p8mul_epi8(_mm512_xor_si512(low, high), diagonal);
    if constexpr (Scatter) {
        destination[plan.route[i]] = block(_mm512_castsi512_si128(u));
        destination[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
        destination[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
        destination[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
    } else _mm512_storeu_si512(destination + outputBegin + 4 * Batch, u);
}

struct Subfield16State {
    __m128i r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128();
    __m128i r2 = _mm_setzero_si128(), r3 = _mm_setzero_si128();
};

constexpr SubfieldIndexCoefficients sf16Pack4(std::array<std::uint8_t, 4> values) {
    SubfieldIndexCoefficients result{};
    for (unsigned lane = 0; lane != 4; ++lane) {
        const auto repeated = 0x0101010101010101ULL * sf16Embed(values[lane]);
        result.bytes[2 * lane] = result.bytes[2 * lane + 1] = repeated;
    }
    return result;
}

inline constexpr auto sf16FftOffsets = sf16Pack4({0, 4, 8, 12});
inline constexpr auto sf16FftSecondOffsets = sf16Pack4({2, 6, 10, 14});
inline constexpr auto sf16FftSubspaces = sf16Pack4({sf16Subspace(1, 0), sf16Subspace(1, 4),
                                                 sf16Subspace(1, 8), sf16Subspace(1, 12)});

inline void sf16TransposeRecords4(__m512i& a, __m512i& b, __m512i& c, __m512i& d) {
    const auto a01 = _mm512_shuffle_i32x4(a, b, 0x44);
    const auto b01 = _mm512_shuffle_i32x4(c, d, 0x44);
    const auto a23 = _mm512_shuffle_i32x4(a, b, 0xee);
    const auto b23 = _mm512_shuffle_i32x4(c, d, 0xee);
    a = _mm512_shuffle_i32x4(a01, b01, 0x88);
    b = _mm512_shuffle_i32x4(a01, b01, 0xdd);
    c = _mm512_shuffle_i32x4(a23, b23, 0x88);
    d = _mm512_shuffle_i32x4(a23, b23, 0xdd);
}

template<unsigned Batch, bool Scatter>
inline void sf16FftOutput(const Subfield16Plan& plan, Block* destination, __m512i value,
                          std::size_t begin, std::size_t outputBegin) {
    const auto i = begin + 4 * Batch;
    const auto u = _mm512_gf2p8mul_epi8(value, sfBroadcastFour(plan.routedDiagonal.data() + i));
    if constexpr (Scatter) {
        destination[plan.route[i]] = block(_mm512_castsi512_si128(u));
        destination[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
        destination[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
        destination[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
    } else _mm512_storeu_si512(destination + outputBegin + 4 * Batch, u);
}

// Same monomial A and C=A^T, factored through basis (1,X,X²+X,X³+X²).
// State-to-basis: b=(r0,r1+r2+r3,r2+r3,r3). The transpose conversion returns
// monomial moments (B0,B1,B1+B2,B1+B2+B3), so neither map changes its basis.
template<bool Scatter>
inline void sf16InnerFft(const Subfield16Plan& plan, const Block* input, Block* destination,
                         std::size_t begin, std::size_t outputBegin, const Subfield16State& state,
                         Subfield16State& moments) {
    auto v0 = _mm512_loadu_si512(input + begin), v1 = _mm512_loadu_si512(input + begin + 4);
    auto v2 = _mm512_loadu_si512(input + begin + 8), v3 = _mm512_loadu_si512(input + begin + 12);
    sf16TransposeRecords4(v0, v1, v2, v3);
    const auto a = _mm512_load_si512(sf16FftOffsets.bytes.data());
    const auto a2 = _mm512_load_si512(sf16FftSecondOffsets.bytes.data());
    const auto s1a = _mm512_load_si512(sf16FftSubspaces.bytes.data());
    const auto six = _mm512_set1_epi8(static_cast<char>(sf16Embed(6)));

    const auto v01 = _mm512_xor_si512(v0, v1), v23 = _mm512_xor_si512(v2, v3);
    const auto w1 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(a, v01), v1);
    const auto w3 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(a2, v23), v3);
    const auto b0 = _mm512_xor_si512(v01, v23), b1 = _mm512_xor_si512(w1, w3);
    const auto b2 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(s1a, b0), _mm512_gf2p8mul_epi8(six, v23));
    const auto b3 = _mm512_xor_si512(_mm512_gf2p8mul_epi8(s1a, b1), _mm512_gf2p8mul_epi8(six, w3));
    const auto moment0 = sfSumRecords4(b0), moment1 = sfSumRecords4(b1);
    const auto moment2 = _mm_xor_si128(moment1, sfSumRecords4(b2));
    const auto moment3 = _mm_xor_si128(moment2, sfSumRecords4(b3));

    const auto c0 = _mm512_broadcast_i32x4(state.r0);
    const auto c2 = _mm512_broadcast_i32x4(_mm_xor_si128(state.r2, state.r3));
    const auto c1 = _mm512_xor_si512(_mm512_broadcast_i32x4(state.r1), c2);
    const auto c3 = _mm512_broadcast_i32x4(state.r3);
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
    sf16FftOutput<0, Scatter>(plan, destination, v0, begin, outputBegin);
    sf16FftOutput<1, Scatter>(plan, destination, v1, begin, outputBegin);
    sf16FftOutput<2, Scatter>(plan, destination, v2, begin, outputBegin);
    sf16FftOutput<3, Scatter>(plan, destination, v3, begin, outputBegin);
    // Avoid a 64-byte aggregate return/copy: GCC with ASAN has emitted an
    // aligned ZMM load from a return temporary that has only 16-byte alignment.
    // Explicit component outputs retain the actual __m128i alignment contract.
    moments.r0 = moment0;
    moments.r1 = moment1;
    moments.r2 = moment2;
    moments.r3 = moment3;
}

template<bool Scatter, bool Fft = false, bool Packed = false>
inline void subfield16InnerRange(const Subfield16Plan& plan, const Block* input, Block* destination,
                                 std::size_t first, std::size_t past, std::size_t outputOffset,
                                 Subfield16State& state) {
    for (std::size_t step = past; step-- > first;) {
        const auto begin = step * 16, outputBegin = begin - outputOffset;
        Subfield16State moments;
        if constexpr (Fft) {
            sf16InnerFft<Scatter>(plan, input, destination, begin, outputBegin, state, moments);
        } else {
            const auto r0 = _mm512_broadcast_i32x4(state.r0), r1 = _mm512_broadcast_i32x4(state.r1);
            const auto r2 = _mm512_broadcast_i32x4(state.r2), r3 = _mm512_broadcast_i32x4(state.r3);
            auto s0 = _mm512_setzero_si512(), s1 = _mm512_setzero_si512();
            auto s2 = _mm512_setzero_si512(), s3 = _mm512_setzero_si512();
            auto packedDiagonal = _mm512_setzero_si512();
            if constexpr (Packed) packedDiagonal = _mm512_broadcast_i32x4(
                _mm_loadu_si128(reinterpret_cast<const __m128i*>(plan.routedDiagonal.data() + begin)));
            sf16InnerBatch<0, Scatter, Packed>(plan, input, destination, begin, outputBegin, r0, r1, r2, r3, s0, s1, s2, s3, packedDiagonal);
            sf16InnerBatch<1, Scatter, Packed>(plan, input, destination, begin, outputBegin, r0, r1, r2, r3, s0, s1, s2, s3, packedDiagonal);
            sf16InnerBatch<2, Scatter, Packed>(plan, input, destination, begin, outputBegin, r0, r1, r2, r3, s0, s1, s2, s3, packedDiagonal);
            sf16InnerBatch<3, Scatter, Packed>(plan, input, destination, begin, outputBegin, r0, r1, r2, r3, s0, s1, s2, s3, packedDiagonal);
            moments.r0 = sfSumRecords4(s0);
            moments.r1 = sfSumRecords4(s1);
            moments.r2 = sfSumRecords4(s2);
            moments.r3 = sfSumRecords4(s3);
        }
        const auto u0 = _mm_xor_si128(state.r0, moments.r0);
        const auto u1 = _mm_xor_si128(state.r1, moments.r1);
        const auto u2 = _mm_xor_si128(state.r2, moments.r2);
        const auto u3 = _mm_xor_si128(state.r3, moments.r3);
        const auto& m = plan.refresh[step].entries;
        __m512i c0, c1, c2, c3;
        if constexpr (Packed) {
            const auto packet = _mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(m.data())));
            c0 = sf16ExpandRow<0>(packet); c1 = sf16ExpandRow<1>(packet);
            c2 = sf16ExpandRow<2>(packet); c3 = sf16ExpandRow<3>(packet);
        } else {
            c0 = sfBroadcastFour(m.data()); c1 = sfBroadcastFour(m.data() + 4);
            c2 = sfBroadcastFour(m.data() + 8); c3 = sfBroadcastFour(m.data() + 12);
        }
        const auto p0 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u0), c0);
        const auto p1 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u1), c1);
        const auto p2 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u2), c2);
        const auto p3 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u3), c3);
        const auto previous = _mm512_xor_si512(_mm512_xor_si512(p0, p1), _mm512_xor_si512(p2, p3));
        state.r0 = _mm512_castsi512_si128(previous);
        state.r1 = _mm512_extracti32x4_epi32(previous, 1);
        state.r2 = _mm512_extracti32x4_epi32(previous, 2);
        state.r3 = _mm512_extracti32x4_epi32(previous, 3);
    }
}

inline __m512i sf16GatherRegion4(const Block* input, const std::uint32_t* index, std::uint32_t offset) {
    auto value = _mm512_castsi128_si512(reg(input[index[0] - offset]));
    value = _mm512_inserti32x4(value, reg(input[index[1] - offset]), 1);
    value = _mm512_inserti32x4(value, reg(input[index[2] - offset]), 2);
    return _mm512_inserti32x4(value, reg(input[index[3] - offset]), 3);
}

template<bool Stream, bool Fft = false, bool Packed = false>
inline void subfield16RegionInner(const Subfield16Plan& plan, const Block* input,
                                  Block* planes, Block* temporary) {
    Subfield16State state;
    for (unsigned region = 16; region--;) {
        const auto offset = region * plan.groups;
        subfield16InnerRange<false, Fft, Packed>(plan, input, temporary, offset / 16,
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

template<unsigned Count, unsigned Offset>
inline void sf16EvalTranspose4(__m512i* values) {
    if constexpr (Count > 1) {
        constexpr unsigned half = Count / 2, degree = std::countr_zero(half);
        sf16EvalTranspose4<half, Offset>(values);
        sf16EvalTranspose4<half, Offset ^ half>(values + half);
        for (unsigned j = 0; j != half; ++j) {
            const auto sum = _mm512_xor_si512(values[j], values[j + half]);
            values[j + half] = _mm512_xor_si512(sfConstant4<sf16Embed(sf16Subspace(degree, Offset))>(sum),
                sfConstant4<sf16Embed(sf16Subspace(degree, half))>(values[j + half]));
            values[j] = sum;
        }
    }
}

template<bool Stream>
inline void sf16Store(Block* output, __m512i value) {
    if constexpr (Stream) _mm512_stream_si512(reinterpret_cast<__m512i*>(output), value);
    else _mm512_storeu_si512(output, value);
}

template<bool Stream>
inline void subfield16OuterTranspose(const Subfield16Plan& plan, const Block* planes, Block* output) {
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

inline std::size_t subfield16TemporaryRecords(const Subfield16Plan& plan, Subfield16Kernel kernel) noexcept {
    return kernel != Subfield16Kernel::Scatter && plan.regionAligned() ? plan.groups : 0;
}

template<bool Fft, bool Packed = false>
inline void subfield16InnerMapDispatch(const Subfield16Plan& plan, const Block* input, Block* planes,
                                       Block* temporary, Subfield16Kernel kernel) {
    if (kernel == Subfield16Kernel::Scatter || !plan.regionAligned()) {
        Subfield16State state;
        subfield16InnerRange<true, Fft, Packed>(plan, input, planes, 0, plan.n / 16, 0, state);
    } else {
        if (!temporary) throw std::invalid_argument("GF16 region schedule requires L temporary records");
        if ((reinterpret_cast<std::uintptr_t>(planes) & 63U) == 0)
            subfield16RegionInner<true, Fft, Packed>(plan, input, planes, temporary);
        else subfield16RegionInner<false, Fft, Packed>(plan, input, planes, temporary);
    }
}

inline void subfield16InnerDispatch(const Subfield16Plan& plan, const Block* input, Block* planes,
                                    Block* temporary, Subfield16Kernel kernel,
                                    Subfield16MapKernel map = Subfield16MapKernel::Direct) {
    if (map == Subfield16MapKernel::AdditiveFft)
        subfield16InnerMapDispatch<true>(plan, input, planes, temporary, kernel);
    else if (map == Subfield16MapKernel::PackedDirect)
        subfield16InnerMapDispatch<false, true>(plan, input, planes, temporary, kernel);
    else subfield16InnerMapDispatch<false>(plan, input, planes, temporary, kernel);
}

inline void subfield16OuterDispatch(const Subfield16Plan& plan, const Block* planes, Block* output,
                                    Subfield16Kernel kernel) {
    if (kernel == Subfield16Kernel::RegionAllNonTemporal &&
        (reinterpret_cast<std::uintptr_t>(output) & 63U) == 0)
        subfield16OuterTranspose<true>(plan, planes, output);
    else subfield16OuterTranspose<false>(plan, planes, output);
}

// Exact in-place alias input==output is supported and leaves [k,n) untouched.
// Plane scratch has n records; temporary has subfield16TemporaryRecords records.
// Scratch buffers do not overlap input, output, or one another. No timed allocation.
inline void subfield16TransposeFast128(const Subfield16Plan& plan, const Block* input, Block* output,
                                      Block* planes, Block* temporary, Subfield16Kernel kernel,
                                      Subfield16MapKernel map = Subfield16MapKernel::Direct) {
    subfield16InnerDispatch(plan, input, planes, temporary, kernel, map);
    subfield16OuterDispatch(plan, planes, output, kernel);
}

} // namespace spin::research::char2
