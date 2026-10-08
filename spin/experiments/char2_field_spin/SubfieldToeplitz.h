#pragma once
#include "Subfield16.h"
#include "SubfieldToeplitzOuterAlternatives.h"
#include <type_traits>
#include <utility>

#if defined(_MSC_VER)
#define SPIN_TOEPLITZ_INLINE __forceinline
#else
#define SPIN_TOEPLITZ_INLINE inline __attribute__((always_inline))
#endif

namespace spin::research::char2 {

// New research family: each message row x in GF16^32 maps to (x,T*x),
// with an independently sampled uniform Toeplitz T. Coefficients below are
// embedded AES-field bytes. T[i,j]=diagonals[31+i-j], including zero entries.
// The setup PRNG supplies reproducible fixtures, not ideal independent setup.
struct SubfieldToeplitzMatrix {
    // Three embedded GF16 coefficients and one zero byte per row. Each
    // four-byte broadcast is bounded; no sixteen-byte load reads past this.
    std::array<std::uint8_t, 12> entries{};
};

inline bool subfieldToeplitzInvertible(std::array<std::uint8_t, 9> matrix) noexcept {
    for (unsigned column = 0; column != 3; ++column) {
        unsigned pivot = column;
        while (pivot != 3 && matrix[3 * pivot + column] == 0) ++pivot;
        if (pivot == 3) return false;
        if (pivot != column)
            for (unsigned j = column; j != 3; ++j)
                std::swap(matrix[3 * pivot + j], matrix[3 * column + j]);
        const auto inverse = sf16Inverse(matrix[4 * column]);
        for (unsigned row = column + 1; row != 3; ++row) {
            const auto factor = sf16Multiply(matrix[3 * row + column], inverse);
            for (unsigned j = column; j != 3; ++j)
                matrix[3 * row + j] ^= sf16Multiply(factor, matrix[3 * column + j]);
        }
    }
    return true;
}

template<class Generator>
inline SubfieldToeplitzMatrix subfieldToeplitzSampleMatrix(Generator& generator) {
    std::array<std::uint8_t, 9> abstract;
    do {
        auto packed = std::uint64_t(generator());
        for (auto& entry : abstract) { entry = std::uint8_t(packed & 15U); packed >>= 4; }
    } while (!subfieldToeplitzInvertible(abstract));
    SubfieldToeplitzMatrix matrix;
    for (unsigned row = 0; row != 3; ++row)
        for (unsigned column = 0; column != 3; ++column)
            matrix.entries[4 * row + column] = sf16Embed(abstract[3 * row + column]);
    return matrix;
}

template<unsigned Size>
inline constexpr unsigned subfieldToeplitzLeaves = Size == 1 ? 1 : 3 * subfieldToeplitzLeaves<Size / 2>;
template<>
inline constexpr unsigned subfieldToeplitzLeaves<1> = 1;
static_assert(subfieldToeplitzLeaves<32> == 243);

// Setup-only coefficient recursion. Diagonals use index Size-1+i-j.
// If T=[[T0,T1],[T2,T0]], the three children are T0,T2+T0,T1+T0.
template<unsigned Size>
inline void subfieldToeplitzPrepare(const std::uint8_t* diagonals, std::uint8_t* leaves,
                                    unsigned stride = 1) {
    if constexpr (Size == 1) *leaves = *diagonals;
    else {
        constexpr unsigned half = Size / 2, count = subfieldToeplitzLeaves<half>;
        std::array<std::uint8_t, 2 * half - 1> low{}, high{};
        for (unsigned j = 0; j != 2 * half - 1; ++j) {
            low[j] = diagonals[2 * half + j] ^ diagonals[half + j];
            high[j] = diagonals[j] ^ diagonals[half + j];
        }
        subfieldToeplitzPrepare<half>(diagonals + half, leaves, stride);
        subfieldToeplitzPrepare<half>(low.data(), leaves + count * stride, stride);
        subfieldToeplitzPrepare<half>(high.data(), leaves + 2 * count * stride, stride);
    }
}

enum class SubfieldToeplitzOuter : std::uint8_t { Karatsuba, Dense, KaratsubaPacked, KaratsubaDense2 };

struct SubfieldToeplitzPlan {
    static constexpr unsigned width = 16, outer = 64, dimension = 32, state = 3;
    static constexpr unsigned leafCount = subfieldToeplitzLeaves<dimension>;
    std::size_t k, n, groups;
    std::uint64_t seed;
    bool withDiagonal;
    std::vector<std::uint32_t> route, inverseRoute;
    std::vector<std::uint8_t> routedDiagonal;
    std::vector<SubfieldToeplitzMatrix> refresh;
    std::vector<std::array<std::uint8_t, 63>> outerDiagonals;
    // Four consecutive source rows per tile; each leaf has four bytes.
    // Last tile is padded with zero rows. No expanded 64-byte coefficients.
    std::vector<std::uint8_t> transposeLeaves;
    std::array<std::array<SubfieldIndexCoefficients, 2>, 4> coefficients;

    SubfieldToeplitzPlan(std::size_t message, std::uint64_t setupSeed, bool inputDiagonal = false)
        : k(message), n(2 * message), groups(message / dimension), seed(setupSeed), withDiagonal(inputDiagonal) {
        if (!k || k % dimension || k > std::numeric_limits<std::uint32_t>::max() / 2)
            throw std::invalid_argument("Toeplitz candidate requires positive K divisible by 32 within uint32 route limits");
        // Distinct streams construct only this profile's setup. All claims
        // about independence refer to the ideal setup, not these seeded PRNGs.
        std::mt19937_64 outerRng(seed ^ 0x07cbe6f912a043d8ULL);
        std::mt19937_64 routeRng(seed ^ 0xdf4a1b3285ec7096ULL);
        std::mt19937_64 matrixRng(seed ^ 0x58a06fc9d137eb24ULL);
        outerDiagonals.resize(groups);
        transposeLeaves.resize(((groups + 3) / 4) * 4 * leafCount);
        for (std::size_t group = 0; group != groups; ++group) {
            auto& diagonal = outerDiagonals[group];
            std::uint64_t packed = 0;
            for (unsigned j = 0; j != 63; ++j) {
                if (j % 16 == 0) packed = std::uint64_t(outerRng());
                diagonal[j] = sf16Embed(std::uint8_t(packed & 15U));
                packed >>= 4;
            }
            std::array<std::uint8_t, 63> reversed;
            std::reverse_copy(diagonal.begin(), diagonal.end(), reversed.begin());
            subfieldToeplitzPrepare<32>(reversed.data(),
                transposeLeaves.data() + (group / 4) * 4 * leafCount + group % 4, 4);
        }
        route.resize(n);
        for (unsigned region = 0; region != outer; ++region) {
            auto begin = route.begin() + region * groups;
            std::iota(begin, begin + groups, std::uint32_t(region * groups));
            for (std::size_t i = groups - 1; i; --i) {
                std::uniform_int_distribution<std::size_t> draw(0, i);
                std::swap(begin[i], begin[draw(routeRng)]);
            }
        }
        inverseRoute.resize(n);
        for (std::size_t i = 0; i != n; ++i)
            inverseRoute[route[i]] = static_cast<std::uint32_t>(i);
        if (withDiagonal) {
            std::mt19937_64 diagonalRng(seed ^ 0x782d043af9165bc1ULL);
            routedDiagonal.resize(n);
            // Sampling in routed order has exactly the same independent law.
            for (auto& entry : routedDiagonal) {
                std::uint8_t coefficient;
                do { coefficient = std::uint8_t(diagonalRng() & 15U); } while (!coefficient);
                entry = sf16Embed(coefficient);
            }
        }
        refresh.resize(n / width);
        for (auto& matrix : refresh) matrix = subfieldToeplitzSampleMatrix(matrixRng);
        for (unsigned batch = 0; batch != 4; ++batch)
            for (unsigned power = 1; power != 3; ++power)
                for (unsigned lane = 0; lane != 4; ++lane) {
                    const auto repeated = 0x0101010101010101ULL * sf16Embed(sf16Monomial(power, 4 * batch + lane));
                    coefficients[batch][power - 1].bytes[2 * lane] = repeated;
                    coefficients[batch][power - 1].bytes[2 * lane + 1] = repeated;
                }
    }

    bool regionAligned() const noexcept { return groups % width == 0; }
    std::size_t setupBytes() const noexcept {
        return (route.size() + inverseRoute.size()) * sizeof(std::uint32_t) + routedDiagonal.size() +
               refresh.size() * sizeof(SubfieldToeplitzMatrix) + outerDiagonals.size() * 63 +
               transposeLeaves.size() + sizeof(coefficients);
    }
};

inline std::vector<std::uint8_t> subfieldToeplitzForwardOracle(const SubfieldToeplitzPlan& plan,
    const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.k * bytes)
        throw std::invalid_argument("Toeplitz forward oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.n * bytes);
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned row = 0; row != 32; ++row)
            for (unsigned component = 0; component != bytes; ++component) {
                planes[(row * plan.groups + group) * bytes + component] = input[(group * 32 + row) * bytes + component];
                std::uint8_t value = 0;
                for (unsigned column = 0; column != 32; ++column)
                    value ^= sfMultiply(input[(group * 32 + column) * bytes + component],
                                        plan.outerDiagonals[group][31 + row - column]);
                planes[((32 + row) * plan.groups + group) * bytes + component] = value;
            }
    std::array<std::array<std::uint8_t, 16>, 3> q{};
    for (std::size_t step = 0; step != plan.refresh.size(); ++step) {
        std::array<std::array<std::uint8_t, 16>, 3> r{}, sum{};
        const auto& m = plan.refresh[step].entries;
        for (unsigned a = 0; a != 3; ++a)
            for (unsigned b = 0; b != 3; ++b)
                for (unsigned component = 0; component != bytes; ++component)
                    r[a][component] ^= sfMultiply(m[4 * a + b], q[b][component]);
        for (unsigned j = 0; j != 16; ++j) {
            const auto i = step * 16 + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto raw = planes[std::size_t(plan.route[i]) * bytes + component];
                const auto x = plan.withDiagonal ? sfMultiply(raw, plan.routedDiagonal[i]) : raw;
                auto y = x;
                for (unsigned a = 0; a != 3; ++a) {
                    const auto coefficient = sf16Embed(sf16Monomial(a, j));
                    y ^= sfMultiply(coefficient, r[a][component]);
                    sum[a][component] ^= sfMultiply(coefficient, x);
                }
                output[i * bytes + component] = y;
            }
        }
        for (unsigned a = 0; a != 3; ++a)
            for (unsigned component = 0; component != bytes; ++component)
                q[a][component] = r[a][component] ^ sum[a][component];
    }
    return output;
}

inline std::vector<std::uint8_t> subfieldToeplitzTransposeOracle(const SubfieldToeplitzPlan& plan,
    const std::vector<std::uint8_t>& input, unsigned bytes) {
    if (!sfValidWidth(bytes) || input.size() != plan.n * bytes)
        throw std::invalid_argument("Toeplitz transpose oracle dimensions");
    std::vector<std::uint8_t> planes(plan.n * bytes), output(plan.k * bytes);
    std::array<std::array<std::uint8_t, 16>, 3> r{};
    for (std::size_t step = plan.refresh.size(); step--;) {
        std::array<std::array<std::uint8_t, 16>, 3> sum{}, previous{};
        for (unsigned j = 0; j != 16; ++j) {
            const auto i = step * 16 + j;
            for (unsigned component = 0; component != bytes; ++component) {
                const auto v = input[i * bytes + component];
                auto u = v;
                for (unsigned a = 0; a != 3; ++a) {
                    const auto coefficient = sf16Embed(sf16Monomial(a, j));
                    u ^= sfMultiply(coefficient, r[a][component]);
                    sum[a][component] ^= sfMultiply(coefficient, v);
                }
                planes[std::size_t(plan.route[i]) * bytes + component] =
                    plan.withDiagonal ? sfMultiply(u, plan.routedDiagonal[i]) : u;
            }
        }
        const auto& m = plan.refresh[step].entries;
        for (unsigned a = 0; a != 3; ++a)
            for (unsigned b = 0; b != 3; ++b)
                for (unsigned component = 0; component != bytes; ++component)
                    previous[a][component] ^= sfMultiply(m[4 * b + a], r[b][component] ^ sum[b][component]);
        r = previous;
    }
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned column = 0; column != 32; ++column)
            for (unsigned component = 0; component != bytes; ++component) {
                auto value = planes[(column * plan.groups + group) * bytes + component];
                for (unsigned row = 0; row != 32; ++row)
                    value ^= sfMultiply(planes[((32 + row) * plan.groups + group) * bytes + component],
                                        plan.outerDiagonals[group][31 + row - column]);
                output[(group * 32 + column) * bytes + component] = value;
            }
    return output;
}

struct SubfieldToeplitzState {
    __m128i r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128(), r2 = _mm_setzero_si128();
};

template<unsigned Batch, bool Scatter, bool Diagonal = false, class Plan>
SPIN_TOEPLITZ_INLINE void subfieldToeplitzInnerBatch(const Plan& plan,
    const Block* input, Block* destination, std::size_t begin, std::size_t outputBegin,
    __m512i r0, __m512i r1, __m512i r2, __m512i& sum0, __m512i& sum1, __m512i& sum2,
    __m512i packedDiagonal) {
    const auto i = begin + 4 * Batch;
    const auto a1 = _mm512_load_si512(plan.coefficients[Batch][0].bytes.data());
    const auto a2 = _mm512_load_si512(plan.coefficients[Batch][1].bytes.data());
    const auto v = _mm512_loadu_si512(input + i);
    sum0 = _mm512_xor_si512(sum0, v);
    sum1 = _mm512_xor_si512(sum1, _mm512_gf2p8mul_epi8(v, a1));
    sum2 = _mm512_xor_si512(sum2, _mm512_gf2p8mul_epi8(v, a2));
    const auto low = _mm512_xor_si512(v, r0);
    const auto high = _mm512_xor_si512(_mm512_gf2p8mul_epi8(r1, a1), _mm512_gf2p8mul_epi8(r2, a2));
    auto u = _mm512_xor_si512(low, high);
    if constexpr (Diagonal) u = _mm512_gf2p8mul_epi8(u, sf16ExpandRow<Batch>(packedDiagonal));
    if constexpr (Scatter) {
        destination[plan.route[i]] = block(_mm512_castsi512_si128(u));
        destination[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
        destination[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
        destination[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
    } else _mm512_storeu_si512(destination + outputBegin + 4 * Batch, u);
}

// Three independent ZMM products compute M^T*u in lanes 0..2. Lane 3
// remains zero because each stored refresh row has a zero padding byte.
SPIN_TOEPLITZ_INLINE __m512i subfieldToeplitzRefresh(const SubfieldToeplitzMatrix& matrix,
    __m128i u0, __m128i u1, __m128i u2) {
    const auto& m = matrix.entries;
    const auto p0 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u0), sfBroadcastFour(m.data()));
    const auto p1 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u1), sfBroadcastFour(m.data() + 4));
    const auto p2 = _mm512_gf2p8mul_epi8(_mm512_broadcast_i32x4(u2), sfBroadcastFour(m.data() + 8));
    return _mm512_xor_si512(_mm512_xor_si512(p0, p1), p2);
}

template<bool Scatter, bool Diagonal = false, class Plan>
inline void subfieldToeplitzInnerRange(const Plan& plan, const Block* input,
    Block* destination, std::size_t firstStep, std::size_t endStep, std::size_t destinationBase,
    SubfieldToeplitzState& state) {
    for (std::size_t step = endStep; step-- > firstStep;) {
        const auto begin = step * 16, outputBegin = begin - destinationBase;
        const auto r0 = _mm512_broadcast_i32x4(state.r0), r1 = _mm512_broadcast_i32x4(state.r1);
        const auto r2 = _mm512_broadcast_i32x4(state.r2);
        auto sum0 = _mm512_setzero_si512(), sum1 = _mm512_setzero_si512(), sum2 = _mm512_setzero_si512();
        auto packed = _mm512_setzero_si512();
        if constexpr (Diagonal) packed = _mm512_broadcast_i32x4(
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(plan.routedDiagonal.data() + begin)));
        subfieldToeplitzInnerBatch<0, Scatter, Diagonal>(plan, input, destination, begin, outputBegin, r0, r1, r2, sum0, sum1, sum2, packed);
        subfieldToeplitzInnerBatch<1, Scatter, Diagonal>(plan, input, destination, begin, outputBegin, r0, r1, r2, sum0, sum1, sum2, packed);
        subfieldToeplitzInnerBatch<2, Scatter, Diagonal>(plan, input, destination, begin, outputBegin, r0, r1, r2, sum0, sum1, sum2, packed);
        subfieldToeplitzInnerBatch<3, Scatter, Diagonal>(plan, input, destination, begin, outputBegin, r0, r1, r2, sum0, sum1, sum2, packed);
        const auto previous = subfieldToeplitzRefresh(plan.refresh[step],
            _mm_xor_si128(state.r0, sfSumRecords4(sum0)),
            _mm_xor_si128(state.r1, sfSumRecords4(sum1)),
            _mm_xor_si128(state.r2, sfSumRecords4(sum2)));
        state.r0 = _mm512_castsi512_si128(previous);
        state.r1 = _mm512_extracti32x4_epi32(previous, 1);
        state.r2 = _mm512_extracti32x4_epi32(previous, 2);
    }
}

template<bool Stream, bool Diagonal = false, class Plan>
inline void subfieldToeplitzRegionInner(const Plan& plan, const Block* input,
                                       Block* planes, Block* temporary) {
    SubfieldToeplitzState state;
    for (unsigned region = Plan::outer; region--;) {
        const auto offset = region * plan.groups;
        subfieldToeplitzInnerRange<false, Diagonal>(plan, input, temporary, offset / 16,
            (offset + plan.groups) / 16, offset, state);
        for (std::size_t j = 0; j != plan.groups; j += 4) {
            const auto value = sf16GatherRegion4(temporary, plan.inverseRoute.data() + offset + j,
                                                static_cast<std::uint32_t>(offset));
            sf16Store<Stream>(planes + offset + j, value);
        }
    }
    if constexpr (Stream) _mm_sfence();
}

// Index packs unroll every fixed-width input XOR, output XOR, and outer
// store. The callable is a compile-time lambda, with no runtime dispatch.
template<class Function, std::size_t... Index>
SPIN_TOEPLITZ_INLINE void subfieldToeplitzUnroll(std::index_sequence<Index...>, Function&& function) {
    (function(std::integral_constant<std::size_t, Index>{}), ...);
}

template<unsigned Size, unsigned Offset = 0>
SPIN_TOEPLITZ_INLINE void subfieldToeplitzApply4(const std::uint8_t* leaves,
                                               const __m512i* input, __m512i* output) {
    if constexpr (Size == 1)
        output[0] = _mm512_gf2p8mul_epi8(input[0], sfBroadcastFour(leaves + 4 * Offset));
    else {
        constexpr unsigned half = Size / 2, count = subfieldToeplitzLeaves<half>;
        __m512i temporary[half];
        subfieldToeplitzUnroll(std::make_index_sequence<half>{}, [&](auto i) {
            temporary[i] = _mm512_xor_si512(input[i], input[i + half]);
        });
        subfieldToeplitzApply4<half, Offset>(leaves, temporary, output);
        subfieldToeplitzApply4<half, Offset + count>(leaves, input, output + half);
        subfieldToeplitzApply4<half, Offset + 2 * count>(leaves, input + half, temporary);
        subfieldToeplitzUnroll(std::make_index_sequence<half>{}, [&](auto i) {
            output[i + half] = _mm512_xor_si512(output[i], output[i + half]);
            output[i] = _mm512_xor_si512(output[i], temporary[i]);
        });
    }
}

template<unsigned Coordinate, bool Stream>
SPIN_TOEPLITZ_INLINE void subfieldToeplitzStore4(const SubfieldToeplitzPlan& plan,
    const Block* planes, Block* output, std::size_t group, const __m512i* parity) {
    auto a = _mm512_xor_si512(parity[Coordinate], _mm512_loadu_si512(planes + Coordinate * plan.groups + group));
    auto b = _mm512_xor_si512(parity[Coordinate + 1], _mm512_loadu_si512(planes + (Coordinate + 1) * plan.groups + group));
    auto c = _mm512_xor_si512(parity[Coordinate + 2], _mm512_loadu_si512(planes + (Coordinate + 2) * plan.groups + group));
    auto d = _mm512_xor_si512(parity[Coordinate + 3], _mm512_loadu_si512(planes + (Coordinate + 3) * plan.groups + group));
    sf16TransposeRecords4(a, b, c, d);
    sf16Store<Stream>(output + 32 * group + Coordinate, a);
    sf16Store<Stream>(output + 32 * (group + 1) + Coordinate, b);
    sf16Store<Stream>(output + 32 * (group + 2) + Coordinate, c);
    sf16Store<Stream>(output + 32 * (group + 3) + Coordinate, d);
}

template<bool Stream, bool Dense, bool Packed = false, bool Dense2 = false>
inline void subfieldToeplitzOuterTranspose(const SubfieldToeplitzPlan& plan,
                                          const Block* planes, Block* output) {
    std::size_t group = 0;
    for (; group + 4 <= plan.groups; group += 4) {
        __m512i input[32], parity[32];
        subfieldToeplitzUnroll(std::make_index_sequence<32>{}, [&](auto i) {
            input[i] = _mm512_loadu_si512(planes + (32 + i) * plan.groups + group);
        });
        if constexpr (Dense) {
            // Independent dense map for A/B and correctness. Only the fast
            // recursive path uses the prepared leaf coefficients.
            for (unsigned column = 0; column != 32; ++column) {
                auto value = _mm512_setzero_si512();
                for (unsigned row = 0; row != 32; ++row) {
                    const unsigned index = 31 + row - column;
                    const std::array<std::uint8_t, 4> coefficients{
                        plan.outerDiagonals[group][index], plan.outerDiagonals[group + 1][index],
                        plan.outerDiagonals[group + 2][index], plan.outerDiagonals[group + 3][index]};
                    value = _mm512_xor_si512(value, _mm512_gf2p8mul_epi8(input[row], sfBroadcastFour(coefficients.data())));
                }
                parity[column] = value;
            }
        } else if constexpr (Packed || Dense2)
            subfieldToeplitzApplyAlternative4<32, Dense2>(plan.transposeLeaves.data() + (group / 4) * 4 * plan.leafCount, input, parity);
        else subfieldToeplitzApply4<32>(plan.transposeLeaves.data() + (group / 4) * 4 * plan.leafCount, input, parity);
        subfieldToeplitzStore4<0, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<4, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<8, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<12, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<16, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<20, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<24, Stream>(plan, planes, output, group, parity);
        subfieldToeplitzStore4<28, Stream>(plan, planes, output, group, parity);
    }
    // At most three natural-size tail rows. All published power-of-two
    // candidates use the four-row SIMD path exclusively.
    for (; group != plan.groups; ++group)
        for (unsigned column = 0; column != 32; ++column) {
            auto result = reg(planes[column * plan.groups + group]);
            for (unsigned row = 0; row != 32; ++row)
                result = _mm_xor_si128(result, _mm_gf2p8mul_epi8(reg(planes[(32 + row) * plan.groups + group]),
                    _mm_set1_epi8(static_cast<char>(plan.outerDiagonals[group][31 + row - column]))));
            output[32 * group + column] = block(result);
        }
    if constexpr (Stream) _mm_sfence();
}

template<class Plan>
inline std::size_t subfieldToeplitzTemporaryRecords(const Plan& plan,
                                                   Subfield16Kernel kernel) noexcept {
    return kernel != Subfield16Kernel::Scatter && plan.regionAligned() ? plan.groups : 0;
}

template<bool Diagonal, class Plan>
inline void subfieldToeplitzInnerDiagonalDispatch(const Plan& plan, const Block* input,
    Block* planes, Block* temporary, Subfield16Kernel kernel) {
    if (kernel == Subfield16Kernel::Scatter || !plan.regionAligned()) {
        SubfieldToeplitzState state;
        subfieldToeplitzInnerRange<true, Diagonal>(plan, input, planes, 0, plan.n / 16, 0, state);
    } else {
        if (!temporary) throw std::invalid_argument("Toeplitz region schedule requires L temporary records");
        if ((reinterpret_cast<std::uintptr_t>(planes) & 63U) == 0)
            subfieldToeplitzRegionInner<true, Diagonal>(plan, input, planes, temporary);
        else subfieldToeplitzRegionInner<false, Diagonal>(plan, input, planes, temporary);
    }
}

template<class Plan>
inline void subfieldToeplitzInnerDispatch(const Plan& plan, const Block* input,
    Block* planes, Block* temporary, Subfield16Kernel kernel) {
    if (plan.withDiagonal) subfieldToeplitzInnerDiagonalDispatch<true>(plan, input, planes, temporary, kernel);
    else subfieldToeplitzInnerDiagonalDispatch<false>(plan, input, planes, temporary, kernel);
}

template<bool Dense, bool Packed = false, bool Dense2 = false>
inline void subfieldToeplitzOuterStoreDispatch(const SubfieldToeplitzPlan& plan, const Block* planes,
                                               Block* output, Subfield16Kernel kernel) {
    if (kernel == Subfield16Kernel::RegionAllNonTemporal && (reinterpret_cast<std::uintptr_t>(output) & 63U) == 0)
        subfieldToeplitzOuterTranspose<true, Dense, Packed, Dense2>(plan, planes, output);
    else subfieldToeplitzOuterTranspose<false, Dense, Packed, Dense2>(plan, planes, output);
}

inline void subfieldToeplitzOuterDispatch(const SubfieldToeplitzPlan& plan, const Block* planes,
    Block* output, Subfield16Kernel kernel, SubfieldToeplitzOuter outer = SubfieldToeplitzOuter::Karatsuba) {
    if (outer == SubfieldToeplitzOuter::Dense) subfieldToeplitzOuterStoreDispatch<true>(plan, planes, output, kernel);
    else if (outer == SubfieldToeplitzOuter::KaratsubaPacked) subfieldToeplitzOuterStoreDispatch<false, true>(plan, planes, output, kernel);
    else if (outer == SubfieldToeplitzOuter::KaratsubaDense2) subfieldToeplitzOuterStoreDispatch<false, false, true>(plan, planes, output, kernel);
    else subfieldToeplitzOuterStoreDispatch<false>(plan, planes, output, kernel);
}

// Exact input==output aliasing is supported, with [k,n) unchanged. Scratch
// contains n plane records plus TemporaryRecords region records, disjoint
// from input/output and each other. Encode has no allocation or setup work.
inline void subfieldToeplitzTransposeFast128(const SubfieldToeplitzPlan& plan,
    const Block* input, Block* output, Block* planes, Block* temporary, Subfield16Kernel kernel,
    SubfieldToeplitzOuter outer = SubfieldToeplitzOuter::Karatsuba) {
    subfieldToeplitzInnerDispatch(plan, input, planes, temporary, kernel);
    subfieldToeplitzOuterDispatch(plan, planes, output, kernel, outer);
}

} // namespace spin::research::char2

#undef SPIN_TOEPLITZ_INLINE
