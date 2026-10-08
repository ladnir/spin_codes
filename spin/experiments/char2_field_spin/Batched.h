#pragma once
#include "Full.h"

namespace spin::research::char2 {

enum class BatchedRefresh : std::uint8_t { PostEmission, PreEmission };
enum class BatchedOuterKernel : std::uint8_t { SingleGroup, Groups4 };
enum class BatchedInnerKernel : std::uint8_t { Baseline, Unrolled, Prefetch128 };

// Research candidate, with no finite distance certificate. An outer group has
// four GF(2^128) inputs and eight GRS evaluations. Coordinate r always goes to
// region r; only the group order inside each region is uniformly shuffled.
// There is deliberately no per-group coordinate shuffle in this candidate.
struct BatchedPlan {
    std::size_t k, n, groups;
    unsigned width;
    std::uint64_t seed;
    BatchedRefresh refresh;
    bool inputDiagonal;
    // inner-position -> outer plane offset r*groups+group. Each region's
    // destinations are confined to its contiguous plane, not the whole array.
    std::vector<std::uint32_t> route;
    std::vector<std::array<std::uint8_t, 3>> shear;
    std::array<Coeff4, 4> indexCoefficients;
    // Independent diagonal coefficients are assigned to original outer outputs,
    // then stored in inner traversal order. This cache changes no distribution
    // and lets four routed stores share a contiguous four-coefficient load.
    std::vector<std::uint8_t> routedDiagonal;

    BatchedPlan(std::size_t message, unsigned t, std::uint64_t setupSeed,
                BatchedRefresh placement = BatchedRefresh::PostEmission,
                bool diagonal = false)
        : k(message), n(2 * message), groups(message / 4), width(t), seed(setupSeed),
          refresh(placement), inputDiagonal(diagonal) {
        if ((t != 8 && t != 16) || !k || k % 4 ||
            k > std::numeric_limits<std::uint32_t>::max() / 2 || n % t ||
            (placement != BatchedRefresh::PostEmission && placement != BatchedRefresh::PreEmission))
            throw std::invalid_argument("batched candidate requires t8/16 and natural size");
        route.resize(n);
        std::mt19937_64 rng(seed);
        for (unsigned region = 0; region != 8; ++region) {
            auto begin = route.begin() + region * groups;
            std::iota(begin, begin + groups, std::uint32_t(region * groups));
            for (std::size_t i = groups - 1; i; --i) {
                std::uniform_int_distribution<std::size_t> draw(0, i);
                std::swap(begin[i], begin[draw(rng)]);
            }
        }
        shear.resize(n / t);
        for (auto& step : shear)
            for (auto& c : step) c = std::uint8_t(1 + (rng() & 127));
        for (unsigned batch = 0; batch != 4; ++batch) {
            const std::array<std::uint8_t, 4> c{
                std::uint8_t(4 * batch), std::uint8_t(4 * batch + 1),
                std::uint8_t(4 * batch + 2), std::uint8_t(4 * batch + 3)};
            indexCoefficients[batch] = prepare4(c.data());
        }
        if (inputDiagonal) {
            // Independent setup domain; legacy post-emission/no-D sampling is
            // byte-for-byte unchanged. The construction descriptor includes
            // seed, width, refresh placement and the diagonal flag.
            std::mt19937_64 diagonalRng(seed ^ 0x34c9d8a072f51be6ULL);
            std::vector<std::uint8_t> original(n);
            for (auto& c : original) c = std::uint8_t(1 + (diagonalRng() & 127));
            routedDiagonal.resize(n);
            for (std::size_t i = 0; i != n; ++i) routedDiagonal[i] = original[route[i]];
        }
    }

    std::size_t setupBytes() const noexcept {
        return route.size() * sizeof(route[0]) + shear.size() * sizeof(shear[0]) +
               sizeof(indexCoefficients) + routedDiagonal.size();
    }
};

// M=U_c L_b U_a; a shear adds a scalar multiple of one state to the other.
// Each shear is invertible for every coefficient, including in characteristic2.
inline std::array<Block, 2> shearForwardOracle(std::array<Block, 2> q,
                                               const std::array<std::uint8_t, 3>& c) {
    q[0] = bxor(q[0], mulSmallOracle(q[1], c[0]));
    q[1] = bxor(q[1], mulSmallOracle(q[0], c[1]));
    q[0] = bxor(q[0], mulSmallOracle(q[1], c[2]));
    return q;
}

inline std::array<Block, 2> shearTransposeOracle(std::array<Block, 2> r,
                                                 const std::array<std::uint8_t, 3>& c) {
    r[1] = bxor(r[1], mulSmallOracle(r[0], c[2]));
    r[0] = bxor(r[0], mulSmallOracle(r[1], c[1]));
    r[1] = bxor(r[1], mulSmallOracle(r[0], c[0]));
    return r;
}

inline std::vector<Block> batchedForwardOracle(const BatchedPlan& plan,
                                              const std::vector<Block>& input) {
    if (input.size() != plan.k) throw std::invalid_argument("batched oracle input size");
    std::vector<Block> planes(plan.n), output(plan.n);
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned coordinate = 0; coordinate != 8; ++coordinate)
            for (unsigned j = 0; j != 4; ++j)
                planes[coordinate * plan.groups + group] = bxor(
                    planes[coordinate * plan.groups + group],
                    mulSmallOracle(input[4 * group + j], basis(j, coordinate)));
    std::array<Block, 2> state{};
    for (std::size_t step = 0; step != plan.shear.size(); ++step) {
        if (plan.refresh == BatchedRefresh::PreEmission)
            state = shearForwardOracle(state, plan.shear[step]);
        std::array<Block, 2> contribution{};
        for (unsigned j = 0; j != plan.width; ++j) {
            const auto i = step * plan.width + j;
            auto x = planes[plan.route[i]];
            if (plan.inputDiagonal) x = mulSmallOracle(x, plan.routedDiagonal[i]);
            output[i] = bxor(x, bxor(state[0], mulSmallOracle(state[1], j)));
            contribution[0] = bxor(contribution[0], x);
            contribution[1] = bxor(contribution[1], mulSmallOracle(x, j));
        }
        if (plan.refresh == BatchedRefresh::PostEmission)
            state = shearForwardOracle(state, plan.shear[step]);
        state[0] = bxor(state[0], contribution[0]);
        state[1] = bxor(state[1], contribution[1]);
    }
    return output;
}

inline std::vector<Block> batchedTransposeOracle(const BatchedPlan& plan,
                                                const std::vector<Block>& input) {
    if (input.size() != plan.n) throw std::invalid_argument("batched transpose input size");
    std::vector<Block> planes(plan.n), output(plan.k);
    std::array<Block, 2> state{};
    for (std::size_t step = plan.shear.size(); step--;) {
        std::array<Block, 2> contribution{};
        for (unsigned j = 0; j != plan.width; ++j) {
            const auto i = step * plan.width + j;
            const auto v = input[i];
            auto u = bxor(v, bxor(state[0], mulSmallOracle(state[1], j)));
            if (plan.inputDiagonal) u = mulSmallOracle(u, plan.routedDiagonal[i]);
            planes[plan.route[i]] = u;
            contribution[0] = bxor(contribution[0], v);
            contribution[1] = bxor(contribution[1], mulSmallOracle(v, j));
        }
        if (plan.refresh == BatchedRefresh::PostEmission)
            state = shearTransposeOracle(state, plan.shear[step]);
        state[0] = bxor(state[0], contribution[0]);
        state[1] = bxor(state[1], contribution[1]);
        if (plan.refresh == BatchedRefresh::PreEmission)
            state = shearTransposeOracle(state, plan.shear[step]);
    }
    for (std::size_t group = 0; group != plan.groups; ++group)
        for (unsigned j = 0; j != 4; ++j)
            for (unsigned coordinate = 0; coordinate != 8; ++coordinate)
                output[4 * group + j] = bxor(output[4 * group + j],
                    mulSmallOracle(planes[coordinate * plan.groups + group], basis(j, coordinate)));
    return output;
}

inline __m128i sumRecords4(__m512i x) noexcept {
    return _mm_xor_si128(
        _mm_xor_si128(_mm512_castsi512_si128(x), _mm512_extracti32x4_epi32(x, 1)),
        _mm_xor_si128(_mm512_extracti32x4_epi32(x, 2), _mm512_extracti32x4_epi32(x, 3)));
}

template<unsigned Width, bool PreEmission = false, bool InputDiagonal = false>
inline void batchedInnerTranspose(const BatchedPlan& plan, const Block* input, Block* planes) {
    static_assert(Width == 8 || Width == 16);
    auto r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128();
    for (std::size_t step = plan.n / Width; step--;) {
        const auto state0 = _mm512_broadcast_i32x4(r0);
        const auto state1 = _mm512_broadcast_i32x4(r1);
        auto sum0 = _mm512_setzero_si512(), sum1 = _mm512_setzero_si512();
        for (unsigned batch = 0; batch != Width / 4; ++batch) {
            const auto i = step * Width + 4 * batch;
            const auto v = _mm512_loadu_si512(input + i);
            sum0 = _mm512_xor_si512(sum0, v);
            sum1 = _mm512_xor_si512(sum1, smallmul4(v, plan.indexCoefficients[batch]));
            auto u = _mm512_xor_si512(v, _mm512_xor_si512(state0,
                smallmul4(state1, plan.indexCoefficients[batch])));
            if constexpr (InputDiagonal) u = smallmul4(u, plan.routedDiagonal.data() + i);
            planes[plan.route[i]] = block(_mm512_castsi512_si128(u));
            planes[plan.route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
            planes[plan.route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
            planes[plan.route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
        }
        if constexpr (PreEmission) {
            r0 = _mm_xor_si128(r0, sumRecords4(sum0));
            r1 = _mm_xor_si128(r1, sumRecords4(sum1));
        }
        const auto& c = plan.shear[step];
        r1 = _mm_xor_si128(r1, smallmul128(r0, c[2]));
        r0 = _mm_xor_si128(r0, smallmul128(r1, c[1]));
        r1 = _mm_xor_si128(r1, smallmul128(r0, c[0]));
        if constexpr (!PreEmission) {
            r0 = _mm_xor_si128(r0, sumRecords4(sum0));
            r1 = _mm_xor_si128(r1, sumRecords4(sum1));
        }
    }
}

// This alternate retains the baseline algebra and store order, but expands
// each four-record batch at compile time. It is separately selectable so an
// apparent improvement does not silently replace the measured baseline.
#if defined(_MSC_VER)
#define SPIN_CHAR2_BATCH_INLINE __forceinline
#else
#define SPIN_CHAR2_BATCH_INLINE inline __attribute__((always_inline))
#endif
template<unsigned Batch, bool InputDiagonal, unsigned Prefetch>
SPIN_CHAR2_BATCH_INLINE void batchedFourRecords(
    const Block* input, Block* planes, const std::uint32_t* route,
    const std::uint8_t* diagonal, const Coeff4* indexCoefficients,
    std::size_t stepBegin, __m512i state0, __m512i state1,
    __m512i& sum0, __m512i& sum1) {
    const auto i = stepBegin + 4 * Batch;
    if constexpr (Prefetch != 0) {
        if (i >= Prefetch) {
            _mm_prefetch(reinterpret_cast<const char*>(planes + route[i - Prefetch]), _MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(planes + route[i - Prefetch + 1]), _MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(planes + route[i - Prefetch + 2]), _MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(planes + route[i - Prefetch + 3]), _MM_HINT_T0);
        }
    }
    const auto v = _mm512_loadu_si512(input + i);
    sum0 = _mm512_xor_si512(sum0, v);
    sum1 = _mm512_xor_si512(sum1, smallmul4(v, indexCoefficients[Batch]));
    auto u = _mm512_xor_si512(v, _mm512_xor_si512(state0,
        smallmul4(state1, indexCoefficients[Batch])));
    if constexpr (InputDiagonal) u = smallmul4(u, diagonal + i);
    planes[route[i]] = block(_mm512_castsi512_si128(u));
    planes[route[i + 1]] = block(_mm512_extracti32x4_epi32(u, 1));
    planes[route[i + 2]] = block(_mm512_extracti32x4_epi32(u, 2));
    planes[route[i + 3]] = block(_mm512_extracti32x4_epi32(u, 3));
}
#undef SPIN_CHAR2_BATCH_INLINE

template<unsigned Width, bool PreEmission, bool InputDiagonal, unsigned Prefetch = 0>
inline void batchedInnerUnrolled(const BatchedPlan& plan, const Block* input, Block* planes) {
    static_assert(Width == 8 || Width == 16);
    auto r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128();
    const auto* route = plan.route.data();
    const auto* diagonal = plan.routedDiagonal.data();
    const auto* maps = plan.indexCoefficients.data();
    for (std::size_t step = plan.n / Width; step--;) {
        const auto begin = step * Width;
        const auto state0 = _mm512_broadcast_i32x4(r0);
        const auto state1 = _mm512_broadcast_i32x4(r1);
        auto sum0 = _mm512_setzero_si512(), sum1 = _mm512_setzero_si512();
        batchedFourRecords<0, InputDiagonal, Prefetch>(input, planes, route, diagonal,
            maps, begin, state0, state1, sum0, sum1);
        batchedFourRecords<1, InputDiagonal, Prefetch>(input, planes, route, diagonal,
            maps, begin, state0, state1, sum0, sum1);
        if constexpr (Width == 16) {
            batchedFourRecords<2, InputDiagonal, Prefetch>(input, planes, route, diagonal,
                maps, begin, state0, state1, sum0, sum1);
            batchedFourRecords<3, InputDiagonal, Prefetch>(input, planes, route, diagonal,
                maps, begin, state0, state1, sum0, sum1);
        }
        if constexpr (PreEmission) {
            r0 = _mm_xor_si128(r0, sumRecords4(sum0));
            r1 = _mm_xor_si128(r1, sumRecords4(sum1));
        }
        const auto& c = plan.shear[step];
        r1 = _mm_xor_si128(r1, smallmul128(r0, c[2]));
        r0 = _mm_xor_si128(r0, smallmul128(r1, c[1]));
        r1 = _mm_xor_si128(r1, smallmul128(r0, c[0]));
        if constexpr (!PreEmission) {
            r0 = _mm_xor_si128(r0, sumRecords4(sum0));
            r1 = _mm_xor_si128(r1, sumRecords4(sum1));
        }
    }
}

inline void batchedOuterTransposeReference(const BatchedPlan& plan, const Block* planes,
                                           Block* output, std::size_t firstGroup = 0) {
    for (std::size_t group = firstGroup; group != plan.groups; ++group) {
        __m128i v[8];
        for (unsigned coordinate = 0; coordinate != 8; ++coordinate)
            v[coordinate] = reg(planes[coordinate * plan.groups + group]);
        evalTranspose<4, 0>(v);
        evalTranspose<4, 4>(v + 4);
        for (unsigned j = 0; j != 4; ++j)
            output[4 * group + j] = block(_mm_xor_si128(v[j], v[j + 4]));
    }
}

template<std::uint32_t C> inline __m512i batchedConstantMul4(__m512i x) {
    static_assert(C < 256); // All GRS[8,4] butterfly constants fit one byte.
    if constexpr (C == 0) return _mm512_setzero_si512();
    else if constexpr (C == 1) return x;
    else return smallmul4(x, Coeff128{affineMatrix(C, false), affineMatrix(C, true)});
}

template<unsigned Count, unsigned Offset>
inline void batchedEvalTranspose4(__m512i* values) {
    if constexpr (Count > 1) {
        constexpr unsigned half = Count / 2, degree = std::countr_zero(half);
        batchedEvalTranspose4<half, Offset>(values);
        batchedEvalTranspose4<half, Offset ^ half>(values + half);
        for (unsigned j = 0; j != half; ++j) {
            const auto sum = _mm512_xor_si512(values[j], values[j + half]);
            values[j + half] = _mm512_xor_si512(
                batchedConstantMul4<subspace(degree, Offset)>(sum),
                batchedConstantMul4<subspace(degree, half)>(values[j + half]));
            values[j] = sum;
        }
    }
}

inline void batchedOuterTransposeFast(const BatchedPlan& plan, const Block* planes, Block* output) {
    std::size_t group = 0;
    for (; group + 4 <= plan.groups; group += 4) {
        __m512i values[8];
        for (unsigned coordinate = 0; coordinate != 8; ++coordinate)
            values[coordinate] = _mm512_loadu_si512(planes + coordinate * plan.groups + group);
        batchedEvalTranspose4<4, 0>(values);
        batchedEvalTranspose4<4, 4>(values + 4);
        for (unsigned j = 0; j != 4; ++j)
            values[j] = _mm512_xor_si512(values[j], values[j + 4]);
        // Values[j] holds message coordinate j for four independent groups.
        // Transpose the 4x4 grid of 128-bit records into contiguous group order.
        const auto a01 = _mm512_shuffle_i32x4(values[0], values[1], 0x44);
        const auto b01 = _mm512_shuffle_i32x4(values[2], values[3], 0x44);
        const auto a23 = _mm512_shuffle_i32x4(values[0], values[1], 0xee);
        const auto b23 = _mm512_shuffle_i32x4(values[2], values[3], 0xee);
        _mm512_storeu_si512(output + 4 * group,      _mm512_shuffle_i32x4(a01, b01, 0x88));
        _mm512_storeu_si512(output + 4 * group + 4,  _mm512_shuffle_i32x4(a01, b01, 0xdd));
        _mm512_storeu_si512(output + 4 * group + 8,  _mm512_shuffle_i32x4(a23, b23, 0x88));
        _mm512_storeu_si512(output + 4 * group + 12, _mm512_shuffle_i32x4(a23, b23, 0xdd));
    }
    batchedOuterTransposeReference(plan, planes, output, group);
}

inline void batchedOuterTranspose(const BatchedPlan& plan, const Block* planes, Block* output) {
    batchedOuterTransposeFast(plan, planes, output);
}

template<unsigned Width, bool PreEmission, bool InputDiagonal>
inline void batchedInnerKernelDispatch(const BatchedPlan& plan, const Block* input,
                                       Block* scratch, BatchedInnerKernel kernel) {
    if (kernel == BatchedInnerKernel::Unrolled)
        batchedInnerUnrolled<Width, PreEmission, InputDiagonal, 0>(plan, input, scratch);
    else if (kernel == BatchedInnerKernel::Prefetch128)
        batchedInnerUnrolled<Width, PreEmission, InputDiagonal, 128>(plan, input, scratch);
    else batchedInnerTranspose<Width, PreEmission, InputDiagonal>(plan, input, scratch);
}

template<unsigned Width>
inline void batchedInnerTransposeDispatch(const BatchedPlan& plan, const Block* input, Block* scratch,
                                          BatchedInnerKernel kernel = BatchedInnerKernel::Baseline) {
    if (plan.refresh == BatchedRefresh::PreEmission) {
        if (plan.inputDiagonal) batchedInnerKernelDispatch<Width, true, true>(plan, input, scratch, kernel);
        else batchedInnerKernelDispatch<Width, true, false>(plan, input, scratch, kernel);
    } else {
        if (plan.inputDiagonal) batchedInnerKernelDispatch<Width, false, true>(plan, input, scratch, kernel);
        else batchedInnerKernelDispatch<Width, false, false>(plan, input, scratch, kernel);
    }
}

inline void batchedTransposeFast(const BatchedPlan& plan, const Block* input,
                                 Block* output, Block* scratch,
                                 BatchedOuterKernel outer = BatchedOuterKernel::Groups4,
                                 BatchedInnerKernel inner = BatchedInnerKernel::Baseline) {
    if (plan.width == 8) batchedInnerTransposeDispatch<8>(plan, input, scratch, inner);
    else batchedInnerTransposeDispatch<16>(plan, input, scratch, inner);
    if (outer == BatchedOuterKernel::Groups4) batchedOuterTransposeFast(plan, scratch, output);
    else batchedOuterTransposeReference(plan, scratch, output);
}

} // namespace spin::research::char2
