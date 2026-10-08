#pragma once
#include "Subfield.h"

namespace spin::research::char2 {

enum class OneRouteFamily : std::uint8_t { Regional, Global };
enum class OneRouteKernel : std::uint8_t {
    Scatter, ScatterPrefetch, OuterGather, GatherPass, BucketPass, RegionGather,
    RegionPrefetch, RegionNonTemporal, RegionHardwareGather, ScatterOuterNonTemporal,
    RegionAllNonTemporal
};

// Execution experiment, not a distance claim. Regional is exactly SubfieldPlan.
// Global replaces its route by one uniform permutation of all n coordinates.
// In both families the diagonal and SL2 refresh distributions stay unchanged.
// Kernels select how the SAME sampled map is evaluated, never its distribution.
// As in SubfieldPlan, mt19937_64 is a reproducible research fixture, not the
// production randomness source required to instantiate the ideal setup law.
struct OneRoutePlan : SubfieldPlan {
    OneRouteFamily family;
    std::vector<std::uint32_t> inverseRoute;
    std::size_t bucketRecords;
    std::vector<std::uint32_t> bucketSources, bucketOffsets;

    OneRoutePlan(std::size_t message, unsigned t, std::uint64_t setupSeed,
                 unsigned b = 16, OneRouteFamily routeFamily = OneRouteFamily::Regional,
                 std::size_t bucketSize = 0)
        : SubfieldPlan(message, t, setupSeed, b), family(routeFamily), bucketRecords(bucketSize) {
        if (family == OneRouteFamily::Global) {
            std::vector<std::uint8_t> originalDiagonal(n);
            for (std::size_t i = 0; i != n; ++i) originalDiagonal[route[i]] = routedDiagonal[i];
            std::mt19937_64 generator(seed ^ 0x5793e18db402c6afULL);
            std::iota(route.begin(), route.end(), std::uint32_t(0));
            for (std::size_t i = n - 1; i; --i) {
                std::uniform_int_distribution<std::size_t> draw(0, i);
                std::swap(route[i], route[draw(generator)]);
            }
            for (std::size_t i = 0; i != n; ++i) routedDiagonal[i] = originalDiagonal[route[i]];
        }
        inverseRoute.resize(n);
        for (std::size_t i = 0; i != n; ++i) inverseRoute[route[i]] = static_cast<std::uint32_t>(i);
        if (bucketRecords) {
            const auto buckets = (n + bucketRecords - 1) / bucketRecords;
            bucketOffsets.assign(buckets + 1, 0);
            for (auto destination : route) ++bucketOffsets[destination / bucketRecords + 1];
            std::partial_sum(bucketOffsets.begin(), bucketOffsets.end(), bucketOffsets.begin());
            auto cursor = bucketOffsets;
            bucketSources.resize(n);
            // Within each destination tile, visit source records in ascending
            // order. Setup does the sorting once; encode does no allocation.
            for (std::size_t source = 0; source != n; ++source)
                bucketSources[cursor[route[source] / bucketRecords]++] = static_cast<std::uint32_t>(source);
        }
    }

    std::size_t setupBytes() const noexcept {
        return SubfieldPlan::setupBytes() + sizeof(std::uint32_t) *
            (inverseRoute.size() + bucketSources.size() + bucketOffsets.size());
    }
};

template<unsigned Batch>
inline void oneRouteLinearBatch(const OneRoutePlan& plan, const Block* input, Block* linear,
                                std::size_t begin, std::size_t outputBegin, __m512i state0, __m512i state1,
                                __m512i& sum0, __m512i& sum1) {
    const auto i = begin + 4 * Batch;
    const auto index = _mm512_load_si512(plan.indexCoefficients[Batch].bytes.data());
    const auto v = _mm512_loadu_si512(input + i);
    sum0 = _mm512_xor_si512(sum0, v);
    sum1 = _mm512_xor_si512(sum1, _mm512_gf2p8mul_epi8(v, index));
    auto u = _mm512_xor_si512(v, _mm512_xor_si512(state0, _mm512_gf2p8mul_epi8(state1, index)));
    u = _mm512_gf2p8mul_epi8(u, sfBroadcastFour(plan.routedDiagonal.data() + i));
    _mm512_storeu_si512(linear + outputBegin + 4 * Batch, u);
}

template<unsigned Width>
inline void oneRouteLinearInnerRange(const OneRoutePlan& plan, const Block* input, Block* linear,
                                     std::size_t firstStep, std::size_t pastStep,
                                     std::size_t outputOffset, __m128i& r0, __m128i& r1) {
    static_assert(Width == 8 || Width == 16);
    for (std::size_t step = pastStep; step-- > firstStep;) {
        const auto begin = step * Width;
        const auto outputBegin = begin - outputOffset;
        const auto state0 = _mm512_broadcast_i32x4(r0), state1 = _mm512_broadcast_i32x4(r1);
        auto sum0 = _mm512_setzero_si512(), sum1 = _mm512_setzero_si512();
        oneRouteLinearBatch<0>(plan, input, linear, begin, outputBegin, state0, state1, sum0, sum1);
        oneRouteLinearBatch<1>(plan, input, linear, begin, outputBegin, state0, state1, sum0, sum1);
        if constexpr (Width == 16) {
            oneRouteLinearBatch<2>(plan, input, linear, begin, outputBegin, state0, state1, sum0, sum1);
            oneRouteLinearBatch<3>(plan, input, linear, begin, outputBegin, state0, state1, sum0, sum1);
        }
        const auto s0 = _mm_xor_si128(r0, sfSumRecords4(sum0));
        const auto s1 = _mm_xor_si128(r1, sfSumRecords4(sum1));
        const auto m = plan.refresh[step];
        r0 = _mm_xor_si128(_mm_gf2p8mul_epi8(s0, _mm_set1_epi8(static_cast<char>(m.a))),
                           _mm_gf2p8mul_epi8(s1, _mm_set1_epi8(static_cast<char>(m.c))));
        r1 = _mm_xor_si128(_mm_gf2p8mul_epi8(s0, _mm_set1_epi8(static_cast<char>(m.b))),
                           _mm_gf2p8mul_epi8(s1, _mm_set1_epi8(static_cast<char>(m.d))));
    }
}

template<unsigned Width>
inline void oneRouteLinearInner(const OneRoutePlan& plan, const Block* input, Block* linear) {
    auto r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128();
    oneRouteLinearInnerRange<Width>(plan, input, linear, 0, plan.n / Width, 0, r0, r1);
}

inline __m512i oneRouteGather4(const Block* input, const std::uint32_t* index,
                              std::uint32_t offset = 0) noexcept {
    auto result = _mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(input + (index[0] - offset))));
    result = _mm512_inserti32x4(result, _mm_loadu_si128(reinterpret_cast<const __m128i*>(input + (index[1] - offset))), 1);
    result = _mm512_inserti32x4(result, _mm_loadu_si128(reinterpret_cast<const __m128i*>(input + (index[2] - offset))), 2);
    return _mm512_inserti32x4(result, _mm_loadu_si128(reinterpret_cast<const __m128i*>(input + (index[3] - offset))), 3);
}

inline __m512i oneRouteHardwareGather4(const Block* input, const std::uint32_t* index,
                                      std::uint32_t offset) noexcept {
    const auto local = _mm_sub_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(index)),
                                    _mm_set1_epi32(static_cast<int>(offset)));
    const auto pairs = _mm256_slli_epi64(_mm256_cvtepu32_epi64(local), 1);
    auto indices = _mm512_permutexvar_epi64(_mm512_setr_epi64(0, 0, 1, 1, 2, 2, 3, 3),
                                           _mm512_castsi256_si512(pairs));
    indices = _mm512_add_epi64(indices, _mm512_setr_epi64(0, 1, 0, 1, 0, 1, 0, 1));
    return _mm512_i64gather_epi64(indices, input, 8);
}

inline bool oneRouteIsRegionKernel(OneRouteKernel kernel) noexcept {
    return kernel == OneRouteKernel::RegionGather || kernel == OneRouteKernel::RegionPrefetch ||
           kernel == OneRouteKernel::RegionNonTemporal || kernel == OneRouteKernel::RegionHardwareGather ||
           kernel == OneRouteKernel::RegionAllNonTemporal;
}

inline bool oneRouteRegionAligned(const OneRoutePlan& plan) noexcept {
    return plan.family == OneRouteFamily::Regional && plan.groups % plan.width == 0;
}

template<unsigned Width, unsigned Prefetch = 0, bool NonTemporal = false, bool HardwareGather = false>
inline void oneRouteRegionalInner(const OneRoutePlan& plan, const Block* input,
                                  Block* planes, Block* temporary) {
    auto r0 = _mm_setzero_si128(), r1 = _mm_setzero_si128();
    for (unsigned region = plan.outer; region--;) {
        const auto offset = region * plan.groups;
        oneRouteLinearInnerRange<Width>(plan, input, temporary, offset / Width,
            (offset + plan.groups) / Width, offset, r0, r1);
        // This source region is still cache-hot. Reorder immediately, then
        // reuse the same temporary for the preceding region. Do not reset r.
        for (std::size_t destination = 0; destination != plan.groups; destination += 4) {
            if constexpr (Prefetch != 0) {
                if (destination + Prefetch < plan.groups) {
                    const auto future = plan.inverseRoute.data() + offset + destination + Prefetch;
                    _mm_prefetch(reinterpret_cast<const char*>(temporary + (future[0] - offset)), _MM_HINT_T0);
                    _mm_prefetch(reinterpret_cast<const char*>(temporary + (future[1] - offset)), _MM_HINT_T0);
                    _mm_prefetch(reinterpret_cast<const char*>(temporary + (future[2] - offset)), _MM_HINT_T0);
                    _mm_prefetch(reinterpret_cast<const char*>(temporary + (future[3] - offset)), _MM_HINT_T0);
                }
            }
            __m512i value;
            if constexpr (HardwareGather)
                value = oneRouteHardwareGather4(temporary, plan.inverseRoute.data() + offset + destination,
                                                static_cast<std::uint32_t>(offset));
            else
                value = oneRouteGather4(temporary, plan.inverseRoute.data() + offset + destination,
                                        static_cast<std::uint32_t>(offset));
            if constexpr (NonTemporal)
                _mm512_stream_si512(reinterpret_cast<__m512i*>(planes + offset + destination), value);
            else _mm512_storeu_si512(planes + offset + destination, value);
        }
    }
    // Every streamed plane must be visible before the outer starts reading.
    if constexpr (NonTemporal) _mm_sfence();
}

template<unsigned Width>
inline void oneRouteRegionalVariant(const OneRoutePlan& plan, const Block* input, Block* planes,
                                    Block* temporary, OneRouteKernel kernel) {
    if (kernel == OneRouteKernel::RegionPrefetch)
        oneRouteRegionalInner<Width, 64>(plan, input, planes, temporary);
    else if (kernel == OneRouteKernel::RegionHardwareGather)
        oneRouteRegionalInner<Width, 0, false, true>(plan, input, planes, temporary);
    else if ((kernel == OneRouteKernel::RegionNonTemporal || kernel == OneRouteKernel::RegionAllNonTemporal) &&
             (reinterpret_cast<std::uintptr_t>(planes) & 63U) == 0)
        oneRouteRegionalInner<Width, 0, true>(plan, input, planes, temporary);
    else oneRouteRegionalInner<Width>(plan, input, planes, temporary);
}

// Destination-linear writes; independent 128-bit loads avoid splitting a
// record into two AVX512 hardware-gather lanes.
inline void oneRouteGatherReorder(const OneRoutePlan& plan, const Block* linear, Block* planes) {
    std::size_t destination = 0;
    for (; destination + 4 <= plan.n; destination += 4)
        _mm512_storeu_si512(planes + destination,
            oneRouteGather4(linear, plan.inverseRoute.data() + destination));
    for (; destination != plan.n; ++destination)
        planes[destination] = linear[plan.inverseRoute[destination]];
}

inline void oneRouteBucketReorder(const OneRoutePlan& plan, const Block* linear, Block* planes) {
    for (std::size_t bucket = 0; bucket + 1 < plan.bucketOffsets.size(); ++bucket) {
        std::size_t cursor = plan.bucketOffsets[bucket];
        const auto end = plan.bucketOffsets[bucket + 1];
        for (; cursor + 4 <= end; cursor += 4) {
            const auto a = plan.bucketSources[cursor], b = plan.bucketSources[cursor + 1];
            const auto c = plan.bucketSources[cursor + 2], d = plan.bucketSources[cursor + 3];
            planes[plan.route[a]] = linear[a];
            planes[plan.route[b]] = linear[b];
            planes[plan.route[c]] = linear[c];
            planes[plan.route[d]] = linear[d];
        }
        for (; cursor != end; ++cursor) {
            const auto source = plan.bucketSources[cursor];
            planes[plan.route[source]] = linear[source];
        }
    }
}

template<unsigned Outer>
inline void oneRouteGatherOuter(const OneRoutePlan& plan, const Block* linear, Block* output) {
    static_assert(Outer == 8 || Outer == 16);
    constexpr unsigned dimension = Outer / 2;
    std::size_t group = 0;
    for (; group + 4 <= plan.groups; group += 4) {
        __m512i values[Outer];
        for (unsigned coordinate = 0; coordinate != Outer; ++coordinate)
            values[coordinate] = oneRouteGather4(linear,
                plan.inverseRoute.data() + coordinate * plan.groups + group);
        sfEvalTranspose4<dimension, 0>(values);
        sfEvalTranspose4<dimension, dimension>(values + dimension);
        for (unsigned j = 0; j != dimension; ++j)
            values[j] = _mm512_xor_si512(values[j], values[j + dimension]);
        for (unsigned j = 0; j != dimension; j += 4) {
            const auto a01 = _mm512_shuffle_i32x4(values[j], values[j + 1], 0x44);
            const auto b01 = _mm512_shuffle_i32x4(values[j + 2], values[j + 3], 0x44);
            const auto a23 = _mm512_shuffle_i32x4(values[j], values[j + 1], 0xee);
            const auto b23 = _mm512_shuffle_i32x4(values[j + 2], values[j + 3], 0xee);
            _mm512_storeu_si512(output + dimension * group + j, _mm512_shuffle_i32x4(a01, b01, 0x88));
            _mm512_storeu_si512(output + dimension * (group + 1) + j, _mm512_shuffle_i32x4(a01, b01, 0xdd));
            _mm512_storeu_si512(output + dimension * (group + 2) + j, _mm512_shuffle_i32x4(a23, b23, 0x88));
            _mm512_storeu_si512(output + dimension * (group + 3) + j, _mm512_shuffle_i32x4(a23, b23, 0xdd));
        }
    }
    for (; group != plan.groups; ++group)
        for (unsigned j = 0; j != dimension; ++j) {
            auto result = _mm_setzero_si128();
            for (unsigned coordinate = 0; coordinate != Outer; ++coordinate)
                result = _mm_xor_si128(result, _mm_gf2p8mul_epi8(
                    reg(linear[plan.inverseRoute[coordinate * plan.groups + group]]),
                    _mm_set1_epi8(static_cast<char>(sfBasis(j, coordinate)))));
            output[dimension * group + j] = block(result);
        }
}

template<bool Stream>
inline void oneRouteOuterStore(Block* output, __m512i value) noexcept {
    if constexpr (Stream) _mm512_stream_si512(reinterpret_cast<__m512i*>(output), value);
    else _mm512_storeu_si512(output, value);
}

// The Subfield.h four-group outer schedule, with only the final store policy
// varied. Stream requires a 64-byte-aligned output. A final partial group uses
// ordinary stores on distinct cache lines, and the fence completes all writes.
template<unsigned Outer, bool Stream>
inline void oneRouteOuterFinish(const OneRoutePlan& plan, const Block* planes, Block* output) {
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
            oneRouteOuterStore<Stream>(output + dimension * group + j,
                _mm512_shuffle_i32x4(a01, b01, 0x88));
            oneRouteOuterStore<Stream>(output + dimension * (group + 1) + j,
                _mm512_shuffle_i32x4(a01, b01, 0xdd));
            oneRouteOuterStore<Stream>(output + dimension * (group + 2) + j,
                _mm512_shuffle_i32x4(a23, b23, 0x88));
            oneRouteOuterStore<Stream>(output + dimension * (group + 3) + j,
                _mm512_shuffle_i32x4(a23, b23, 0xdd));
        }
    }
    for (; group != plan.groups; ++group)
        for (unsigned j = 0; j != dimension; ++j) {
            auto result = _mm_setzero_si128();
            for (unsigned coordinate = 0; coordinate != Outer; ++coordinate)
                result = _mm_xor_si128(result, _mm_gf2p8mul_epi8(reg(planes[coordinate * plan.groups + group]),
                    _mm_set1_epi8(static_cast<char>(sfBasis(j, coordinate)))));
            output[dimension * group + j] = block(result);
        }
    if constexpr (Stream) _mm_sfence();
}

template<bool Stream>
inline void oneRouteOuterFinishDispatch(const OneRoutePlan& plan, const Block* planes, Block* output) {
    if (plan.outer == 8) oneRouteOuterFinish<8, Stream>(plan, planes, output);
    else oneRouteOuterFinish<16, Stream>(plan, planes, output);
}

inline bool oneRouteNeedsSecondScratch(OneRouteKernel kernel) noexcept {
    return kernel == OneRouteKernel::GatherPass || kernel == OneRouteKernel::BucketPass ||
           oneRouteIsRegionKernel(kernel);
}

inline std::size_t oneRouteSecondScratchRecords(const OneRoutePlan& plan, OneRouteKernel kernel) noexcept {
    if (!oneRouteNeedsSecondScratch(kernel)) return 0;
    return oneRouteIsRegionKernel(kernel) && oneRouteRegionAligned(plan) ? plan.groups : plan.n;
}

inline void oneRouteInnerDispatch(const OneRoutePlan& plan, const Block* input, Block* scratch,
                                  OneRouteKernel kernel, Block* secondScratch = nullptr) {
    if (oneRouteIsRegionKernel(kernel) && oneRouteRegionAligned(plan)) {
        if (!secondScratch) throw std::invalid_argument("one-route region gather needs L temporary records");
        if (plan.width == 8) oneRouteRegionalVariant<8>(plan, input, scratch, secondScratch, kernel);
        else oneRouteRegionalVariant<16>(plan, input, scratch, secondScratch, kernel);
        return;
    }
    if (kernel == OneRouteKernel::Scatter || kernel == OneRouteKernel::ScatterPrefetch ||
        kernel == OneRouteKernel::ScatterOuterNonTemporal) {
        subfieldInnerTransposeDispatch(plan, input, scratch,
            kernel == OneRouteKernel::ScatterPrefetch ? SubfieldInnerKernel::Prefetch128 : SubfieldInnerKernel::Baseline);
    } else if (plan.width == 8) oneRouteLinearInner<8>(plan, input, scratch);
    else oneRouteLinearInner<16>(plan, input, scratch);
}

inline const Block* oneRouteReorderDispatch(const OneRoutePlan& plan, const Block* scratch,
                                           Block* secondScratch, OneRouteKernel kernel) {
    if (oneRouteIsRegionKernel(kernel) && oneRouteRegionAligned(plan)) return scratch;
    if (oneRouteNeedsSecondScratch(kernel)) {
        if (!secondScratch) throw std::invalid_argument("one-route reorder needs n records of second scratch");
        if (kernel == OneRouteKernel::GatherPass || oneRouteIsRegionKernel(kernel))
            oneRouteGatherReorder(plan, scratch, secondScratch);
        else {
            if (plan.bucketSources.size() != plan.n)
                throw std::invalid_argument("one-route bucket metadata was not prepared");
            oneRouteBucketReorder(plan, scratch, secondScratch);
        }
        return secondScratch;
    }
    return scratch;
}

inline void oneRouteOuterDispatch(const OneRoutePlan& plan, const Block* scratch, Block* output,
                                  OneRouteKernel kernel) {
    if (kernel == OneRouteKernel::OuterGather) {
        if (plan.outer == 8) oneRouteGatherOuter<8>(plan, scratch, output);
        else oneRouteGatherOuter<16>(plan, scratch, output);
    } else if (kernel == OneRouteKernel::ScatterOuterNonTemporal || kernel == OneRouteKernel::RegionAllNonTemporal) {
        if ((reinterpret_cast<std::uintptr_t>(output) & 63U) == 0)
            oneRouteOuterFinishDispatch<true>(plan, scratch, output);
        else oneRouteOuterFinishDispatch<false>(plan, scratch, output);
    } else subfieldOuterTransposeDispatch(plan, scratch, output);
}

// input has n records; output receives k records and may equal input exactly.
// Its suffix [k,n) stays unchanged. Scratch has n records. The second buffer
// has oneRouteSecondScratchRecords(plan,kernel) records: L for region-gather
// when regional boundaries align with inner steps, n for its exact full-pass
// fallback and other pass kernels, or zero otherwise. Scratch buffers must not
// overlap one another, input, or output. Setup and allocation belong to the caller.
inline void oneRouteTransposeFast128(const OneRoutePlan& plan, const Block* input, Block* output,
                                     Block* scratch, Block* secondScratch, OneRouteKernel kernel) {
    oneRouteInnerDispatch(plan, input, scratch, kernel, secondScratch);
    const auto routed = oneRouteReorderDispatch(plan, scratch, secondScratch, kernel);
    oneRouteOuterDispatch(plan, routed, output, kernel);
}

} // namespace spin::research::char2
