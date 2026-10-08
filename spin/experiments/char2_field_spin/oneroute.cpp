#include "OneRoute.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace c2 = spin::research::char2;
using Clock = std::chrono::steady_clock;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static c2::OneRouteKernel parseKernel(const std::string& name) {
    if (name == "scatter") return c2::OneRouteKernel::Scatter;
    if (name == "prefetch") return c2::OneRouteKernel::ScatterPrefetch;
    if (name == "outer-gather") return c2::OneRouteKernel::OuterGather;
    if (name == "gather-pass") return c2::OneRouteKernel::GatherPass;
    if (name == "bucket-pass") return c2::OneRouteKernel::BucketPass;
    if (name == "region-gather") return c2::OneRouteKernel::RegionGather;
    if (name == "region-prefetch") return c2::OneRouteKernel::RegionPrefetch;
    if (name == "region-nt") return c2::OneRouteKernel::RegionNonTemporal;
    if (name == "region-hwgather") return c2::OneRouteKernel::RegionHardwareGather;
    if (name == "scatter-outer-nt") return c2::OneRouteKernel::ScatterOuterNonTemporal;
    if (name == "region-nt-all") return c2::OneRouteKernel::RegionAllNonTemporal;
    throw std::invalid_argument("kernel must be scatter, prefetch, outer-gather, gather-pass, bucket-pass, region-gather, region-prefetch, region-nt, region-hwgather, scatter-outer-nt or region-nt-all");
}

static std::uint8_t byteDot(const std::vector<std::uint8_t>& a,
                            const std::vector<std::uint8_t>& b) {
    require(a.size() == b.size(), "dot dimensions");
    std::uint8_t value = 0;
    for (std::size_t i = 0; i != a.size(); ++i) value ^= c2::sfMultiply(a[i], b[i]);
    return value;
}

static void check() {
    std::mt19937_64 rng(20261006);
    const c2::Block guard{0x712938cf5b83d714ULL, 0x823f971bad04ea57ULL};
    for (auto family : {c2::OneRouteFamily::Regional, c2::OneRouteFamily::Global})
      for (unsigned outer : {8, 16})
       for (unsigned width : {8, 16})
        for (std::size_t k : {4, 8, 12, 16, 24, 28, 32, 64, 256, 1024}) {
            if (k % (outer / 2) || 2 * k % width) continue;
            for (std::uint64_t seed : {1, 91}) {
                c2::OneRoutePlan plan(k, width, seed, outer, family, 7);
                c2::SubfieldPlan baseline(k, width, seed, outer);
                std::vector<unsigned> seen(plan.n);
                for (std::size_t i = 0; i != plan.n; ++i) {
                    require(plan.route[i] < plan.n, "route index outside domain");
                    require(plan.inverseRoute[plan.route[i]] == i, "inverse route mismatch");
                    ++seen[plan.route[i]];
                    if (family == c2::OneRouteFamily::Regional)
                        require(plan.route[i] / plan.groups == i / plan.groups, "regional domain mismatch");
                }
                require(std::all_of(seen.begin(), seen.end(), [](auto count) { return count == 1; }),
                        "route is not a permutation");
                if (family == c2::OneRouteFamily::Regional) {
                    require(plan.route == baseline.route, "regional realization changed");
                    require(plan.routedDiagonal == baseline.routedDiagonal, "regional diagonal changed");
                }
                std::vector<std::uint8_t> diagonal(plan.n), referenceDiagonal(plan.n);
                for (std::size_t i = 0; i != plan.n; ++i) {
                    diagonal[plan.route[i]] = plan.routedDiagonal[i];
                    referenceDiagonal[baseline.route[i]] = baseline.routedDiagonal[i];
                }
                require(diagonal == referenceDiagonal, "route change modified original-coordinate diagonal");
                require(plan.refresh.size() == baseline.refresh.size() &&
                        std::memcmp(plan.refresh.data(), baseline.refresh.data(),
                                    plan.refresh.size() * sizeof(c2::SubfieldMatrix)) == 0,
                        "route change modified refresh matrices");
                for (std::size_t bucket = 0; bucket + 1 < plan.bucketOffsets.size(); ++bucket) {
                    const auto begin = plan.bucketOffsets[bucket], end = plan.bucketOffsets[bucket + 1];
                    require(std::is_sorted(plan.bucketSources.begin() + begin, plan.bucketSources.begin() + end),
                            "bucket sources are not ascending");
                    for (auto p = begin; p != end; ++p)
                        require(plan.route[plan.bucketSources[p]] / plan.bucketRecords == bucket,
                                "bucket contains foreign destination");
                }
                for (unsigned bytes : {1, 4, 8, 16}) {
                    std::vector<std::uint8_t> x(plan.k * bytes), v(plan.n * bytes);
                    for (auto& value : x) value = std::uint8_t(rng());
                    for (auto& value : v) value = std::uint8_t(rng());
                    const auto y = c2::subfieldForwardOracle(plan, x, bytes);
                    const auto expected = c2::subfieldTransposeOracle(plan, v, bytes);
                    require(byteDot(y, v) == byteDot(x, expected), "scalar GF256 adjoint mismatch");
                    if (bytes != 16) continue;
                    std::vector<c2::Block> input(plan.n), fast(plan.k + 2), scratch(plan.n + 2);
                    std::memcpy(input.data(), v.data(), v.size());
                    for (auto kernel : {c2::OneRouteKernel::Scatter, c2::OneRouteKernel::ScatterPrefetch,
                                        c2::OneRouteKernel::OuterGather, c2::OneRouteKernel::GatherPass,
                                        c2::OneRouteKernel::BucketPass, c2::OneRouteKernel::RegionGather,
                                        c2::OneRouteKernel::RegionPrefetch, c2::OneRouteKernel::RegionNonTemporal,
                                        c2::OneRouteKernel::RegionHardwareGather, c2::OneRouteKernel::ScatterOuterNonTemporal,
                                        c2::OneRouteKernel::RegionAllNonTemporal}) {
                        const auto secondCount = c2::oneRouteSecondScratchRecords(plan, kernel);
                        std::vector<c2::Block> second(secondCount + 2);
                        for (unsigned offset : {0, 1}) {
                            std::fill(fast.begin(), fast.end(), guard);
                            std::fill(scratch.begin(), scratch.end(), guard);
                            std::fill(second.begin(), second.end(), guard);
                            c2::oneRouteTransposeFast128(plan, input.data(), fast.data() + offset,
                                scratch.data() + offset, second.data() + offset, kernel);
                            require(std::memcmp(fast.data() + offset, expected.data(), expected.size()) == 0,
                                    "SIMD execution schedule changed the encoded map");
                            if (offset) require(fast[0] == guard && scratch[0] == guard && second[0] == guard,
                                                "prefix guard changed");
                            require(fast[plan.k + offset] == guard && scratch[plan.n + offset] == guard &&
                                    second[secondCount + offset] == guard, "suffix guard changed");
                        }
                        auto alias = input;
                        c2::oneRouteTransposeFast128(plan, alias.data(), alias.data(), scratch.data(), second.data(), kernel);
                        require(std::memcmp(alias.data(), expected.data(), expected.size()) == 0, "in-place prefix mismatch");
                        require(std::equal(alias.begin() + plan.k, alias.end(), input.begin() + plan.k),
                                "in-place suffix changed");
                    }
                    // Guarantee the streaming-store branches are exercised:
                    // std::vector<Block> need not have 64-byte alignment.
                    for (auto kernel : {c2::OneRouteKernel::RegionNonTemporal,
                                        c2::OneRouteKernel::ScatterOuterNonTemporal,
                                        c2::OneRouteKernel::RegionAllNonTemporal}) {
                        c2::BenchPayload alignedInput(plan.n, false), alignedScratch(plan.n, false);
                        const auto temporarySize = c2::oneRouteSecondScratchRecords(plan, kernel);
                        c2::BenchPayload temporary(std::max<std::size_t>(1, temporarySize), false);
                        std::memcpy(alignedInput.p, input.data(), plan.n * sizeof(c2::Block));
                        c2::oneRouteTransposeFast128(plan, alignedInput.p, alignedInput.p,
                            alignedScratch.p, temporary.p, kernel);
                        require(std::memcmp(alignedInput.p, expected.data(), expected.size()) == 0,
                                "aligned non-temporal in-place schedule mismatch");
                        require(std::equal(alignedInput.p + plan.k, alignedInput.p + plan.n,
                                           input.data() + plan.k), "aligned non-temporal suffix changed");
                    }
                }
            }
        }
    std::cout << "one-route: regional/global bijections, identical diagonals/refresh, scalar adjoints, "
                 "eleven exact SIMD schedules, natural sizes, guards and in-place suffix checks pass\n";
}

static double micros(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::micro>(end - begin).count();
}

static void benchmark(std::size_t k, unsigned width, unsigned outer, std::uint64_t seed,
                      unsigned calls, bool huge, bool phases, const std::string& familyName,
                      const std::string& kernelName, std::size_t bucketRecords) {
    const auto kernel = parseKernel(kernelName);
    const auto family = familyName == "regional" ? c2::OneRouteFamily::Regional : c2::OneRouteFamily::Global;
    const auto setupStart = Clock::now();
    c2::OneRoutePlan plan(k, width, seed, outer, family,
        kernel == c2::OneRouteKernel::BucketPass ? bucketRecords : 0);
    const auto setupMs = micros(setupStart, Clock::now()) / 1000;
    c2::BenchPayload input(plan.n, huge), scratch(plan.n, huge);
    std::unique_ptr<c2::BenchPayload> second;
    const auto secondCount = c2::oneRouteSecondScratchRecords(plan, kernel);
    if (secondCount) second = std::make_unique<c2::BenchPayload>(secondCount, huge);
    auto secondPointer = second ? second->p : nullptr;
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
    auto encode = [&] {
        c2::oneRouteTransposeFast128(plan, input.p, input.p, scratch.p, secondPointer, kernel);
    };
    for (unsigned i = 0; i != 5; ++i) encode();
    std::vector<double> total(calls), inner(calls), route(calls), finish(calls);
    for (unsigned sample = 0; sample != calls; ++sample) {
        const auto start = Clock::now();
        if (phases) {
            c2::oneRouteInnerDispatch(plan, input.p, scratch.p, kernel, secondPointer);
            const auto innerEnd = Clock::now();
            const auto routed = c2::oneRouteReorderDispatch(plan, scratch.p, secondPointer, kernel);
            const auto routeEnd = Clock::now();
            c2::oneRouteOuterDispatch(plan, routed, input.p, kernel);
            const auto end = Clock::now();
            inner[sample] = micros(start, innerEnd);
            route[sample] = micros(innerEnd, routeEnd);
            finish[sample] = micros(routeEnd, end);
            total[sample] = micros(start, end);
        } else {
            encode();
            total[sample] = micros(start, Clock::now());
        }
    }
    std::sort(total.begin(), total.end());
    std::uint64_t checksum = 0;
    for (std::size_t i = 0; i != k; ++i) {
        checksum = (checksum ^ input.p[i].lo) * 0x100000001b3ULL;
        checksum = (checksum ^ input.p[i].hi) * 0x100000001b3ULL;
    }
    std::cout << "candidate,K,t,outer,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,mode,kernel,bucket_records\n"
              << "gf256-sl2-pre-d-" << familyName << ',' << k << ',' << width << ',' << outer << ','
              << seed << ',' << calls << ',' << std::fixed << std::setprecision(4) << total[calls / 2] << ','
              << total[calls / 10] << ',' << total[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes()
              << ',' << (plan.n + secondCount) * sizeof(c2::Block) << ',' << std::hex << checksum << std::dec
              << ',' << (huge ? "huge-preferred" : "normal") << ',' << (phases ? "phases" : "full") << ','
              << kernelName << ',' << plan.bucketRecords << '\n';
    if (phases) {
        std::sort(inner.begin(), inner.end());
        std::sort(route.begin(), route.end());
        std::sort(finish.begin(), finish.end());
        std::cout << "inner_us,reorder_us,outer_us,total_us\n" << inner[calls / 2] << ','
                  << route[calls / 2] << ',' << finish[calls / 2] << ',' << total[calls / 2] << '\n';
        // Scatter includes routing in inner_us; outer-gather includes routing
        // in outer_us. Region-gather includes each regional reorder in inner_us
        // when aligned, or uses the explicit reorder phase in its fallback.
    }
}

int main(int argc, char** argv) {
    try {
        if (!c2::gfni512Available()) return 77;
        if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
        if (argc < 3 || argc > 11)
            throw std::invalid_argument("usage: oneroute K t [outer=16] [seed=113] [calls=51] [huge|normal] [full|phases] [regional|global] [scatter|prefetch|outer-gather|gather-pass|bucket-pass|region-gather|region-prefetch|region-nt|region-hwgather|scatter-outer-nt|region-nt-all] [bucket_records=1024]");
        const auto k = std::stoull(argv[1]);
        const auto width = std::stoul(argv[2]);
        const auto outer = argc > 3 ? std::stoul(argv[3]) : 16;
        const auto seed = argc > 4 ? std::stoull(argv[4]) : 113;
        const auto calls = argc > 5 ? std::stoul(argv[5]) : 51;
        const std::string pages = argc > 6 ? argv[6] : "huge";
        const std::string mode = argc > 7 ? argv[7] : "full";
        const std::string family = argc > 8 ? argv[8] : "regional";
        const std::string kernel = argc > 9 ? argv[9] : "scatter";
        const auto bucket = argc > 10 ? std::stoull(argv[10]) : 1024;
        if (calls < 1 || calls > 1000000 || (pages != "huge" && pages != "normal") ||
            (mode != "full" && mode != "phases") || (family != "regional" && family != "global") ||
            !bucket || bucket > std::numeric_limits<std::uint32_t>::max())
            throw std::invalid_argument("one-route benchmark options");
        benchmark(k, static_cast<unsigned>(width), static_cast<unsigned>(outer), seed,
                  static_cast<unsigned>(calls), pages == "huge", mode == "phases", family, kernel, bucket);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
