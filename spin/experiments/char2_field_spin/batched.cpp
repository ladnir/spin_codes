#include "Batched.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

namespace c2 = spin::research::char2;
using Clock = std::chrono::steady_clock;

static void fill(std::vector<c2::Block>& x, std::mt19937_64& rng) {
    for (auto& value : x) value = {rng(), rng()};
}

static void check() {
    std::mt19937_64 rng(20261004);
    for (const auto refresh : {c2::BatchedRefresh::PostEmission, c2::BatchedRefresh::PreEmission})
      for (const bool diagonal : {false, true})
    for (const unsigned width : {8, 16})
        for (const std::size_t k : {4, 8, 12, 16, 24, 28, 32, 64, 256}) {
            if ((2 * k) % width) continue;
            for (const unsigned seed : {1, 91}) {
                c2::BatchedPlan plan(k, width, seed, refresh, diagonal);
                std::vector<c2::Block> x(k), v(2 * k), fast(k), scratch(2 * k);
                fill(x, rng);
                fill(v, rng);
                const auto expected = c2::batchedTransposeOracle(plan, v);
                c2::batchedTransposeFast(plan, v.data(), fast.data(), scratch.data());
                if (fast != expected) throw std::runtime_error("batched scalar/SIMD mismatch");
                c2::batchedTransposeFast(plan, v.data(), fast.data(), scratch.data(),
                                        c2::BatchedOuterKernel::SingleGroup);
                if (fast != expected) throw std::runtime_error("batched single-group outer mismatch");
                for (const auto inner : {c2::BatchedInnerKernel::Unrolled, c2::BatchedInnerKernel::Prefetch128}) {
                    c2::batchedTransposeFast(plan, v.data(), fast.data(), scratch.data(),
                                            c2::BatchedOuterKernel::Groups4, inner);
                    if (fast != expected) throw std::runtime_error("batched alternate inner mismatch");
                }
                const auto y = c2::batchedForwardOracle(plan, x);
                c2::Block lhs{}, rhs{};
                for (std::size_t i = 0; i != 2 * k; ++i)
                    lhs = c2::bxor(lhs, c2::mulOracle(y[i], v[i]));
                for (std::size_t i = 0; i != k; ++i)
                    rhs = c2::bxor(rhs, c2::mulOracle(x[i], fast[i]));
                if (lhs != rhs) throw std::runtime_error("batched field adjoint mismatch");
                auto aliased = v;
                c2::batchedTransposeFast(plan, aliased.data(), aliased.data(), scratch.data());
                if (!std::equal(expected.begin(), expected.end(), aliased.begin()) ||
                    !std::equal(v.begin() + k, v.end(), aliased.begin() + k))
                    throw std::runtime_error("batched alias or suffix mismatch");
                const c2::Block scalar{rng(), rng()};
                auto scaled = x;
                for (auto& value : scaled) value = c2::mulOracle(value, scalar);
                const auto ys = c2::batchedForwardOracle(plan, scaled);
                for (std::size_t i = 0; i != 2 * k; ++i)
                    if (ys[i] != c2::mulOracle(y[i], scalar))
                        throw std::runtime_error("batched full-field linearity mismatch");
                // Each region is a bijection onto its own contiguous plane.
                for (unsigned region = 0; region != 8; ++region) {
                    std::vector<std::uint32_t> seen(plan.route.begin() + region * plan.groups,
                                                  plan.route.begin() + (region + 1) * plan.groups);
                    std::sort(seen.begin(), seen.end());
                    for (std::size_t j = 0; j != plan.groups; ++j)
                        if (seen[j] != region * plan.groups + j)
                            throw std::runtime_error("batched region permutation mismatch");
                }
            }
        }
    std::cout << "batched candidate: both placements, both diagonals, both SIMD outer kernels, "
                 "scalar, field adjoint, field linearity, aliases, routes pass\n";
}

static void benchmark(std::size_t k, unsigned width, std::uint64_t seed, unsigned calls,
                      bool preferHuge, const std::string& mode, const std::string& kernel,
                      const std::string& innerName) {
    const auto refresh = mode == "pre" || mode == "pre-d"
        ? c2::BatchedRefresh::PreEmission : c2::BatchedRefresh::PostEmission;
    const bool diagonal = mode == "pre-d" || mode == "post-d";
    const auto outer = kernel == "scalar" ? c2::BatchedOuterKernel::SingleGroup
                                          : c2::BatchedOuterKernel::Groups4;
    const auto innerKernel = innerName == "unroll" ? c2::BatchedInnerKernel::Unrolled :
                             innerName == "prefetch" ? c2::BatchedInnerKernel::Prefetch128 :
                             c2::BatchedInnerKernel::Baseline;
    const auto setupStart = Clock::now();
    c2::BatchedPlan plan(k, width, seed, refresh, diagonal);
    const auto setupMs = std::chrono::duration<double, std::milli>(Clock::now() - setupStart).count();
    c2::BenchPayload input(plan.n, preferHuge), scratch(plan.n, preferHuge);
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
    auto encode = [&] { c2::batchedTransposeFast(plan, input.p, input.p, scratch.p, outer, innerKernel); };
    for (unsigned i = 0; i != 5; ++i) encode();
    std::vector<double> times(calls), inner(calls), finish(calls);
    for (unsigned sample = 0; sample != calls; ++sample) {
        const auto start = Clock::now();
        if (kernel == "phases") {
            if (width == 8) c2::batchedInnerTransposeDispatch<8>(plan, input.p, scratch.p, innerKernel);
            else c2::batchedInnerTransposeDispatch<16>(plan, input.p, scratch.p, innerKernel);
            const auto split = Clock::now();
            c2::batchedOuterTransposeFast(plan, scratch.p, input.p);
            const auto end = Clock::now();
            inner[sample] = std::chrono::duration<double, std::micro>(split - start).count();
            finish[sample] = std::chrono::duration<double, std::micro>(end - split).count();
            times[sample] = std::chrono::duration<double, std::micro>(end - start).count();
        } else {
            encode();
            times[sample] = std::chrono::duration<double, std::micro>(Clock::now() - start).count();
        }
    }
    std::sort(times.begin(), times.end());
    std::uint64_t hash = 0;
    for (std::size_t i = 0; i != k; ++i) {
        hash = (hash ^ input.p[i].lo) * 0x100000001b3ULL;
        hash = (hash ^ input.p[i].hi) * 0x100000001b3ULL;
    }
    std::cout << "candidate,K,t,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,mode,kernel,inner\n"
              << "fixed-column-regions," << k << ',' << width << ',' << seed << ',' << calls << ','
              << std::fixed << std::setprecision(4) << times[calls / 2] << ',' << times[calls / 10] << ','
              << times[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes() << ','
              << plan.n * sizeof(c2::Block) << ',' << std::hex << hash << std::dec << ','
              << (preferHuge ? "huge-preferred" : "normal") << ',' << mode << ',' << kernel << ',' << innerName << '\n';
    if (kernel == "phases") {
        std::sort(inner.begin(), inner.end());
        std::sort(finish.begin(), finish.end());
        std::cout << "inner_route_us,outer_us,total_us\n" << inner[calls / 2] << ','
                  << finish[calls / 2] << ',' << times[calls / 2] << '\n';
    }
}

int main(int argc, char** argv) {
    try {
        if (!c2::gfni512Available()) return 77;
        if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
        if (argc < 3 || argc > 9)
            throw std::invalid_argument("usage: batched K t [seed=113] [calls=51] [huge|normal] "
                                        "[post|pre|post-d|pre-d] [wide|scalar|phases] [base|unroll|prefetch]");
        const auto k = std::stoull(argv[1]);
        const auto width = std::stoul(argv[2]);
        const auto seed = argc > 3 ? std::stoull(argv[3]) : 113;
        const auto calls = argc > 4 ? std::stoul(argv[4]) : 51;
        const std::string pages = argc > 5 ? argv[5] : "huge";
        const std::string mode = argc > 6 ? argv[6] : "post";
        const std::string kernel = argc > 7 ? argv[7] : "wide";
        const std::string innerName = argc > 8 ? argv[8] : "base";
        if (pages != "huge" && pages != "normal") throw std::invalid_argument("page policy");
        if (mode != "post" && mode != "pre" && mode != "post-d" && mode != "pre-d")
            throw std::invalid_argument("construction mode");
        if (kernel != "wide" && kernel != "scalar" && kernel != "phases")
            throw std::invalid_argument("kernel mode");
        if (innerName != "base" && innerName != "unroll" && innerName != "prefetch")
            throw std::invalid_argument("inner kernel mode");
        if (calls < 1 || calls > 1000000) throw std::invalid_argument("call count");
        benchmark(k, static_cast<unsigned>(width), seed, static_cast<unsigned>(calls), pages == "huge", mode, kernel, innerName);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
