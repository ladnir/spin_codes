#include "Subfield.h"
#include "BenchPayload.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>

namespace c2 = spin::research::char2;
using Clock = std::chrono::steady_clock;
using Element = std::array<std::uint8_t, 16>;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

// Tests only: independent schoolbook multiplication in the declared extension
// basis. Nothing in the encoding hot path calls this general multiplication.
static Element extensionMultiply(const std::uint8_t* left, const std::uint8_t* right, unsigned bytes) {
    Element result{};
    if (bytes == 1) { result[0] = c2::sfMultiply(left[0], right[0]); return result; }
    require(bytes == 4 || bytes == 8 || bytes == 16, "extension width");
    const std::uint8_t constant = bytes == 4 ? 0x3b : bytes == 8 ? 0x57 : 0x63;
    const std::uint8_t linear = bytes == 4 ? 0xd1 : bytes == 8 ? 0x48 : 0xa5;
    const std::uint8_t penultimate = bytes == 4 ? 0xd4 : bytes == 8 ? 0x0f : 0x41;
    std::array<std::uint8_t, 31> product{};
    for (unsigned i = 0; i != bytes; ++i)
        for (unsigned j = 0; j != bytes; ++j)
            product[i + j] ^= c2::sfMultiply(left[i], right[j]);
    for (unsigned degree = 2 * bytes - 1; degree-- > bytes;) {
        const auto coefficient = product[degree];
        product[degree] = 0;
        product[degree - bytes] ^= c2::sfMultiply(coefficient, constant);
        product[degree - bytes + 1] ^= c2::sfMultiply(coefficient, linear);
        product[degree - 1] ^= c2::sfMultiply(coefficient, penultimate);
    }
    std::copy_n(product.begin(), bytes, result.begin());
    return result;
}

static Element dot(const std::vector<std::uint8_t>& left, const std::vector<std::uint8_t>& right,
                    unsigned bytes) {
    require(left.size() == right.size() && left.size() % bytes == 0, "extension dot dimensions");
    Element result{};
    for (std::size_t i = 0; i != left.size(); i += bytes) {
        const auto product = extensionMultiply(left.data() + i, right.data() + i, bytes);
        for (unsigned component = 0; component != bytes; ++component) result[component] ^= product[component];
    }
    return result;
}

static std::vector<std::uint8_t> scaled(const std::vector<std::uint8_t>& input,
                                       const Element& scalar, unsigned bytes) {
    auto result = input;
    for (std::size_t i = 0; i != input.size(); i += bytes) {
        const auto product = extensionMultiply(input.data() + i, scalar.data(), bytes);
        std::copy_n(product.begin(), bytes, result.begin() + i);
    }
    return result;
}

static void samplerChecks() {
    require(c2::sfMultiply(0x57, 0x83) == 0xc1, "AES field product");
    for (unsigned value = 1; value != 256; ++value)
        require(c2::sfMultiply(std::uint8_t(value), c2::sfInverse(std::uint8_t(value))) == 1,
                "AES field inverse");
    // Recover the parameterizing first column and free byte for every nonzero
    // column, with representative free values. The algebraic bijection covers
    // all256 free values; no empirical-distribution argument is being used.
    for (unsigned column = 1; column != 65536; ++column)
        for (const unsigned free : {0, 1, 127, 255}) {
            const auto matrix = c2::sfMatrixFromColumn(std::uint8_t(column),
                std::uint8_t(column >> 8), std::uint8_t(free));
            require((c2::sfMultiply(matrix.a, matrix.d) ^ c2::sfMultiply(matrix.b, matrix.c)) == 1,
                    "SL2 parameterization determinant");
            require((matrix.a ? matrix.b : matrix.d) == free, "SL2 free-byte recovery");
        }
    std::mt19937_64 rng(0x93805e91);
    for (unsigned i = 0; i != 10000; ++i) {
        const auto matrix = c2::sfSampleMatrix(rng);
        require((c2::sfMultiply(matrix.a, matrix.d) ^ c2::sfMultiply(matrix.b, matrix.c)) == 1,
                "SL2 sampled determinant");
        require(c2::sfSampleNonzero(rng) != 0, "diagonal has zero coefficient");
    }
}

static void extensionChecks() {
    std::mt19937_64 rng(0x1029384756);
    for (const unsigned bytes : {1, 4, 8, 16})
        for (unsigned trial = 0; trial != 64; ++trial) {
            Element a{}, b{}, c{}, bPlusC{};
            for (unsigned component = 0; component != bytes; ++component) {
                a[component] = std::uint8_t(rng());
                b[component] = std::uint8_t(rng());
                c[component] = std::uint8_t(rng());
                bPlusC[component] = b[component] ^ c[component];
            }
            const auto ab = extensionMultiply(a.data(), b.data(), bytes);
            const auto ba = extensionMultiply(b.data(), a.data(), bytes);
            require(ab == ba, "extension multiplication commutativity");
            const auto bc = extensionMultiply(b.data(), c.data(), bytes);
            require(extensionMultiply(ab.data(), c.data(), bytes) ==
                    extensionMultiply(a.data(), bc.data(), bytes), "extension multiplication associativity");
            const auto ac = extensionMultiply(a.data(), c.data(), bytes);
            auto sum = ab;
            for (unsigned component = 0; component != bytes; ++component) sum[component] ^= ac[component];
            require(sum == extensionMultiply(a.data(), bPlusC.data(), bytes), "extension distributivity");
            Element coefficient{};
            coefficient[0] = std::uint8_t(rng());
            const auto product = extensionMultiply(a.data(), coefficient.data(), bytes);
            for (unsigned component = 0; component != bytes; ++component)
                require(product[component] == c2::sfMultiply(a[component], coefficient[0]),
                        "subfield action is not componentwise");
        }
}

static void check() {
    samplerChecks();
    extensionChecks();
    std::mt19937_64 rng(20261005);
    for (const unsigned outer : {8, 16})
      for (const unsigned width : {8, 16})
        for (const std::size_t k : {4, 8, 12, 16, 24, 28, 32, 64, 256}) {
            if (k % (outer / 2) || (2 * k) % width) continue;
            for (const unsigned seed : {1, 91}) {
                c2::SubfieldPlan plan(k, width, seed, outer);
                for (unsigned region = 0; region != outer; ++region) {
                    std::vector<std::uint32_t> seen(plan.route.begin() + region * plan.groups,
                                                  plan.route.begin() + (region + 1) * plan.groups);
                    std::sort(seen.begin(), seen.end());
                    for (std::size_t j = 0; j != plan.groups; ++j)
                        require(seen[j] == region * plan.groups + j, "subfield route bijection");
                }
                for (const unsigned bytes : {1, 4, 8, 16}) {
                    std::vector<std::uint8_t> x(k * bytes), v(2 * k * bytes);
                    for (auto& component : x) component = std::uint8_t(rng());
                    for (auto& component : v) component = std::uint8_t(rng());
                    const auto expected = c2::subfieldTransposeOracle(plan, v, bytes);
                    const auto y = c2::subfieldForwardOracle(plan, x, bytes);
                    require(dot(y, v, bytes) == dot(x, expected, bytes), "full-extension adjoint mismatch");
                    Element scalar{};
                    for (unsigned component = 0; component != bytes; ++component)
                        scalar[component] = std::uint8_t(rng());
                    require(c2::subfieldForwardOracle(plan, scaled(x, scalar, bytes), bytes) == scaled(y, scalar, bytes),
                            "full-extension forward linearity mismatch");
                    require(c2::subfieldTransposeOracle(plan, scaled(v, scalar, bytes), bytes) == scaled(expected, scalar, bytes),
                            "full-extension transpose linearity mismatch");
                    if (bytes == 16) {
                        std::vector<c2::Block> input(2 * k), fast(k + 2), scratch(2 * k + 2);
                        std::memcpy(input.data(), v.data(), v.size());
                        const c2::Block guard{0x173859012934efa1ULL, 0x831afb36917502ceULL};
                        for (const auto inner : {c2::SubfieldInnerKernel::Baseline, c2::SubfieldInnerKernel::Prefetch128})
                          for (const unsigned offset : {0, 1}) {
                            std::fill(fast.begin(), fast.end(), guard);
                            std::fill(scratch.begin(), scratch.end(), guard);
                            c2::subfieldTransposeFast128(plan, input.data(), fast.data() + offset, scratch.data() + offset, inner);
                            require(std::memcmp(fast.data() + offset, expected.data(), expected.size()) == 0,
                                    "subfield SIMD mismatch");
                            for (std::size_t i = 0; i != offset; ++i)
                                require(fast[i] == guard && scratch[i] == guard, "subfield prefix guard");
                            require(fast[k + offset] == guard && scratch[2 * k + offset] == guard,
                                    "subfield suffix guard");
                        }
                        for (const auto inner : {c2::SubfieldInnerKernel::Baseline, c2::SubfieldInnerKernel::Prefetch128}) {
                            auto alias = input;
                            c2::subfieldTransposeFast128(plan, alias.data(), alias.data(), scratch.data(), inner);
                            require(std::memcmp(alias.data(), expected.data(), expected.size()) == 0,
                                    "subfield alias prefix");
                            require(std::equal(input.begin() + k, input.end(), alias.begin() + k),
                                    "subfield alias suffix changed");
                        }
                    }
                }
            }
        }
    for (const auto bad : {0U, 1U, 7U, 9U}) {
        bool rejected = false;
        try { c2::SubfieldPlan invalid(bad, 16, 1, 16); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid subfield natural size accepted");
    }
    std::cout << "subfield: AES arithmetic, SL2 sampling, GF256-extension adjoints/linearity, "
                 "GRS8/16, t8/16, baseline/prefetch SIMD, routes, guards, aliases and natural sizes pass\n";
}

static void benchmark(std::size_t k, unsigned width, unsigned outer, std::uint64_t seed,
                      unsigned calls, bool preferHuge, bool phases, const std::string& innerName) {
    const auto innerKernel = innerName == "prefetch" ? c2::SubfieldInnerKernel::Prefetch128
                                                    : c2::SubfieldInnerKernel::Baseline;
    const auto setupStart = Clock::now();
    c2::SubfieldPlan plan(k, width, seed, outer);
    const auto setupMs = std::chrono::duration<double, std::milli>(Clock::now() - setupStart).count();
    c2::BenchPayload input(plan.n, preferHuge), scratch(plan.n, preferHuge);
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
    auto encode = [&] { c2::subfieldTransposeFast128(plan, input.p, input.p, scratch.p, innerKernel); };
    for (unsigned i = 0; i != 5; ++i) encode();
    std::vector<double> total(calls), inner(calls), finish(calls);
    for (unsigned sample = 0; sample != calls; ++sample) {
        const auto start = Clock::now();
        if (phases) {
            c2::subfieldInnerTransposeDispatch(plan, input.p, scratch.p, innerKernel);
            const auto split = Clock::now();
            c2::subfieldOuterTransposeDispatch(plan, scratch.p, input.p);
            const auto end = Clock::now();
            inner[sample] = std::chrono::duration<double, std::micro>(split - start).count();
            finish[sample] = std::chrono::duration<double, std::micro>(end - split).count();
            total[sample] = std::chrono::duration<double, std::micro>(end - start).count();
        } else {
            encode();
            total[sample] = std::chrono::duration<double, std::micro>(Clock::now() - start).count();
        }
    }
    std::sort(total.begin(), total.end());
    std::uint64_t hash = 0;
    for (std::size_t i = 0; i != k; ++i) {
        hash = (hash ^ input.p[i].lo) * 0x100000001b3ULL;
        hash = (hash ^ input.p[i].hi) * 0x100000001b3ULL;
    }
    std::cout << "candidate,K,t,outer,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,kernel,inner\n"
              << "gf256-sl2-pre-d," << k << ',' << width << ',' << outer << ',' << seed << ',' << calls << ','
              << std::fixed << std::setprecision(4) << total[calls / 2] << ',' << total[calls / 10] << ','
              << total[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes() << ','
              << plan.n * sizeof(c2::Block) << ',' << std::hex << hash << std::dec << ','
              << (preferHuge ? "huge-preferred" : "normal") << ',' << (phases ? "phases" : "full")
              << ',' << innerName << '\n';
    if (phases) {
        std::sort(inner.begin(), inner.end());
        std::sort(finish.begin(), finish.end());
        std::cout << "inner_route_us,outer_us,total_us\n" << inner[calls / 2] << ','
                  << finish[calls / 2] << ',' << total[calls / 2] << '\n';
    }
}

int main(int argc, char** argv) {
    try {
        if (!c2::gfni512Available()) return 77;
        if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
        if (argc < 3 || argc > 9)
            throw std::invalid_argument("usage: subfield K t [outer=16] [seed=113] [calls=51] [huge|normal] [full|phases] [base|prefetch]");
        const auto k = std::stoull(argv[1]);
        const auto width = std::stoul(argv[2]);
        const auto outer = argc > 3 ? std::stoul(argv[3]) : 16;
        const auto seed = argc > 4 ? std::stoull(argv[4]) : 113;
        const auto calls = argc > 5 ? std::stoul(argv[5]) : 51;
        const std::string pages = argc > 6 ? argv[6] : "huge";
        const std::string kernel = argc > 7 ? argv[7] : "full";
        const std::string innerName = argc > 8 ? argv[8] : "base";
        if (calls < 1 || calls > 1000000 || (pages != "huge" && pages != "normal") ||
            (kernel != "full" && kernel != "phases") ||
            (innerName != "base" && innerName != "prefetch")) throw std::invalid_argument("benchmark options");
        benchmark(k, static_cast<unsigned>(width), static_cast<unsigned>(outer), seed,
                  static_cast<unsigned>(calls), pages == "huge", kernel == "phases", innerName);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
