#include "SubfieldToeplitz.h"
#include "BenchPayload.h"
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace c2 = spin::research::char2;
using Clock = std::chrono::steady_clock;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static c2::Subfield16Kernel parseKernel(const std::string& name) {
    if (name == "scatter") return c2::Subfield16Kernel::Scatter;
    if (name == "region-nt") return c2::Subfield16Kernel::RegionNonTemporal;
    if (name == "region-nt-all") return c2::Subfield16Kernel::RegionAllNonTemporal;
    throw std::invalid_argument("Toeplitz kernel must be scatter, region-nt or region-nt-all");
}

static c2::SubfieldToeplitzOuter parseOuter(const std::string& name) {
    if (name == "kara") return c2::SubfieldToeplitzOuter::Karatsuba;
    if (name == "dense") return c2::SubfieldToeplitzOuter::Dense;
    if (name == "kara-packed") return c2::SubfieldToeplitzOuter::KaratsubaPacked;
    if (name == "kara-dense2") return c2::SubfieldToeplitzOuter::KaratsubaDense2;
    throw std::invalid_argument("Toeplitz outer must be kara, dense, kara-packed or kara-dense2");
}

static std::uint8_t dot(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b) {
    require(a.size() == b.size(), "GF256 dot dimensions");
    std::uint8_t result = 0;
    for (std::size_t i = 0; i != a.size(); ++i) result ^= c2::sfMultiply(a[i], b[i]);
    return result;
}

// Independent of Gaussian elimination: all six determinant terms, with
// signs omitted because the field has characteristic two.
static std::uint8_t determinant(const c2::SubfieldToeplitzMatrix& matrix) {
    std::array<unsigned, 3> permutation{0, 1, 2};
    std::uint8_t result = 0;
    do {
        std::uint8_t term = 1;
        for (unsigned row = 0; row != 3; ++row)
            term = c2::sfMultiply(term, matrix.entries[4 * row + permutation[row]]);
        result ^= term;
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

static void matrixChecks() {
    std::mt19937_64 rng(0x82fd390127ULL);
    std::array<bool, 256> inSubfield{};
    for (unsigned j = 0; j != 16; ++j) inSubfield[c2::sf16Embed(std::uint8_t(j))] = true;
    for (unsigned trial = 0; trial != 512; ++trial) {
        std::array<std::uint8_t, 9> abstract;
        c2::SubfieldToeplitzMatrix embedded;
        for (unsigned row = 0; row != 3; ++row)
            for (unsigned column = 0; column != 3; ++column) {
                const auto coefficient = std::uint8_t(rng() & 15U);
                abstract[3 * row + column] = coefficient;
                embedded.entries[4 * row + column] = c2::sf16Embed(coefficient);
            }
        require(c2::subfieldToeplitzInvertible(abstract) == (determinant(embedded) != 0),
                "GL3 elimination disagrees with determinant");
        const auto matrix = c2::subfieldToeplitzSampleMatrix(rng);
        require(determinant(matrix) != 0, "sampled GL3 matrix is singular");
        for (unsigned row = 0; row != 3; ++row) {
            require(matrix.entries[4 * row + 3] == 0, "GL3 padding is not zero");
            for (unsigned column = 0; column != 3; ++column)
                require(inSubfield[matrix.entries[4 * row + column]], "GL3 coefficient outside embedded GF16");
        }
        for (unsigned basis = 0; basis != 3; ++basis) {
            const auto zero = _mm_setzero_si128(), one = c2::reg(c2::Block{0x57, 0});
            const auto value = c2::subfieldToeplitzRefresh(matrix, basis == 0 ? one : zero,
                basis == 1 ? one : zero, basis == 2 ? one : zero);
            std::array<c2::Block, 4> actual;
            _mm512_storeu_si512(actual.data(), value);
            for (unsigned j = 0; j != 4; ++j)
                require(actual[j] == c2::Block{c2::sfMultiply(0x57, matrix.entries[4 * basis + j]), 0},
                        "packed GL3 transpose product mismatch");
        }
    }
}

template<unsigned Size>
static void outerRecursionChecks(std::mt19937_64& rng) {
    constexpr auto leaves = c2::subfieldToeplitzLeaves<Size>;
    std::array<std::array<std::uint8_t, 2 * Size - 1>, 4> diagonals;
    std::array<std::uint8_t, 4 * leaves> coefficients{};
    alignas(64) c2::Block inputs[Size][4], outputs[Size][4];
    __m512i input[Size], output[Size];
    // Zero, identity, asymmetric shift and random diagonals. Independent
    // matrices in the four lanes catch erroneous coefficient broadcasts.
    for (unsigned trial = 0; trial != 12; ++trial) {
        for (unsigned lane = 0; lane != 4; ++lane) {
            for (unsigned j = 0; j != 2 * Size - 1; ++j) {
                auto abstract = std::uint8_t(rng() & 15U);
                if (trial == 0) abstract = 0;
                if (trial == 1) abstract = std::uint8_t(j == Size - 1);
                if (trial == 2) abstract = std::uint8_t(j == Size - 1 + (Size > 1));
                diagonals[lane][j] = c2::sf16Embed(abstract);
            }
            c2::subfieldToeplitzPrepare<Size>(diagonals[lane].data(), coefficients.data() + lane, 4);
        }
        for (unsigned j = 0; j != Size; ++j) {
            for (unsigned lane = 0; lane != 4; ++lane) inputs[j][lane] = {rng(), rng()};
            input[j] = _mm512_load_si512(inputs[j]);
        }
        for (unsigned kernel = 0; kernel != 3; ++kernel) {
        if (kernel == 0) c2::subfieldToeplitzApply4<Size>(coefficients.data(), input, output);
        else if (kernel == 1) c2::subfieldToeplitzApplyAlternative4<Size, false>(coefficients.data(), input, output);
        else c2::subfieldToeplitzApplyAlternative4<Size, true>(coefficients.data(), input, output);
        for (unsigned j = 0; j != Size; ++j) _mm512_store_si512(outputs[j], output[j]);
        for (unsigned row = 0; row != Size; ++row)
            for (unsigned lane = 0; lane != 4; ++lane) {
                std::array<std::uint8_t, 16> expected{};
                for (unsigned column = 0; column != Size; ++column) {
                    std::array<std::uint8_t, 16> value;
                    std::memcpy(value.data(), &inputs[column][lane], 16);
                    for (unsigned byte = 0; byte != 16; ++byte)
                        expected[byte] ^= c2::sfMultiply(value[byte], diagonals[lane][Size - 1 + row - column]);
                }
                require(std::memcmp(&outputs[row][lane], expected.data(), 16) == 0,
                        "Toeplitz recursion disagrees with independent dense multiplication");
            }
        }
    }
}

template<bool Diagonal>
static void innerBasisChecks() {
    c2::SubfieldToeplitzPlan plan(32, 13, Diagonal);
    for (unsigned basis = 0; basis != 19; ++basis) {
        std::array<c2::Block, 16> input{}, output{};
        c2::SubfieldToeplitzState state;
        if (basis < 16) input[basis].lo = 1;
        else {
            const auto one = c2::reg(c2::Block{1, 0});
            if (basis == 16) state.r0 = one;
            if (basis == 17) state.r1 = one;
            if (basis == 18) state.r2 = one;
        }
        const auto r0 = _mm512_broadcast_i32x4(state.r0), r1 = _mm512_broadcast_i32x4(state.r1);
        const auto r2 = _mm512_broadcast_i32x4(state.r2);
        auto sum0 = _mm512_setzero_si512(), sum1 = _mm512_setzero_si512(), sum2 = _mm512_setzero_si512();
        auto packed = _mm512_setzero_si512();
        if constexpr (Diagonal) packed = _mm512_broadcast_i32x4(_mm_loadu_si128(
            reinterpret_cast<const __m128i*>(plan.routedDiagonal.data())));
        c2::subfieldToeplitzInnerBatch<0, false, Diagonal>(plan, input.data(), output.data(), 0, 0, r0, r1, r2, sum0, sum1, sum2, packed);
        c2::subfieldToeplitzInnerBatch<1, false, Diagonal>(plan, input.data(), output.data(), 0, 0, r0, r1, r2, sum0, sum1, sum2, packed);
        c2::subfieldToeplitzInnerBatch<2, false, Diagonal>(plan, input.data(), output.data(), 0, 0, r0, r1, r2, sum0, sum1, sum2, packed);
        c2::subfieldToeplitzInnerBatch<3, false, Diagonal>(plan, input.data(), output.data(), 0, 0, r0, r1, r2, sum0, sum1, sum2, packed);
        for (unsigned j = 0; j != 16; ++j) {
            auto expected = std::uint8_t(basis == j);
            if (basis >= 16) expected ^= c2::sf16Embed(c2::sf16Monomial(basis - 16, j));
            if constexpr (Diagonal) expected = c2::sfMultiply(expected, plan.routedDiagonal[j]);
            require(output[j] == c2::Block{expected, 0}, "Toeplitz inner basis evaluation differs from monomial map");
        }
        const std::array<c2::Block, 3> moments{c2::block(c2::sfSumRecords4(sum0)),
            c2::block(c2::sfSumRecords4(sum1)), c2::block(c2::sfSumRecords4(sum2))};
        for (unsigned power = 0; power != 3; ++power) {
            const auto expected = basis < 16 ? c2::sf16Embed(c2::sf16Monomial(power, basis)) : 0;
            require(moments[power] == c2::Block{static_cast<std::uint64_t>(expected), 0},
                    "Toeplitz inner transpose basis mismatch");
        }
    }
}

static void check() {
    matrixChecks();
    innerBasisChecks<false>();
    innerBasisChecks<true>();
    std::mt19937_64 rng(20261008);
    outerRecursionChecks<1>(rng);
    outerRecursionChecks<2>(rng);
    outerRecursionChecks<4>(rng);
    outerRecursionChecks<8>(rng);
    outerRecursionChecks<16>(rng);
    outerRecursionChecks<32>(rng);
    const c2::Block guard{0x32bcf571ed348a09ULL, 0x92017bed523ac816ULL};
    for (std::size_t k : {32, 64, 96, 128, 160, 224, 256, 480, 512, 544, 640, 1024})
      for (std::uint64_t seed : {1, 91})
       for (bool diagonal : {false, true}) {
        c2::SubfieldToeplitzPlan plan(k, seed, diagonal);
        require(plan.routedDiagonal.size() == (diagonal ? plan.n : 0), "Toeplitz diagonal storage disagrees with profile");
        std::array<bool, 256> inSubfield{};
        for (unsigned j = 0; j != 16; ++j) inSubfield[c2::sf16Embed(std::uint8_t(j))] = true;
        for (const auto& diagonal : plan.outerDiagonals)
            for (auto value : diagonal) require(inSubfield[value], "Toeplitz outer coefficient outside embedded GF16");
        for (std::size_t i = 0; i != plan.n; ++i) {
            require(plan.route[i] < plan.n && plan.route[i] / plan.groups == i / plan.groups,
                    "Toeplitz route left its region");
            require(plan.inverseRoute[plan.route[i]] == i, "Toeplitz route inverse mismatch");
            if (diagonal) require(plan.routedDiagonal[i] != 0 && inSubfield[plan.routedDiagonal[i]],
                                  "Toeplitz input diagonal outside GF16 star");
        }
        for (unsigned j = 0; j != 16; ++j)
            for (unsigned power = 1; power != 3; ++power)
                require(plan.coefficients[j / 4][power - 1].bytes[2 * (j % 4)] ==
                        0x0101010101010101ULL * c2::sf16Embed(c2::sf16Monomial(power, j)),
                        "Toeplitz inner map differs from monomial specification");
        for (unsigned bytes : {1, 4, 8, 16}) {
            std::vector<std::uint8_t> message(k * bytes), inputBytes(plan.n * bytes);
            for (auto& value : message) value = std::uint8_t(rng());
            for (auto& value : inputBytes) value = std::uint8_t(rng());
            const auto encoded = c2::subfieldToeplitzForwardOracle(plan, message, bytes);
            const auto expected = c2::subfieldToeplitzTransposeOracle(plan, inputBytes, bytes);
            require(dot(encoded, inputBytes) == dot(message, expected), "Toeplitz scalar adjoint mismatch");
            const auto scalar = std::uint8_t(rng());
            auto scaledMessage = message, scaledEncoded = encoded;
            for (auto& value : scaledMessage) value = c2::sfMultiply(value, scalar);
            for (auto& value : scaledEncoded) value = c2::sfMultiply(value, scalar);
            require(c2::subfieldToeplitzForwardOracle(plan, scaledMessage, bytes) == scaledEncoded,
                    "Toeplitz GF256 scalar extension mismatch");
            if (bytes != 16) continue;
            c2::BenchPayload input(plan.n + 4, false), output(plan.k + 4, false), planes(plan.n + 4, false);
            for (auto kernel : {c2::Subfield16Kernel::Scatter, c2::Subfield16Kernel::RegionNonTemporal,
                                c2::Subfield16Kernel::RegionAllNonTemporal})
              for (auto outer : {c2::SubfieldToeplitzOuter::Karatsuba, c2::SubfieldToeplitzOuter::Dense,
                                 c2::SubfieldToeplitzOuter::KaratsubaPacked, c2::SubfieldToeplitzOuter::KaratsubaDense2}) {
                const auto temporarySize = c2::subfieldToeplitzTemporaryRecords(plan, kernel);
                c2::BenchPayload temporary(temporarySize + 4, false);
                for (unsigned offset : {0, 1, 2, 3}) {
                    std::fill(input.p, input.p + plan.n + 4, guard);
                    std::fill(output.p, output.p + plan.k + 4, guard);
                    std::fill(planes.p, planes.p + plan.n + 4, guard);
                    std::fill(temporary.p, temporary.p + temporarySize + 4, guard);
                    std::memcpy(input.p + offset, inputBytes.data(), inputBytes.size());
                    c2::subfieldToeplitzTransposeFast128(plan, input.p + offset, output.p + offset,
                        planes.p + offset, temporary.p + offset, kernel, outer);
                    require(std::memcmp(output.p + offset, expected.data(), expected.size()) == 0,
                            "Toeplitz SIMD transpose differs from scalar oracle");
                    for (unsigned i = 0; i != offset; ++i)
                        require(input.p[i] == guard && output.p[i] == guard && planes.p[i] == guard && temporary.p[i] == guard,
                                "Toeplitz prefix guard changed");
                    require(input.p[plan.n + offset] == guard && output.p[plan.k + offset] == guard &&
                            planes.p[plan.n + offset] == guard && temporary.p[temporarySize + offset] == guard,
                            "Toeplitz suffix guard changed");
                    std::vector<std::uint8_t> fastBytes(expected.size());
                    std::memcpy(fastBytes.data(), output.p + offset, expected.size());
                    require(dot(encoded, inputBytes) == dot(message, fastBytes), "Toeplitz SIMD adjoint mismatch");
                    require(std::memcmp(input.p + offset, inputBytes.data(), inputBytes.size()) == 0,
                            "Toeplitz out-of-place input changed");
                    c2::subfieldToeplitzTransposeFast128(plan, input.p + offset, input.p + offset,
                        planes.p + offset, temporary.p + offset, kernel, outer);
                    require(std::memcmp(input.p + offset, expected.data(), expected.size()) == 0,
                            "Toeplitz in-place prefix mismatch");
                    require(std::memcmp(input.p + offset + k, inputBytes.data() + k * 16, k * 16) == 0,
                            "Toeplitz in-place suffix changed");
                }
            }
        }
        if (plan.regionAligned()) {
            c2::BenchPayload input(plan.n, false), planes(plan.n, false);
            bool rejected = false;
            try { c2::subfieldToeplitzInnerDispatch(plan, input.p, planes.p, nullptr, c2::Subfield16Kernel::RegionAllNonTemporal); }
            catch (const std::invalid_argument&) { rejected = true; }
            require(rejected, "Toeplitz region schedule accepted missing temporary");
        }
    }
    for (std::size_t invalid : {0, 1, 31, 33}) {
        bool rejected = false;
        try { c2::SubfieldToeplitzPlan bad(invalid, 1); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid Toeplitz natural size accepted");
    }
    if (sizeof(std::size_t) >= 8) {
        bool rejected = false;
        try { c2::SubfieldToeplitzPlan bad(std::size_t(1) << 31, 1); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "Toeplitz uint32 route overflow accepted");
    }
    std::cout << "GF16 Toeplitz[64,32] s3: GL3 determinant/sampler and packed update, dense/recursive outer, "
                 "inner basis maps, scalar/SIMD adjoints, diagonal/no-diagonal, region schedules, natural sizes, all alignments, guards and aliases pass\n";
}

// Optional scale regression, without timing. The dense SIMD outer has no
// dependency on the prepared recursive leaves; all supported large geometries
// must agree for both routing schedules and preserve the in-place suffix.
static void checkLarge() {
    std::mt19937_64 rng(713);
    for (std::size_t k : {std::size_t(1) << 16, std::size_t(1) << 18, std::size_t(1) << 20}) {
        c2::SubfieldToeplitzPlan plan(k, 91);
        c2::BenchPayload input(plan.n, false), planes(plan.n, false), temporary(plan.groups, false);
        c2::BenchPayload reference(plan.k, false), actual(plan.k, false);
        for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
        c2::subfieldToeplitzTransposeFast128(plan, input.p, reference.p, planes.p, nullptr,
            c2::Subfield16Kernel::Scatter, c2::SubfieldToeplitzOuter::Dense);
        for (auto kernel : {c2::Subfield16Kernel::Scatter, c2::Subfield16Kernel::RegionAllNonTemporal})
          for (auto outer : {c2::SubfieldToeplitzOuter::Karatsuba, c2::SubfieldToeplitzOuter::KaratsubaPacked,
                             c2::SubfieldToeplitzOuter::KaratsubaDense2}) {
            c2::subfieldToeplitzTransposeFast128(plan, input.p, actual.p, planes.p, temporary.p, kernel, outer);
            require(std::equal(actual.p, actual.p + plan.k, reference.p), "large Toeplitz dense/Karatsuba mismatch");
        }
        std::vector<c2::Block> prefix(input.p, input.p + plan.k), suffix(input.p + plan.k, input.p + plan.n);
        for (auto outer : {c2::SubfieldToeplitzOuter::Karatsuba, c2::SubfieldToeplitzOuter::KaratsubaPacked,
                           c2::SubfieldToeplitzOuter::KaratsubaDense2}) {
            std::copy(prefix.begin(), prefix.end(), input.p);
            c2::subfieldToeplitzTransposeFast128(plan, input.p, input.p, planes.p, temporary.p,
                c2::Subfield16Kernel::RegionAllNonTemporal, outer);
            require(std::equal(input.p, input.p + plan.k, reference.p), "large Toeplitz in-place mismatch");
            require(std::equal(input.p + plan.k, input.p + plan.n, suffix.begin()), "large Toeplitz suffix changed");
        }
    }
    std::cout << "Toeplitz K16/K18/K20 dense/recursive, scatter/region and in-place checks pass\n";
}

static double micros(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::micro>(end - begin).count();
}

static void benchmark(std::size_t k, std::uint64_t seed, unsigned calls, bool huge,
                      bool phases, const std::string& kernelName, const std::string& outerName, bool diagonal) {
    const auto kernel = parseKernel(kernelName);
    const auto outerKernel = parseOuter(outerName);
    const auto setupStart = Clock::now();
    c2::SubfieldToeplitzPlan plan(k, seed, diagonal);
    const auto setupMs = micros(setupStart, Clock::now()) / 1000;
    c2::BenchPayload input(plan.n, huge), planes(plan.n, huge);
    const auto temporarySize = c2::subfieldToeplitzTemporaryRecords(plan, kernel);
    std::unique_ptr<c2::BenchPayload> temporary;
    if (temporarySize) temporary = std::make_unique<c2::BenchPayload>(temporarySize, huge);
    auto temporaryPointer = temporary ? temporary->p : nullptr;
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
    auto encode = [&] { c2::subfieldToeplitzTransposeFast128(plan, input.p, input.p, planes.p, temporaryPointer, kernel, outerKernel); };
    for (unsigned i = 0; i != 5; ++i) encode();
    std::vector<double> total(calls), inner(calls), outer(calls);
    for (unsigned sample = 0; sample != calls; ++sample) {
        const auto begin = Clock::now();
        if (phases) {
            c2::subfieldToeplitzInnerDispatch(plan, input.p, planes.p, temporaryPointer, kernel);
            const auto split = Clock::now();
            c2::subfieldToeplitzOuterDispatch(plan, planes.p, input.p, kernel, outerKernel);
            const auto end = Clock::now();
            inner[sample] = micros(begin, split);
            outer[sample] = micros(split, end);
            total[sample] = micros(begin, end);
        } else {
            encode();
            total[sample] = micros(begin, Clock::now());
        }
    }
    std::sort(total.begin(), total.end());
    std::uint64_t checksum = 0;
    for (std::size_t i = 0; i != plan.k; ++i) {
        checksum = (checksum ^ input.p[i].lo) * 0x100000001b3ULL;
        checksum = (checksum ^ input.p[i].hi) * 0x100000001b3ULL;
    }
    std::cout << "candidate,K,t,s,outer,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,mode,kernel,outer_kernel,diagonal\n"
              << "gf16-toeplitz64-gl3-pre," << k << ",16,3,64," << seed << ',' << calls << ',' << std::fixed
              << std::setprecision(4) << total[calls / 2] << ',' << total[calls / 10] << ','
              << total[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes() << ','
              << (plan.n + temporarySize) * sizeof(c2::Block) << ',' << std::hex << checksum << std::dec << ','
              << (huge ? "huge-preferred" : "normal") << ',' << (phases ? "phases" : "full") << ','
              << kernelName << ',' << outerName << ',' << (diagonal ? "diag" : "nod") << '\n';
    if (phases) {
        std::sort(inner.begin(), inner.end());
        std::sort(outer.begin(), outer.end());
        std::cout << "inner_route_us,outer_us,total_us\n" << inner[calls / 2] << ','
                  << outer[calls / 2] << ',' << total[calls / 2] << '\n';
    }
}

// Compare full encodes in one process, with identical plan and payload addresses.
// Each small batch starts from the same input; restores and checksums are outside
// timed calls. The six orders balance each variant's position across passes.
static void benchmarkCompare(std::size_t k, std::uint64_t seed, unsigned calls, bool huge,
                             const std::string& kernelName) {
    const auto kernel = parseKernel(kernelName);
    constexpr std::array<c2::SubfieldToeplitzOuter, 3> variants{
        c2::SubfieldToeplitzOuter::Karatsuba, c2::SubfieldToeplitzOuter::KaratsubaPacked,
        c2::SubfieldToeplitzOuter::KaratsubaDense2};
    constexpr std::array<const char*, 3> names{"kara", "kara-packed", "kara-dense2"};
    constexpr std::array<std::array<unsigned, 3>, 6> orders{{
        {{0, 1, 2}}, {{2, 1, 0}}, {{1, 2, 0}},
        {{0, 2, 1}}, {{2, 0, 1}}, {{1, 0, 2}}}};
    constexpr unsigned batchCalls = 5, warmupCalls = 5;
    const auto setupStart = Clock::now();
    c2::SubfieldToeplitzPlan plan(k, seed);
    const auto setupMs = micros(setupStart, Clock::now()) / 1000;
    c2::BenchPayload input(plan.n, huge), planes(plan.n, huge), initial(plan.n, huge);
    const auto temporarySize = c2::subfieldToeplitzTemporaryRecords(plan, kernel);
    std::unique_ptr<c2::BenchPayload> temporary;
    if (temporarySize) temporary = std::make_unique<c2::BenchPayload>(temporarySize, huge);
    auto temporaryPointer = temporary ? temporary->p : nullptr;
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) initial.p[i] = {rng(), rng()};
    auto restore = [&] { std::memcpy(input.p, initial.p, plan.n * sizeof(c2::Block)); };
    auto encode = [&](c2::SubfieldToeplitzOuter outer) {
        c2::subfieldToeplitzTransposeFast128(plan, input.p, input.p, planes.p, temporaryPointer, kernel, outer);
    };
    auto checksum = [&] {
        std::uint64_t result = 0;
        for (std::size_t i = 0; i != plan.k; ++i) {
            result = (result ^ input.p[i].lo) * 0x100000001b3ULL;
            result = (result ^ input.p[i].hi) * 0x100000001b3ULL;
        }
        return result;
    };
    std::array<std::vector<double>, 3> samples;
    for (auto& values : samples) values.resize(calls);
    std::array<unsigned, 3> measured{};
    std::array<std::uint64_t, 3> checksums{};
    for (auto outer : variants) {
        restore();
        for (unsigned i = 0; i != warmupCalls; ++i) encode(outer);
    }
    unsigned passes = 0;
    for (unsigned completed = 0; completed != calls; ++passes) {
        const auto count = std::min(batchCalls, calls - completed);
        std::array<std::uint64_t, 3> batchChecksums{};
        for (auto variant : orders[passes % orders.size()]) {
            restore();
            const auto outer = variants[variant];
            for (unsigned i = 0; i != count; ++i) {
                const auto begin = Clock::now();
                encode(outer);
                const auto end = Clock::now();
                samples[variant][measured[variant]++] = micros(begin, end);
            }
            batchChecksums[variant] = checksum();
            // Include every batch's result in the one reported checksum.
            checksums[variant] = (checksums[variant] ^ batchChecksums[variant]) * 0x100000001b3ULL;
        }
        require(batchChecksums[0] == batchChecksums[1] && batchChecksums[0] == batchChecksums[2],
                "Toeplitz comparison batch checksums disagree");
        completed += count;
    }
    require(measured[0] == calls && measured[1] == calls && measured[2] == calls,
            "Toeplitz comparison call counts disagree");
    require(checksums[0] == checksums[1] && checksums[0] == checksums[2],
            "Toeplitz comparison checksums disagree");
    std::cout << "candidate,K,t,s,outer,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,mode,kernel,outer_kernel,diagonal,batch_calls,batches,warmup_calls,input_address,planes_address,temporary_address\n";
    for (unsigned variant = 0; variant != variants.size(); ++variant) {
        auto& total = samples[variant];
        std::sort(total.begin(), total.end());
        std::cout << "gf16-toeplitz64-gl3-pre," << k << ",16,3,64," << seed << ',' << measured[variant] << ','
                  << std::fixed << std::setprecision(4) << total[calls / 2] << ',' << total[calls / 10] << ','
                  << total[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes() << ','
                  << (plan.n + temporarySize) * sizeof(c2::Block) << ',' << std::hex << checksums[variant]
                  << std::dec << ',' << (huge ? "huge-preferred" : "normal") << ",compare-full," << kernelName
                  << ',' << names[variant] << ",nod," << batchCalls << ',' << passes << ',' << warmupCalls << ','
                  << static_cast<const void*>(input.p) << ',' << static_cast<const void*>(planes.p) << ','
                  << static_cast<const void*>(temporaryPointer) << '\n';
    }
}

int main(int argc, char** argv) {
    try {
        if (!c2::gfni512Available()) return 77;
        if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
        if (argc == 2 && std::string(argv[1]) == "--check-large") { checkLarge(); return 0; }
        if (argc > 1 && std::string(argv[1]) == "--compare") {
            if (argc < 3 || argc > 7)
                throw std::invalid_argument("usage: toeplitz --compare K [seed=113] [calls=101] [huge|normal] [scatter|region-nt|region-nt-all]");
            const auto k = std::stoull(argv[2]);
            const auto seed = argc > 3 ? std::stoull(argv[3]) : 113;
            const auto calls = argc > 4 ? std::stoul(argv[4]) : 101;
            const std::string pages = argc > 5 ? argv[5] : "huge";
            const std::string kernel = argc > 6 ? argv[6] : "scatter";
            if (calls < 1 || calls > 1000000 || (pages != "huge" && pages != "normal"))
                throw std::invalid_argument("Toeplitz comparison options");
            benchmarkCompare(k, seed, static_cast<unsigned>(calls), pages == "huge", kernel);
            return 0;
        }
        if (argc < 2 || argc > 9)
            throw std::invalid_argument("usage: toeplitz K [seed=113] [calls=51] [huge|normal] [full|phases] [scatter|region-nt|region-nt-all] [kara|dense|kara-packed|kara-dense2] [nod|diag]");
        const auto k = std::stoull(argv[1]);
        const auto seed = argc > 2 ? std::stoull(argv[2]) : 113;
        const auto calls = argc > 3 ? std::stoul(argv[3]) : 51;
        const std::string pages = argc > 4 ? argv[4] : "huge";
        const std::string mode = argc > 5 ? argv[5] : "full";
        const std::string kernel = argc > 6 ? argv[6] : "scatter";
        const std::string outer = argc > 7 ? argv[7] : "kara";
        const std::string diagonal = argc > 8 ? argv[8] : "nod";
        if (calls < 1 || calls > 1000000 || (pages != "huge" && pages != "normal") ||
            (mode != "full" && mode != "phases") || (diagonal != "nod" && diagonal != "diag"))
            throw std::invalid_argument("Toeplitz benchmark options");
        benchmark(k, seed, static_cast<unsigned>(calls), pages == "huge", mode == "phases", kernel, outer, diagonal == "diag");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
