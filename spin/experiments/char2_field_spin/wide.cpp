#include "SubfieldWide.h"
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
    throw std::invalid_argument("wide kernel must be scatter, region-nt or region-nt-all");
}

static c2::SubfieldWideMapKernel parseMap(const std::string& name) {
    if (name == "direct") return c2::SubfieldWideMapKernel::Direct;
    if (name == "fft") return c2::SubfieldWideMapKernel::AdditiveFft;
    throw std::invalid_argument("wide map must be direct or fft");
}

static std::uint8_t dot(const std::vector<std::uint8_t>& a, const std::vector<std::uint8_t>& b) {
    require(a.size() == b.size(), "GF256 dot dimensions");
    std::uint8_t result = 0;
    for (std::size_t i = 0; i != a.size(); ++i) result ^= c2::sfMultiply(a[i], b[i]);
    return result;
}

// Independent of the sampler's Gaussian elimination: sum all 720 determinant
// terms, whose signs disappear in characteristic two.
static std::uint8_t determinant(const c2::SubfieldWideMatrix& matrix) {
    std::array<unsigned, 6> permutation{0, 1, 2, 3, 4, 5};
    std::uint8_t result = 0;
    do {
        std::uint8_t term = 1;
        for (unsigned row = 0; row != 6; ++row)
            term = c2::sfMultiply(term, matrix.entries[8 * row + permutation[row]]);
        result ^= term;
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

static void matrixChecks() {
    std::mt19937_64 rng(0x82fd390127ULL);
    std::array<bool, 256> inSubfield{};
    for (unsigned j = 0; j != 16; ++j) inSubfield[c2::sf16Embed(std::uint8_t(j))] = true;
    for (unsigned trial = 0; trial != 128; ++trial) {
        std::array<std::uint8_t, 36> abstract;
        c2::SubfieldWideMatrix embedded;
        for (unsigned row = 0; row != 6; ++row)
            for (unsigned column = 0; column != 6; ++column) {
                const auto coefficient = std::uint8_t(rng() & 15U);
                abstract[6 * row + column] = coefficient;
                embedded.entries[8 * row + column] = c2::sf16Embed(coefficient);
            }
        require(c2::subfieldWideInvertible(abstract) == (determinant(embedded) != 0),
                "GL6 elimination disagrees with independent determinant");
        const auto matrix = c2::subfieldWideSampleMatrix(rng);
        require(determinant(matrix) != 0, "sampled GL6 matrix is singular");
        for (unsigned row = 0; row != 6; ++row) {
            require(matrix.entries[8 * row + 6] == 0 && matrix.entries[8 * row + 7] == 0,
                    "GL6 matrix row padding is not zero");
            for (unsigned column = 0; column != 6; ++column)
                require(inSubfield[matrix.entries[8 * row + column]], "GL6 coefficient outside GF16");
        }
        // Each input-state basis vector must map to one matrix row under M^T.
        for (unsigned basis = 0; basis != 6; ++basis) {
            const auto zero = _mm_setzero_si128(), one = c2::reg(c2::Block{0x57, 0});
            auto low = _mm512_setzero_si512();
            auto high = _mm256_setzero_si256();
            c2::subfieldWideRefreshAdd<0>(matrix, basis == 0 ? one : zero, low, high);
            c2::subfieldWideRefreshAdd<1>(matrix, basis == 1 ? one : zero, low, high);
            c2::subfieldWideRefreshAdd<2>(matrix, basis == 2 ? one : zero, low, high);
            c2::subfieldWideRefreshAdd<3>(matrix, basis == 3 ? one : zero, low, high);
            c2::subfieldWideRefreshAdd<4>(matrix, basis == 4 ? one : zero, low, high);
            c2::subfieldWideRefreshAdd<5>(matrix, basis == 5 ? one : zero, low, high);
            std::array<c2::Block, 6> actual;
            _mm512_storeu_si512(actual.data(), low);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(actual.data() + 4), high);
            for (unsigned j = 0; j != 6; ++j)
                require(actual[j] == c2::Block{c2::sfMultiply(0x57, matrix.entries[8 * basis + j]), 0},
                        "packed GL6 transpose product mismatch");
        }
    }
}

static void fftMapChecks() {
    c2::SubfieldWidePlan plan(8, 13);
    for (unsigned j = 0; j != 16; ++j)
        require(c2::sf16Subspace(2, j) == (c2::sf16Monomial(4, j) ^
                    c2::sf16Multiply(7, c2::sf16Monomial(2, j)) ^ c2::sf16Multiply(6, std::uint8_t(j))),
                "GF16 S2 polynomial coefficients differ from FFT factorization");
    // Isolate all sixteen input and six state basis coordinates. The output
    // object deliberately has only 16-byte alignment, not incidental ZMM alignment.
    for (unsigned basis = 0; basis != 22; ++basis) {
        std::array<c2::Block, 16> input{}, output{};
        c2::SubfieldWideState state;
        if (basis < 16) input[basis].lo = 1;
        else {
            const auto one = c2::reg(c2::Block{1, 0});
            if (basis == 16) state.r0 = one;
            if (basis == 17) state.r1 = one;
            if (basis == 18) state.r2 = one;
            if (basis == 19) state.r3 = one;
            if (basis == 20) state.r4 = one;
            if (basis == 21) state.r5 = one;
        }
        static_assert(alignof(c2::SubfieldWideState) == 16);
        alignas(64) std::array<std::byte, sizeof(c2::SubfieldWideState) + 64> storage{};
        auto& moments = *std::construct_at(reinterpret_cast<c2::SubfieldWideState*>(storage.data() + 16));
        c2::subfieldWideInnerFft<false>(plan, input.data(), output.data(), 0, 0, state, moments);
        for (unsigned j = 0; j != 16; ++j) {
            auto expected = std::uint8_t(basis == j);
            if (basis >= 16) expected ^= c2::sf16Embed(c2::sf16Monomial(basis - 16, j));
            expected = c2::sfMultiply(expected, plan.routedDiagonal[j]);
            require(output[j] == c2::Block{expected, 0}, "wide FFT evaluation differs from monomial basis map");
        }
        const std::array<c2::Block, 6> actual{c2::block(moments.r0), c2::block(moments.r1),
            c2::block(moments.r2), c2::block(moments.r3), c2::block(moments.r4), c2::block(moments.r5)};
        for (unsigned power = 0; power != 6; ++power) {
            const auto expected = basis < 16 ? c2::sf16Embed(c2::sf16Monomial(power, basis)) : 0;
            require(actual[power] == c2::Block{static_cast<std::uint64_t>(expected), 0},
                    "wide FFT transpose differs from monomial basis map");
        }
        std::destroy_at(&moments);
    }
}

static void check() {
    matrixChecks();
    fftMapChecks();
    std::mt19937_64 rng(20261008);
    const c2::Block guard{0x32bcf571ed348a09ULL, 0x92017bed523ac816ULL};
    for (std::size_t k : {8, 16, 24, 32, 56, 64, 128, 136, 256})
      for (std::uint64_t seed : {1, 91}) {
        c2::SubfieldWidePlan plan(k, seed);
        for (std::size_t i = 0; i != plan.n; ++i) {
            require(plan.route[i] < plan.n && plan.route[i] / plan.groups == i / plan.groups,
                    "wide route left its region");
            require(plan.inverseRoute[plan.route[i]] == i, "wide route inverse mismatch");
            require(plan.routedDiagonal[i] != 0, "wide diagonal contains zero");
        }
        for (unsigned j = 0; j != 16; ++j)
            for (unsigned power = 1; power != 6; ++power)
                require(plan.coefficients[j / 4][power - 1].bytes[2 * (j % 4)] ==
                        0x0101010101010101ULL * c2::sf16Embed(c2::sf16Monomial(power, j)),
                        "wide map differs from monomial specification");
        for (unsigned bytes : {1, 4, 8, 16}) {
            std::vector<std::uint8_t> message(k * bytes), inputBytes(plan.n * bytes);
            for (auto& value : message) value = std::uint8_t(rng());
            for (auto& value : inputBytes) value = std::uint8_t(rng());
            const auto encoded = c2::subfieldWideForwardOracle(plan, message, bytes);
            const auto expected = c2::subfieldWideTransposeOracle(plan, inputBytes, bytes);
            require(dot(encoded, inputBytes) == dot(message, expected), "wide scalar adjoint mismatch");
            const auto scalar = std::uint8_t(rng());
            auto scaledMessage = message, scaledEncoded = encoded;
            for (auto& value : scaledMessage) value = c2::sfMultiply(value, scalar);
            for (auto& value : scaledEncoded) value = c2::sfMultiply(value, scalar);
            require(c2::subfieldWideForwardOracle(plan, scaledMessage, bytes) == scaledEncoded,
                    "wide GF256 scalar extension mismatch");
            if (bytes != 16) continue;
            std::vector<c2::Block> input(plan.n), output(plan.k + 2), planes(plan.n + 2);
            std::memcpy(input.data(), inputBytes.data(), inputBytes.size());
            for (auto kernel : {c2::Subfield16Kernel::Scatter, c2::Subfield16Kernel::RegionNonTemporal,
                                c2::Subfield16Kernel::RegionAllNonTemporal})
              for (auto map : {c2::SubfieldWideMapKernel::Direct, c2::SubfieldWideMapKernel::AdditiveFft}) {
                const auto temporarySize = c2::subfieldWideTemporaryRecords(plan, kernel);
                std::vector<c2::Block> temporary(temporarySize + 2);
                for (unsigned offset : {0, 1}) {
                    std::fill(output.begin(), output.end(), guard);
                    std::fill(planes.begin(), planes.end(), guard);
                    std::fill(temporary.begin(), temporary.end(), guard);
                    c2::subfieldWideTransposeFast128(plan, input.data(), output.data() + offset,
                        planes.data() + offset, temporary.data() + offset, kernel, map);
                    require(std::memcmp(output.data() + offset, expected.data(), expected.size()) == 0,
                            "wide SIMD transpose differs from scalar oracle");
                    if (offset) require(output[0] == guard && planes[0] == guard && temporary[0] == guard,
                                        "wide prefix guard changed");
                    require(output[plan.k + offset] == guard && planes[plan.n + offset] == guard &&
                            temporary[temporarySize + offset] == guard, "wide suffix guard changed");
                    std::vector<std::uint8_t> fastBytes(expected.size());
                    std::memcpy(fastBytes.data(), output.data() + offset, expected.size());
                    require(dot(encoded, inputBytes) == dot(message, fastBytes), "wide SIMD adjoint mismatch");
                }
                auto alias = input;
                c2::subfieldWideTransposeFast128(plan, alias.data(), alias.data(), planes.data(), temporary.data(), kernel, map);
                require(std::memcmp(alias.data(), expected.data(), expected.size()) == 0, "wide in-place output mismatch");
                require(std::equal(alias.begin() + k, alias.end(), input.begin() + k), "wide in-place suffix changed");
                c2::BenchPayload alignedInput(plan.n, false), alignedPlanes(plan.n, false);
                c2::BenchPayload alignedTemporary(std::max<std::size_t>(1, temporarySize), false);
                std::memcpy(alignedInput.p, input.data(), plan.n * sizeof(c2::Block));
                c2::subfieldWideTransposeFast128(plan, alignedInput.p, alignedInput.p,
                    alignedPlanes.p, alignedTemporary.p, kernel, map);
                require(std::memcmp(alignedInput.p, expected.data(), expected.size()) == 0, "wide aligned SIMD mismatch");
                require(std::equal(alignedInput.p + k, alignedInput.p + plan.n, input.data() + k),
                        "wide aligned suffix changed");
            }
        }
      }
    for (std::size_t invalid : {0, 1, 7, 9}) {
        bool rejected = false;
        try { c2::SubfieldWidePlan bad(invalid, 1); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid wide natural size accepted");
    }
    std::cout << "GF16 s6: GL6 determinant/sampler and packed update, monomial A/C, scalar and SIMD adjoints, "
                 "direct/FFT basis maps, three schedules, natural sizes, guards and aliases pass\n";
}

static double micros(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::micro>(end - begin).count();
}

static void benchmark(std::size_t k, std::uint64_t seed, unsigned calls, bool huge,
                      bool phases, const std::string& kernelName, const std::string& mapName) {
    const auto kernel = parseKernel(kernelName);
    const auto map = parseMap(mapName);
    const auto setupStart = Clock::now();
    c2::SubfieldWidePlan plan(k, seed);
    const auto setupMs = micros(setupStart, Clock::now()) / 1000;
    c2::BenchPayload input(plan.n, huge), planes(plan.n, huge);
    const auto temporarySize = c2::subfieldWideTemporaryRecords(plan, kernel);
    std::unique_ptr<c2::BenchPayload> temporary;
    if (temporarySize) temporary = std::make_unique<c2::BenchPayload>(temporarySize, huge);
    auto temporaryPointer = temporary ? temporary->p : nullptr;
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
    auto encode = [&] { c2::subfieldWideTransposeFast128(plan, input.p, input.p, planes.p, temporaryPointer, kernel, map); };
    for (unsigned i = 0; i != 5; ++i) encode();
    std::vector<double> total(calls), inner(calls), outer(calls);
    for (unsigned sample = 0; sample != calls; ++sample) {
        const auto begin = Clock::now();
        if (phases) {
            c2::subfieldWideInnerDispatch(plan, input.p, planes.p, temporaryPointer, kernel, map);
            const auto split = Clock::now();
            c2::subfieldWideOuterDispatch(plan, planes.p, input.p, kernel);
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
    std::cout << "candidate,K,t,s,outer,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,mode,kernel,map\n"
              << "gf16-gl6-pre-d," << k << ",16,6,16," << seed << ',' << calls << ',' << std::fixed
              << std::setprecision(4) << total[calls / 2] << ',' << total[calls / 10] << ','
              << total[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes() << ','
              << (plan.n + temporarySize) * sizeof(c2::Block) << ',' << std::hex << checksum << std::dec << ','
              << (huge ? "huge-preferred" : "normal") << ',' << (phases ? "phases" : "full") << ','
              << kernelName << ',' << mapName << '\n';
    if (phases) {
        std::sort(inner.begin(), inner.end());
        std::sort(outer.begin(), outer.end());
        std::cout << "inner_route_us,outer_us,total_us\n" << inner[calls / 2] << ','
                  << outer[calls / 2] << ',' << total[calls / 2] << '\n';
    }
}

int main(int argc, char** argv) {
    try {
        if (!c2::gfni512Available()) return 77;
        if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
        if (argc < 2 || argc > 8)
            throw std::invalid_argument("usage: wide K [seed=113] [calls=51] [huge|normal] [full|phases] [scatter|region-nt|region-nt-all] [direct|fft]");
        const auto k = std::stoull(argv[1]);
        const auto seed = argc > 2 ? std::stoull(argv[2]) : 113;
        const auto calls = argc > 3 ? std::stoul(argv[3]) : 51;
        const std::string pages = argc > 4 ? argv[4] : "huge";
        const std::string mode = argc > 5 ? argv[5] : "full";
        const std::string kernel = argc > 6 ? argv[6] : "scatter";
        const std::string map = argc > 7 ? argv[7] : "direct";
        if (calls < 1 || calls > 1000000 || (pages != "huge" && pages != "normal") ||
            (mode != "full" && mode != "phases")) throw std::invalid_argument("wide benchmark options");
        benchmark(k, seed, static_cast<unsigned>(calls), pages == "huge", mode == "phases", kernel, map);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
