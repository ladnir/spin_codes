#include "Subfield16.h"
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
    throw std::invalid_argument("GF16 kernel must be scatter, region-nt or region-nt-all");
}

static c2::Subfield16MapKernel parseMap(const std::string& name) {
    if (name == "direct") return c2::Subfield16MapKernel::Direct;
    if (name == "fft") return c2::Subfield16MapKernel::AdditiveFft;
    if (name == "packed") return c2::Subfield16MapKernel::PackedDirect;
    throw std::invalid_argument("GF16 map must be direct, fft or packed");
}

static std::uint8_t dot(const std::vector<std::uint8_t>& left,
                         const std::vector<std::uint8_t>& right) {
    require(left.size() == right.size(), "GF256 component dot dimensions");
    std::uint8_t result = 0;
    for (std::size_t i = 0; i != left.size(); ++i) result ^= c2::sfMultiply(left[i], right[i]);
    return result;
}

// Determinant via all 24 permutations, independent of the setup elimination.
// Characteristic two removes the permutation signs.
static std::uint8_t determinant(const std::array<std::uint8_t, 16>& matrix) {
    std::array<unsigned, 4> permutation{0, 1, 2, 3};
    std::uint8_t result = 0;
    do {
        std::uint8_t term = 1;
        for (unsigned row = 0; row != 4; ++row)
            term = c2::sfMultiply(term, matrix[4 * row + permutation[row]]);
        result ^= term;
    } while (std::next_permutation(permutation.begin(), permutation.end()));
    return result;
}

static void arithmeticChecks() {
    std::array<bool, 256> seen{};
    for (unsigned a = 0; a != 16; ++a) {
        const auto image = c2::sf16Embed(std::uint8_t(a));
        require(!seen[image], "GF16 embedding is not injective");
        seen[image] = true;
        auto frobenius = image;
        for (unsigned power = 0; power != 4; ++power) frobenius = c2::sfMultiply(frobenius, frobenius);
        require(frobenius == image, "embedded scalar not in GF16");
        if (a) require(c2::sf16Multiply(std::uint8_t(a), c2::sf16Inverse(std::uint8_t(a))) == 1,
                       "GF16 inverse failed");
        for (unsigned b = 0; b != 16; ++b) {
            require((image ^ c2::sf16Embed(std::uint8_t(b))) == c2::sf16Embed(std::uint8_t(a ^ b)),
                    "GF16 embedding does not preserve addition");
            require(c2::sfMultiply(image, c2::sf16Embed(std::uint8_t(b))) ==
                    c2::sf16Embed(c2::sf16Multiply(std::uint8_t(a), std::uint8_t(b))),
                    "GF16 embedding does not preserve multiplication");
        }
    }
    require(c2::sf16Basis(2, 2) != c2::sf16Monomial(2, 2), "test must distinguish additive and monomial bases");
    std::mt19937_64 rng(0x471234678abULL);
    for (unsigned trial = 0; trial != 2048; ++trial) {
        std::array<std::uint8_t, 16> abstract, embedded;
        for (unsigned j = 0; j != 16; ++j) {
            abstract[j] = std::uint8_t(rng() & 15);
            embedded[j] = c2::sf16Embed(abstract[j]);
        }
        require(c2::sf16Invertible(abstract) == (determinant(embedded) != 0),
                "GL4 elimination disagrees with independent determinant");
        const auto sampled = c2::sf16SampleMatrix(rng);
        require(determinant(sampled.entries) != 0, "sampled GL4 matrix is singular");
        for (auto entry : sampled.entries) require(seen[entry], "GL4 entry outside embedded GF16");
    }
}

static void fftMapChecks() {
    c2::Subfield16Plan plan(8, 13);
    // Basis vectors isolate all sixteen input coordinates and four state
    // coordinates; these are tests of the monomial maps, before refresh/outer.
    for (unsigned basis = 0; basis != 20; ++basis) {
        std::array<c2::Block, 16> input{}, output{};
        c2::Subfield16State state;
        if (basis < 16) input[basis].lo = 1;
        else {
            const auto one = c2::reg(c2::Block{1, 0});
            if (basis == 16) state.r0 = one;
            if (basis == 17) state.r1 = one;
            if (basis == 18) state.r2 = one;
            if (basis == 19) state.r3 = one;
        }
        // Deliberately provide only the state's required 16-byte alignment,
        // not incidental 32/64-byte stack alignment, for the FFT output.
        static_assert(alignof(c2::Subfield16State) == 16);
        alignas(64) std::array<std::byte, sizeof(c2::Subfield16State) + 64> storage{};
        auto& moments = *std::construct_at(reinterpret_cast<c2::Subfield16State*>(storage.data() + 16));
        c2::sf16InnerFft<false>(plan, input.data(), output.data(), 0, 0, state, moments);
        for (unsigned j = 0; j != 16; ++j) {
            auto expected = std::uint8_t(basis == j);
            if (basis >= 16) expected ^= c2::sf16Embed(c2::sf16Monomial(basis - 16, j));
            expected = c2::sfMultiply(expected, plan.routedDiagonal[j]);
            require(output[j] == c2::Block{expected, 0}, "FFT evaluation differs from monomial basis map");
        }
        const std::array<c2::Block, 4> actual{c2::block(moments.r0), c2::block(moments.r1),
                                             c2::block(moments.r2), c2::block(moments.r3)};
        for (unsigned power = 0; power != 4; ++power) {
            const auto expected = basis < 16 ? c2::sf16Embed(c2::sf16Monomial(power, basis)) : 0;
            require(actual[power] == c2::Block{static_cast<std::uint64_t>(expected), 0},
                    "FFT transpose differs from monomial basis map");
        }
        std::destroy_at(&moments);
    }
}

static void check() {
    arithmeticChecks();
    fftMapChecks();
    std::mt19937_64 rng(20261007);
    const c2::Block guard{0xc51af32071de5694ULL, 0x59fe7823ac16094bULL};
    for (std::size_t k : {8, 16, 24, 32, 56, 64, 128, 136, 256, 1024})
      for (std::uint64_t seed : {1, 91}) {
        c2::Subfield16Plan plan(k, seed);
        for (std::size_t i = 0; i != plan.n; ++i) {
            require(plan.route[i] < plan.n && plan.route[i] / plan.groups == i / plan.groups,
                    "GF16 route left its region");
            require(plan.inverseRoute[plan.route[i]] == i, "GF16 route is not a bijection");
            require(plan.routedDiagonal[i] != 0, "GF16 diagonal contains zero");
        }
        for (unsigned j = 0; j != 16; ++j)
            for (unsigned power = 1; power != 4; ++power) {
                const auto coefficient = c2::sf16Embed(c2::sf16Monomial(power, j));
                require(plan.coefficients[j / 4][power - 1].bytes[2 * (j % 4)] ==
                        coefficient * 0x0101010101010101ULL, "inner map is not the specified monomial map");
            }
        for (unsigned bytes : {1, 4, 8, 16}) {
            std::vector<std::uint8_t> message(plan.k * bytes), inputBytes(plan.n * bytes);
            for (auto& value : message) value = std::uint8_t(rng());
            for (auto& value : inputBytes) value = std::uint8_t(rng());
            const auto encoded = c2::subfield16ForwardOracle(plan, message, bytes);
            const auto expected = c2::subfield16TransposeOracle(plan, inputBytes, bytes);
            require(dot(encoded, inputBytes) == dot(message, expected), "GF16 scalar adjoint mismatch");
            // The component matrix is over the embedded GF16 subfield. Check
            // linearity under arbitrary GF256 scalars, not only nibble scalars.
            const auto scalar = std::uint8_t(rng());
            auto scaledMessage = message, scaledExpected = encoded;
            for (auto& value : scaledMessage) value = c2::sfMultiply(value, scalar);
            for (auto& value : scaledExpected) value = c2::sfMultiply(value, scalar);
            require(c2::subfield16ForwardOracle(plan, scaledMessage, bytes) == scaledExpected,
                    "GF256 scalar extension linearity mismatch");
            if (bytes != 16) continue;
            std::vector<c2::Block> input(plan.n), output(plan.k + 2), planes(plan.n + 2);
            std::memcpy(input.data(), inputBytes.data(), inputBytes.size());
            for (auto kernel : {c2::Subfield16Kernel::Scatter, c2::Subfield16Kernel::RegionNonTemporal,
                                c2::Subfield16Kernel::RegionAllNonTemporal})
              for (auto map : {c2::Subfield16MapKernel::Direct, c2::Subfield16MapKernel::AdditiveFft,
                               c2::Subfield16MapKernel::PackedDirect}) {
                const auto temporarySize = c2::subfield16TemporaryRecords(plan, kernel);
                std::vector<c2::Block> temporary(temporarySize + 2);
                for (unsigned offset : {0, 1}) {
                    std::fill(output.begin(), output.end(), guard);
                    std::fill(planes.begin(), planes.end(), guard);
                    std::fill(temporary.begin(), temporary.end(), guard);
                    c2::subfield16TransposeFast128(plan, input.data(), output.data() + offset,
                        planes.data() + offset, temporary.data() + offset, kernel, map);
                    require(std::memcmp(output.data() + offset, expected.data(), expected.size()) == 0,
                            "GF16 SIMD differs from scalar transpose");
                    if (offset) require(output[0] == guard && planes[0] == guard && temporary[0] == guard,
                                        "GF16 prefix guard changed");
                    require(output[plan.k + offset] == guard && planes[plan.n + offset] == guard &&
                            temporary[temporarySize + offset] == guard, "GF16 suffix guard changed");
                    std::vector<std::uint8_t> fastBytes(expected.size());
                    std::memcpy(fastBytes.data(), output.data() + offset, fastBytes.size());
                    require(dot(encoded, inputBytes) == dot(message, fastBytes), "GF16 SIMD adjoint mismatch");
                }
                auto alias = input;
                c2::subfield16TransposeFast128(plan, alias.data(), alias.data(), planes.data(), temporary.data(), kernel, map);
                require(std::memcmp(alias.data(), expected.data(), expected.size()) == 0, "GF16 alias output mismatch");
                require(std::equal(alias.begin() + plan.k, alias.end(), input.begin() + plan.k), "GF16 alias suffix changed");
                // Exercise streaming branches even if vector allocation was
                // not 64-byte-aligned; also verify exact in-place alias there.
                c2::BenchPayload alignedInput(plan.n, false), alignedPlanes(plan.n, false);
                c2::BenchPayload alignedTemporary(std::max<std::size_t>(1, temporarySize), false);
                std::memcpy(alignedInput.p, input.data(), plan.n * sizeof(c2::Block));
                c2::subfield16TransposeFast128(plan, alignedInput.p, alignedInput.p, alignedPlanes.p,
                    alignedTemporary.p, kernel, map);
                require(std::memcmp(alignedInput.p, expected.data(), expected.size()) == 0, "GF16 aligned output mismatch");
                require(std::equal(alignedInput.p + plan.k, alignedInput.p + plan.n, input.data() + plan.k),
                        "GF16 aligned alias suffix changed");
            }
        }
      }
    for (std::size_t invalid : {0, 1, 7, 9, 15}) {
        bool rejected = false;
        try { c2::Subfield16Plan bad(invalid, 1); }
        catch (const std::invalid_argument&) { rejected = true; }
        require(rejected, "invalid GF16 message length accepted");
    }
    std::cout << "GF16: beta=" << unsigned(c2::sf16Beta)
              << "; embedding, GL4 sampler/determinants, monomial A/C, scalar and SIMD adjoints, "
                 "GRS16 outer, natural sizes, three exact schedules, direct/FFT maps, guards and aliases pass\n";
}

static double micros(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::micro>(end - begin).count();
}

static void benchmark(std::size_t k, std::uint64_t seed, unsigned calls, bool huge,
                      bool phases, const std::string& kernelName, const std::string& mapName) {
    const auto kernel = parseKernel(kernelName);
    const auto map = parseMap(mapName);
    const auto setupStart = Clock::now();
    c2::Subfield16Plan plan(k, seed);
    const auto setupMs = micros(setupStart, Clock::now()) / 1000;
    c2::BenchPayload input(plan.n, huge), planes(plan.n, huge);
    const auto temporarySize = c2::subfield16TemporaryRecords(plan, kernel);
    std::unique_ptr<c2::BenchPayload> temporary;
    if (temporarySize) temporary = std::make_unique<c2::BenchPayload>(temporarySize, huge);
    auto temporaryPointer = temporary ? temporary->p : nullptr;
    std::mt19937_64 rng(913);
    for (std::size_t i = 0; i != plan.n; ++i) input.p[i] = {rng(), rng()};
    auto encode = [&] {
        c2::subfield16TransposeFast128(plan, input.p, input.p, planes.p, temporaryPointer, kernel, map);
    };
    for (unsigned i = 0; i != 5; ++i) encode();
    std::vector<double> total(calls), inner(calls), outer(calls);
    for (unsigned sample = 0; sample != calls; ++sample) {
        const auto begin = Clock::now();
        if (phases) {
            c2::subfield16InnerDispatch(plan, input.p, planes.p, temporaryPointer, kernel, map);
            const auto split = Clock::now();
            c2::subfield16OuterDispatch(plan, planes.p, input.p, kernel);
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
    std::cout << "candidate,K,t,s,outer,seed,calls,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum,page_policy,mode,kernel,map,beta\n"
              << "gf16-gl4-pre-d," << k << ",16,4,16," << seed << ',' << calls << ',' << std::fixed
              << std::setprecision(4) << total[calls / 2] << ',' << total[calls / 10] << ','
              << total[9 * calls / 10] << ',' << setupMs << ',' << plan.setupBytes() << ','
              << (plan.n + temporarySize) * sizeof(c2::Block) << ',' << std::hex << checksum << std::dec << ','
              << (huge ? "huge-preferred" : "normal") << ',' << (phases ? "phases" : "full") << ','
              << kernelName << ',' << mapName << ',' << unsigned(c2::sf16Beta) << '\n';
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
            throw std::invalid_argument("usage: subfield16 K [seed=113] [calls=51] [huge|normal] [full|phases] [scatter|region-nt|region-nt-all] [direct|fft|packed]");
        const auto k = std::stoull(argv[1]);
        const auto seed = argc > 2 ? std::stoull(argv[2]) : 113;
        const auto calls = argc > 3 ? std::stoul(argv[3]) : 51;
        const std::string pages = argc > 4 ? argv[4] : "huge";
        const std::string mode = argc > 5 ? argv[5] : "full";
        const std::string kernel = argc > 6 ? argv[6] : "scatter";
        const std::string map = argc > 7 ? argv[7] : "direct";
        if (calls < 1 || calls > 1000000 || (pages != "huge" && pages != "normal") ||
            (mode != "full" && mode != "phases")) throw std::invalid_argument("GF16 benchmark options");
        benchmark(k, seed, static_cast<unsigned>(calls), pages == "huge", mode == "phases", kernel, map);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
