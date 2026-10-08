#include "SmallField.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace c2 = spin::research::char2;

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

static std::vector<c2::Coeff4> prepare(const std::vector<std::uint8_t>& coefficients) {
    std::vector<c2::Coeff4> result((coefficients.size() + 3) / 4);
    for (std::size_t i = 0; i != result.size(); ++i) {
        std::array<std::uint8_t, 4> c{1,1,1,1};
        for (unsigned j = 0; j != 4 && 4 * i + j < coefficients.size(); ++j)
            c[j] = coefficients[4 * i + j];
        result[i] = c2::prepare4(c.data());
    }
    return result;
}

static void checkCase(unsigned bits, const std::vector<std::uint8_t>& input,
                      const std::vector<std::uint8_t>& coefficients, unsigned offset) {
    const auto n = coefficients.size();
    const auto bytes = bits / 8;
    require(input.size() == n * bytes, "test input shape");
    std::vector<std::uint8_t> reference(input.size());
    for (std::size_t i = 0; i != n; ++i)
        c2::multiplyOracle(input.data() + bytes * i, coefficients[i],
                           reference.data() + bytes * i, bits);
    const auto start = std::size_t(17 + offset);
    std::vector<std::uint8_t> storage(start + input.size() + 33, 0xa7);
    std::vector<std::uint8_t> output(storage.size(), 0xa7);
    std::copy(input.begin(), input.end(), storage.begin() + start);
    const auto original = storage;
    auto verify = [&] {
        require(std::equal(reference.begin(), reference.end(), output.begin() + start),
                "coefficient kernel differs from polynomial oracle");
        require(std::all_of(output.begin(), output.begin() + start,
                            [](auto x) { return x == 0xa7; }), "prefix guard overwritten");
        require(std::all_of(output.begin() + start + input.size(), output.end(),
                            [](auto x) { return x == 0xa7; }), "suffix guard overwritten");
        require(storage == original, "out-of-place kernel changed input");
    };
    c2::multiplyPortable(storage.data() + start, coefficients.data(),
                          output.data() + start, n, bits);
    verify();
    if (bits == 128 && c2::gfni128Available()) {
        std::fill(output.begin(), output.end(), 0xa7);
        c2::multiplyGfni128(storage.data() + start, coefficients.data(), output.data() + start, n);
        verify();
    }
    if (c2::gfni512Available()) {
        std::fill(output.begin(), output.end(), 0xa7);
        c2::multiplyGfni512(storage.data() + start, coefficients.data(),
                           output.data() + start, n, bits);
        verify();
        if (bits == 128) {
            const auto prepared = prepare(coefficients);
            std::fill(output.begin(), output.end(), 0xa7);
            c2::multiplyGfniPrepared(storage.data() + start, prepared.data(), output.data() + start, n);
            verify();
        }
    }
    output = storage;
    if (c2::gfni512Available())
        c2::multiplyGfni512(output.data() + start, coefficients.data(), output.data() + start, n, bits);
    else c2::multiplyPortable(output.data() + start, coefficients.data(), output.data() + start, n, bits);
    verify();
}

static void checks() {
    std::mt19937_64 rng(20261003);
    for (const unsigned bits : {8,32,64,128}) {
        const auto bytes = bits / 8;
        // Every input basis vector and every coefficient, including zero as an
        // arithmetic edge case. The intended sampling alphabet excludes zero.
        std::vector<std::uint8_t> basis(bits * 256 * bytes, 0), coefficients(bits * 256);
        for (unsigned bit = 0; bit != bits; ++bit)
            for (unsigned coefficient = 0; coefficient != 256; ++coefficient) {
                const auto record = bit * 256 + coefficient;
                basis[record * bytes + bit / 8] = std::uint8_t(1U << (bit % 8));
                coefficients[record] = std::uint8_t(coefficient);
            }
        checkCase(bits, basis, coefficients, 0);
        for (const std::size_t n : {0,1,2,3,4,5,7,8,9,31,63,64,65,131})
            for (const unsigned offset : {0,1,7,15}) {
                std::vector<std::uint8_t> values(n * bytes), c(n);
                for (auto& x : values) x = std::uint8_t(rng());
                for (auto& x : c) x = std::uint8_t(rng());
                checkCase(bits, values, c, offset);
            }
        // The c=X operator crosses bytes; it is not bytewise GF256 scaling.
        std::array<std::uint8_t, 16> input{}, output{};
        input[0] = 0x80;
        c2::multiplyOracle(input.data(), 2, output.data(), bits);
        if (bits > 8) require(output[0] == 0 && output[1] == 1, "cross-byte shift missing");
        input.fill(0);
        input[bytes - 1] = 0x80;
        c2::multiplyOracle(input.data(), 2, output.data(), bits);
        require(output[0] == c2::reduction(bits), "native-field reduction mismatch");
        std::cout << "PASS bits=" << bits << " all basis/coefficient pairs, tails, guards, aliases\n";
    }
    if (c2::gfni128Available()) {
        std::array<std::uint8_t, 16> expected{}, actual{}, preparedActual{};
        std::vector<std::uint8_t> coefficients(257);
        std::vector<c2::Coeff128> prepared(coefficients.size());
        for (auto& x : expected) x = std::uint8_t(rng());
        actual = expected;
        preparedActual = expected;
        for (std::size_t i = 0; i != coefficients.size(); ++i) {
            coefficients[i] = std::uint8_t(1 + rng() % 255);
            prepared[i] = c2::coefficient128(coefficients[i]);
            c2::multiplyOracle(expected.data(), coefficients[i], expected.data(), 128);
        }
        c2::chainGfni128(actual.data(), coefficients.data(), coefficients.size());
        c2::chainGfni128Prepared(preparedActual.data(), prepared.data(), prepared.size());
        require(expected == actual && expected == preparedActual, "latency-chain oracle mismatch");
    }
    std::cout << "PASS gfni128=" << c2::gfni128Available()
              << " gfni512=" << c2::gfni512Available() << '\n';
}

static std::uint64_t checksum(const std::vector<std::uint8_t>& data) {
    std::uint64_t value = 0xcbf29ce484222325ULL;
    for (auto byte : data) value = (value ^ byte) * 0x100000001b3ULL;
    return value;
}

static int benchmark(unsigned bits, std::size_t count, unsigned iterations,
                      std::uint64_t seed, const std::string& kernel) {
    if (!c2::validBits(bits) || count == 0 || iterations == 0 ||
        count > std::numeric_limits<std::size_t>::max() / (bits / 8))
        throw std::invalid_argument("invalid benchmark dimensions");
    const bool chain = kernel == "chain" || kernel == "chain-prepared";
    const bool scalar = kernel == "gfni128" || chain;
    const bool wide = kernel == "gfni512" || kernel == "prepared";
    if (!scalar && !wide && kernel != "portable") throw std::invalid_argument("unknown kernel");
    if ((scalar || kernel == "prepared") && bits != 128)
        throw std::invalid_argument("this kernel requires bits=128");
    if ((scalar && !c2::gfni128Available()) || (wide && !c2::gfni512Available())) {
        std::cout << "SKIP requested ISA unavailable\n";
        return 77;
    }
    std::mt19937_64 rng(seed);
    std::vector<std::uint8_t> input((bits / 8) * count), output(input.size()), coefficients(count);
    for (auto& x : input) x = std::uint8_t(rng());
    for (auto& x : coefficients) x = std::uint8_t(1 + rng() % 255);
    const auto prepared = prepare(coefficients);
    std::vector<c2::Coeff128> scalarPrepared(count);
    for (std::size_t i = 0; i != count; ++i) scalarPrepared[i] = c2::coefficient128(coefficients[i]);
    if (chain) {
        output.resize(16);
        std::copy_n(input.begin(), 16, output.begin());
    }
    auto invoke = [&] {
        if (kernel == "portable") c2::multiplyPortable(input.data(), coefficients.data(), output.data(), count, bits);
        else if (kernel == "gfni128") c2::multiplyGfni128(input.data(), coefficients.data(), output.data(), count);
        else if (kernel == "gfni512") c2::multiplyGfni512(input.data(), coefficients.data(), output.data(), count, bits);
        else if (kernel == "prepared") c2::multiplyGfniPrepared(input.data(), prepared.data(), output.data(), count);
        else if (kernel == "chain") c2::chainGfni128(output.data(), coefficients.data(), count);
        else c2::chainGfni128Prepared(output.data(), scalarPrepared.data(), count);
    };
    for (unsigned i = 0; i != 5; ++i) invoke();
    std::vector<double> timings;
    timings.reserve(iterations);
    for (unsigned i = 0; i != iterations; ++i) {
        const auto start = std::chrono::steady_clock::now();
        invoke();
        const auto finish = std::chrono::steady_clock::now();
        timings.push_back(std::chrono::duration<double, std::nano>(finish - start).count());
    }
    std::sort(timings.begin(), timings.end());
    const auto median = timings[timings.size() / 2];
    std::cout << std::setprecision(9) << "kernel=" << kernel << " bits=" << bits
              << " count=" << count << " iterations=" << iterations << " seed=" << seed
              << " median_us=" << median / 1000 << " ns_per_record=" << median / count
              << " min_us=" << timings.front() / 1000 << " checksum=" << std::hex << checksum(output) << '\n';
    return 0;
}

int main(int argc, char** argv) {
    try {
        bool bench = false;
        unsigned bits = 128, iterations = 101;
        std::size_t count = 65536;
        std::uint64_t seed = 17;
        std::string kernel = "gfni512";
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--check") continue;
            if (arg == "--bench") { bench = true; continue; }
            if (i + 1 == argc) throw std::invalid_argument("missing option value");
            const std::string value = argv[++i];
            if (arg == "--bits") bits = std::stoul(value);
            else if (arg == "--count") count = std::stoull(value);
            else if (arg == "--iterations") iterations = std::stoul(value);
            else if (arg == "--seed") seed = std::stoull(value);
            else if (arg == "--kernel") kernel = value;
            else throw std::invalid_argument("unknown option");
        }
        if (bench) return benchmark(bits, count, iterations, seed, kernel);
        checks();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
