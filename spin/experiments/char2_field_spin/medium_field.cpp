#include "MediumField32.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace mf = spin::research::char2::medium32;
using Clock = std::chrono::steady_clock;

static bool available() noexcept {
#if defined(_MSC_VER)
    int cpu[4]{};
    __cpuid(cpu, 0);
    if (cpu[0] < 7) return false;
    __cpuidex(cpu, 1, 0);
    if (!(cpu[2] & (1 << 27)) || (_xgetbv(0) & 0xe6) != 0xe6) return false;
    __cpuidex(cpu, 7, 0);
    return (cpu[1] & (1 << 16)) && (cpu[1] & (1 << 30)) && (cpu[1] & (1U << 31)) &&
           (cpu[2] & (1 << 8)) && (cpu[2] & (1 << 10));
#else
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") &&
           __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("gfni") &&
           __builtin_cpu_supports("vpclmulqdq");
#endif
}

static void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<mf::Kernel K>
static std::uint32_t oracle(std::uint32_t x, std::uint32_t c32, std::uint8_t c8) {
    if constexpr (K == mf::Kernel::Poly32) return mf::multiplyOracle(x, c32);
    else if constexpr (K == mf::Kernel::Tower32) return mf::towerMultiplyOracle(x, c32);
    else if constexpr (K == mf::Kernel::Poly32Small) return mf::multiplyOracle(x, c8);
    else return std::uint32_t(mf::byteOracle(std::uint8_t(x), c8)) |
        (std::uint32_t(mf::byteOracle(std::uint8_t(x >> 8), c8)) << 8) |
        (std::uint32_t(mf::byteOracle(std::uint8_t(x >> 16), c8)) << 16) |
        (std::uint32_t(mf::byteOracle(std::uint8_t(x >> 24), c8)) << 24);
}

static void arithmeticChecks() {
    require(mf::frobenius(32) == 2 &&
        mf::polynomialGcd(mf::frobenius(16) ^ 2U, mf::polynomial) == 1,
        "polynomial GF32 irreducibility criterion failed");
    require(mf::byteTrace(mf::delta) == 1, "tower quadratic trace criterion failed");
    for (unsigned a = 0; a != 256; ++a)
        for (unsigned b = 0; b != 256; ++b)
            require(mf::byteMultiply(std::uint8_t(a), std::uint8_t(b)) ==
                    mf::byteOracle(std::uint8_t(a), std::uint8_t(b)), "AES byte product mismatch");

    alignas(64) std::array<std::uint32_t, 16> a{}, b{}, p{}, t{};
    for (unsigned left = 0; left != 32; ++left)
        for (unsigned firstRight = 0; firstRight != 32; firstRight += 16) {
            for (unsigned lane = 0; lane != 16; ++lane) {
                a[lane] = 1U << left;
                b[lane] = 1U << (firstRight + lane);
            }
            _mm512_store_si512(p.data(), mf::multiply512(_mm512_load_si512(a.data()), _mm512_load_si512(b.data())));
            _mm512_store_si512(t.data(), mf::towerMultiply512(_mm512_load_si512(a.data()), _mm512_load_si512(b.data())));
            for (unsigned lane = 0; lane != 16; ++lane) {
                require(p[lane] == mf::multiplyOracle(a[lane], b[lane]), "polynomial SIMD basis mismatch");
                require(p[lane] == mf::multiplyScalar(a[lane], b[lane]), "polynomial scalar basis mismatch");
                require(t[lane] == mf::towerMultiplyOracle(a[lane], b[lane]), "tower SIMD basis mismatch");
                require(t[lane] == mf::towerMultiplyScalar(a[lane], b[lane]), "tower scalar basis mismatch");
            }
        }
    std::mt19937_64 rng(20261005);
    for (unsigned sample = 0; sample != 1024; ++sample) {
        for (unsigned lane = 0; lane != 16; ++lane) { a[lane] = std::uint32_t(rng()); b[lane] = std::uint32_t(rng()); }
        _mm512_store_si512(p.data(), mf::multiply512(_mm512_load_si512(a.data()), _mm512_load_si512(b.data())));
        _mm512_store_si512(t.data(), mf::towerMultiply512(_mm512_load_si512(a.data()), _mm512_load_si512(b.data())));
        for (unsigned lane = 0; lane != 16; ++lane) {
            require(p[lane] == mf::multiplyOracle(a[lane], b[lane]), "polynomial SIMD random mismatch");
            require(p[lane] == mf::multiplyScalar(a[lane], b[lane]), "polynomial scalar random mismatch");
            require(t[lane] == mf::towerMultiplyOracle(a[lane], b[lane]), "tower SIMD random mismatch");
            require(t[lane] == mf::towerMultiplyScalar(a[lane], b[lane]), "tower scalar random mismatch");
        }
    }
    for (unsigned coefficient = 0; coefficient != 256; ++coefficient) {
        alignas(16) std::array<std::uint8_t, 4> c{};
        c.fill(std::uint8_t(coefficient));
        for (auto& lane : a) lane = std::uint32_t(rng());
        _mm512_store_si512(p.data(), mf::smallMultiply512(_mm512_load_si512(a.data()), c.data()));
        const auto embedded = _mm512_set1_epi32(coefficient);
        _mm512_store_si512(t.data(), mf::towerMultiply512(_mm512_load_si512(a.data()), embedded));
        for (unsigned lane = 0; lane != 16; ++lane) {
            require(p[lane] == mf::multiplyOracle(a[lane], coefficient), "polynomial small coefficient mismatch");
            require(t[lane] == oracle<mf::Kernel::Gf8>(a[lane], 0, std::uint8_t(coefficient)),
                    "tower embedded GF256 multiplication is not bytewise GFNI");
        }
    }
    // A high tower component must mix basis components: this catches the
    // accidental replacement of GF32 multiplication by four byte multiplies.
    require(mf::towerMultiplyScalar(0x100, 0x100) == (0x100U ^ mf::delta), "u^2 relation failed");
    require(mf::towerMultiplyScalar(0x10000, 0x10000) == (0x10000U ^ (unsigned(mf::delta) << 8)),
            "v^2 relation failed");
}

template<mf::Kernel K>
static void diagonalChecks() {
    std::mt19937_64 rng(518731);
    const mf::Record128 guard{{0x736baa91U, 0xf152793dU, 0xa0407b1eU, 0x103abc47U}};
    for (std::size_t count : {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 63, 64, 65, 257}) {
        std::vector<mf::Record128> input(count), expected(count), output(count + 2, guard);
        std::vector<std::uint32_t> c32(count);
        std::vector<std::uint8_t> c8(count);
        for (std::size_t i = 0; i != count; ++i) {
            c32[i] = std::uint32_t(rng()); c8[i] = std::uint8_t(rng());
            for (unsigned j = 0; j != 4; ++j) {
                input[i].lane[j] = std::uint32_t(rng());
                expected[i].lane[j] = oracle<K>(input[i].lane[j], c32[i], c8[i]);
            }
        }
        for (unsigned offset = 0; offset != 2; ++offset) {
            std::fill(output.begin(), output.end(), guard);
            mf::diagonal<K>(input.data(), c32.data(), c8.data(), output.data() + offset, count);
            require(std::equal(expected.begin(), expected.end(), output.begin() + offset), "array kernel/oracle mismatch");
            require(output[count + offset] == guard, "array suffix guard changed");
            require(!offset || output[0] == guard, "array prefix guard changed");
        }
        auto alias = input;
        mf::diagonal<K>(alias.data(), c32.data(), c8.data(), alias.data(), count);
        require(alias == expected, "array in-place alias mismatch");
    }
}

static void check() {
    arithmeticChecks();
    diagonalChecks<mf::Kernel::Gf8>();
    diagonalChecks<mf::Kernel::Poly32>();
    diagonalChecks<mf::Kernel::Tower32>();
    diagonalChecks<mf::Kernel::Poly32Small>();
    std::cout << "GF32 polynomial=0x10000008d: Rabin irreducibility pass; tower delta=0x20: trace=1; "
                 "AES products, all 1024 binary basis products in both representations, independent quartic oracle, "
                 "random 16-lane SIMD products, every byte coefficient, embedded GF256 action, tails, guards, aliases pass\n";
}

template<mf::Kernel K>
static constexpr const char* kernelName() {
    if constexpr (K == mf::Kernel::Gf8) return "gf8";
    else if constexpr (K == mf::Kernel::Poly32) return "poly32";
    else if constexpr (K == mf::Kernel::Tower32) return "tower32";
    else return "poly32-small";
}

template<mf::Kernel K>
static void benchmark(std::size_t count, unsigned calls, unsigned passes, std::uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::vector<mf::Record128> input(count), output(count);
    std::vector<std::uint32_t> c32(count);
    std::vector<std::uint8_t> c8(count);
    for (std::size_t i = 0; i != count; ++i) {
        for (auto& lane : input[i].lane) lane = std::uint32_t(rng());
        // Rejection gives uniform nonzero coefficients under an ideal uniform
        // word source. This deterministic benchmark PRNG is only a fixture.
        do { c32[i] = std::uint32_t(rng()); } while (!c32[i]);
        do { c8[i] = std::uint8_t(rng()); } while (!c8[i]);
    }
    const auto run = [&] { mf::diagonal<K>(input.data(), c32.data(), c8.data(), output.data(), count); };
    for (unsigned i = 0; i != 3; ++i) run();
    std::vector<double> timings(calls);
    for (auto& timing : timings) {
        const auto begin = Clock::now();
        for (unsigned pass = 0; pass != passes; ++pass) run();
        timing = std::chrono::duration<double, std::micro>(Clock::now() - begin).count() / passes;
    }
    std::sort(timings.begin(), timings.end());
    std::uint64_t checksum = 0;
    for (const auto& record : output)
        for (const auto lane : record.lane) checksum = (checksum ^ lane) * 0x100000001b3ULL;
    constexpr bool wide = K == mf::Kernel::Poly32 || K == mf::Kernel::Tower32;
    std::cout << kernelName<K>() << ',' << count << ",128," << (K == mf::Kernel::Gf8 ? 16 : 4) << ','
              << (wide ? 32 : 8) << ',' << calls << ',' << passes << ',' << seed << ',' << std::fixed
              << std::setprecision(6) << timings[calls / 2] << ',' << timings[calls / 10] << ','
              << timings[9 * calls / 10] << ',' << 1000 * timings[calls / 2] / count << ','
              << (wide ? 4 : 1) * count << ',' << 32 * count << ',' << std::hex << checksum << std::dec << '\n';
}

int main(int argc, char** argv) {
    try {
        if (!available()) {
            std::cerr << "medium_field requires OS-enabled AVX512F/BW/VL, GFNI, and VPCLMULQDQ\n";
            return 77;
        }
        if (argc == 2 && std::string(argv[1]) == "--check") { check(); return 0; }
        if (argc < 2 || argc > 6)
            throw std::invalid_argument("usage: medium_field --check | records [calls=51] [all|gf8|poly32|tower32|poly32-small] [passes=1] [seed=113]");
        const auto count = std::stoull(argv[1]);
        const auto calls = argc > 2 ? std::stoul(argv[2]) : 51;
        const std::string kernel = argc > 3 ? argv[3] : "all";
        const auto passes = argc > 4 ? std::stoul(argv[4]) : 1;
        const auto seed = argc > 5 ? std::stoull(argv[5]) : 113;
        if (!count || count > (1ULL << 28) || !calls || calls > 1000000 || !passes || passes > 1000000 ||
            (kernel != "all" && kernel != "gf8" && kernel != "poly32" && kernel != "tower32" && kernel != "poly32-small"))
            throw std::invalid_argument("invalid medium-field benchmark dimensions or kernel");
        std::cout << "kernel,records,payload_bits,field_lanes,coefficient_bits,calls,passes,seed,median_us,p10_us,p90_us,ns_per_record,coefficient_bytes,payload_io_bytes,checksum\n";
        if (kernel == "all" || kernel == "gf8") benchmark<mf::Kernel::Gf8>(count, unsigned(calls), unsigned(passes), seed);
        if (kernel == "all" || kernel == "poly32") benchmark<mf::Kernel::Poly32>(count, unsigned(calls), unsigned(passes), seed);
        if (kernel == "all" || kernel == "tower32") benchmark<mf::Kernel::Tower32>(count, unsigned(calls), unsigned(passes), seed);
        if (kernel == "all" || kernel == "poly32-small") benchmark<mf::Kernel::Poly32Small>(count, unsigned(calls), unsigned(passes), seed);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
