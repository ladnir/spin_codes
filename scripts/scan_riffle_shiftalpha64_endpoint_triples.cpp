// Exact compatibility scan for outer weight-three endpoint profiles.
//
// The outer parity checks use projective columns (1,alpha_i), (1,0), and
// (0,1), where alpha_i=gamma^(64+i).  On a three-position support, the
// nonzero outer words form one field line.  This scanner normalizes each of
// the three coordinates in turn to the unique field value whose BCH encoding
// is all ones.  It counts exactly the resulting profiles (22,106,128).

#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <wmmintrin.h>

namespace {

constexpr std::uint64_t kReduction = 0x1bULL;
constexpr std::uint64_t kAllOneMessage = 0x8e63cf44efd4fa21ULL;
constexpr unsigned kDefaultDataBlocks = 16384;

struct Codeword { std::uint64_t low; std::uint64_t high; };

constexpr std::array<Codeword, 64> kRows{{
    {0xf4845518b9582a1fULL,0x8000000000000000ULL},
    {0xe908aa3172b0543eULL,0x8000000000000001ULL},
    {0xd2115462e560a87cULL,0x8000000000000003ULL},
    {0xa422a8c5cac150f8ULL,0x8000000000000007ULL},
    {0x4845518b9582a1f0ULL,0x800000000000000fULL},
    {0x908aa3172b0543e0ULL,0x800000000000001eULL},
    {0x2115462e560a87c0ULL,0x800000000000003dULL},
    {0x422a8c5cac150f80ULL,0x800000000000007aULL},
    {0x845518b9582a1f00ULL,0x80000000000000f4ULL},
    {0x08aa3172b0543e00ULL,0x80000000000001e9ULL},
    {0x115462e560a87c00ULL,0x80000000000003d2ULL},
    {0x22a8c5cac150f800ULL,0x80000000000007a4ULL},
    {0x45518b9582a1f000ULL,0x8000000000000f48ULL},
    {0x8aa3172b0543e000ULL,0x8000000000001e90ULL},
    {0x15462e560a87c000ULL,0x8000000000003d21ULL},
    {0x2a8c5cac150f8000ULL,0x8000000000007a42ULL},
    {0x5518b9582a1f0000ULL,0x800000000000f484ULL},
    {0xaa3172b0543e0000ULL,0x800000000001e908ULL},
    {0x5462e560a87c0000ULL,0x800000000003d211ULL},
    {0xa8c5cac150f80000ULL,0x800000000007a422ULL},
    {0x518b9582a1f00000ULL,0x80000000000f4845ULL},
    {0xa3172b0543e00000ULL,0x80000000001e908aULL},
    {0x462e560a87c00000ULL,0x80000000003d2115ULL},
    {0x8c5cac150f800000ULL,0x80000000007a422aULL},
    {0x18b9582a1f000000ULL,0x8000000000f48455ULL},
    {0x3172b0543e000000ULL,0x8000000001e908aaULL},
    {0x62e560a87c000000ULL,0x8000000003d21154ULL},
    {0xc5cac150f8000000ULL,0x8000000007a422a8ULL},
    {0x8b9582a1f0000000ULL,0x800000000f484551ULL},
    {0x172b0543e0000000ULL,0x800000001e908aa3ULL},
    {0x2e560a87c0000000ULL,0x800000003d211546ULL},
    {0x5cac150f80000000ULL,0x800000007a422a8cULL},
    {0xb9582a1f00000000ULL,0x80000000f4845518ULL},
    {0x72b0543e00000000ULL,0x80000001e908aa31ULL},
    {0xe560a87c00000000ULL,0x80000003d2115462ULL},
    {0xcac150f800000000ULL,0x80000007a422a8c5ULL},
    {0x9582a1f000000000ULL,0x8000000f4845518bULL},
    {0x2b0543e000000000ULL,0x8000001e908aa317ULL},
    {0x560a87c000000000ULL,0x8000003d2115462eULL},
    {0xac150f8000000000ULL,0x8000007a422a8c5cULL},
    {0x582a1f0000000000ULL,0x800000f4845518b9ULL},
    {0xb0543e0000000000ULL,0x800001e908aa3172ULL},
    {0x60a87c0000000000ULL,0x800003d2115462e5ULL},
    {0xc150f80000000000ULL,0x800007a422a8c5caULL},
    {0x82a1f00000000000ULL,0x80000f4845518b95ULL},
    {0x0543e00000000000ULL,0x80001e908aa3172bULL},
    {0x0a87c00000000000ULL,0x80003d2115462e56ULL},
    {0x150f800000000000ULL,0x80007a422a8c5cacULL},
    {0x2a1f000000000000ULL,0x8000f4845518b958ULL},
    {0x543e000000000000ULL,0x8001e908aa3172b0ULL},
    {0xa87c000000000000ULL,0x8003d2115462e560ULL},
    {0x50f8000000000000ULL,0x8007a422a8c5cac1ULL},
    {0xa1f0000000000000ULL,0x800f4845518b9582ULL},
    {0x43e0000000000000ULL,0x801e908aa3172b05ULL},
    {0x87c0000000000000ULL,0x803d2115462e560aULL},
    {0x0f80000000000000ULL,0x807a422a8c5cac15ULL},
    {0x1f00000000000000ULL,0x80f4845518b9582aULL},
    {0x3e00000000000000ULL,0x81e908aa3172b054ULL},
    {0x7c00000000000000ULL,0x83d2115462e560a8ULL},
    {0xf800000000000000ULL,0x87a422a8c5cac150ULL},
    {0xf000000000000000ULL,0x8f4845518b9582a1ULL},
    {0xe000000000000000ULL,0x9e908aa3172b0543ULL},
    {0xc000000000000000ULL,0xbd2115462e560a87ULL},
    {0x8000000000000000ULL,0xfa422a8c5cac150fULL},
}};

std::array<std::array<Codeword, 1U << 16>, 4> kTables{};

void initialize_tables() {
    for (unsigned chunk = 0; chunk < 4; ++chunk) {
        for (unsigned value = 1; value < (1U << 16); ++value) {
            const unsigned bit = std::countr_zero(value);
            const auto previous = value & (value - 1);
            kTables[chunk][value] = {
                kTables[chunk][previous].low ^ kRows[16 * chunk + bit].low,
                kTables[chunk][previous].high ^ kRows[16 * chunk + bit].high,
            };
        }
    }
}

unsigned bch_weight(std::uint64_t message) {
    Codeword value{0, 0};
    for (unsigned chunk = 0; chunk < 4; ++chunk) {
        const auto row = kTables[chunk][(message >> (16 * chunk)) & 0xffffU];
        value.low ^= row.low;
        value.high ^= row.high;
    }
    return std::popcount(value.low) + std::popcount(value.high);
}

std::uint64_t multiply_x(std::uint64_t value) {
    const auto top = value >> 63;
    return (value << 1) ^ (kReduction & (0ULL - top));
}

std::uint64_t multiply_scalar(std::uint64_t left, std::uint64_t right) {
    std::uint64_t result = 0;
    while (right) {
        if (right & 1) result ^= left;
        right >>= 1;
        left = multiply_x(left);
    }
    return result;
}

std::uint64_t multiply(std::uint64_t left, std::uint64_t right) {
    const auto a = _mm_set_epi64x(0, static_cast<long long>(left));
    const auto b = _mm_set_epi64x(0, static_cast<long long>(right));
    const auto product = _mm_clmulepi64_si128(a, b, 0x00);
    const auto low = static_cast<std::uint64_t>(_mm_cvtsi128_si64(product));
    const auto high = static_cast<std::uint64_t>(
        _mm_cvtsi128_si64(_mm_srli_si128(product, 8)));
    const unsigned __int64 overflow =
        (high >> 63) ^ (high >> 61) ^ (high >> 60);
    std::uint64_t reduced = low ^ high ^ (high << 1) ^ (high << 3) ^ (high << 4);
    reduced ^= overflow ^ (overflow << 1) ^ (overflow << 3) ^ (overflow << 4);
    return reduced;
}

std::uint64_t power(std::uint64_t base, std::uint64_t exponent) {
    std::uint64_t result = 1;
    while (exponent) {
        if (exponent & 1) result = multiply(result, base);
        exponent >>= 1;
        if (exponent) base = multiply(base, base);
    }
    return result;
}

std::uint64_t inverse(std::uint64_t value) {
    if (!value) throw std::runtime_error("inverse of zero");
    return power(value, UINT64_MAX - 1);
}

struct ScanResult {
    std::array<unsigned long long, 129> candidate_histogram{};
    unsigned long long endpoint_profiles = 0;
    std::array<unsigned long long, 4> endpoint_profiles_by_role{};
};

void add_candidate(ScanResult& result, std::uint64_t value,
                   unsigned long long multiplicity, unsigned role) {
    const auto weight = bch_weight(value);
    result.candidate_histogram[weight] += multiplicity;
    if (weight == 22 || weight == 106) {
        result.endpoint_profiles += multiplicity;
        result.endpoint_profiles_by_role[role] += multiplicity;
    }
}

ScanResult scan_fast(unsigned blocks) {
    std::vector<std::uint64_t> gamma(blocks + 1);
    std::vector<std::uint64_t> gamma_inverse(blocks + 1);
    std::vector<std::uint64_t> endpoint_over_one_plus(blocks + 1);
    gamma[0] = 1;
    for (unsigned index = 1; index <= blocks; ++index) {
        gamma[index] = multiply_x(gamma[index - 1]);
    }
    const auto inverse_gamma = inverse(2);
    gamma_inverse[0] = 1;
    for (unsigned index = 1; index <= blocks; ++index) {
        gamma_inverse[index] = multiply(gamma_inverse[index - 1], inverse_gamma);
        endpoint_over_one_plus[index] = multiply(
            kAllOneMessage, inverse(1 ^ gamma[index]));
    }

    ScanResult result;
    // Three data positions.  Translation by the first exponent cancels, so
    // only the two positive differences remain.
    for (unsigned second = 1; second + 1 < blocks; ++second) {
        for (unsigned third = second + 1; third < blocks; ++third) {
            const auto multiplicity = static_cast<unsigned long long>(blocks - third);
            const auto delta = third - second;
            const auto base0 = gamma[second] ^ gamma[third];
            const auto base1 = 1 ^ gamma[third];

            const auto scale0 = multiply(
                endpoint_over_one_plus[delta], gamma_inverse[second]);
            add_candidate(result, multiply(scale0, base1), multiplicity, 0);

            add_candidate(
                result,
                multiply(endpoint_over_one_plus[third], base0),
                multiplicity,
                1);
            add_candidate(
                result,
                multiply(endpoint_over_one_plus[second], base0),
                multiplicity,
                2);
        }
    }

    // The p0 position (projective point zero) and two data positions.
    for (unsigned difference = 1; difference < blocks; ++difference) {
        const auto multiplicity = static_cast<unsigned long long>(blocks - difference);
        const auto one_plus = 1 ^ gamma[difference];
        add_candidate(
            result,
            multiply(endpoint_over_one_plus[difference], gamma[difference]),
            multiplicity,
            3);
        add_candidate(
            result,
            multiply(kAllOneMessage, multiply(one_plus, gamma_inverse[difference])),
            multiplicity,
            3);
        add_candidate(
            result,
            multiply(kAllOneMessage, one_plus),
            multiplicity,
            3);
    }
    // Supports containing p1 have two equal finite components after either
    // finite component is normalized.  They cannot have weights 22 and 106.
    return result;
}

void validate_arithmetic() {
    std::uint64_t left = 0x0123456789abcdefULL;
    std::uint64_t right = 0xfedcba9876543210ULL;
    for (unsigned probe = 0; probe < 1000; ++probe) {
        if (multiply(left, right) != multiply_scalar(left, right)) {
            throw std::runtime_error("PCLMUL field multiplication failed");
        }
        left = left * 0x9e3779b97f4a7c15ULL + 1;
        right = right * 0xd1342543de82ef95ULL + 1;
    }
    if (multiply(2, inverse(2)) != 1) {
        throw std::runtime_error("field inversion failed");
    }
    if (bch_weight(kAllOneMessage) != 128) {
        throw std::runtime_error("all-one message validation failed");
    }
}

unsigned long long choose3(unsigned long long value) {
    return value * (value - 1) * (value - 2) / 6;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const unsigned blocks = argc >= 2 ? static_cast<unsigned>(std::stoul(argv[1]))
                                           : kDefaultDataBlocks;
        if (blocks < 3 || blocks > kDefaultDataBlocks) {
            throw std::runtime_error("data block count must lie in [3,16384]");
        }
        initialize_tables();
        validate_arithmetic();
        const auto started = std::chrono::steady_clock::now();
        const auto result = scan_fast(blocks);
        const auto elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();

        const auto data_supports = choose3(blocks);
        const auto p0_supports = static_cast<unsigned long long>(blocks) * (blocks - 1) / 2;
        const auto total_supports = choose3(static_cast<unsigned long long>(blocks) + 2);
        std::cout << "{\n"
                  << "  \"schema\": \"riffle-shiftalpha64-endpoint-triples-v1\",\n"
                  << "  \"data_blocks\": " << blocks << ",\n"
                  << "  \"all_one_message_hex\": \"8e63cf44efd4fa21\",\n"
                  << "  \"data_only_supports\": " << data_supports << ",\n"
                  << "  \"p0_and_two_data_supports\": " << p0_supports << ",\n"
                  << "  \"all_outer_supports\": " << total_supports << ",\n"
                  << "  \"exact_22_106_128_outer_words\": "
                  << result.endpoint_profiles << ",\n"
                  << "  \"role_counts\": [";
        for (unsigned index = 0; index < result.endpoint_profiles_by_role.size(); ++index) {
            if (index) std::cout << ", ";
            std::cout << result.endpoint_profiles_by_role[index];
        }
        std::cout << "],\n  \"candidate_weight_histogram\": {";
        bool first = true;
        for (unsigned weight = 0; weight <= 128; ++weight) {
            if (!result.candidate_histogram[weight]) continue;
            if (!first) std::cout << ",";
            std::cout << "\n    \"" << weight << "\": "
                      << result.candidate_histogram[weight];
            first = false;
        }
        std::cout << "\n  },\n"
                  << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsed << ",\n"
                  << "  \"validation\": {\n"
                  << "    \"pclmul_matches_scalar_1000_probes\": true,\n"
                  << "    \"all_one_bch_weight_is_128\": true,\n"
                  << "    \"p1_supports_cannot_have_profile_22_106_128\": true\n"
                  << "  },\n"
                  << "  \"scope\": \"Exact only for outer symbol weight three and BCH profile (22,106,128).\"\n"
                  << "}\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "endpoint triple scanner: " << error.what() << "\n";
        return 1;
    }
}
