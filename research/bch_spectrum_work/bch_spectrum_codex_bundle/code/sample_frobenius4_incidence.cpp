#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <immintrin.h>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace {

constexpr unsigned kBatch = 32;
constexpr unsigned kFieldSize = 256;
constexpr std::array<std::array<std::uint8_t, 4>, 4> kMul4{{
    {{0, 0, 0, 0}}, {{0, 1, 2, 3}}, {{0, 2, 3, 1}}, {{0, 3, 1, 2}}
}};

using MulTable = std::array<std::array<std::uint8_t, kFieldSize>, kFieldSize>;

struct Factor {
    unsigned degree{};
    std::array<std::uint8_t, 5> coefficient{};
    std::vector<std::uint8_t> points;
};

struct CandidateType {
    unsigned singles;
    unsigned pairs;
    unsigned quads;
    std::uint64_t population;
};

std::uint8_t multiply_slow(unsigned left, unsigned right) noexcept {
    unsigned result = 0;
    while (right != 0) {
        if (right & 1U) result ^= left;
        right >>= 1U;
        left <<= 1U;
        if (left & 0x100U) left ^= 0x14dU;
    }
    return static_cast<std::uint8_t>(result);
}

std::uint64_t choose(unsigned n, unsigned k) noexcept {
    if (k > n) return 0;
    if (k > n - k) k = n - k;
    std::uint64_t result = 1;
    for (unsigned i = 1; i <= k; ++i) result = result * (n - k + i) / i;
    return result;
}

template <std::size_t N>
std::array<unsigned, N> choose_distinct(unsigned universe,
                                        std::mt19937_64& rng) noexcept {
    std::array<unsigned, N> result{};
    std::uniform_int_distribution<unsigned> distribution(0, universe - 1);
    for (unsigned i = 0; i < N; ++i) {
        unsigned value;
        bool duplicate;
        do {
            value = distribution(rng);
            duplicate = false;
            for (unsigned j = 0; j < i; ++j) duplicate |= result[j] == value;
        } while (duplicate);
        result[i] = value;
    }
    return result;
}

std::array<unsigned, 6> unrank_combination(unsigned universe, unsigned size,
                                           std::uint64_t rank) {
    std::array<unsigned, 6> result{};
    unsigned minimum = 0;
    for (unsigned position = 0; position < size; ++position) {
        bool selected = false;
        const unsigned remaining = size - position - 1;
        for (unsigned value = minimum; value + remaining < universe; ++value) {
            const auto suffixes = choose(universe - value - 1, remaining);
            if (rank < suffixes) {
                result[position] = value;
                minimum = value + 1;
                selected = true;
                break;
            }
            rank -= suffixes;
        }
        if (!selected) throw std::runtime_error("combination rank out of range");
    }
    if (rank != 0) throw std::runtime_error("combination rank residue");
    return result;
}

bool scalar_incidence(const std::array<std::uint8_t, 25>& elementary) noexcept {
    constexpr std::array<std::uint8_t, 4> square{{0, 1, 3, 2}};
    std::array<std::uint8_t, 128> homogeneous{};
    homogeneous[0] = 1;
    const auto compute_h = [&](unsigned n) {
        std::uint8_t value = 0;
        for (unsigned i = n & 1U; i <= std::min(24U, n); i += 2)
            value ^= kMul4[elementary[i]][square[homogeneous[(n - i) / 2]]];
        homogeneous[n] = value;
    };
    for (unsigned n = 1; n <= 31; ++n) compute_h(n);
    for (unsigned n = 49; n <= 63; ++n) compute_h(n);
    for (unsigned n = 122; n <= 127; ++n) compute_h(n);
    if (kMul4[elementary[24]][homogeneous[122]] != 1) return false;
    for (unsigned n = 123; n <= 127; ++n)
        if (homogeneous[n] != 0) return false;
    return true;
}

void multiply_factor(std::array<std::uint8_t, 25>& polynomial,
                     unsigned& degree, const Factor& factor) noexcept {
    std::array<std::uint8_t, 25> next{};
    for (unsigned i = 0; i <= degree; ++i) {
        const auto left = polynomial[i];
        for (unsigned j = 0; j <= factor.degree; ++j)
            next[i + j] ^= kMul4[left][factor.coefficient[j]];
    }
    degree += factor.degree;
    polynomial = next;
}

std::array<std::uint8_t, 19> monomial_remainder(
    const std::array<std::uint8_t, 25>& elementary) {
    // P(x)=x^24+e_1*x^23+...+e_24, with elementary stored as e_i.
    std::array<std::uint8_t, 24> modulus{};
    for (unsigned i = 0; i < 24; ++i) modulus[i] = elementary[24 - i];
    std::array<std::uint8_t, 24> remainder{};
    remainder[0] = 1;
    for (unsigned exponent = 0; exponent < 146; ++exponent) {
        const auto leading = remainder[23];
        for (int i = 23; i > 0; --i) remainder[i] = remainder[i - 1];
        remainder[0] = 0;
        if (leading)
            for (unsigned i = 0; i < 24; ++i)
                remainder[i] ^= kMul4[leading][modulus[i]];
    }
    if (remainder[0] != 1)
        throw std::runtime_error("incidence test/remainder constant mismatch");
    for (unsigned i = 19; i < 24; ++i)
        if (remainder[i] != 0)
            throw std::runtime_error("incidence test/remainder degree mismatch");
    std::array<std::uint8_t, 19> result{};
    std::copy_n(remainder.begin(), 19, result.begin());
    return result;
}

std::string coefficient_key(const std::array<std::uint8_t, 19>& coefficients,
                            const std::array<std::uint8_t, 4>& lift) {
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const auto coefficient : coefficients)
        output << std::setw(2) << static_cast<unsigned>(lift[coefficient]);
    return output.str();
}

unsigned agreement_count(const std::array<std::uint8_t, 19>& coefficients,
                         const std::array<std::uint8_t, 4>& lift,
                         const MulTable& multiply,
                         const std::array<std::uint8_t, 256>& target) noexcept {
    unsigned agreements = 0;
    for (unsigned point = 1; point < 256; ++point) {
        std::uint8_t value = 0;
        for (int i = 18; i >= 0; --i)
            value = multiply[value][point] ^ lift[coefficients[i]];
        agreements += value == target[point];
    }
    return agreements;
}

}  // namespace

int main(int argc, char** argv) try {
    if (argc != 4) {
        std::cerr << "usage: sample_frobenius4_incidence SAMPLES|all SEED OUTPUT.json\n";
        return 2;
    }
    const bool exhaustive = std::string(argv[1]) == "all";
    std::uint64_t requested_samples = exhaustive ? 0 : std::stoull(argv[1]);
    const std::uint64_t seed = std::stoull(argv[2], nullptr, 0);
    const std::string output_path = argv[3];
    if (!exhaustive && requested_samples == 0)
        throw std::runtime_error("SAMPLES must be positive");

    MulTable multiply{};
    for (unsigned left = 0; left < 256; ++left)
        for (unsigned right = 0; right < 256; ++right)
            multiply[left][right] = multiply_slow(left, right);

    std::array<std::uint8_t, 4> lift{};
    lift[0] = 0;
    lift[1] = 1;
    for (unsigned value = 2; value < 256; ++value) {
        const auto square = multiply[value][value];
        const auto fourth = multiply[square][square];
        if (fourth == value) {
            lift[2] = static_cast<std::uint8_t>(value);
            lift[3] = static_cast<std::uint8_t>(value ^ 1U);
            break;
        }
    }
    if (lift[2] == 0 || multiply[lift[2]][lift[2]] != lift[3])
        throw std::runtime_error("failed to identify GF(4) subfield");
    std::array<std::uint8_t, 256> project{};
    project.fill(0xff);
    for (unsigned code = 0; code < 4; ++code) project[lift[code]] = code;
    std::array<std::uint8_t, 256> target{};
    for (unsigned point = 1; point < 256; ++point) {
        std::uint8_t value = 1;
        std::uint8_t base = static_cast<std::uint8_t>(point);
        unsigned exponent = 146;
        while (exponent) {
            if (exponent & 1U) value = multiply[value][base];
            base = multiply[base][base];
            exponent >>= 1U;
        }
        target[point] = value;
    }

    std::array<bool, 256> seen{};
    std::array<std::vector<Factor>, 5> factors;
    for (unsigned start = 1; start < 256; ++start) {
        if (seen[start]) continue;
        Factor factor;
        unsigned current = start;
        do {
            factor.points.push_back(static_cast<std::uint8_t>(current));
            seen[current] = true;
            const auto square = multiply[current][current];
            current = multiply[square][square];
        } while (current != start);
        factor.degree = static_cast<unsigned>(factor.points.size());
        std::array<std::uint8_t, 5> over256{};
        over256[0] = 1;
        unsigned degree = 0;
        for (const auto point : factor.points) {
            for (int i = static_cast<int>(degree); i >= 0; --i)
                over256[i + 1] ^= multiply[over256[i]][point];
            ++degree;
        }
        for (unsigned i = 0; i <= degree; ++i) {
            if (project[over256[i]] == 0xff)
                throw std::runtime_error("orbit factor coefficient not in GF(4)");
            factor.coefficient[i] = project[over256[i]];
        }
        factors[factor.degree].push_back(factor);
    }
    if (factors[1].size() != 3 || factors[2].size() != 6 ||
        factors[4].size() != 60)
        throw std::runtime_error("unexpected x -> x^4 orbit structure");

    const std::array<CandidateType, 7> types{{
        {0, 0, 6, choose(60, 6)},
        {0, 2, 5, choose(6, 2) * choose(60, 5)},
        {0, 4, 4, choose(6, 4) * choose(60, 4)},
        {0, 6, 3, choose(60, 3)},
        {2, 1, 5, choose(3, 2) * choose(6, 1) * choose(60, 5)},
        {2, 3, 4, choose(3, 2) * choose(6, 3) * choose(60, 4)},
        {2, 5, 3, choose(3, 2) * choose(6, 5) * choose(60, 3)},
    }};
    std::uint64_t total_population = 0;
    for (const auto& type : types) total_population += type.population;
    if (total_population != 267516561ULL)
        throw std::runtime_error("candidate population mismatch");
    if (exhaustive) requested_samples = total_population;

    alignas(32) std::array<std::array<std::uint8_t, kBatch>, 25> elementary{};
    std::array<std::array<std::uint8_t, 25>, kBatch> lane_elementary{};
    std::array<__m256i, 128> homogeneous{};
    const __m256i one = _mm256_set1_epi8(1);
    const __m256i zero = _mm256_setzero_si256();
    const __m256i square_lookup = _mm256_setr_epi8(
        0, 1, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    const __m256i multiply_lookup = _mm256_setr_epi8(
        0, 0, 0, 0, 0, 1, 2, 3, 0, 2, 3, 1, 0, 3, 1, 2,
        0, 0, 0, 0, 0, 1, 2, 3, 0, 2, 3, 1, 0, 3, 1, 2);

    std::mt19937_64 rng(seed);
    std::array<std::uint64_t, 7> type_samples{};
    std::array<std::uint64_t, 7> type_hits{};
    std::array<std::uint64_t, 256> agreement_histogram{};
    struct PolynomialRecord { unsigned agreements; std::uint64_t multiplicity; };
    std::unordered_map<std::string, PolynomialRecord> polynomials;
    std::uint64_t samples = 0;
    std::uint64_t hits = 0;
    const auto start_time = std::chrono::steady_clock::now();

    while (samples < requested_samples) {
        const unsigned lanes = static_cast<unsigned>(
            std::min<std::uint64_t>(kBatch, requested_samples - samples));
        std::array<unsigned, kBatch> lane_type{};
        for (unsigned lane = 0; lane < lanes; ++lane) {
            const auto draw = exhaustive
                ? samples + lane
                : std::uniform_int_distribution<std::uint64_t>(
                      0, total_population - 1)(rng);
            std::uint64_t cumulative = 0;
            unsigned type_index = 0;
            while (type_index + 1 < types.size() &&
                   draw >= cumulative + types[type_index].population) {
                cumulative += types[type_index].population;
                ++type_index;
            }
            lane_type[lane] = type_index;
            ++type_samples[type_index];
            const auto& type = types[type_index];
            std::array<std::uint8_t, 25> poly{};
            poly[0] = 1;
            unsigned degree = 0;
            std::array<unsigned, 6> selected_singles{};
            std::array<unsigned, 6> selected_pairs{};
            std::array<unsigned, 6> selected_quads{};
            if (exhaustive) {
                auto local_rank = draw - cumulative;
                const auto quad_count = choose(60, type.quads);
                const auto pair_count = choose(6, type.pairs);
                const auto quad_rank = local_rank % quad_count;
                local_rank /= quad_count;
                const auto pair_rank = local_rank % pair_count;
                const auto single_rank = local_rank / pair_count;
                selected_singles = unrank_combination(3, type.singles, single_rank);
                selected_pairs = unrank_combination(6, type.pairs, pair_rank);
                selected_quads = unrank_combination(60, type.quads, quad_rank);
            } else {
                const auto random_singles = choose_distinct<2>(3, rng);
                const auto random_pairs = choose_distinct<6>(6, rng);
                const auto random_quads = choose_distinct<6>(60, rng);
                std::copy(random_singles.begin(), random_singles.end(),
                          selected_singles.begin());
                selected_pairs = random_pairs;
                selected_quads = random_quads;
            }
            for (unsigned i = 0; i < type.singles; ++i)
                multiply_factor(poly, degree, factors[1][selected_singles[i]]);
            for (unsigned i = 0; i < type.pairs; ++i)
                multiply_factor(poly, degree, factors[2][selected_pairs[i]]);
            for (unsigned i = 0; i < type.quads; ++i)
                multiply_factor(poly, degree, factors[4][selected_quads[i]]);
            if (degree != 24) throw std::runtime_error("candidate degree mismatch");
            lane_elementary[lane] = poly;
            for (unsigned i = 0; i <= 24; ++i) elementary[i][lane] = poly[i];
        }
        for (unsigned lane = lanes; lane < kBatch; ++lane)
            for (unsigned i = 0; i <= 24; ++i) elementary[i][lane] = 0;

        homogeneous.fill(zero);
        homogeneous[0] = one;
        const auto compute_h = [&](unsigned n) {
            __m256i accumulator = zero;
            for (unsigned i = n & 1U; i <= std::min(24U, n); i += 2) {
                const __m256i e = _mm256_load_si256(
                    reinterpret_cast<const __m256i*>(elementary[i].data()));
                const __m256i squared = _mm256_shuffle_epi8(
                    square_lookup, homogeneous[(n - i) / 2]);
                const __m256i index = _mm256_or_si256(e, _mm256_slli_epi16(squared, 2));
                accumulator = _mm256_xor_si256(
                    accumulator, _mm256_shuffle_epi8(multiply_lookup, index));
            }
            homogeneous[n] = accumulator;
        };
        for (unsigned n = 1; n <= 31; ++n) compute_h(n);
        for (unsigned n = 49; n <= 63; ++n) compute_h(n);
        for (unsigned n = 122; n <= 127; ++n) compute_h(n);

        const __m256i e24 = _mm256_load_si256(
            reinterpret_cast<const __m256i*>(elementary[24].data()));
        const __m256i product_index = _mm256_or_si256(
            e24, _mm256_slli_epi16(homogeneous[122], 2));
        __m256i valid = _mm256_cmpeq_epi8(
            _mm256_shuffle_epi8(multiply_lookup, product_index), one);
        for (unsigned n = 123; n <= 127; ++n)
            valid = _mm256_and_si256(valid, _mm256_cmpeq_epi8(homogeneous[n], zero));
        unsigned valid_mask = static_cast<unsigned>(_mm256_movemask_epi8(valid));
        if (lanes < kBatch) valid_mask &= (1U << lanes) - 1U;
        if (samples < 1024) {
            for (unsigned lane = 0; lane < lanes; ++lane) {
                const bool vector_result = ((valid_mask >> lane) & 1U) != 0;
                if (vector_result != scalar_incidence(lane_elementary[lane]))
                    throw std::runtime_error("scalar/AVX2 incidence mismatch");
            }
        }

        while (valid_mask != 0) {
            unsigned lane;
#if defined(_MSC_VER)
            unsigned long bit;
            _BitScanForward(&bit, valid_mask);
            lane = static_cast<unsigned>(bit);
#else
            lane = static_cast<unsigned>(__builtin_ctz(valid_mask));
#endif
            valid_mask &= valid_mask - 1;
            ++hits;
            ++type_hits[lane_type[lane]];
            const auto remainder = monomial_remainder(lane_elementary[lane]);
            const auto agreements = agreement_count(remainder, lift, multiply, target);
            if (agreements < 24 || agreements > 37)
                throw std::runtime_error("agreement count outside RS range");
            ++agreement_histogram[agreements];
            const auto key = coefficient_key(remainder, lift);
            auto [position, inserted] = polynomials.emplace(
                key, PolynomialRecord{agreements, 0});
            if (!inserted && position->second.agreements != agreements)
                throw std::runtime_error("polynomial agreement mismatch");
            ++position->second.multiplicity;
        }
        samples += lanes;
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time).count();

    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("failed to open output file");
    output << "{\n"
           << "  \"classification\": \""
           << (exhaustive ? "exact exhaustive x -> x^4 invariant 24-set computation"
                          : "Monte Carlo over uniform x -> x^4 invariant 24-sets")
           << "\",\n"
           << "  \"exhaustive\": " << (exhaustive ? "true" : "false") << ",\n"
           << "  \"samples\": " << samples << ",\n"
           << "  \"seed\": " << seed << ",\n"
           << "  \"seconds\": " << std::setprecision(12) << seconds << ",\n"
           << "  \"invariant_24_set_population\": " << total_population << ",\n"
           << "  \"gf4_lift_hex\": [\"00\", \"01\", \""
           << std::hex << std::setw(2) << std::setfill('0') << static_cast<unsigned>(lift[2])
           << "\", \"" << std::setw(2) << static_cast<unsigned>(lift[3])
           << "\"],\n" << std::dec << std::setfill(' ')
           << "  \"admissible_hits\": " << hits << ",\n"
           << "  \"candidate_types\": [\n";
    for (unsigned i = 0; i < types.size(); ++i) {
        const auto& type = types[i];
        output << "    {\"singles\": " << type.singles
               << ", \"pairs\": " << type.pairs
               << ", \"quads\": " << type.quads
               << ", \"population\": " << type.population
               << ", \"samples\": " << type_samples[i]
               << ", \"hits\": " << type_hits[i] << "}"
               << (i + 1 == types.size() ? "\n" : ",\n");
    }
    output << "  ],\n  \"agreement_count_histogram\": {";
    bool first = true;
    for (unsigned n = 0; n < agreement_histogram.size(); ++n) {
        if (agreement_histogram[n] == 0) continue;
        output << (first ? "\n" : ",\n") << "    \"" << n << "\": "
               << agreement_histogram[n];
        first = false;
    }
    if (!first) output << '\n';
    output << "  },\n  \"distinct_polynomials\": [";
    first = true;
    std::vector<std::pair<std::string, PolynomialRecord>> ordered(
        polynomials.begin(), polynomials.end());
    std::sort(ordered.begin(), ordered.end(),
              [](const auto& left, const auto& right) {
                  return left.first < right.first;
              });
    for (const auto& [key, record] : ordered) {
        output << (first ? "\n" : ",\n")
               << "    {\"coefficients_hex_ascending\": \"" << key
               << "\", \"agreements\": " << record.agreements
               << ", \"sample_multiplicity\": " << record.multiplicity << "}";
        first = false;
    }
    if (!first) output << '\n';
    output << "  ]\n}\n";
    std::cout << "sampled " << samples << " candidates; " << hits
              << " admissible; " << polynomials.size() << " distinct; "
              << seconds << " seconds\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
