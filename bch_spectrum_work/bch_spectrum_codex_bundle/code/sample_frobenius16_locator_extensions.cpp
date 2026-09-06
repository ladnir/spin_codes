#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

constexpr unsigned kFieldSize = 256;
using MulTable = std::array<std::array<std::uint8_t, kFieldSize>, kFieldSize>;

struct Orbit {
    std::array<std::uint8_t, 2> point{};
    unsigned size{};
};

struct BaseType {
    unsigned singles;
    unsigned pairs;
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

std::uint8_t power(std::uint8_t value, unsigned exponent,
                   const MulTable& multiply) noexcept {
    std::uint8_t result = 1;
    while (exponent != 0) {
        if (exponent & 1U) result = multiply[result][value];
        value = multiply[value][value];
        exponent >>= 1U;
    }
    return result;
}

std::uint64_t choose(unsigned n, unsigned k) noexcept {
    if (k > n) return 0;
    if (k > n - k) k = n - k;
    std::uint64_t result = 1;
    for (unsigned i = 1; i <= k; ++i) result = result * (n - k + i) / i;
    return result;
}

std::uint8_t evaluate(const std::array<std::uint8_t, 19>& coefficients,
                      std::uint8_t point, const MulTable& multiply) noexcept {
    std::uint8_t value = 0;
    for (int index = 18; index >= 0; --index)
        value = multiply[value][point] ^ coefficients[index];
    return value;
}

std::array<std::uint8_t, 19> interpolate(
    const std::array<std::uint8_t, 19>& nodes,
    const std::array<std::uint8_t, 19>& values,
    const std::array<std::uint8_t, 256>& inverse,
    const MulTable& multiply) noexcept {
    std::array<std::uint8_t, 19> divided = values;
    for (unsigned order = 1; order <= 18; ++order) {
        for (int index = 18; index >= static_cast<int>(order); --index) {
            const auto denominator = nodes[index] ^ nodes[index - order];
            divided[index] = multiply[divided[index] ^ divided[index - 1]]
                                     [inverse[denominator]];
        }
    }
    std::array<std::uint8_t, 19> polynomial{};
    polynomial[0] = divided[18];
    unsigned degree = 0;
    for (int index = 17; index >= 0; --index) {
        std::array<std::uint8_t, 19> updated{};
        for (unsigned j = 0; j <= degree; ++j) {
            updated[j] ^= multiply[polynomial[j]][nodes[index]];
            updated[j + 1] ^= polynomial[j];
        }
        updated[0] ^= divided[index];
        polynomial = updated;
        ++degree;
    }
    return polynomial;
}

template <std::size_t N>
void select_without_replacement(std::array<unsigned, N>& pool, unsigned count,
                                std::array<unsigned, N>& swaps,
                                std::mt19937_64& rng) noexcept {
    for (unsigned i = 0; i < count; ++i) {
        const unsigned selected = std::uniform_int_distribution<unsigned>(i, N - 1)(rng);
        swaps[i] = selected;
        std::swap(pool[i], pool[selected]);
    }
}

template <std::size_t N>
void restore_pool(std::array<unsigned, N>& pool, unsigned count,
                  const std::array<unsigned, N>& swaps) noexcept {
    for (int i = static_cast<int>(count) - 1; i >= 0; --i)
        std::swap(pool[i], pool[swaps[i]]);
}

std::uint64_t invariant_subset_count(unsigned singles, unsigned pairs,
                                     unsigned size) noexcept {
    std::uint64_t result = 0;
    for (unsigned selected_singles = size & 1U;
         selected_singles <= std::min(singles, size); selected_singles += 2) {
        const unsigned remaining = size - selected_singles;
        if (remaining / 2 <= pairs)
            result += choose(singles, selected_singles) *
                      choose(pairs, remaining / 2);
    }
    return result;
}

std::string coefficient_key(const std::array<std::uint8_t, 19>& coefficients) {
    std::ostringstream output;
    output << std::hex << std::setfill('0');
    for (const auto coefficient : coefficients)
        output << std::setw(2) << static_cast<unsigned>(coefficient);
    return output.str();
}

}  // namespace

int main(int argc, char** argv) try {
    if (argc != 4) {
        std::cerr << "usage: sample_frobenius16_locator_extensions SAMPLES SEED OUTPUT.json\n";
        return 2;
    }
    const std::uint64_t samples = std::stoull(argv[1]);
    const std::uint64_t seed = std::stoull(argv[2], nullptr, 0);
    const std::string output_path = argv[3];
    if (samples == 0) throw std::runtime_error("SAMPLES must be positive");

    MulTable multiply{};
    for (unsigned left = 0; left < 256; ++left)
        for (unsigned right = 0; right < 256; ++right)
            multiply[left][right] = multiply_slow(left, right);
    std::array<std::uint8_t, 256> inverse{};
    std::array<std::uint8_t, 256> target{};
    std::array<bool, 256> in_gf16{};
    for (unsigned value = 0; value < 256; ++value) {
        if (value) inverse[value] = power(static_cast<std::uint8_t>(value), 254, multiply);
        target[value] = value ? power(static_cast<std::uint8_t>(value), 146, multiply) : 0;
        in_gf16[value] = power(static_cast<std::uint8_t>(value), 16, multiply) == value;
    }

    std::array<bool, 256> seen{};
    std::vector<Orbit> singles;
    std::vector<Orbit> pairs;
    for (unsigned start = 1; start < 256; ++start) {
        if (seen[start]) continue;
        Orbit orbit;
        orbit.point[0] = static_cast<std::uint8_t>(start);
        orbit.point[1] = power(static_cast<std::uint8_t>(start), 16, multiply);
        orbit.size = orbit.point[1] == orbit.point[0] ? 1 : 2;
        seen[orbit.point[0]] = true;
        seen[orbit.point[1]] = true;
        (orbit.size == 1 ? singles : pairs).push_back(orbit);
    }
    if (singles.size() != 15 || pairs.size() != 120)
        throw std::runtime_error("unexpected x -> x^16 orbit structure");

    std::array<BaseType, 8> base_types{};
    std::uint64_t base_population = 0;
    for (unsigned index = 0; index < base_types.size(); ++index) {
        const unsigned singleton_count = 2 * index;
        const unsigned pair_count = (18 - singleton_count) / 2;
        const auto population = choose(15, singleton_count) * choose(120, pair_count);
        base_types[index] = {singleton_count, pair_count, population};
        base_population += population;
    }
    std::uint64_t target_population = 0;
    for (unsigned singleton_count = 0; singleton_count <= 14; singleton_count += 2)
        target_population += choose(15, singleton_count) *
                             choose(120, (24 - singleton_count) / 2);
    if (base_population != 199417781755185ULL ||
        target_population != 348770148711057905ULL)
        throw std::runtime_error("invariant-set population mismatch");

    std::array<unsigned, 15> singleton_pool{};
    std::array<unsigned, 120> pair_pool{};
    for (unsigned i = 0; i < singleton_pool.size(); ++i) singleton_pool[i] = i;
    for (unsigned i = 0; i < pair_pool.size(); ++i) pair_pool[i] = i;
    std::array<unsigned, 15> singleton_swaps{};
    std::array<unsigned, 120> pair_swaps{};
    std::mt19937_64 rng(seed);

    std::array<std::uint64_t, 8> type_samples{};
    std::array<std::uint64_t, 8> type_endpoint_hits{};
    std::array<std::uint64_t, 20> extra_histogram{};
    long double weighted_sum = 0;
    long double weighted_second_sum = 0;
    struct EndpointRecord { std::uint64_t multiplicity; unsigned fixed_points; };
    std::unordered_map<std::string, EndpointRecord> endpoints;

    const auto start_time = std::chrono::steady_clock::now();
    for (std::uint64_t sample = 0; sample < samples; ++sample) {
        const auto draw = std::uniform_int_distribution<std::uint64_t>(
            0, base_population - 1)(rng);
        std::uint64_t cumulative = 0;
        unsigned type_index = 0;
        while (type_index + 1 < base_types.size() &&
               draw >= cumulative + base_types[type_index].population) {
            cumulative += base_types[type_index].population;
            ++type_index;
        }
        const auto& type = base_types[type_index];
        ++type_samples[type_index];
        select_without_replacement(singleton_pool, type.singles, singleton_swaps, rng);
        select_without_replacement(pair_pool, type.pairs, pair_swaps, rng);

        std::array<bool, 15> selected_single{};
        std::array<bool, 120> selected_pair{};
        std::array<std::uint8_t, 19> nodes{};
        std::array<std::uint8_t, 19> values{};
        values[0] = 1;
        unsigned node_count = 1;
        for (unsigned i = 0; i < type.singles; ++i) {
            const auto index = singleton_pool[i];
            selected_single[index] = true;
            nodes[node_count] = singles[index].point[0];
            values[node_count] = target[nodes[node_count]];
            ++node_count;
        }
        for (unsigned i = 0; i < type.pairs; ++i) {
            const auto index = pair_pool[i];
            selected_pair[index] = true;
            for (unsigned member = 0; member < 2; ++member) {
                nodes[node_count] = pairs[index].point[member];
                values[node_count] = target[nodes[node_count]];
                ++node_count;
            }
        }
        if (node_count != 19) throw std::runtime_error("base size mismatch");
        const auto polynomial = interpolate(nodes, values, inverse, multiply);
        if (!std::all_of(polynomial.begin(), polynomial.end(),
                         [&](std::uint8_t coefficient) { return in_gf16[coefficient]; }))
            throw std::runtime_error("invariant interpolation left GF(16)");

        unsigned agreements = 0;
        unsigned remaining_singles = 0;
        unsigned remaining_pairs = 0;
        for (unsigned index = 0; index < singles.size(); ++index) {
            const bool agrees = evaluate(polynomial, singles[index].point[0], multiply) ==
                                target[singles[index].point[0]];
            agreements += agrees;
            remaining_singles += agrees && !selected_single[index];
        }
        for (unsigned index = 0; index < pairs.size(); ++index) {
            const bool first_agrees =
                evaluate(polynomial, pairs[index].point[0], multiply) ==
                target[pairs[index].point[0]];
            const bool second_agrees =
                evaluate(polynomial, pairs[index].point[1], multiply) ==
                target[pairs[index].point[1]];
            if (first_agrees != second_agrees)
                throw std::runtime_error("agreement set is not x -> x^16 invariant");
            agreements += 2 * first_agrees;
            remaining_pairs += first_agrees && !selected_pair[index];
        }
        if (agreements < 18 || agreements > 37)
            throw std::runtime_error("agreement count outside locator range");
        ++extra_histogram[agreements - 18];

        long double statistic = 0;
        for (unsigned added_singles = 0; added_singles <= 6; added_singles += 2) {
            const unsigned added_pairs = (6 - added_singles) / 2;
            if (added_singles > remaining_singles || added_pairs > remaining_pairs)
                continue;
            const unsigned target_singles = type.singles + added_singles;
            const unsigned target_pairs = type.pairs + added_pairs;
            const auto extensions = choose(remaining_singles, added_singles) *
                                    choose(remaining_pairs, added_pairs);
            const auto bases_per_target =
                invariant_subset_count(target_singles, target_pairs, 18);
            statistic += static_cast<long double>(extensions) / bases_per_target;
        }
        weighted_sum += statistic;
        weighted_second_sum += statistic * statistic;

        if (agreements == 37) {
            ++type_endpoint_hits[type_index];
            const auto key = coefficient_key(polynomial);
            auto [position, inserted] = endpoints.emplace(
                key, EndpointRecord{0, 0});
            ++position->second.multiplicity;
            if (inserted) {
                position->second.fixed_points = 0;
                for (const auto& orbit : singles)
                    position->second.fixed_points +=
                        evaluate(polynomial, orbit.point[0], multiply) ==
                        target[orbit.point[0]];
            }
        }

        restore_pool(singleton_pool, type.singles, singleton_swaps);
        restore_pool(pair_pool, type.pairs, pair_swaps);
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time).count();
    const long double mean = weighted_sum / samples;
    const long double second = weighted_second_sum / samples;
    const long double standard_error = std::sqrt(std::max(
        0.0L, (second - mean * mean) / samples));
    const long double estimated_sets = mean * base_population;
    const long double reference_sets =
        static_cast<long double>(target_population) / std::pow(16.0L, 6);
    const long double ratio = estimated_sets / reference_sets;
    const long double ratio_se = standard_error * base_population / reference_sets;

    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("cannot open output path");
    output << std::setprecision(17)
           << "{\n"
           << "  \"classification\": \"diagnostic invariant-base Monte Carlo; not a proof or upper bound\",\n"
           << "  \"samples\": " << samples << ",\n"
           << "  \"seed\": " << seed << ",\n"
           << "  \"seconds\": " << seconds << ",\n"
           << "  \"action\": \"x -> x^16 on GF(256)^*\",\n"
           << "  \"field_orbit_size_histogram\": {\"1\": 15, \"2\": 120},\n"
           << "  \"invariant_18_set_population\": " << base_population << ",\n"
           << "  \"invariant_24_set_population\": " << target_population << ",\n"
           << "  \"base_type_records\": [\n";
    for (unsigned i = 0; i < base_types.size(); ++i) {
        output << "    {\"singletons\": " << base_types[i].singles
               << ", \"pairs\": " << base_types[i].pairs
               << ", \"population\": " << base_types[i].population
               << ", \"samples\": " << type_samples[i]
               << ", \"endpoint_hits\": " << type_endpoint_hits[i] << "}"
               << (i + 1 == base_types.size() ? "\n" : ",\n");
    }
    output << "  ],\n  \"extra_agreement_histogram\": {";
    bool first = true;
    for (unsigned extra = 0; extra < extra_histogram.size(); ++extra) {
        if (extra_histogram[extra] == 0) continue;
        output << (first ? "\n" : ",\n") << "    \"" << extra << "\": "
               << extra_histogram[extra];
        first = false;
    }
    if (!first) output << '\n';
    output << "  },\n"
           << "  \"weighted_invariant_extension_estimator\": {\n"
           << "    \"sample_mean\": " << static_cast<double>(mean) << ",\n"
           << "    \"sample_standard_error\": " << static_cast<double>(standard_error) << ",\n"
           << "    \"estimated_admissible_24_sets\": " << static_cast<double>(estimated_sets) << ",\n"
           << "    \"independent_rank_reference_sets\": " << static_cast<double>(reference_sets) << ",\n"
           << "    \"ratio\": " << static_cast<double>(ratio) << ",\n"
           << "    \"ratio_standard_error\": " << static_cast<double>(ratio_se) << ",\n"
           << "    \"normal_95_upper_ratio\": " << static_cast<double>(ratio + 1.96L * ratio_se) << "\n"
           << "  },\n"
           << "  \"endpoint_hits\": "
           << std::accumulate(type_endpoint_hits.begin(), type_endpoint_hits.end(), std::uint64_t{0})
           << ",\n  \"distinct_endpoint_polynomials\": [";
    std::vector<std::pair<std::string, EndpointRecord>> ordered(endpoints.begin(), endpoints.end());
    std::sort(ordered.begin(), ordered.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    first = true;
    for (const auto& [key, record] : ordered) {
        output << (first ? "\n" : ",\n")
               << "    {\"coefficients_hex_ascending\": \"" << key
               << "\", \"sample_multiplicity\": " << record.multiplicity
               << ", \"agreement_fixed_points_in_GF16_star\": "
               << record.fixed_points << "}";
        first = false;
    }
    if (!first) output << '\n';
    output << "  ]\n}\n";
    std::cout << "samples=" << samples << " seconds=" << seconds
              << " ratio=" << static_cast<double>(ratio)
              << " endpoint_hits="
              << std::accumulate(type_endpoint_hits.begin(), type_endpoint_hits.end(), std::uint64_t{0})
              << " distinct_endpoints=" << endpoints.size() << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
