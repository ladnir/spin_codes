#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

namespace {

constexpr unsigned kFieldSize = 256;

using MulTable = std::array<std::array<std::uint8_t, kFieldSize>, kFieldSize>;

std::uint8_t multiply_slow(unsigned left, unsigned right, unsigned field_size,
                           unsigned modulus) noexcept {
    unsigned result = 0;
    while (right != 0) {
        if (right & 1U) result ^= left;
        right >>= 1U;
        left <<= 1U;
        if (left & field_size) left ^= modulus;
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
    // Keep this fixed-width loop: it is the dominant kernel, and the compiler
    // fully unrolls it for the q=256 experiment.  Lower-degree families leave
    // their unused leading coefficients equal to zero.
    for (int index = 18; index >= 0; --index)
        value = multiply[value][point] ^ coefficients[index];
    return value;
}

std::array<std::uint8_t, 19> interpolate(
    const std::array<std::uint8_t, 19>& nodes,
    const std::array<std::uint8_t, 19>& values,
    const std::array<std::uint8_t, 256>& inverse,
    const MulTable& multiply, unsigned base_agreements) noexcept {
    // Newton divided differences, followed by a conversion to monomials.
    std::array<std::uint8_t, 19> divided = values;
    for (unsigned order = 1; order <= base_agreements; ++order) {
        for (int index = static_cast<int>(base_agreements);
             index >= static_cast<int>(order); --index) {
            const auto denominator = nodes[index] ^ nodes[index - order];
            divided[index] = multiply[divided[index] ^ divided[index - 1]]
                                     [inverse[denominator]];
        }
    }

    std::array<std::uint8_t, 19> polynomial{};
    polynomial[0] = divided[base_agreements];
    unsigned degree = 0;
    for (int index = static_cast<int>(base_agreements) - 1; index >= 0; --index) {
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

}  // namespace

int main(int argc, char** argv) try {
    if (argc != 3 && argc != 9) {
        std::cerr
            << "usage: sample_locator_extensions SAMPLES OUTPUT.json "
               "[FIELD_BITS MODULUS BASE_AGREEMENTS MAX_AGREEMENTS "
               "TARGET_EXPONENT SEED]\n";
        return 2;
    }
    const std::uint64_t samples = std::stoull(argv[1]);
    const std::string output_path = argv[2];
    const unsigned field_bits = argc == 9 ? std::stoul(argv[3], nullptr, 0) : 8;
    const unsigned modulus = argc == 9 ? std::stoul(argv[4], nullptr, 0) : 0x14d;
    const unsigned base_agreements = argc == 9 ? std::stoul(argv[5], nullptr, 0) : 18;
    const unsigned maximum_agreements = argc == 9 ? std::stoul(argv[6], nullptr, 0) : 37;
    const unsigned target_exponent = argc == 9 ? std::stoul(argv[7], nullptr, 0) : 146;
    const std::uint64_t seed = argc == 9 ? std::stoull(argv[8], nullptr, 0)
                                        : 0xb37c0deULL;
    const unsigned field_size = 1U << field_bits;
    if (samples == 0) throw std::runtime_error("SAMPLES must be positive");
    if (field_bits < 2 || field_bits > 8 || base_agreements > 18 ||
        maximum_agreements != 2 * base_agreements + 1 ||
        maximum_agreements >= field_size)
        throw std::runtime_error("invalid locator-family parameters");
    const unsigned expected_exponent =
        (maximum_agreements * (field_size / 2)) % (field_size - 1);
    if (target_exponent != expected_exponent)
        throw std::runtime_error("target exponent does not match locator substitution");
    const unsigned maximum_moment =
        std::min(8U, maximum_agreements - base_agreements);
    const unsigned remaining_points = field_size - 1 - base_agreements;

    MulTable multiply{};
    for (unsigned left = 0; left < field_size; ++left)
        for (unsigned right = 0; right < field_size; ++right)
            multiply[left][right] =
                multiply_slow(left, right, field_size, modulus);

    std::array<std::uint8_t, 256> inverse{};
    for (unsigned value = 1; value < field_size; ++value)
        inverse[value] =
            power(static_cast<std::uint8_t>(value), field_size - 2, multiply);

    std::array<std::uint8_t, 256> target{};
    for (unsigned value = 1; value < field_size; ++value)
        target[value] =
            power(static_cast<std::uint8_t>(value), target_exponent, multiply);
    const bool has_middle_subfield = field_bits % 2 == 0;
    const unsigned middle_subfield_size =
        has_middle_subfield ? (1U << (field_bits / 2)) : 0;
    std::array<bool, 256> in_middle_subfield{};
    if (has_middle_subfield) {
        for (unsigned value = 0; value < field_size; ++value)
            in_middle_subfield[value] =
                power(static_cast<std::uint8_t>(value), middle_subfield_size,
                      multiply) == value;
    }
    const bool has_affine_ridge = field_size == 256 && target_exponent == 146;
    std::array<bool, 256> in_affine_ridge_image{};
    if (has_affine_ridge)
        for (unsigned value = 1; value < field_size; ++value)
            in_affine_ridge_image[power(static_cast<std::uint8_t>(value), 130,
                                        multiply)] = true;

    std::array<std::uint8_t, 255> pool{};
    for (unsigned i = 0; i < field_size - 1; ++i)
        pool[i] = static_cast<std::uint8_t>(i + 1);
    std::array<unsigned, 18> swap_indices{};
    std::mt19937_64 random(seed);

    std::array<std::uint64_t, 20> extra_histogram{};
    std::array<long double, 9> empirical_moments{};
    std::array<long double, 9> empirical_second_moments{};
    std::array<std::uint64_t, 19> subfield_stratum_samples{};
    std::array<long double, 19> subfield_stratum_sixth_moment{};
    std::array<std::uint64_t, 19> degree_stratum_samples{};
    std::array<long double, 19> degree_stratum_sixth_moment{};
    std::array<long double, 19> degree_stratum_sixth_second_moment{};
    std::array<std::uint64_t, 2> ridge_image_stratum_samples{};
    std::array<long double, 2> ridge_image_stratum_sixth_moment{};
    std::array<long double, 2> ridge_image_stratum_sixth_second_moment{};
    const auto start = std::chrono::steady_clock::now();
    for (std::uint64_t sample = 0; sample < samples; ++sample) {
        std::array<std::uint8_t, 19> nodes{};
        std::array<std::uint8_t, 19> values{};
        values[0] = 1;  // g(0)=1.
        unsigned subfield_points = 0;
        for (unsigned i = 0; i < base_agreements; ++i) {
            std::uniform_int_distribution<unsigned> distribution(
                i, field_size - 2);
            const unsigned selected = distribution(random);
            swap_indices[i] = selected;
            std::swap(pool[i], pool[selected]);
            nodes[i + 1] = pool[i];
            values[i + 1] = target[pool[i]];
            subfield_points += in_middle_subfield[pool[i]];
        }

        const auto polynomial =
            interpolate(nodes, values, inverse, multiply, base_agreements);
        unsigned agreements = 0;
        for (unsigned point = 1; point < field_size; ++point)
            agreements += evaluate(polynomial, static_cast<std::uint8_t>(point),
                                   multiply) == target[point];
        if (agreements < base_agreements ||
            agreements > maximum_agreements)
            throw std::runtime_error("agreement count violates locator-degree bound");
        const unsigned extra = agreements - base_agreements;
        const long double sixth_moment = choose(extra, 6);
        ++extra_histogram[extra];
        unsigned polynomial_degree = base_agreements;
        while (polynomial_degree != 0 && polynomial[polynomial_degree] == 0)
            --polynomial_degree;
        ++degree_stratum_samples[polynomial_degree];
        degree_stratum_sixth_moment[polynomial_degree] += sixth_moment;
        degree_stratum_sixth_second_moment[polynomial_degree] +=
            sixth_moment * sixth_moment;
        if (has_affine_ridge) {
            const unsigned ridge_image = in_affine_ridge_image[polynomial[1]];
            ++ridge_image_stratum_samples[ridge_image];
            ridge_image_stratum_sixth_moment[ridge_image] += sixth_moment;
            ridge_image_stratum_sixth_second_moment[ridge_image] +=
                sixth_moment * sixth_moment;
        }
        for (unsigned order = 0; order <= maximum_moment; ++order) {
            const long double moment = choose(extra, order);
            empirical_moments[order] += moment;
            empirical_second_moments[order] += moment * moment;
        }
        if (has_middle_subfield) {
            ++subfield_stratum_samples[subfield_points];
            subfield_stratum_sixth_moment[subfield_points] += choose(extra, 6);
        }

        for (int i = static_cast<int>(base_agreements) - 1; i >= 0; --i)
            std::swap(pool[i], pool[swap_indices[i]]);
    }
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();

    std::ofstream output(output_path);
    if (!output) throw std::runtime_error("cannot open output path");
    output << std::setprecision(17);
    output << "{\n"
           << "  \"classification\": \"diagnostic Monte Carlo; not a proof or upper bound\",\n"
           << "  \"samples\": " << samples << ",\n"
           << "  \"seed\": " << seed << ",\n"
           << "  \"field_bits\": " << field_bits << ",\n"
           << "  \"field_size\": " << field_size << ",\n"
           << "  \"field_modulus\": " << modulus << ",\n"
           << "  \"target_exponent\": " << target_exponent << ",\n"
           << "  \"base_points_per_interpolant\": " << base_agreements << ",\n"
           << "  \"remaining_evaluation_points\": " << remaining_points << ",\n"
           << "  \"agreement_count_upper\": " << maximum_agreements << ",\n"
           << "  \"seconds\": " << seconds << ",\n"
           << "  \"extra_agreement_histogram\": {\n";
    bool first = true;
    for (unsigned extra = 0; extra < extra_histogram.size(); ++extra) {
        if (extra_histogram[extra] == 0) continue;
        if (!first) output << ",\n";
        first = false;
        output << "    \"" << extra << "\": " << extra_histogram[extra];
    }
    output << "\n  },\n  \"binomial_moment_comparison\": [\n";
    for (unsigned order = 0; order <= maximum_moment; ++order) {
        const long double empirical = empirical_moments[order] / samples;
        const long double empirical_second =
            empirical_second_moments[order] / samples;
        const long double standard_error = std::sqrt(std::max(
            0.0L, (empirical_second - empirical * empirical) / samples));
        const long double reference =
            static_cast<long double>(choose(remaining_points, order)) /
            std::pow(static_cast<long double>(field_size),
                     static_cast<int>(order));
        if (order != 0) output << ",\n";
        output << "    {\"order\": " << order
               << ", \"empirical\": " << static_cast<double>(empirical)
               << ", \"standard_error\": " << static_cast<double>(standard_error)
               << ", \"binomial_reference\": " << static_cast<double>(reference)
               << ", \"ratio\": " << static_cast<double>(empirical / reference)
               << ", \"normal_95_upper_ratio\": "
               << static_cast<double>((empirical + 1.96L * standard_error) / reference)
               << "}";
    }
    output << "\n  ]";
    if (has_middle_subfield) {
        const long double sixth_reference =
            static_cast<long double>(choose(remaining_points, 6)) /
            std::pow(static_cast<long double>(field_size), 6);
        output << ",\n  \"middle_subfield\": {\n"
               << "    \"field_size\": " << middle_subfield_size << ",\n"
               << "    \"strata\": [\n";
        bool first_stratum = true;
        for (unsigned count = 0; count <= base_agreements; ++count) {
            if (subfield_stratum_samples[count] == 0) continue;
            if (!first_stratum) output << ",\n";
            first_stratum = false;
            const long double conditional =
                subfield_stratum_sixth_moment[count] /
                subfield_stratum_samples[count];
            const long double contribution_ratio =
                (subfield_stratum_sixth_moment[count] / samples) /
                sixth_reference;
            output << "      {\"base_points_in_subfield\": " << count
                   << ", \"samples\": " << subfield_stratum_samples[count]
                   << ", \"conditional_sixth_moment\": "
                   << static_cast<double>(conditional)
                   << ", \"global_reference_contribution_ratio\": "
                   << static_cast<double>(contribution_ratio) << "}";
        }
        output << "\n    ]\n  }";
    }
    {
        const long double sixth_reference =
            static_cast<long double>(choose(remaining_points, 6)) /
            std::pow(static_cast<long double>(field_size), 6);
        output << ",\n  \"interpolant_degree\": {\n"
               << "    \"strata\": [\n";
        bool first_stratum = true;
        for (unsigned degree = 0; degree <= base_agreements; ++degree) {
            if (degree_stratum_samples[degree] == 0) continue;
            if (!first_stratum) output << ",\n";
            first_stratum = false;
            const long double conditional =
                degree_stratum_sixth_moment[degree] /
                degree_stratum_samples[degree];
            const long double conditional_second =
                degree_stratum_sixth_second_moment[degree] /
                degree_stratum_samples[degree];
            const long double standard_error = std::sqrt(std::max(
                0.0L, (conditional_second - conditional * conditional) /
                          degree_stratum_samples[degree]));
            const long double contribution_ratio =
                (degree_stratum_sixth_moment[degree] / samples) /
                sixth_reference;
            output << "      {\"degree\": " << degree
                   << ", \"samples\": " << degree_stratum_samples[degree]
                   << ", \"conditional_sixth_moment\": "
                   << static_cast<double>(conditional)
                   << ", \"conditional_standard_error\": "
                   << static_cast<double>(standard_error)
                   << ", \"conditional_ratio\": "
                   << static_cast<double>(conditional / sixth_reference)
                   << ", \"normal_95_upper_conditional_ratio\": "
                   << static_cast<double>(
                          (conditional + 1.96L * standard_error) / sixth_reference)
                   << ", \"global_reference_contribution_ratio\": "
                   << static_cast<double>(contribution_ratio) << "}";
        }
        output << "\n    ]\n  }";
    }
    if (has_affine_ridge) {
        const long double sixth_reference =
            static_cast<long double>(choose(remaining_points, 6)) /
            std::pow(static_cast<long double>(field_size), 6);
        output << ",\n  \"affine_ridge_image_by_linear_coefficient\": {\n"
               << "    \"ridge_condition\": \"c_1 is a^130 for some nonzero a\",\n"
               << "    \"strata\": [\n";
        for (unsigned in_image = 0; in_image < 2; ++in_image) {
            if (in_image) output << ",\n";
            const long double conditional =
                ridge_image_stratum_sixth_moment[in_image] /
                ridge_image_stratum_samples[in_image];
            const long double conditional_second =
                ridge_image_stratum_sixth_second_moment[in_image] /
                ridge_image_stratum_samples[in_image];
            const long double standard_error = std::sqrt(std::max(
                0.0L, (conditional_second - conditional * conditional) /
                          ridge_image_stratum_samples[in_image]));
            const long double contribution_ratio =
                (ridge_image_stratum_sixth_moment[in_image] / samples) /
                sixth_reference;
            output << "      {\"in_ridge_image\": "
                   << (in_image ? "true" : "false")
                   << ", \"samples\": " << ridge_image_stratum_samples[in_image]
                   << ", \"conditional_sixth_moment\": "
                   << static_cast<double>(conditional)
                   << ", \"conditional_standard_error\": "
                   << static_cast<double>(standard_error)
                   << ", \"conditional_ratio\": "
                   << static_cast<double>(conditional / sixth_reference)
                   << ", \"normal_95_upper_conditional_ratio\": "
                   << static_cast<double>(
                          (conditional + 1.96L * standard_error) / sixth_reference)
                   << ", \"global_reference_contribution_ratio\": "
                   << static_cast<double>(contribution_ratio) << "}";
        }
        output << "\n    ]\n  }";
    }
    output << "\n}\n";
    std::cout << "samples=" << samples << " seconds=" << seconds
              << " max_extra=";
    for (int extra = 19; extra >= 0; --extra) {
        if (extra_histogram[extra] != 0) {
            std::cout << extra << '\n';
            break;
        }
    }
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
