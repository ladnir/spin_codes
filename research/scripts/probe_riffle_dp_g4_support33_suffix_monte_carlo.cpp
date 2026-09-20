#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned PacketSlots = 16;
    constexpr unsigned PacketSupport = 33;
    constexpr unsigned SampleSpace = 524'352;
    constexpr unsigned DistanceThreshold = 188'766;

    struct Word128
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
    };

    std::uint64_t carrylessMultiplyLow(std::uint64_t left, std::uint64_t right)
    {
        std::uint64_t result = 0;
        while (right)
        {
            const auto bit = std::countr_zero(right);
            result ^= left << bit;
            right &= right - 1;
        }
        return result;
    }

    std::uint64_t generatorInverse()
    {
        std::uint64_t inverse = 1;
        for (unsigned bit = 1; bit < 64; ++bit)
        {
            if ((carrylessMultiplyLow(Generator, inverse) >> bit) & 1)
                inverse |= std::uint64_t{1} << bit;
        }
        if (carrylessMultiplyLow(Generator, inverse) != 1)
            throw std::runtime_error("suffix probe: generator inverse failed");
        return inverse;
    }

    Word128 encode(std::uint64_t message)
    {
        Word128 result{};
        while (message)
        {
            const auto bit = std::countr_zero(message);
            result.low ^= Generator << bit;
            result.high ^=
                (bit == 0 ? 0 : Generator >> (64 - bit)) |
                (std::uint64_t{1} << 63);
            message &= message - 1;
        }
        return result;
    }

    std::array<std::uint64_t, 64> systematicRightColumns()
    {
        const auto inverse = generatorInverse();
        std::array<std::uint64_t, 64> columns{};
        for (unsigned column = 0; column < 64; ++column)
        {
            const auto message = carrylessMultiplyLow(std::uint64_t{1} << column, inverse);
            const auto word = encode(message);
            if (word.low != (std::uint64_t{1} << column))
                throw std::runtime_error("suffix probe: systematic construction failed");
            columns[column] = word.high;
        }
        return columns;
    }

    std::uint64_t accumulate(std::uint64_t value)
    {
        value ^= value << 1;
        value ^= value << 2;
        value ^= value << 4;
        value ^= value << 8;
        value ^= value << 16;
        value ^= value << 32;
        return value;
    }

    class SystematicRight
    {
    public:
        SystematicRight()
        {
            const auto columns = systematicRightColumns();
            for (unsigned byteIndex = 0; byteIndex < 8; ++byteIndex)
            {
                for (unsigned byte = 0; byte < 256; ++byte)
                {
                    std::uint64_t value = 0;
                    for (unsigned bit = 0; bit < 8; ++bit)
                    {
                        if ((byte >> bit) & 1)
                            value ^= columns[8 * byteIndex + bit];
                    }
                    mTables[byteIndex][byte] = value;
                }
            }
        }

        std::uint64_t operator()(std::uint64_t output) const
        {
            return
                mTables[0][output & 0xff] ^
                mTables[1][(output >> 8) & 0xff] ^
                mTables[2][(output >> 16) & 0xff] ^
                mTables[3][(output >> 24) & 0xff] ^
                mTables[4][(output >> 32) & 0xff] ^
                mTables[5][(output >> 40) & 0xff] ^
                mTables[6][(output >> 48) & 0xff] ^
                mTables[7][output >> 56];
        }

    private:
        std::array<std::array<std::uint64_t, 256>, 8> mTables{};
    };

    struct Options
    {
        unsigned length = 5'940;
        unsigned trials = 100'000;
        std::uint64_t seed = 20'260'818;
        std::string output = "explorations/riffle_dp_g4_g2_support33_suffix_monte_carlo.json";
    };

    Options parseOptions(int argc, char** argv)
    {
        Options options;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            if (index + 1 >= argc)
                throw std::invalid_argument("suffix probe: missing option value");
            if (argument == "--length")
                options.length = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--trials")
                options.trials = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--seed")
                options.seed = std::stoull(argv[++index]);
            else if (argument == "--output")
                options.output = argv[++index];
            else
                throw std::invalid_argument("suffix probe: unknown option " + argument);
        }
        if (options.length == 0 || options.length * PacketSlots > SampleSpace || options.trials == 0)
            throw std::invalid_argument("suffix probe: invalid dimensions");
        return options;
    }
}

int main(int argc, char** argv)
{
    try
    {
        const auto options = parseOptions(argc, argv);
        const SystematicRight systematicRight;
        std::mt19937_64 generator(options.seed);
        std::uniform_int_distribution<unsigned> locationDistribution(
            0,
            options.length * PacketSlots - 1);
        std::array<unsigned, PacketSupport> locations{};
        std::array<unsigned, PacketSupport> values{
            4,4,4,4,4,4,
            5,5,5,5,5,5,
            6,6,6,
            7,7,7,7,7,7,7,7,7,
            14,14,14,14,14,14,14,14,14};
        std::vector<std::uint64_t> drives(options.length);
        std::uint64_t failureCount = 0;
        std::uint64_t totalWeight = 0;
        unsigned minimumWeight = std::numeric_limits<unsigned>::max();
        unsigned maximumWeight = 0;

        for (unsigned trial = 0; trial < options.trials; ++trial)
        {
            std::fill(drives.begin(), drives.end(), 0);
            for (unsigned packet = 0; packet < PacketSupport; ++packet)
            {
                while (true)
                {
                    const auto candidate = locationDistribution(generator);
                    bool duplicate = false;
                    for (unsigned previous = 0; previous < packet; ++previous)
                        duplicate |= locations[previous] == candidate;
                    if (!duplicate)
                    {
                        locations[packet] = candidate;
                        break;
                    }
                }
            }
            std::shuffle(values.begin(), values.end(), generator);
            for (unsigned packet = 0; packet < PacketSupport; ++packet)
            {
                const auto node = locations[packet] / PacketSlots;
                const auto slot = locations[packet] % PacketSlots;
                drives[node] |= std::uint64_t{values[packet]} << (4 * slot);
            }

            std::uint64_t state = 0;
            unsigned weight = 0;
            for (unsigned node = 0; node < options.length; ++node)
            {
                const auto output = accumulate(state ^ drives[node]);
                weight += std::popcount(output);
                state = systematicRight(output);
            }
            failureCount += weight <= DistanceThreshold;
            totalWeight += weight;
            minimumWeight = std::min(minimumWeight, weight);
            maximumWeight = std::max(maximumWeight, weight);
        }

        long double log2WindowProbability = 0;
        const auto windowCells = options.length * PacketSlots;
        for (unsigned packet = 0; packet < PacketSupport; ++packet)
        {
            log2WindowProbability += std::log2(
                static_cast<long double>(windowCells - packet) /
                static_cast<long double>(SampleSpace - packet));
        }
        const auto conditionalFailure =
            static_cast<long double>(failureCount) / options.trials;
        const auto log2SingleContribution = failureCount == 0
            ? -std::numeric_limits<long double>::infinity()
            : log2WindowProbability + std::log2(conditionalFailure);
        const auto log2ThreeFamilyContribution =
            log2SingleContribution + std::log2(static_cast<long double>(3));

        std::ofstream output(options.output);
        output << std::setprecision(18)
            << "{\n"
            << "  \"schema\": \"riffle-dp-g4-g2-support33-suffix-monte-carlo-v1\",\n"
            << "  \"candidate\": \"Riffle DP g=4@g0-v1\",\n"
            << "  \"evidence_label\": \"DIAGNOSTIC\",\n"
            << "  \"window_length\": " << options.length << ",\n"
            << "  \"window_cells\": " << windowCells << ",\n"
            << "  \"trials\": " << options.trials << ",\n"
            << "  \"seed\": " << options.seed << ",\n"
            << "  \"failure_count\": " << failureCount << ",\n"
            << "  \"conditional_failure_estimate\": " << conditionalFailure << ",\n"
            << "  \"minimum_weight\": " << minimumWeight << ",\n"
            << "  \"maximum_weight\": " << maximumWeight << ",\n"
            << "  \"mean_weight\": "
            << static_cast<long double>(totalWeight) / options.trials << ",\n"
            << "  \"log2_window_event_probability\": " << log2WindowProbability << ",\n"
            << "  \"log2_single_family_contribution_estimate\": " << log2SingleContribution << ",\n"
            << "  \"log2_three_family_contribution_estimate\": " << log2ThreeFamilyContribution << ",\n"
            << "  \"scope_limitation\": \"Pseudorandom Monte Carlo estimates the conditional failure probability. It is not a rigorous lower or upper bound.\"\n"
            << "}\n";
        if (!output)
            throw std::runtime_error("suffix probe: output write failed");

        std::cout << "candidate=Riffle DP g=4@g0-v1\n";
        std::cout << "length=" << options.length << '\n';
        std::cout << "trials=" << options.trials << '\n';
        std::cout << "failures=" << failureCount << '\n';
        std::cout << "conditional_failure=" << static_cast<double>(conditionalFailure) << '\n';
        std::cout << "minimum_weight=" << minimumWeight << '\n';
        std::cout << "mean_weight=" << static_cast<double>(static_cast<long double>(totalWeight) / options.trials) << '\n';
        std::cout << "log2_three_family_contribution=" << static_cast<double>(log2ThreeFamilyContribution) << '\n';
        std::cout << "output=" << options.output << '\n';
        std::cout << "status=DIAGNOSTIC_SUPPORT33_SUFFIX_MONTE_CARLO\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
