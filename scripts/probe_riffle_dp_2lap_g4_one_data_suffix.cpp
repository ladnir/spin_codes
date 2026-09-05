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
    using u8 = std::uint8_t;
    using u16 = std::uint16_t;
    using u64 = std::uint64_t;

    constexpr u64 FieldReduction = 0x1b;
    constexpr u64 LocalGenerator = 0xF4845518B9582A1FULL;
    constexpr unsigned PacketSlots = 16;
    constexpr unsigned SampleSpace = 524'352;
    constexpr unsigned TotalNodes = 32'772;
    constexpr unsigned DistanceThreshold = 188'766;
    constexpr unsigned MaximumSupport = 39;

    struct Word128
    {
        u64 low = 0;
        u64 high = 0;
    };

    u64 multiplyX(u64 value)
    {
        const auto high = value >> 63;
        return (value << 1) ^ (FieldReduction & (u64{0} - high));
    }

    Word128 encodeLocal(u64 message)
    {
        Word128 result{};
        while (message)
        {
            const auto bit = std::countr_zero(message);
            result.low ^= LocalGenerator << bit;
            result.high ^=
                (bit == 0 ? 0 : LocalGenerator >> (64 - bit)) |
                (u64{1} << 63);
            message &= message - 1;
        }
        return result;
    }

    void appendNibbles(
        const Word128& word,
        std::array<u8, MaximumSupport>& values,
        unsigned& size)
    {
        for (unsigned slot = 0; slot < 16; ++slot)
        {
            const auto value = static_cast<u8>((word.low >> (4 * slot)) & 0xf);
            if (value)
                values[size++] = value;
        }
        for (unsigned slot = 0; slot < 16; ++slot)
        {
            const auto value = static_cast<u8>((word.high >> (4 * slot)) & 0xf);
            if (value)
                values[size++] = value;
        }
    }

    u64 carrylessMultiplyLow(u64 left, u64 right)
    {
        u64 result = 0;
        while (right)
        {
            const auto bit = std::countr_zero(right);
            result ^= left << bit;
            right &= right - 1;
        }
        return result;
    }

    u64 generatorInverse()
    {
        u64 inverse = 1;
        for (unsigned bit = 1; bit < 64; ++bit)
        {
            if ((carrylessMultiplyLow(LocalGenerator, inverse) >> bit) & 1)
                inverse |= u64{1} << bit;
        }
        if (carrylessMultiplyLow(LocalGenerator, inverse) != 1)
            throw std::runtime_error("DP-2Lap one-data probe: inverse failed");
        return inverse;
    }

    std::array<u64, 64> systematicRightColumns()
    {
        const auto inverse = generatorInverse();
        std::array<u64, 64> columns{};
        for (unsigned column = 0; column < 64; ++column)
        {
            const auto message = carrylessMultiplyLow(u64{1} << column, inverse);
            const auto word = encodeLocal(message);
            if (word.low != (u64{1} << column))
                throw std::runtime_error(
                    "DP-2Lap one-data probe: systematic construction failed");
            columns[column] = word.high;
        }
        return columns;
    }

    u64 accumulate(u64 value)
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
                    for (unsigned bit = 0; bit < 8; ++bit)
                    {
                        if ((byte >> bit) & 1)
                            mTables[byteIndex][byte] ^= columns[8 * byteIndex + bit];
                    }
                }
            }
        }

        u64 operator()(u64 output) const
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
        std::array<std::array<u64, 256>, 8> mTables{};
    };

    struct Profile
    {
        u8 support = 0;
        std::array<u8, MaximumSupport> values{};
    };

    struct Options
    {
        std::string records = "explorations/riffle_dp_g4_one_data_support39.bin";
        std::string output =
            "explorations/riffle_dp_2lap_g4_one_data_suffix_probe.json";
        unsigned length = 5'940;
        unsigned trialsPerSupport = 100'000;
        u64 seed = 20'260'820;
    };

    Options parseOptions(int argc, char** argv)
    {
        Options options;
        for (int index = 1; index < argc; ++index)
        {
            if (index + 1 >= argc)
                throw std::invalid_argument("DP-2Lap one-data probe: missing option value");
            const std::string argument = argv[index];
            if (argument == "--records")
                options.records = argv[++index];
            else if (argument == "--output")
                options.output = argv[++index];
            else if (argument == "--length")
                options.length = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--trials-per-support")
                options.trialsPerSupport = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--seed")
                options.seed = std::stoull(argv[++index]);
            else
                throw std::invalid_argument(
                    "DP-2Lap one-data probe: unknown option " + argument);
        }
        if (!options.length || options.length > TotalNodes ||
            options.length * PacketSlots > SampleSpace || !options.trialsPerSupport)
            throw std::invalid_argument("DP-2Lap one-data probe: invalid dimensions");
        return options;
    }

    struct TrialResult
    {
        bool terminalZero = false;
        bool firstFailure = false;
        bool finalFailure = false;
        bool prefixCrossed = false;
        unsigned crossingNode = 0;
        unsigned firstWeight = 0;
        unsigned evaluatedFinalWeight = 0;
    };

    TrialResult runTrial(
        std::vector<u64>& drives,
        unsigned prefixLength,
        const SystematicRight& systematicRight)
    {
        TrialResult result;
        u64 state = 0;
        for (auto& drive : drives)
        {
            const auto output = accumulate(state ^ drive);
            drive = output;
            result.firstWeight += std::popcount(output);
            state = systematicRight(output);
        }
        result.terminalZero = state == 0;
        result.firstFailure = result.firstWeight < DistanceThreshold;

        unsigned weight = 0;
        for (unsigned node = 0; node < prefixLength; ++node)
        {
            const auto output = accumulate(state);
            weight += std::popcount(output);
            state = systematicRight(output);
            if (weight >= DistanceThreshold)
            {
                result.prefixCrossed = true;
                result.crossingNode = node + 1;
                result.evaluatedFinalWeight = weight;
                return result;
            }
        }
        for (const auto input : drives)
        {
            const auto output = accumulate(state ^ input);
            weight += std::popcount(output);
            state = systematicRight(output);
        }
        result.evaluatedFinalWeight = weight;
        result.finalFailure = weight < DistanceThreshold;
        return result;
    }
}

int main(int argc, char** argv)
{
    try
    {
        const auto options = parseOptions(argc, argv);
        const auto prefixLength = TotalNodes - options.length;
        std::ifstream records(options.records, std::ios::binary);
        if (!records)
            throw std::runtime_error("DP-2Lap one-data probe: could not open records");

        std::array<std::vector<Profile>, MaximumSupport + 1> profiles;
        while (true)
        {
            u16 exponent = 0;
            u64 parameter = 0;
            u8 declaredSupport = 0;
            records.read(reinterpret_cast<char*>(&exponent), sizeof(exponent));
            if (!records)
            {
                if (records.eof())
                    break;
                throw std::runtime_error("DP-2Lap one-data probe: record read failed");
            }
            records.read(reinterpret_cast<char*>(&parameter), sizeof(parameter));
            records.read(reinterpret_cast<char*>(&declaredSupport), sizeof(declaredSupport));
            if (!records || declaredSupport > MaximumSupport)
                throw std::runtime_error("DP-2Lap one-data probe: invalid record");

            auto weighted = parameter;
            for (unsigned index = 0; index < exponent; ++index)
                weighted = multiplyX(weighted);
            Profile profile;
            unsigned actualSupport = 0;
            appendNibbles(encodeLocal(parameter), profile.values, actualSupport);
            appendNibbles(encodeLocal(parameter), profile.values, actualSupport);
            appendNibbles(encodeLocal(weighted), profile.values, actualSupport);
            if (actualSupport != declaredSupport)
                throw std::runtime_error(
                    "DP-2Lap one-data probe: profile support mismatch");
            profile.support = declaredSupport;
            profiles[profile.support].push_back(profile);
        }

        const SystematicRight systematicRight;
        std::mt19937_64 generator(options.seed);
        const auto windowCells = options.length * PacketSlots;
        std::uniform_int_distribution<unsigned> locationDistribution(0, windowCells - 1);
        std::vector<u64> drives(options.length);
        std::array<unsigned, MaximumSupport> locations{};
        const auto nonemptySupportCount = static_cast<unsigned>(std::count_if(
            profiles.begin(), profiles.end(),
            [](const auto& supportProfiles) { return !supportProfiles.empty(); }));
        const auto simultaneousAlphaPerSupport = 0.05L / nonemptySupportCount;
        long double aggregateSimultaneousUpper95 = 0;
        std::size_t totalPopulation = 0;

        std::ofstream output(options.output);
        output << std::setprecision(18)
            << "{\n"
            << "  \"schema\": \"riffle-dp-2lap-g4-one-data-suffix-probe-v1\",\n"
            << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n"
            << "  \"evidence_label\": \"DIAGNOSTIC\",\n"
            << "  \"records\": \"" << options.records << "\",\n"
            << "  \"total_nodes\": " << TotalNodes << ",\n"
            << "  \"window_length\": " << options.length << ",\n"
            << "  \"prefix_length\": " << prefixLength << ",\n"
            << "  \"window_cells\": " << windowCells << ",\n"
            << "  \"trials_per_nonempty_support\": " << options.trialsPerSupport << ",\n"
            << "  \"seed\": " << options.seed << ",\n"
            << "  \"support_rows\": [\n";

        bool firstRow = true;
        for (unsigned support = 0; support <= MaximumSupport; ++support)
        {
            if (profiles[support].empty())
                continue;
            totalPopulation += profiles[support].size();
            std::uniform_int_distribution<std::size_t> profileDistribution(
                0, profiles[support].size() - 1);
            u64 terminalZeroCount = 0;
            u64 firstFailureCount = 0;
            u64 finalFailureCount = 0;
            u64 prefixCrossingCount = 0;
            u64 totalCrossingNode = 0;
            unsigned earliestCrossingNode = std::numeric_limits<unsigned>::max();
            unsigned latestCrossingNode = 0;
            unsigned minimumFirstWeight = std::numeric_limits<unsigned>::max();
            unsigned minimumEvaluatedFinalWeight = std::numeric_limits<unsigned>::max();

            for (unsigned trial = 0; trial < options.trialsPerSupport; ++trial)
            {
                std::fill(drives.begin(), drives.end(), 0);
                const auto& profile = profiles[support][profileDistribution(generator)];
                for (unsigned packet = 0; packet < support; ++packet)
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
                auto values = profile.values;
                std::shuffle(values.begin(), values.begin() + support, generator);
                for (unsigned packet = 0; packet < support; ++packet)
                {
                    const auto node = locations[packet] / PacketSlots;
                    const auto slot = locations[packet] % PacketSlots;
                    drives[node] |= u64{values[packet]} << (4 * slot);
                }

                const auto result = runTrial(drives, prefixLength, systematicRight);
                terminalZeroCount += result.terminalZero;
                firstFailureCount += result.firstFailure;
                finalFailureCount += result.finalFailure;
                minimumFirstWeight = std::min(minimumFirstWeight, result.firstWeight);
                minimumEvaluatedFinalWeight = std::min(
                    minimumEvaluatedFinalWeight, result.evaluatedFinalWeight);
                if (result.prefixCrossed)
                {
                    ++prefixCrossingCount;
                    totalCrossingNode += result.crossingNode;
                    earliestCrossingNode = std::min(
                        earliestCrossingNode, result.crossingNode);
                    latestCrossingNode = std::max(latestCrossingNode, result.crossingNode);
                }
            }

            long double log2Window = 0;
            for (unsigned packet = 0; packet < support; ++packet)
            {
                log2Window += std::log2(
                    static_cast<long double>(windowCells - packet) /
                    static_cast<long double>(SampleSpace - packet));
            }
            const auto failureRate =
                static_cast<long double>(finalFailureCount) / options.trialsPerSupport;
            const auto zeroFailureSimultaneousUpper95 = finalFailureCount == 0
                ? 1 - std::pow(
                    simultaneousAlphaPerSupport,
                    1.0L / options.trialsPerSupport)
                : std::numeric_limits<long double>::quiet_NaN();
            const auto contribution = static_cast<long double>(profiles[support].size())
                * std::exp2(log2Window) * failureRate;
            const auto upperContribution95 = finalFailureCount == 0
                ? static_cast<long double>(profiles[support].size())
                    * std::exp2(log2Window) * zeroFailureSimultaneousUpper95
                : std::numeric_limits<long double>::quiet_NaN();
            if (finalFailureCount == 0)
                aggregateSimultaneousUpper95 += upperContribution95;

            if (!firstRow)
                output << ",\n";
            firstRow = false;
            output << "    {\"support\": " << support
                << ", \"population\": " << profiles[support].size()
                << ", \"terminal_zero_count\": " << terminalZeroCount
                << ", \"first_lap_failures_below_distance\": " << firstFailureCount
                << ", \"final_failures_below_distance\": " << finalFailureCount
                << ", \"conditional_failure_estimate\": " << failureRate
                << ", \"prefix_crossing_count\": " << prefixCrossingCount
                << ", \"earliest_prefix_crossing_node\": " << earliestCrossingNode
                << ", \"latest_prefix_crossing_node\": " << latestCrossingNode
                << ", \"mean_prefix_crossing_node\": "
                << static_cast<long double>(totalCrossingNode) / prefixCrossingCount
                << ", \"minimum_first_lap_weight\": " << minimumFirstWeight
                << ", \"minimum_evaluated_final_weight\": "
                << minimumEvaluatedFinalWeight
                << ", \"log2_window_event_probability\": " << log2Window;
            if (contribution > 0)
                output << ", \"log2_contribution_estimate\": "
                    << std::log2(contribution);
            else
                output << ", \"log2_contribution_estimate\": null";
            if (finalFailureCount == 0)
                output << ", \"zero_failure_binomial_upper_simultaneous_95\": "
                    << zeroFailureSimultaneousUpper95
                    << ", \"log2_contribution_upper_simultaneous_95\": "
                    << std::log2(upperContribution95);
            else
                output << ", \"zero_failure_binomial_upper_simultaneous_95\": null"
                    << ", \"log2_contribution_upper_simultaneous_95\": null";
            output << "}";

            std::cout << "support=" << support
                << " population=" << profiles[support].size()
                << " terminal_zero=" << terminalZeroCount
                << " first_failures=" << firstFailureCount
                << " final_failures=" << finalFailureCount
                << " latest_crossing=" << latestCrossingNode
                << '\n';
        }
        output << "\n  ],\n"
            << "  \"total_outer_population\": " << totalPopulation << ",\n"
            << "  \"log2_aggregate_contribution_estimate\": null,\n"
            << "  \"log2_aggregate_contribution_upper_simultaneous_95\": "
            << std::log2(aggregateSimultaneousUpper95) << ",\n"
            << "  \"scope_limitation\": \"The one-data outer population and packet profiles are exact. Pseudorandom sampling estimates the conditional inner behavior. The simultaneous 95% value is a Bonferroni binomial confidence bound, not a theorem bound for the construction. Prefix early exit preserves exact failure detection.\"\n"
            << "}\n";
        if (!output)
            throw std::runtime_error("DP-2Lap one-data probe: output write failed");

        std::cout << "total_outer_population=" << totalPopulation << '\n';
        std::cout << "aggregate_log2_contribution_upper_simultaneous_95="
            << static_cast<double>(std::log2(aggregateSimultaneousUpper95)) << '\n';
        std::cout << "output=" << options.output << '\n';
        std::cout << "status=DIAGNOSTIC_DP_2LAP_EXACT_ONE_DATA_POPULATION\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
