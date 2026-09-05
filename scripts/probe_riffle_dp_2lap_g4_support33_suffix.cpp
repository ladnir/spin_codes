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
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned PacketSlots = 16;
    constexpr unsigned PacketSupport = 33;
    constexpr unsigned SampleSpace = 524'352;
    constexpr unsigned TotalNodes = 32'772;
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
            throw std::runtime_error("DP-2Lap probe: generator inverse failed");
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
                throw std::runtime_error("DP-2Lap probe: systematic construction failed");
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
        std::uint64_t seed = 20'260'819;
        std::string output =
            "explorations/riffle_dp_2lap_g4_support33_suffix_probe.json";
    };

    Options parseOptions(int argc, char** argv)
    {
        Options options;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            if (index + 1 >= argc)
                throw std::invalid_argument("DP-2Lap probe: missing option value");
            if (argument == "--length")
                options.length = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--trials")
                options.trials = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--seed")
                options.seed = std::stoull(argv[++index]);
            else if (argument == "--output")
                options.output = argv[++index];
            else
                throw std::invalid_argument("DP-2Lap probe: unknown option " + argument);
        }
        if (options.length == 0 || options.length > TotalNodes ||
            options.length * PacketSlots > SampleSpace || options.trials == 0)
            throw std::invalid_argument("DP-2Lap probe: invalid dimensions");
        return options;
    }

    struct LapResult
    {
        std::vector<std::uint64_t> output;
        std::uint64_t terminal = 0;
    };

    LapResult runLap(
        const std::vector<std::uint64_t>& input,
        std::uint64_t initialState,
        const SystematicRight& systematicRight)
    {
        LapResult result;
        result.output.resize(input.size());
        auto state = initialState;
        for (std::size_t node = 0; node < input.size(); ++node)
        {
            const auto output = accumulate(state ^ input[node]);
            result.output[node] = output;
            state = systematicRight(output);
        }
        result.terminal = state;
        return result;
    }

    unsigned verifyDecomposition(const SystematicRight& systematicRight)
    {
        std::mt19937_64 generator(0x4450324c41504d41ULL);
        constexpr std::array<unsigned, 7> Lengths{1, 2, 3, 7, 19, 64, 127};
        unsigned cases = 0;
        for (const auto length : Lengths)
        {
            std::vector<std::uint64_t> input(length);
            std::vector<std::uint64_t> zeros(length);
            for (unsigned trial = 0; trial < 64; ++trial)
            {
                for (auto& word : input)
                    word = generator();
                const auto first = runLap(input, 0, systematicRight);
                const auto literal = runLap(first.output, first.terminal, systematicRight);
                const auto zeroStart = runLap(first.output, 0, systematicRight);
                const auto autonomous = runLap(zeros, first.terminal, systematicRight);
                for (unsigned node = 0; node < length; ++node)
                {
                    if (literal.output[node] !=
                        (zeroStart.output[node] ^ autonomous.output[node]))
                        throw std::runtime_error(
                            "DP-2Lap probe: G=F^2+JL verification failed");
                }
                ++cases;
            }
        }
        return cases;
    }

    struct TwoLapMeasurement
    {
        unsigned firstWeight = 0;
        unsigned prefixWeight = 0;
        unsigned suffixWeight = 0;
        std::uint64_t firstTerminal = 0;
        std::uint64_t finalTerminal = 0;

        unsigned totalWeight() const { return prefixWeight + suffixWeight; }
    };

    TwoLapMeasurement measureSuffix(
        const std::vector<std::uint64_t>& drives,
        unsigned prefixLength,
        std::vector<std::uint64_t>& firstOutput,
        const SystematicRight& systematicRight)
    {
        TwoLapMeasurement result;
        auto state = std::uint64_t{0};
        for (std::size_t node = 0; node < drives.size(); ++node)
        {
            const auto output = accumulate(state ^ drives[node]);
            firstOutput[node] = output;
            result.firstWeight += std::popcount(output);
            state = systematicRight(output);
        }
        result.firstTerminal = state;

        for (unsigned node = 0; node < prefixLength; ++node)
        {
            const auto output = accumulate(state);
            result.prefixWeight += std::popcount(output);
            state = systematicRight(output);
        }
        for (const auto input : firstOutput)
        {
            const auto output = accumulate(state ^ input);
            result.suffixWeight += std::popcount(output);
            state = systematicRight(output);
        }
        result.finalTerminal = state;
        return result;
    }

    std::string hex64(std::uint64_t value)
    {
        std::ostringstream output;
        output << "0x" << std::hex << std::setw(16) << std::setfill('0') << value;
        return output.str();
    }

    unsigned autonomousLowerBound(unsigned length)
    {
        if (length == 0)
            return 0;
        return 9 * length - 8 - 8 * ((length - 1) / 3);
    }
}

int main(int argc, char** argv)
{
    try
    {
        const auto options = parseOptions(argc, argv);
        const auto prefixLength = TotalNodes - options.length;
        const SystematicRight systematicRight;
        const auto decompositionCases = verifyDecomposition(systematicRight);

        std::vector<std::uint64_t> clusterDrives(options.length);
        if (options.length >= 3)
        {
            clusterDrives[0] = 0x777507e004ee7e70ULL;
            clusterDrives[1] = 0xe507060705e0e4e0ULL;
            clusterDrives[2] = 0x0640e46704050554ULL;
        }
        std::vector<std::uint64_t> firstOutput(options.length);
        const auto cluster = measureSuffix(
            clusterDrives, prefixLength, firstOutput, systematicRight);

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

        std::uint64_t terminalZeroCount = 0;
        std::uint64_t firstFailureCount = 0;
        std::uint64_t finalFailureCount = 0;
        std::uint64_t finalAtDistanceCount = 0;
        std::uint64_t totalFirstWeight = 0;
        std::uint64_t totalPrefixWeight = 0;
        std::uint64_t totalSuffixWeight = 0;
        std::uint64_t totalFinalWeight = 0;
        unsigned minimumFirstWeight = std::numeric_limits<unsigned>::max();
        unsigned maximumFirstWeight = 0;
        unsigned minimumPrefixWeight = std::numeric_limits<unsigned>::max();
        unsigned maximumPrefixWeight = 0;
        unsigned minimumSuffixWeight = std::numeric_limits<unsigned>::max();
        unsigned maximumSuffixWeight = 0;
        unsigned minimumFinalWeight = std::numeric_limits<unsigned>::max();
        unsigned maximumFinalWeight = 0;

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

            const auto measurement = measureSuffix(
                drives, prefixLength, firstOutput, systematicRight);
            const auto finalWeight = measurement.totalWeight();
            terminalZeroCount += measurement.firstTerminal == 0;
            firstFailureCount += measurement.firstWeight < DistanceThreshold;
            finalFailureCount += finalWeight < DistanceThreshold;
            finalAtDistanceCount += finalWeight == DistanceThreshold;
            totalFirstWeight += measurement.firstWeight;
            totalPrefixWeight += measurement.prefixWeight;
            totalSuffixWeight += measurement.suffixWeight;
            totalFinalWeight += finalWeight;
            minimumFirstWeight = std::min(minimumFirstWeight, measurement.firstWeight);
            maximumFirstWeight = std::max(maximumFirstWeight, measurement.firstWeight);
            minimumPrefixWeight = std::min(minimumPrefixWeight, measurement.prefixWeight);
            maximumPrefixWeight = std::max(maximumPrefixWeight, measurement.prefixWeight);
            minimumSuffixWeight = std::min(minimumSuffixWeight, measurement.suffixWeight);
            maximumSuffixWeight = std::max(maximumSuffixWeight, measurement.suffixWeight);
            minimumFinalWeight = std::min(minimumFinalWeight, finalWeight);
            maximumFinalWeight = std::max(maximumFinalWeight, finalWeight);
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
            static_cast<long double>(finalFailureCount) / options.trials;
        const auto zeroFailureUpper95 = finalFailureCount == 0
            ? 1 - std::pow(0.05L, 1.0L / options.trials)
            : std::numeric_limits<long double>::quiet_NaN();

        std::ofstream output(options.output);
        output << std::setprecision(18)
            << "{\n"
            << "  \"schema\": \"riffle-dp-2lap-g4-support33-suffix-probe-v1\",\n"
            << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n"
            << "  \"evidence_label\": \"DIAGNOSTIC\",\n"
            << "  \"model\": \"G(x)=F(F(x)) xor J(L(x)); the second lap retains the first-lap terminal state and reads the overwritten first-lap word\",\n"
            << "  \"decomposition_verification_cases\": " << decompositionCases << ",\n"
            << "  \"total_nodes\": " << TotalNodes << ",\n"
            << "  \"window_length\": " << options.length << ",\n"
            << "  \"prefix_length\": " << prefixLength << ",\n"
            << "  \"window_cells\": " << windowCells << ",\n"
            << "  \"trials\": " << options.trials << ",\n"
            << "  \"seed\": " << options.seed << ",\n"
            << "  \"distance_threshold\": " << DistanceThreshold << ",\n"
            << "  \"terminal_zero_count\": " << terminalZeroCount << ",\n"
            << "  \"first_lap_failure_count_below_distance\": " << firstFailureCount << ",\n"
            << "  \"final_failure_count_below_distance\": " << finalFailureCount << ",\n"
            << "  \"final_at_distance_count\": " << finalAtDistanceCount << ",\n"
            << "  \"conditional_failure_estimate\": " << conditionalFailure << ",\n";
        if (finalFailureCount == 0)
            output << "  \"zero_failure_binomial_upper_95\": " << zeroFailureUpper95 << ",\n";
        else
            output << "  \"zero_failure_binomial_upper_95\": null,\n";
        output
            << "  \"certified_autonomous_prefix_lower_bound_for_nonzero_terminal\": "
            << autonomousLowerBound(prefixLength) << ",\n"
            << "  \"first_lap_weight\": {\"minimum\": " << minimumFirstWeight
            << ", \"maximum\": " << maximumFirstWeight
            << ", \"mean\": "
            << static_cast<long double>(totalFirstWeight) / options.trials << "},\n"
            << "  \"second_lap_prefix_weight\": {\"minimum\": " << minimumPrefixWeight
            << ", \"maximum\": " << maximumPrefixWeight
            << ", \"mean\": "
            << static_cast<long double>(totalPrefixWeight) / options.trials << "},\n"
            << "  \"second_lap_suffix_weight\": {\"minimum\": " << minimumSuffixWeight
            << ", \"maximum\": " << maximumSuffixWeight
            << ", \"mean\": "
            << static_cast<long double>(totalSuffixWeight) / options.trials << "},\n"
            << "  \"final_weight\": {\"minimum\": " << minimumFinalWeight
            << ", \"maximum\": " << maximumFinalWeight
            << ", \"mean\": "
            << static_cast<long double>(totalFinalWeight) / options.trials << "},\n"
            << "  \"known_cluster_replay\": {\n"
            << "    \"first_lap_weight\": " << cluster.firstWeight << ",\n"
            << "    \"first_terminal_state\": \"" << hex64(cluster.firstTerminal) << "\",\n"
            << "    \"second_lap_prefix_weight\": " << cluster.prefixWeight << ",\n"
            << "    \"second_lap_suffix_weight\": " << cluster.suffixWeight << ",\n"
            << "    \"final_weight\": " << cluster.totalWeight() << ",\n"
            << "    \"final_terminal_state\": \"" << hex64(cluster.finalTerminal) << "\"\n"
            << "  },\n"
            << "  \"log2_window_event_probability\": " << log2WindowProbability << ",\n"
            << "  \"scope_limitation\": \"The recurrence and cluster replay are exact. The pseudorandom Monte Carlo sample is diagnostic, not a rigorous probability bound.\"\n"
            << "}\n";
        if (!output)
            throw std::runtime_error("DP-2Lap probe: output write failed");

        std::cout << "candidate=Riffle DP-2Lap g=4\n";
        std::cout << "decomposition_cases=" << decompositionCases << '\n';
        std::cout << "length=" << options.length << '\n';
        std::cout << "prefix_length=" << prefixLength << '\n';
        std::cout << "trials=" << options.trials << '\n';
        std::cout << "terminal_zero_count=" << terminalZeroCount << '\n';
        std::cout << "first_lap_failures=" << firstFailureCount << '\n';
        std::cout << "final_failures=" << finalFailureCount << '\n';
        std::cout << "minimum_final_weight=" << minimumFinalWeight << '\n';
        std::cout << "mean_final_weight="
            << static_cast<double>(static_cast<long double>(totalFinalWeight) / options.trials)
            << '\n';
        std::cout << "cluster_first_weight=" << cluster.firstWeight << '\n';
        std::cout << "cluster_final_weight=" << cluster.totalWeight() << '\n';
        std::cout << "output=" << options.output << '\n';
        std::cout << "status=DIAGNOSTIC_DP_2LAP_SUPPORT33_SUFFIX\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
