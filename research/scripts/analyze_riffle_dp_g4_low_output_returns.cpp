#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned InnerNodes = 32'772;
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
            throw std::runtime_error("generator inverse construction failed");
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
            const auto message = carrylessMultiplyLow(
                std::uint64_t{1} << column,
                inverse);
            const auto word = encode(message);
            if (word.low != (std::uint64_t{1} << column))
                throw std::runtime_error("systematic column construction failed");
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

    class OutputMap
    {
    public:
        OutputMap()
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
            const auto parity =
                mTables[0][output & 0xff] ^
                mTables[1][(output >> 8) & 0xff] ^
                mTables[2][(output >> 16) & 0xff] ^
                mTables[3][(output >> 24) & 0xff] ^
                mTables[4][(output >> 32) & 0xff] ^
                mTables[5][(output >> 40) & 0xff] ^
                mTables[6][(output >> 48) & 0xff] ^
                mTables[7][output >> 56];
            return accumulate(parity);
        }

    private:
        std::array<std::array<std::uint64_t, 256>, 8> mTables{};
    };

    std::uint64_t nextCombination(std::uint64_t value)
    {
        const auto low = value & (~value + 1);
        const auto ripple = value + low;
        if (ripple == 0)
            return 0;
        return ripple | (((ripple ^ value) >> 2) / low);
    }
}

int main(int argc, char** argv)
{
    unsigned lowWeightMaximum = 5;
    unsigned maximumGap = 20;
    std::string outputPath =
        "explorations/riffle_dp_g4_g2_low_output_returns_gap20.json";
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        if (argument == "--low-weight" && index + 1 < argc)
            lowWeightMaximum = static_cast<unsigned>(std::stoul(argv[++index]));
        else if (argument == "--max-gap" && index + 1 < argc)
            maximumGap = static_cast<unsigned>(std::stoul(argv[++index]));
        else if (argument == "--output" && index + 1 < argc)
            outputPath = argv[++index];
        else
            throw std::invalid_argument("unknown or incomplete argument: " + argument);
    }
    if (maximumGap == 0 || maximumGap > 32771)
        throw std::invalid_argument("maximum gap is out of range");
    if (lowWeightMaximum == 0 || lowWeightMaximum > 10)
        throw std::invalid_argument("low weight maximum is out of range");

    const OutputMap step;
    std::vector<std::uint64_t> returnCounts(maximumGap + 1);
    std::vector<std::uint64_t> firstReturnCounts(maximumGap + 1);
    std::vector<unsigned> minimumFirstSegmentWeights(
        maximumGap + 1,
        std::numeric_limits<unsigned>::max());
    struct Example
    {
        std::uint64_t start;
        std::uint64_t finish;
        unsigned gap;
    };
    std::vector<Example> examples;
    std::uint64_t stateCount = 0;
    std::uint64_t noFirstReturnThroughBound = 0;
    std::int64_t minimumFirstSegmentExcessNumerator =
        std::numeric_limits<std::int64_t>::max();
    Example minimumFirstSegmentExample{};
    unsigned minimumFirstSegmentWeight = 0;
    for (unsigned weight = 1; weight <= lowWeightMaximum; ++weight)
    {
        auto start = (std::uint64_t{1} << weight) - 1;
        while (start)
        {
            ++stateCount;
            auto current = start;
            unsigned segmentWeight = std::popcount(start);
            bool foundFirstReturn = false;
            for (unsigned gap = 1; gap <= maximumGap; ++gap)
            {
                current = step(current);
                const auto currentWeight = std::popcount(current);
                if (currentWeight <= lowWeightMaximum)
                {
                    ++returnCounts[gap];
                    if (examples.size() < 100)
                        examples.push_back(Example{start, current, gap});
                    if (!foundFirstReturn)
                    {
                        foundFirstReturn = true;
                        ++firstReturnCounts[gap];
                        minimumFirstSegmentWeights[gap] = std::min(
                            minimumFirstSegmentWeights[gap],
                            segmentWeight);
                        const auto excess =
                            static_cast<std::int64_t>(segmentWeight) * InnerNodes -
                            static_cast<std::int64_t>(DistanceThreshold) * gap;
                        if (excess < minimumFirstSegmentExcessNumerator)
                        {
                            minimumFirstSegmentExcessNumerator = excess;
                            minimumFirstSegmentExample = Example{start, current, gap};
                            minimumFirstSegmentWeight = segmentWeight;
                        }
                    }
                }
                else if (!foundFirstReturn)
                    segmentWeight += currentWeight;
            }
            if (!foundFirstReturn)
                ++noFirstReturnThroughBound;
            start = nextCombination(start);
        }
    }

    unsigned minimumReturnGap = 0;
    for (unsigned gap = 1; gap <= maximumGap; ++gap)
    {
        if (returnCounts[gap])
        {
            minimumReturnGap = gap;
            break;
        }
    }
    std::ofstream output(outputPath, std::ios::trunc);
    if (!output)
        throw std::runtime_error("failed to open output path");
    output << "{\n"
        << "  \"schema\": \"riffle-dp-g4-g2-low-output-returns-v1\",\n"
        << "  \"candidate\": \"Riffle DP g=4@g0-v1\",\n"
        << "  \"evidence_label\": \"EXACT\",\n"
        << "  \"low_output_weight_maximum\": " << lowWeightMaximum << ",\n"
        << "  \"high_output_weight_minimum\": " << lowWeightMaximum + 1 << ",\n"
        << "  \"maximum_gap\": " << maximumGap << ",\n"
        << "  \"enumerated_low_states\": " << stateCount << ",\n"
        << "  \"minimum_return_gap\": " << minimumReturnGap << ",\n"
        << "  \"no_first_return_through_bound\": " << noFirstReturnThroughBound << ",\n"
        << "  \"minimum_first_return_segment\": {\"start_hex\":\"0x"
        << std::hex << minimumFirstSegmentExample.start
        << "\",\"next_low_hex\":\"0x" << minimumFirstSegmentExample.finish
        << std::dec << "\",\"gap\":" << minimumFirstSegmentExample.gap
        << ",\"segment_weight\":" << minimumFirstSegmentWeight
        << ",\"excess_numerator_over_d_over_N\":"
        << minimumFirstSegmentExcessNumerator << "},\n"
        << "  \"return_counts_by_gap\": [";
    for (unsigned gap = 1; gap <= maximumGap; ++gap)
        output << (gap == 1 ? "" : ",") << returnCounts[gap];
    output << "],\n  \"first_return_counts_by_gap\": [";
    for (unsigned gap = 1; gap <= maximumGap; ++gap)
        output << (gap == 1 ? "" : ",") << firstReturnCounts[gap];
    output << "],\n  \"minimum_first_return_segment_weight_by_gap\": [";
    for (unsigned gap = 1; gap <= maximumGap; ++gap)
    {
        if (gap != 1)
            output << ',';
        if (minimumFirstSegmentWeights[gap] == std::numeric_limits<unsigned>::max())
            output << "null";
        else
            output << minimumFirstSegmentWeights[gap];
    }
    output << "],\n  \"minimum_first_return_excess_over_high_baseline_by_gap\": [";
    for (unsigned gap = 1; gap <= maximumGap; ++gap)
    {
        if (gap != 1)
            output << ',';
        if (minimumFirstSegmentWeights[gap] == std::numeric_limits<unsigned>::max())
            output << "null";
        else
            output << static_cast<int>(minimumFirstSegmentWeights[gap]) -
                static_cast<int>((lowWeightMaximum + 1) * gap);
    }
    output << "],\n  \"examples\": [";
    for (std::size_t index = 0; index < examples.size(); ++index)
    {
        const auto& example = examples[index];
        output << (index ? ",\n    {" : "\n    {")
            << "\"start_hex\":\"0x" << std::hex << example.start
            << "\",\"finish_hex\":\"0x" << example.finish << std::dec
            << "\",\"gap\":" << example.gap << "}";
    }
    if (!examples.empty())
        output << '\n';
    output << "  ]\n}\n";
    output.close();

    std::cout << "candidate=Riffle DP g=4@g0-v1\n";
    std::cout << "low_output_weight_maximum=" << lowWeightMaximum << '\n';
    std::cout << "enumerated_low_states=" << stateCount << '\n';
    std::cout << "maximum_gap=" << maximumGap << '\n';
    std::cout << "minimum_return_gap=" << minimumReturnGap << '\n';
    std::cout << "minimum_first_return_segment_gap="
        << minimumFirstSegmentExample.gap << '\n';
    std::cout << "minimum_first_return_segment_weight="
        << minimumFirstSegmentWeight << '\n';
    std::cout << "minimum_first_return_excess_numerator="
        << minimumFirstSegmentExcessNumerator << '\n';
    std::cout << "output=" << outputPath << '\n';
    std::cout << "status=EXACT_G2_LOW_OUTPUT_RETURN_SEARCH\n";
}
