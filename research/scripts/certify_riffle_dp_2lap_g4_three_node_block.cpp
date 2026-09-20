#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned LowOutputMaximum = 8;
    constexpr unsigned RejectedTripleMaximum = 24;
    constexpr std::uint64_t ExpectedLowOutputs = 5'130'659'560ULL;
    constexpr unsigned PrefixLength = 26'832;
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
            throw std::runtime_error("three-node certificate: generator inverse failed");
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
                throw std::runtime_error(
                    "three-node certificate: systematic construction failed");
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

    std::uint64_t applyColumns(
        const std::array<std::uint64_t, 64>& columns,
        std::uint64_t value)
    {
        std::uint64_t result = 0;
        while (value)
        {
            const auto bit = std::countr_zero(value);
            result ^= columns[bit];
            value &= value - 1;
        }
        return result;
    }

    std::array<std::uint64_t, 64> inverseColumns(
        const std::array<std::uint64_t, 64>& columns)
    {
        struct Row
        {
            std::uint64_t left;
            std::uint64_t right;
        };
        std::array<Row, 64> rows{};
        for (unsigned output = 0; output < 64; ++output)
        {
            std::uint64_t left = 0;
            for (unsigned input = 0; input < 64; ++input)
                left |= ((columns[input] >> output) & 1) << input;
            rows[output] = Row{left, std::uint64_t{1} << output};
        }
        for (unsigned column = 0; column < 64; ++column)
        {
            unsigned pivot = column;
            while (pivot < 64 && ((rows[pivot].left >> column) & 1) == 0)
                ++pivot;
            if (pivot == 64)
                throw std::runtime_error("three-node certificate: output map is singular");
            std::swap(rows[column], rows[pivot]);
            for (unsigned row = 0; row < 64; ++row)
            {
                if (row != column && ((rows[row].left >> column) & 1))
                {
                    rows[row].left ^= rows[column].left;
                    rows[row].right ^= rows[column].right;
                }
            }
        }
        std::array<std::uint64_t, 64> inverse{};
        for (unsigned input = 0; input < 64; ++input)
        {
            for (unsigned output = 0; output < 64; ++output)
                inverse[input] |= ((rows[output].right >> input) & 1) << output;
        }
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            if (applyColumns(columns, inverse[bit]) != (std::uint64_t{1} << bit))
                throw std::runtime_error(
                    "three-node certificate: inverse verification failed");
        }
        return inverse;
    }

    std::uint64_t nextCombination(std::uint64_t value)
    {
        const auto low = value & (~value + 1);
        const auto ripple = value + low;
        if (ripple == 0)
            return 0;
        return ripple | (((ripple ^ value) >> 2) / low);
    }

    struct Minimum
    {
        unsigned weight = std::numeric_limits<unsigned>::max();
        std::uint64_t low = 0;
        std::uint64_t left2 = 0;
        std::uint64_t left1 = 0;
        std::uint64_t right1 = 0;
        std::uint64_t right2 = 0;
        const char* alignment = "";
    };

    void updateMinimum(
        Minimum& minimum,
        unsigned weight,
        std::uint64_t low,
        std::uint64_t left2,
        std::uint64_t left1,
        std::uint64_t right1,
        std::uint64_t right2,
        const char* alignment)
    {
        if (weight < minimum.weight)
        {
            minimum = Minimum{
                weight, low, left2, left1, right1, right2, alignment};
        }
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::string outputPath = argc > 1
            ? argv[1]
            : "explorations/riffle_dp_2lap_g4_three_node_block_certificate.json";

        const auto parityColumns = systematicRightColumns();
        std::array<std::uint64_t, 64> forwardColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            forwardColumns[bit] = accumulate(parityColumns[bit]);
        const auto backwardColumns = inverseColumns(forwardColumns);

        std::array<std::uint64_t, 64> forward2Columns{};
        std::array<std::uint64_t, 64> backward2Columns{};
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            forward2Columns[bit] = applyColumns(
                forwardColumns, forwardColumns[bit]);
            backward2Columns[bit] = applyColumns(
                backwardColumns, backwardColumns[bit]);
        }

        std::uint64_t enumerated = 0;
        std::uint64_t rejectedTripleCount = 0;
        Minimum minimum;
        for (unsigned lowWeight = 1; lowWeight <= LowOutputMaximum; ++lowWeight)
        {
            auto low = (std::uint64_t{1} << lowWeight) - 1;
            auto left2 = applyColumns(backward2Columns, low);
            auto left1 = applyColumns(backwardColumns, low);
            auto right1 = applyColumns(forwardColumns, low);
            auto right2 = applyColumns(forward2Columns, low);
            while (low)
            {
                ++enumerated;
                const auto endingWeight = static_cast<unsigned>(
                    std::popcount(left2) + std::popcount(left1) + lowWeight);
                const auto centeredWeight = static_cast<unsigned>(
                    std::popcount(left1) + lowWeight + std::popcount(right1));
                const auto startingWeight = static_cast<unsigned>(
                    lowWeight + std::popcount(right1) + std::popcount(right2));
                rejectedTripleCount += endingWeight <= RejectedTripleMaximum;
                rejectedTripleCount += centeredWeight <= RejectedTripleMaximum;
                rejectedTripleCount += startingWeight <= RejectedTripleMaximum;
                updateMinimum(
                    minimum, endingWeight, low, left2, left1, right1, right2, "ending");
                updateMinimum(
                    minimum, centeredWeight, low, left2, left1, right1, right2, "centered");
                updateMinimum(
                    minimum, startingWeight, low, left2, left1, right1, right2, "starting");

                const auto next = nextCombination(low);
                if (next == 0)
                    break;
                auto difference = low ^ next;
                while (difference)
                {
                    const auto bit = std::countr_zero(difference);
                    left2 ^= backward2Columns[bit];
                    left1 ^= backwardColumns[bit];
                    right1 ^= forwardColumns[bit];
                    right2 ^= forward2Columns[bit];
                    difference &= difference - 1;
                }
                low = next;
            }
        }

        if (enumerated != ExpectedLowOutputs)
            throw std::runtime_error("three-node certificate: enumeration count mismatch");
        if (rejectedTripleCount != 0)
            throw std::runtime_error("three-node certificate: found a rejected triple");
        if (PrefixLength % 3 != 0)
            throw std::runtime_error("three-node certificate: prefix partition mismatch");
        constexpr unsigned TripleLowerBound = RejectedTripleMaximum + 1;
        constexpr unsigned TripleCount = PrefixLength / 3;
        constexpr unsigned PrefixLowerBound = TripleCount * TripleLowerBound;
        static_assert(PrefixLowerBound == 223'600);
        static_assert(PrefixLowerBound > DistanceThreshold);

        std::ofstream output(outputPath, std::ios::trunc);
        if (!output)
            throw std::runtime_error("three-node certificate: output open failed");
        output << "{\n"
            << "  \"schema\": \"riffle-dp-2lap-g4-three-node-block-certificate-v2\",\n"
            << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n"
            << "  \"evidence_label\": \"EXACT\",\n"
            << "  \"autonomous_output_map\": \"U(o)=Acc(P(o))\",\n"
            << "  \"enumerated_output_weight_maximum\": " << LowOutputMaximum << ",\n"
            << "  \"enumerated_nonzero_outputs\": " << enumerated << ",\n"
            << "  \"rejected_three_node_weight_maximum\": "
            << RejectedTripleMaximum << ",\n"
            << "  \"rejected_triples_found\": " << rejectedTripleCount << ",\n"
            << "  \"minimum_weight_among_triples_containing_an_enumerated_output\": "
            << minimum.weight << ",\n"
            << "  \"minimum_witness\": {\n"
            << "    \"alignment\": \"" << minimum.alignment << "\",\n"
            << "    \"low_output_hex\": \"0x" << std::hex << minimum.low << "\",\n"
            << "    \"left2_hex\": \"0x" << minimum.left2 << "\",\n"
            << "    \"left1_hex\": \"0x" << minimum.left1 << "\",\n"
            << "    \"right1_hex\": \"0x" << minimum.right1 << "\",\n"
            << "    \"right2_hex\": \"0x" << minimum.right2 << "\"\n"
            << std::dec
            << "  },\n"
            << "  \"three_node_weight_lower_bound\": " << TripleLowerBound << ",\n"
            << "  \"prefix_length\": " << PrefixLength << ",\n"
            << "  \"three_node_blocks\": " << TripleCount << ",\n"
            << "  \"certified_nonzero_terminal_prefix_weight\": "
            << PrefixLowerBound << ",\n"
            << "  \"distance_threshold\": " << DistanceThreshold << ",\n"
            << "  \"margin_above_distance\": "
            << PrefixLowerBound - DistanceThreshold << ",\n"
            << "  \"proof\": \"A triple of total weight at most 24 has an output of weight at most 8. The exhaustive enumeration checks all three alignments containing every such nonzero output. The autonomous map is invertible, so each block entered from a nonzero terminal state is nonzero.\",\n"
            << "  \"scope_limitation\": \"The certificate proves sufficient wrapped-prefix weight for every nonzero first-lap terminal state. It gives no weight when that terminal state is zero.\"\n"
            << "}\n";
        if (!output)
            throw std::runtime_error("three-node certificate: output write failed");

        std::cout << "candidate=Riffle DP-2Lap g=4\n";
        std::cout << "enumerated_nonzero_outputs=" << enumerated << '\n';
        std::cout << "rejected_triples_found=" << rejectedTripleCount << '\n';
        std::cout << "minimum_enumerated_triple_weight=" << minimum.weight << '\n';
        std::cout << "three_node_weight_lower_bound=" << TripleLowerBound << '\n';
        std::cout << "certified_prefix_weight=" << PrefixLowerBound << '\n';
        std::cout << "margin_above_distance=" << PrefixLowerBound - DistanceThreshold << '\n';
        std::cout << "output=" << outputPath << '\n';
        std::cout << "status=EXACT_DP_2LAP_THREE_NODE_BLOCK_CERTIFICATE\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
