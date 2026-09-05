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
    constexpr unsigned RejectedWindowMaximum = 35;
    constexpr std::uint64_t ExpectedLowOutputs = 5'130'659'560ULL;
    constexpr unsigned BoundaryPrefixNodes = 20'976;
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
            throw std::runtime_error("four-node certificate: generator inverse failed");
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
                    "four-node certificate: systematic construction failed");
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

    std::array<std::uint64_t, 64> composeColumns(
        const std::array<std::uint64_t, 64>& outer,
        const std::array<std::uint64_t, 64>& inner)
    {
        std::array<std::uint64_t, 64> result{};
        for (unsigned bit = 0; bit < 64; ++bit)
            result[bit] = applyColumns(outer, inner[bit]);
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
                throw std::runtime_error("four-node certificate: autonomous map singular");
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
                throw std::runtime_error("four-node certificate: inverse check failed");
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
        std::uint64_t left3 = 0;
        std::uint64_t left2 = 0;
        std::uint64_t left1 = 0;
        std::uint64_t right1 = 0;
        std::uint64_t right2 = 0;
        std::uint64_t right3 = 0;
        unsigned alignment = 0;
    };

    void updateMinimum(
        Minimum& minimum,
        unsigned weight,
        std::uint64_t low,
        std::uint64_t left3,
        std::uint64_t left2,
        std::uint64_t left1,
        std::uint64_t right1,
        std::uint64_t right2,
        std::uint64_t right3,
        unsigned alignment)
    {
        if (weight < minimum.weight)
        {
            minimum = Minimum{
                weight,
                low,
                left3,
                left2,
                left1,
                right1,
                right2,
                right3,
                alignment};
        }
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::string outputPath = argc > 1
            ? argv[1]
            : "constructions/riffle_packetmul_2lap_g4/receipts/goal04_four_node_primary.json";

        const auto parityColumns = systematicRightColumns();
        std::array<std::uint64_t, 64> forward1{};
        for (unsigned bit = 0; bit < 64; ++bit)
            forward1[bit] = accumulate(parityColumns[bit]);
        const auto backward1 = inverseColumns(forward1);
        const auto forward2 = composeColumns(forward1, forward1);
        const auto forward3 = composeColumns(forward1, forward2);
        const auto backward2 = composeColumns(backward1, backward1);
        const auto backward3 = composeColumns(backward1, backward2);

        std::uint64_t enumerated = 0;
        std::uint64_t rejected = 0;
        Minimum minimum;
        for (unsigned lowWeight = 1; lowWeight <= LowOutputMaximum; ++lowWeight)
        {
            auto low = (std::uint64_t{1} << lowWeight) - 1;
            auto left3 = applyColumns(backward3, low);
            auto left2 = applyColumns(backward2, low);
            auto left1 = applyColumns(backward1, low);
            auto right1 = applyColumns(forward1, low);
            auto right2 = applyColumns(forward2, low);
            auto right3 = applyColumns(forward3, low);
            while (low)
            {
                ++enumerated;
                const auto w0 = static_cast<unsigned>(
                    std::popcount(left3) + std::popcount(left2) +
                    std::popcount(left1) + lowWeight);
                const auto w1 = static_cast<unsigned>(
                    std::popcount(left2) + std::popcount(left1) +
                    lowWeight + std::popcount(right1));
                const auto w2 = static_cast<unsigned>(
                    std::popcount(left1) + lowWeight +
                    std::popcount(right1) + std::popcount(right2));
                const auto w3 = static_cast<unsigned>(
                    lowWeight + std::popcount(right1) +
                    std::popcount(right2) + std::popcount(right3));
                rejected += w0 <= RejectedWindowMaximum;
                rejected += w1 <= RejectedWindowMaximum;
                rejected += w2 <= RejectedWindowMaximum;
                rejected += w3 <= RejectedWindowMaximum;
                updateMinimum(minimum, w0, low, left3, left2, left1, right1, right2, right3, 0);
                updateMinimum(minimum, w1, low, left3, left2, left1, right1, right2, right3, 1);
                updateMinimum(minimum, w2, low, left3, left2, left1, right1, right2, right3, 2);
                updateMinimum(minimum, w3, low, left3, left2, left1, right1, right2, right3, 3);

                const auto next = nextCombination(low);
                if (next == 0)
                    break;
                auto difference = low ^ next;
                while (difference)
                {
                    const auto bit = std::countr_zero(difference);
                    left3 ^= backward3[bit];
                    left2 ^= backward2[bit];
                    left1 ^= backward1[bit];
                    right1 ^= forward1[bit];
                    right2 ^= forward2[bit];
                    right3 ^= forward3[bit];
                    difference &= difference - 1;
                }
                low = next;
            }
        }

        if (enumerated != ExpectedLowOutputs)
            throw std::runtime_error("four-node certificate: enumeration count mismatch");

        constexpr unsigned ProposedLowerBound = RejectedWindowMaximum + 1;
        constexpr unsigned CompleteBlocks = BoundaryPrefixNodes / 4;
        constexpr unsigned PrefixLowerBound = CompleteBlocks * ProposedLowerBound;
        static_assert(BoundaryPrefixNodes % 4 == 0);
        static_assert(PrefixLowerBound == 188'784);
        static_assert(PrefixLowerBound > DistanceThreshold);
        const bool proved = rejected == 0;

        std::ofstream output(outputPath, std::ios::trunc);
        if (!output)
            throw std::runtime_error("four-node certificate: output open failed");
        output << "{\n"
            << "  \"schema\": \"riffle-packetmul-2lap-g4-goal04-four-node-primary-v1\",\n"
            << "  \"candidate\": \"Riffle PacketMul-2Lap g=4\",\n"
            << "  \"evidence_label\": \"" << (proved ? "EXACT" : "EXACT_REFUTATION") << "\",\n"
            << "  \"autonomous_output_map\": \"U(o)=Acc(P(o))\",\n"
            << "  \"enumerated_output_weight_maximum\": " << LowOutputMaximum << ",\n"
            << "  \"enumerated_nonzero_outputs\": " << enumerated << ",\n"
            << "  \"rejected_four_node_weight_maximum\": " << RejectedWindowMaximum << ",\n"
            << "  \"rejected_aligned_windows_found\": " << rejected << ",\n"
            << "  \"minimum_weight_among_windows_containing_an_enumerated_output\": "
            << minimum.weight << ",\n"
            << "  \"minimum_witness\": {\n"
            << "    \"anchor_alignment\": " << minimum.alignment << ",\n"
            << "    \"anchor_hex\": \"0x" << std::hex << minimum.low << "\",\n"
            << "    \"left3_hex\": \"0x" << minimum.left3 << "\",\n"
            << "    \"left2_hex\": \"0x" << minimum.left2 << "\",\n"
            << "    \"left1_hex\": \"0x" << minimum.left1 << "\",\n"
            << "    \"right1_hex\": \"0x" << minimum.right1 << "\",\n"
            << "    \"right2_hex\": \"0x" << minimum.right2 << "\",\n"
            << "    \"right3_hex\": \"0x" << minimum.right3 << "\"\n"
            << std::dec
            << "  },\n"
            << "  \"four_node_weight_lower_bound\": " << (proved ? ProposedLowerBound : 0) << ",\n"
            << "  \"boundary_prefix_nodes\": " << BoundaryPrefixNodes << ",\n"
            << "  \"complete_four_node_blocks\": " << CompleteBlocks << ",\n"
            << "  \"certified_prefix_weight\": " << (proved ? PrefixLowerBound : 0) << ",\n"
            << "  \"distance_threshold\": " << DistanceThreshold << ",\n"
            << "  \"target_first_node_strata\": [22650, 22651, 22652],\n"
            << "  \"target_strata_closed\": " << (proved ? "true" : "false") << ",\n"
            << "  \"proof\": \"A four-node window of total weight at most 35 contains an output of weight at most 8. The exhaustive search checks all four alignments containing every such nonzero output.\",\n"
            << "  \"scope_limitation\": \"The certificate concerns only a zero-input autonomous prefix entered from a nonzero state. It does not control terminal-zero events.\"\n"
            << "}\n";
        if (!output)
            throw std::runtime_error("four-node certificate: output write failed");

        std::cout << "candidate=Riffle PacketMul-2Lap g=4\n";
        std::cout << "enumerated_nonzero_outputs=" << enumerated << '\n';
        std::cout << "rejected_aligned_windows_found=" << rejected << '\n';
        std::cout << "minimum_enumerated_four_node_weight=" << minimum.weight << '\n';
        std::cout << "target_strata_closed=" << (proved ? "true" : "false") << '\n';
        std::cout << "output=" << outputPath << '\n';
        std::cout << "status=" << (proved
            ? "EXACT_PACKETMUL_2LAP_FOUR_NODE_CERTIFICATE"
            : "EXACT_PACKETMUL_2LAP_FOUR_NODE_REFUTATION") << '\n';
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
