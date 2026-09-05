#include <algorithm>
#include <array>
#include <bit>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr std::uint64_t Mask64 = ~std::uint64_t{0};
    constexpr unsigned LowOutputMaximum = 8;
    constexpr unsigned RejectedWindowMaximum = 35;
    constexpr std::uint64_t ExpectedLowOutputs = 5'130'659'560ULL;
    constexpr unsigned ExpectedMinimumCoveredWindow = 42;
    constexpr unsigned ExpectedPrefixNodes = 20'976;
    constexpr unsigned ExpectedPrefixWeight = 188'784;

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
                throw std::runtime_error("four-node verifier: singular matrix");
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
        return inverse;
    }

    std::array<std::uint64_t, 64> systematicRightColumns()
    {
        std::array<std::uint64_t, 64> left{};
        std::array<std::uint64_t, 64> right{};
        for (unsigned messageBit = 0; messageBit < 64; ++messageBit)
        {
            left[messageBit] = Generator << messageBit;
            right[messageBit] =
                (messageBit == 0 ? 0 : Generator >> (64 - messageBit)) |
                (std::uint64_t{1} << 63);
        }
        const auto inverseLeft = inverseColumns(left);
        std::array<std::uint64_t, 64> result{};
        for (unsigned column = 0; column < 64; ++column)
        {
            result[column] = applyColumns(right, inverseLeft[column]);
            if (applyColumns(left, inverseLeft[column]) != (std::uint64_t{1} << column))
                throw std::runtime_error(
                    "four-node verifier: systematic reconstruction failed");
        }
        return result;
    }

    std::uint64_t accumulateByBits(std::uint64_t value)
    {
        std::uint64_t running = 0;
        std::uint64_t result = 0;
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            running ^= (value >> bit) & 1;
            result |= running << bit;
        }
        return result & Mask64;
    }

    class ByteMap
    {
    public:
        explicit ByteMap(const std::array<std::uint64_t, 64>& columns)
        {
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

        __forceinline std::uint64_t operator()(std::uint64_t value) const
        {
            return
                mTables[0][value & 0xff] ^
                mTables[1][(value >> 8) & 0xff] ^
                mTables[2][(value >> 16) & 0xff] ^
                mTables[3][(value >> 24) & 0xff] ^
                mTables[4][(value >> 32) & 0xff] ^
                mTables[5][(value >> 40) & 0xff] ^
                mTables[6][(value >> 48) & 0xff] ^
                mTables[7][value >> 56];
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

    std::uint64_t extractUnsigned(const std::string& text, const std::string& key)
    {
        const auto location = text.find("\"" + key + "\"");
        if (location == std::string::npos)
            throw std::runtime_error("four-node verifier: missing key " + key);
        auto cursor = text.find(':', location);
        if (cursor == std::string::npos)
            throw std::runtime_error("four-node verifier: malformed key " + key);
        do
            ++cursor;
        while (cursor < text.size() && std::isspace(static_cast<unsigned char>(text[cursor])));
        std::uint64_t value = 0;
        if (cursor == text.size() || !std::isdigit(static_cast<unsigned char>(text[cursor])))
            throw std::runtime_error("four-node verifier: nonnumeric key " + key);
        while (cursor < text.size() && std::isdigit(static_cast<unsigned char>(text[cursor])))
        {
            value = 10 * value + static_cast<unsigned>(text[cursor] - '0');
            ++cursor;
        }
        return value;
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::string receiptPath = argc > 1
            ? argv[1]
            : "constructions/riffle_packetmul_2lap_g4/receipts/goal04_four_node_primary.json";
        const std::string outputPath = argc > 2
            ? argv[2]
            : "constructions/riffle_packetmul_2lap_g4/receipts/goal04_four_node_independent.json";
        std::ifstream receiptInput(receiptPath);
        if (!receiptInput)
            throw std::runtime_error("four-node verifier: receipt open failed");
        const std::string receipt(
            (std::istreambuf_iterator<char>(receiptInput)),
            std::istreambuf_iterator<char>());

        const auto parity = systematicRightColumns();
        std::array<std::uint64_t, 64> forwardColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            forwardColumns[bit] = accumulateByBits(parity[bit]);
        const auto backwardColumns = inverseColumns(forwardColumns);
        const ByteMap forward(forwardColumns);
        const ByteMap backward(backwardColumns);
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            const auto basis = std::uint64_t{1} << bit;
            if (backward(forward(basis)) != basis || forward(backward(basis)) != basis)
                throw std::runtime_error("four-node verifier: inverse map check failed");
        }

        std::uint64_t enumerated = 0;
        std::uint64_t rejected = 0;
        unsigned minimum = std::numeric_limits<unsigned>::max();
        for (unsigned lowWeight = 1; lowWeight <= LowOutputMaximum; ++lowWeight)
        {
            auto low = (std::uint64_t{1} << lowWeight) - 1;
            while (low)
            {
                ++enumerated;
                const auto left1 = backward(low);
                const auto left2 = backward(left1);
                const auto left3 = backward(left2);
                const auto right1 = forward(low);
                const auto right2 = forward(right1);
                const auto right3 = forward(right2);
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
                minimum = std::min({minimum, w0, w1, w2, w3});
                low = nextCombination(low);
            }
        }

        if (enumerated != ExpectedLowOutputs || rejected != 0 ||
            minimum != ExpectedMinimumCoveredWindow)
            throw std::runtime_error("four-node verifier: exhaustive result mismatch");
        if (extractUnsigned(receipt, "enumerated_nonzero_outputs") != enumerated ||
            extractUnsigned(receipt, "rejected_aligned_windows_found") != rejected ||
            extractUnsigned(receipt, "minimum_weight_among_windows_containing_an_enumerated_output") != minimum ||
            extractUnsigned(receipt, "four_node_weight_lower_bound") != 36 ||
            extractUnsigned(receipt, "boundary_prefix_nodes") != ExpectedPrefixNodes ||
            extractUnsigned(receipt, "certified_prefix_weight") != ExpectedPrefixWeight)
            throw std::runtime_error("four-node verifier: primary receipt mismatch");

        std::ofstream output(outputPath, std::ios::trunc);
        if (!output)
            throw std::runtime_error("four-node verifier: output open failed");
        output << "{\n"
            << "  \"schema\": \"riffle-packetmul-2lap-g4-goal04-four-node-independent-v1\",\n"
            << "  \"candidate\": \"Riffle PacketMul-2Lap g=4\",\n"
            << "  \"evidence_label\": \"EXACT_INDEPENDENT_VERIFICATION\",\n"
            << "  \"enumerated_nonzero_outputs\": " << enumerated << ",\n"
            << "  \"rejected_aligned_windows_found\": " << rejected << ",\n"
            << "  \"minimum_covered_four_node_weight\": " << minimum << ",\n"
            << "  \"verified_four_node_weight_lower_bound\": 36,\n"
            << "  \"verified_boundary_prefix_nodes\": " << ExpectedPrefixNodes << ",\n"
            << "  \"verified_complete_four_node_blocks\": 5244,\n"
            << "  \"verified_prefix_weight\": " << ExpectedPrefixWeight << ",\n"
            << "  \"target_strata_closed\": true,\n"
            << "  \"method\": \"Reconstruct the systematic matrix by binary Gaussian elimination, evaluate accumulation bit by bit, and recompute each neighboring output through independent byte tables.\"\n"
            << "}\n";
        if (!output)
            throw std::runtime_error("four-node verifier: output write failed");

        std::cout << "candidate=Riffle PacketMul-2Lap g=4\n";
        std::cout << "enumerated_nonzero_outputs=" << enumerated << '\n';
        std::cout << "rejected_aligned_windows_found=" << rejected << '\n';
        std::cout << "minimum_covered_four_node_weight=" << minimum << '\n';
        std::cout << "verified_prefix_weight=" << ExpectedPrefixWeight << '\n';
        std::cout << "output=" << outputPath << '\n';
        std::cout << "status=EXACT_PACKETMUL_2LAP_FOUR_NODE_INDEPENDENTLY_VERIFIED\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
