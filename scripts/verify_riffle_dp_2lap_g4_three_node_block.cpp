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
    constexpr unsigned LowOutputMaximum = 8;
    constexpr unsigned RejectedTripleMaximum = 24;
    constexpr std::uint64_t ExpectedLowOutputs = 5'130'659'560ULL;
    constexpr unsigned ExpectedMinimumEnumeratedTriple = 25;
    constexpr unsigned ExpectedPrefixWeight = 223'600;

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
            throw std::runtime_error("three-node verifier: generator inverse failed");
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
                    "three-node verifier: systematic construction failed");
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
                throw std::runtime_error("three-node verifier: output map is singular");
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

        std::uint64_t operator()(std::uint64_t value) const
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
            throw std::runtime_error("three-node verifier: missing receipt key " + key);
        auto cursor = text.find(':', location);
        if (cursor == std::string::npos)
            throw std::runtime_error("three-node verifier: malformed receipt key " + key);
        do
            ++cursor;
        while (cursor < text.size() && std::isspace(static_cast<unsigned char>(text[cursor])));
        std::uint64_t value = 0;
        if (cursor == text.size() || !std::isdigit(static_cast<unsigned char>(text[cursor])))
            throw std::runtime_error("three-node verifier: nonnumeric receipt key " + key);
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
            : "explorations/riffle_dp_2lap_g4_three_node_block_certificate.json";
        std::ifstream receiptInput(receiptPath);
        if (!receiptInput)
            throw std::runtime_error("three-node verifier: receipt open failed");
        const std::string receipt(
            (std::istreambuf_iterator<char>(receiptInput)),
            std::istreambuf_iterator<char>());

        const auto parityColumns = systematicRightColumns();
        std::array<std::uint64_t, 64> forwardColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            forwardColumns[bit] = accumulate(parityColumns[bit]);
        const auto backwardColumns = inverseColumns(forwardColumns);
        const ByteMap forward(forwardColumns);
        const ByteMap backward(backwardColumns);
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            const auto basis = std::uint64_t{1} << bit;
            if (backward(forward(basis)) != basis || forward(backward(basis)) != basis)
                throw std::runtime_error("three-node verifier: inverse map check failed");
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
                const auto right1 = forward(low);
                const auto right2 = forward(right1);
                const auto ending = static_cast<unsigned>(
                    std::popcount(left2) + std::popcount(left1) + lowWeight);
                const auto centered = static_cast<unsigned>(
                    std::popcount(left1) + lowWeight + std::popcount(right1));
                const auto starting = static_cast<unsigned>(
                    lowWeight + std::popcount(right1) + std::popcount(right2));
                rejected += ending <= RejectedTripleMaximum;
                rejected += centered <= RejectedTripleMaximum;
                rejected += starting <= RejectedTripleMaximum;
                minimum = std::min({minimum, ending, centered, starting});
                low = nextCombination(low);
            }
        }

        if (enumerated != ExpectedLowOutputs || rejected != 0 ||
            minimum != ExpectedMinimumEnumeratedTriple)
            throw std::runtime_error("three-node verifier: exhaustive result mismatch");
        if (extractUnsigned(receipt, "enumerated_nonzero_outputs") != enumerated ||
            extractUnsigned(receipt, "rejected_triples_found") != rejected ||
            extractUnsigned(receipt, "minimum_weight_among_triples_containing_an_enumerated_output") != minimum ||
            extractUnsigned(receipt, "certified_nonzero_terminal_prefix_weight") != ExpectedPrefixWeight)
            throw std::runtime_error("three-node verifier: receipt mismatch");

        std::cout << "candidate=Riffle DP-2Lap g=4\n";
        std::cout << "enumerated_nonzero_outputs=" << enumerated << '\n';
        std::cout << "rejected_triples_found=" << rejected << '\n';
        std::cout << "minimum_enumerated_triple_weight=" << minimum << '\n';
        std::cout << "certified_prefix_weight=" << ExpectedPrefixWeight << '\n';
        std::cout << "status=EXACT_DP_2LAP_THREE_NODE_BLOCK_INDEPENDENTLY_VERIFIED\n";
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
