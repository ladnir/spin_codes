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
    constexpr unsigned LowMaximum = 8;
    constexpr std::uint64_t ExpectedStates = 5'130'659'560ULL;

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
                throw std::runtime_error("output map is singular");
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
                throw std::runtime_error("output-map inverse check failed");
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
            throw std::runtime_error("missing receipt key: " + key);
        auto cursor = text.find(':', location);
        if (cursor == std::string::npos)
            throw std::runtime_error("malformed receipt key: " + key);
        do
            ++cursor;
        while (cursor < text.size() && std::isspace(static_cast<unsigned char>(text[cursor])));
        std::uint64_t value = 0;
        if (cursor == text.size() || !std::isdigit(static_cast<unsigned char>(text[cursor])))
            throw std::runtime_error("receipt value is not unsigned: " + key);
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
    const std::string receiptPath = argc > 1
        ? argv[1]
        : "explorations/riffle_dp_g4_g2_low_output_returns_w8_gap2.json";
    std::ifstream receiptInput(receiptPath);
    if (!receiptInput)
        throw std::runtime_error("failed to open forward receipt");
    const std::string receipt(
        (std::istreambuf_iterator<char>(receiptInput)),
        std::istreambuf_iterator<char>());

    const auto parityColumns = systematicRightColumns();
    std::array<std::uint64_t, 64> forwardColumns{};
    for (unsigned bit = 0; bit < 64; ++bit)
        forwardColumns[bit] = accumulate(parityColumns[bit]);
    const ByteMap inverse(inverseColumns(forwardColumns));

    std::uint64_t stateCount = 0;
    std::uint64_t gapOneReturns = 0;
    std::uint64_t gapTwoReturns = 0;
    unsigned minimumGapTwoSegmentWeight = std::numeric_limits<unsigned>::max();
    for (unsigned weight = 1; weight <= LowMaximum; ++weight)
    {
        auto finish = (std::uint64_t{1} << weight) - 1;
        while (finish)
        {
            ++stateCount;
            const auto previous = inverse(finish);
            if (std::popcount(previous) <= LowMaximum)
                ++gapOneReturns;
            const auto start = inverse(previous);
            if (std::popcount(start) <= LowMaximum)
            {
                ++gapTwoReturns;
                minimumGapTwoSegmentWeight = std::min(
                    minimumGapTwoSegmentWeight,
                    static_cast<unsigned>(
                        std::popcount(start) + std::popcount(previous)));
            }
            finish = nextCombination(finish);
        }
    }
    if (stateCount != ExpectedStates)
        throw std::runtime_error("reverse verifier state count mismatch");
    if (gapOneReturns != 0 || gapTwoReturns != 11'039'336 ||
        minimumGapTwoSegmentWeight != 19)
    {
        throw std::runtime_error("reverse verifier return statistics changed");
    }
    if (extractUnsigned(receipt, "enumerated_low_states") != stateCount ||
        extractUnsigned(receipt, "minimum_return_gap") != 2)
    {
        throw std::runtime_error("forward receipt scalar mismatch");
    }
    std::cout << "candidate=Riffle DP g=4@g0-v1\n";
    std::cout << "enumerated_low_endpoints=" << stateCount << '\n';
    std::cout << "reverse_gap_one_returns=" << gapOneReturns << '\n';
    std::cout << "reverse_gap_two_returns=" << gapTwoReturns << '\n';
    std::cout << "reverse_minimum_gap_two_segment_weight="
        << minimumGapTwoSegmentWeight << '\n';
    std::cout << "status=EXACT_G2_LOW_OUTPUT_RETURNS_REVERSE_VERIFIED\n";
}
