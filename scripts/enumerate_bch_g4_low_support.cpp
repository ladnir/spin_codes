#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned CoordinateCount = 128;
    constexpr unsigned Dimension = 64;
    constexpr unsigned PacketBits = 4;
    constexpr unsigned PacketCount = CoordinateCount / PacketBits;

    struct Word128
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
    };

    std::array<std::uint64_t, CoordinateCount> parityCheckColumns()
    {
        std::array<Word128, Dimension> rows{};
        for (unsigned row = 0; row < Dimension; ++row)
        {
            rows[row].low = Generator << row;
            rows[row].high = row == 0 ? 0 : Generator >> (Dimension - row);
            rows[row].high |= std::uint64_t{1} << 63;
        }

        for (unsigned column = 0; column < Dimension; ++column)
        {
            unsigned pivot = column;
            while (pivot < Dimension && ((rows[pivot].low >> column) & 1) == 0)
                ++pivot;
            if (pivot == Dimension)
                throw std::runtime_error("singular systematic half");
            std::swap(rows[column], rows[pivot]);
            for (unsigned row = 0; row < Dimension; ++row)
            {
                if (row != column && ((rows[row].low >> column) & 1))
                {
                    rows[row].low ^= rows[column].low;
                    rows[row].high ^= rows[column].high;
                }
            }
        }

        std::array<std::uint64_t, CoordinateCount> columns{};
        for (unsigned coordinate = 0; coordinate < Dimension; ++coordinate)
        {
            columns[coordinate] = rows[coordinate].high;
            columns[Dimension + coordinate] = std::uint64_t{1} << coordinate;
        }

        for (unsigned row = 0; row < Dimension; ++row)
        {
            std::uint64_t syndrome = 0;
            const auto encodedLow = Generator << row;
            const auto encodedHigh =
                (row == 0 ? 0 : Generator >> (Dimension - row)) |
                (std::uint64_t{1} << 63);
            for (unsigned coordinate = 0; coordinate < Dimension; ++coordinate)
            {
                if ((encodedLow >> coordinate) & 1)
                    syndrome ^= columns[coordinate];
                if ((encodedHigh >> coordinate) & 1)
                    syndrome ^= columns[Dimension + coordinate];
            }
            if (syndrome != 0)
                throw std::runtime_error("parity-check construction failed");
        }
        return columns;
    }

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
        for (unsigned bit = 1; bit < Dimension; ++bit)
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
                (bit == 0 ? 0 : Generator >> (Dimension - bit)) |
                (std::uint64_t{1} << 63);
            message &= message - 1;
        }
        return result;
    }

    unsigned nibbleSupport(const Word128& word)
    {
        unsigned support = 0;
        for (unsigned packet = 0; packet < 16; ++packet)
        {
            support += ((word.low >> (PacketBits * packet)) & 0xF) != 0;
            support += ((word.high >> (PacketBits * packet)) & 0xF) != 0;
        }
        return support;
    }

    class Enumerator
    {
    public:
        Enumerator(unsigned maximumSupport, const std::string& outputPath)
            : mMaximumSupport(maximumSupport),
              mColumns(parityCheckColumns()),
              mGeneratorInverse(generatorInverse())
        {
            if (maximumSupport > 16)
                throw std::invalid_argument("maximum support must not exceed 16");
            if (!outputPath.empty())
            {
                mOutput.open(outputPath, std::ios::binary | std::ios::trunc);
                if (!mOutput)
                    throw std::runtime_error("failed to open output file");
            }
        }

        void run()
        {
            visit(0, 0);
            if (mOutput)
                mOutput.flush();
        }

        const std::array<std::uint64_t, 17>& counts() const { return mCounts; }
        std::uint64_t visitedSubsets() const { return mVisitedSubsets; }
        unsigned maximumKernelDimension() const { return mMaximumKernelDimension; }

    private:
        unsigned mMaximumSupport;
        std::array<std::uint64_t, CoordinateCount> mColumns;
        std::uint64_t mGeneratorInverse;
        std::array<std::uint64_t, Dimension> mPivotValues{};
        std::array<std::uint64_t, Dimension> mPivotRepresentations{};
        std::array<std::uint64_t, Dimension> mKernelBasis{};
        std::array<unsigned, PacketCount> mSelectedPackets{};
        std::array<std::uint64_t, 17> mCounts{};
        std::uint64_t mVisitedSubsets = 0;
        unsigned mKernelDimension = 0;
        unsigned mMaximumKernelDimension = 0;
        std::ofstream mOutput;

        bool insertColumn(
            std::uint64_t value,
            std::uint64_t representation,
            std::array<unsigned, PacketBits>& addedPivots,
            unsigned& addedPivotCount)
        {
            while (value)
            {
                const auto pivot = std::bit_width(value) - 1;
                if (mPivotValues[pivot] == 0)
                {
                    mPivotValues[pivot] = value;
                    mPivotRepresentations[pivot] = representation;
                    addedPivots[addedPivotCount++] = pivot;
                    return true;
                }
                value ^= mPivotValues[pivot];
                representation ^= mPivotRepresentations[pivot];
            }
            mKernelBasis[mKernelDimension++] = representation;
            return false;
        }

        void evaluate(unsigned support)
        {
            if (mKernelDimension == 0)
                return;
            mMaximumKernelDimension = std::max(
                mMaximumKernelDimension, mKernelDimension);
            const auto selectorLimit = std::uint64_t{1} << mKernelDimension;
            for (std::uint64_t selector = 1; selector < selectorLimit; ++selector)
            {
                std::uint64_t localCoordinates = 0;
                for (unsigned basis = 0; basis < mKernelDimension; ++basis)
                {
                    if ((selector >> basis) & 1)
                        localCoordinates ^= mKernelBasis[basis];
                }

                Word128 codeword{};
                for (unsigned local = 0; local < PacketBits * support; ++local)
                {
                    if (((localCoordinates >> local) & 1) == 0)
                        continue;
                    const auto coordinate =
                        PacketBits * mSelectedPackets[local / PacketBits] +
                        local % PacketBits;
                    if (coordinate < Dimension)
                        codeword.low |= std::uint64_t{1} << coordinate;
                    else
                        codeword.high |= std::uint64_t{1} << (coordinate - Dimension);
                }
                if (nibbleSupport(codeword) != support)
                    continue;

                const auto message =
                    carrylessMultiplyLow(codeword.low, mGeneratorInverse);
                const auto check = encode(message);
                if (check.low != codeword.low || check.high != codeword.high)
                    throw std::runtime_error("enumerated codeword failed encoding check");
                ++mCounts[support];
                if (mOutput)
                {
                    const auto supportByte = static_cast<std::uint8_t>(support);
                    mOutput.write(
                        reinterpret_cast<const char*>(&supportByte),
                        sizeof(supportByte));
                    mOutput.write(
                        reinterpret_cast<const char*>(&message), sizeof(message));
                }
            }
        }

        void visit(unsigned nextPacket, unsigned depth)
        {
            if (depth == mMaximumSupport)
                return;
            for (unsigned packet = nextPacket; packet < PacketCount; ++packet)
            {
                ++mVisitedSubsets;
                mSelectedPackets[depth] = packet;
                const auto oldKernelDimension = mKernelDimension;
                std::array<unsigned, PacketBits> addedPivots{};
                unsigned addedPivotCount = 0;
                for (unsigned bit = 0; bit < PacketBits; ++bit)
                {
                    const auto local = PacketBits * depth + bit;
                    insertColumn(
                        mColumns[PacketBits * packet + bit],
                        std::uint64_t{1} << local,
                        addedPivots,
                        addedPivotCount);
                }

                const auto support = depth + 1;
                evaluate(support);
                visit(packet + 1, support);

                mKernelDimension = oldKernelDimension;
                while (addedPivotCount)
                {
                    const auto pivot = addedPivots[--addedPivotCount];
                    mPivotValues[pivot] = 0;
                    mPivotRepresentations[pivot] = 0;
                }
            }
        }
    };
}

int main(int argc, char** argv)
{
    try
    {
        unsigned maximumSupport = 11;
        std::string outputPath;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            if (argument == "--max-support" && index + 1 < argc)
                maximumSupport = static_cast<unsigned>(std::stoul(argv[++index]));
            else if (argument == "--output" && index + 1 < argc)
                outputPath = argv[++index];
            else
                throw std::invalid_argument("unknown or incomplete argument: " + argument);
        }

        Enumerator enumerator(maximumSupport, outputPath);
        enumerator.run();
        std::cout << "exact BCH g=4 low-support enumeration\n";
        std::cout << "maximum_support=" << maximumSupport << '\n';
        std::cout << "visited_packet_subsets=" << enumerator.visitedSubsets() << '\n';
        std::cout << "maximum_shortened_dimension="
                  << enumerator.maximumKernelDimension() << '\n';
        for (unsigned support = 1; support <= maximumSupport; ++support)
            std::cout << "support_" << support << "_codewords="
                      << enumerator.counts()[support] << '\n';
        if (!outputPath.empty())
            std::cout << "output=" << outputPath << '\n';
        std::cout << "status=EXACT_VERIFIED\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
