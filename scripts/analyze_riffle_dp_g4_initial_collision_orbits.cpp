#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned DistanceThreshold = 188'766;
    constexpr unsigned InnerNodes = 32'772;

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

    class AutonomousMap
    {
    public:
        AutonomousMap()
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

        std::uint64_t operator()(std::uint64_t state) const
        {
            const auto output = accumulate(state);
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

    template<typename Callback>
    void enumerateDrives(unsigned support, Callback&& callback)
    {
        if (support == 2)
        {
            for (unsigned a = 0; a < 16; ++a)
            for (unsigned b = a + 1; b < 16; ++b)
            for (unsigned x = 1; x < 16; ++x)
            for (unsigned y = 1; y < 16; ++y)
                callback((std::uint64_t{x} << (4 * a)) | (std::uint64_t{y} << (4 * b)));
        }
        else if (support == 3)
        {
            for (unsigned a = 0; a < 16; ++a)
            for (unsigned b = a + 1; b < 16; ++b)
            for (unsigned c = b + 1; c < 16; ++c)
            for (unsigned x = 1; x < 16; ++x)
            for (unsigned y = 1; y < 16; ++y)
            for (unsigned z = 1; z < 16; ++z)
                callback(
                    (std::uint64_t{x} << (4 * a)) |
                    (std::uint64_t{y} << (4 * b)) |
                    (std::uint64_t{z} << (4 * c)));
        }
        else
            throw std::invalid_argument("support must be two or three");
    }
}

int main(int argc, char** argv)
{
    unsigned support = 2;
    std::string outputPath =
        "explorations/riffle_dp_g4_g2_initial_collision_orbits_s2.json";
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        if (argument == "--support" && index + 1 < argc)
            support = static_cast<unsigned>(std::stoul(argv[++index]));
        else if (argument == "--output" && index + 1 < argc)
            outputPath = argv[++index];
        else
            throw std::invalid_argument("unknown or incomplete argument: " + argument);
    }
    if (support != 2 && support != 3)
        throw std::invalid_argument("support must be two or three");

    const AutonomousMap step;
    std::vector<unsigned> minimumPrefix(InnerNodes + 1, std::numeric_limits<unsigned>::max());
    std::vector<std::uint64_t> crossingCounts(InnerNodes + 2);
    std::uint64_t driveCount = 0;
    unsigned minimumCrossing = InnerNodes + 1;
    unsigned maximumCrossing = 0;
    std::uint64_t worstDrive = 0;
    enumerateDrives(support, [&](std::uint64_t drive)
    {
        ++driveCount;
        auto state = drive;
        unsigned cumulative = 0;
        unsigned crossing = InnerNodes + 1;
        for (unsigned length = 1; length <= InnerNodes; ++length)
        {
            cumulative += std::popcount(accumulate(state));
            minimumPrefix[length] = std::min(minimumPrefix[length], cumulative);
            if (cumulative > DistanceThreshold)
            {
                crossing = length;
                break;
            }
            state = step(state);
        }
        ++crossingCounts[crossing];
        minimumCrossing = std::min(minimumCrossing, crossing);
        if (crossing > maximumCrossing)
        {
            maximumCrossing = crossing;
            worstDrive = drive;
        }
    });
    if (maximumCrossing > InnerNodes)
        throw std::runtime_error("an initial collision orbit did not cross distance");
    for (unsigned length = 1; length < maximumCrossing; ++length)
    {
        if (minimumPrefix[length] == std::numeric_limits<unsigned>::max())
            throw std::runtime_error("minimum prefix table has a gap");
    }

    std::ofstream output(outputPath, std::ios::trunc);
    if (!output)
        throw std::runtime_error("failed to open output path");
    output << "{\n"
        << "  \"schema\": \"riffle-dp-g4-g2-initial-collision-orbits-v1\",\n"
        << "  \"candidate\": \"Riffle DP g=4@g0-v1\",\n"
        << "  \"evidence_label\": \"EXACT\",\n"
        << "  \"initial_packet_support\": " << support << ",\n"
        << "  \"drive_count\": " << driveCount << ",\n"
        << "  \"first_step_above_distance_min\": " << minimumCrossing << ",\n"
        << "  \"first_step_above_distance_max\": " << maximumCrossing << ",\n"
        << "  \"worst_drive_hex\": \"0x" << std::hex << worstDrive << std::dec << "\",\n"
        << "  \"crossing_counts\": {";
    bool first = true;
    for (unsigned crossing = minimumCrossing; crossing <= maximumCrossing; ++crossing)
    {
        if (!crossingCounts[crossing])
            continue;
        output << (first ? "" : ",") << "\"" << crossing << "\":"
            << crossingCounts[crossing];
        first = false;
    }
    output << "},\n  \"minimum_prefix_weights_through_last_safe_length\": [";
    for (unsigned length = 1; length < maximumCrossing; ++length)
        output << (length == 1 ? "" : ",") << minimumPrefix[length];
    output << "]\n}\n";
    output.close();

    std::cout << "candidate=Riffle DP g=4@g0-v1\n";
    std::cout << "initial_packet_support=" << support << '\n';
    std::cout << "drive_count=" << driveCount << '\n';
    std::cout << "crossing_step_min=" << minimumCrossing << '\n';
    std::cout << "crossing_step_max=" << maximumCrossing << '\n';
    std::cout << "worst_drive_hex=0x" << std::hex << worstDrive << std::dec << '\n';
    std::cout << "output=" << outputPath << '\n';
    std::cout << "status=EXACT_G2_INITIAL_COLLISION_ORBITS\n";
}
