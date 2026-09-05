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

    struct Witness
    {
        unsigned firstValue = 0;
        unsigned secondValue = 0;
        unsigned initialFirstCount = 0;
        std::uint64_t initialDrive = 0;
        std::uint64_t resetDrive = 0;
        unsigned resetNode = 0;
        unsigned firstWeight = 0;
        unsigned secondWeight = 0;
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

        std::uint64_t operator()(std::uint64_t value) const
        {
            const auto output = accumulate(value);
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

    std::vector<std::uint64_t> threePacketDrives(
        unsigned firstValue,
        unsigned secondValue,
        unsigned firstCount)
    {
        std::vector<std::uint64_t> result;
        for (unsigned firstSlot = 0; firstSlot < 14; ++firstSlot)
        {
            for (unsigned secondSlot = firstSlot + 1; secondSlot < 15; ++secondSlot)
            {
                for (unsigned thirdSlot = secondSlot + 1; thirdSlot < 16; ++thirdSlot)
                {
                    const std::array<unsigned, 3> slots{firstSlot, secondSlot, thirdSlot};
                    for (unsigned mask = 0; mask < 8; ++mask)
                    {
                        if (static_cast<unsigned>(std::popcount(mask)) != firstCount)
                            continue;
                        std::uint64_t drive = 0;
                        for (unsigned index = 0; index < 3; ++index)
                        {
                            const auto value = ((mask >> index) & 1)
                                ? firstValue
                                : secondValue;
                            drive |= std::uint64_t{value} << (4 * slots[index]);
                        }
                        result.push_back(drive);
                    }
                }
            }
        }
        return result;
    }

    bool hasValueCounts(
        std::uint64_t drive,
        unsigned firstValue,
        unsigned secondValue,
        unsigned firstCount)
    {
        unsigned observedFirst = 0;
        unsigned observedSecond = 0;
        for (unsigned slot = 0; slot < 16; ++slot)
        {
            const auto value = static_cast<unsigned>((drive >> (4 * slot)) & 0xf);
            observedFirst += value == firstValue;
            observedSecond += value == secondValue;
            if (value != 0 && value != firstValue && value != secondValue)
                return false;
        }
        return observedFirst == firstCount && observedSecond == 3 - firstCount;
    }

    void writeNibbles(std::ostream& output, std::uint64_t drive)
    {
        output << '[';
        bool first = true;
        for (unsigned slot = 0; slot < 16; ++slot)
        {
            const auto value = static_cast<unsigned>((drive >> (4 * slot)) & 0xf);
            if (!value)
                continue;
            output << (first ? "" : ",")
                << "{\"slot\":" << slot << ",\"value\":" << value << '}';
            first = false;
        }
        output << ']';
    }
}

int main(int argc, char** argv)
{
    const std::string outputPath = argc > 1
        ? argv[1]
        : "explorations/riffle_dp_2lap_g4_six_packet_3x3_turnoffs.json";
    const AutonomousMap step;
    const std::array<std::array<unsigned, 2>, 2> profiles{{
        {10, 13},
        {5, 13},
    }};
    std::uint64_t testedOrbits = 0;
    std::vector<Witness> witnesses;
    std::array<std::array<std::uint64_t, 4>, 2> initialCounts{};
    std::array<std::uint64_t, 2> targetFirstTurnoffs{};
    std::array<std::uint64_t, 2> targetDoubleTurnoffs{};

    for (unsigned profileIndex = 0; profileIndex < profiles.size(); ++profileIndex)
    {
        const auto firstValue = profiles[profileIndex][0];
        const auto secondValue = profiles[profileIndex][1];
        for (unsigned initialFirstCount = 0; initialFirstCount <= 3; ++initialFirstCount)
        {
            const auto initialDrives = threePacketDrives(
                firstValue,
                secondValue,
                initialFirstCount);
            initialCounts[profileIndex][initialFirstCount] = initialDrives.size();
            for (const auto initialDrive : initialDrives)
            {
                auto firstState = initialDrive;
                std::uint64_t secondState = 0;
                unsigned firstWeight = 0;
                unsigned secondWeight = 0;
                for (unsigned resetNode = 1; resetNode < InnerNodes; ++resetNode)
                {
                    const auto firstOutput = accumulate(firstState);
                    firstWeight += std::popcount(firstOutput);
                    const auto secondOutput = accumulate(secondState ^ firstOutput);
                    secondWeight += std::popcount(secondOutput);
                    secondState = step(secondState ^ firstOutput);
                    firstState = step(firstState);
                    ++testedOrbits;

                    const auto resetFirstCount = 3 - initialFirstCount;
                    if (!hasValueCounts(
                            firstState,
                            firstValue,
                            secondValue,
                            resetFirstCount))
                        continue;
                    ++targetFirstTurnoffs[profileIndex];
                    if (secondState != 0)
                        continue;
                    ++targetDoubleTurnoffs[profileIndex];
                    if (secondWeight < DistanceThreshold)
                    {
                        witnesses.push_back(Witness{
                            firstValue,
                            secondValue,
                            initialFirstCount,
                            initialDrive,
                            firstState,
                            resetNode,
                            firstWeight,
                            secondWeight,
                        });
                    }
                }
            }
        }
    }

    std::ofstream output(outputPath, std::ios::trunc);
    if (!output)
        throw std::runtime_error("failed to open output path");
    output << "{\n"
        << "  \"schema\": \"riffle-dp-2lap-g4-six-packet-3x3-turnoffs-v1\",\n"
        << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n"
        << "  \"evidence_label\": \"EXACT_BOUNDED_SEARCH\",\n"
        << "  \"coverage\": \"Two active nodes, each containing three packets; value multisets 10^3 13^3 and 5^3 13^3; every reset-node offset through 32771\",\n"
        << "  \"maximum_reset_node\": " << InnerNodes - 1 << ",\n"
        << "  \"distance_threshold\": " << DistanceThreshold << ",\n"
        << "  \"initial_drive_counts\": {";
    for (unsigned profileIndex = 0; profileIndex < profiles.size(); ++profileIndex)
    {
        output << (profileIndex ? "," : "") << "\""
            << profiles[profileIndex][0] << "_" << profiles[profileIndex][1] << "\":{";
        for (unsigned count = 0; count <= 3; ++count)
            output << (count ? "," : "") << "\"" << count << "\":"
                << initialCounts[profileIndex][count];
        output << '}';
    }
    output << "},\n"
        << "  \"tested_orbit_steps\": " << testedOrbits << ",\n"
        << "  \"target_first_lap_turnoffs\": {\"10_13\":"
        << targetFirstTurnoffs[0] << ",\"5_13\":" << targetFirstTurnoffs[1]
        << "},\n"
        << "  \"target_double_turnoffs\": {\"10_13\":"
        << targetDoubleTurnoffs[0] << ",\"5_13\":" << targetDoubleTurnoffs[1]
        << "},\n"
        << "  \"below_threshold_double_turnoffs\": " << witnesses.size() << ",\n"
        << "  \"witnesses\": [";
    for (std::size_t index = 0; index < witnesses.size(); ++index)
    {
        const auto& witness = witnesses[index];
        output << (index ? ",\n    {" : "\n    {")
            << "\"packet_value_counts\":{\"" << witness.firstValue
            << "\":3,\"" << witness.secondValue << "\":3},"
            << "\"initial_first_value_count\":" << witness.initialFirstCount
            << ",\"initial_nibbles\":";
        writeNibbles(output, witness.initialDrive);
        output << ",\"reset_nibbles\":";
        writeNibbles(output, witness.resetDrive);
        output << ",\"reset_node\":" << witness.resetNode
            << ",\"first_lap_weight\":" << witness.firstWeight
            << ",\"two_lap_weight\":" << witness.secondWeight
            << '}';
    }
    if (!witnesses.empty())
        output << '\n';
    output << "  ],\n"
        << "  \"scope_limitation\": \"The search excludes six-packet episodes spread across three or more active nodes and excludes other six-packet value multisets.\"\n"
        << "}\n";

    std::cout << "candidate=Riffle DP-2Lap g=4\n";
    std::cout << "tested_orbit_steps=" << testedOrbits << '\n';
    std::cout << "target_first_lap_turnoffs="
        << targetFirstTurnoffs[0] + targetFirstTurnoffs[1] << '\n';
    std::cout << "target_double_turnoffs="
        << targetDoubleTurnoffs[0] + targetDoubleTurnoffs[1] << '\n';
    std::cout << "below_threshold_double_turnoffs=" << witnesses.size() << '\n';
    std::cout << "output=" << outputPath << '\n';
    std::cout << "status=EXACT_BOUNDED_SIX_PACKET_3X3_SEARCH\n";
}
