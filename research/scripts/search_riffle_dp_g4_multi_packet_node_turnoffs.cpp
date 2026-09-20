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
        std::uint64_t initialDrive = 0;
        std::uint64_t resetDrive = 0;
        unsigned initialPackets = 0;
        unsigned resetPackets = 0;
        unsigned resetNode = 0;
        unsigned emittedWeight = 0;
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

    unsigned nibbleSupport(std::uint64_t value)
    {
        unsigned support = 0;
        for (unsigned slot = 0; slot < 16; ++slot)
            support += ((value >> (4 * slot)) & 0xf) != 0;
        return support;
    }

    std::vector<std::uint64_t> drives(unsigned support)
    {
        std::vector<std::uint64_t> result;
        if (support == 1)
        {
            result.reserve(240);
            for (unsigned slot = 0; slot < 16; ++slot)
            {
                for (unsigned value = 1; value < 16; ++value)
                    result.push_back(std::uint64_t{value} << (4 * slot));
            }
        }
        else if (support == 2)
        {
            result.reserve(27'000);
            for (unsigned firstSlot = 0; firstSlot < 16; ++firstSlot)
            {
                for (unsigned secondSlot = firstSlot + 1; secondSlot < 16; ++secondSlot)
                {
                    for (unsigned firstValue = 1; firstValue < 16; ++firstValue)
                    {
                        for (unsigned secondValue = 1; secondValue < 16; ++secondValue)
                        {
                            result.push_back(
                                (std::uint64_t{firstValue} << (4 * firstSlot)) |
                                (std::uint64_t{secondValue} << (4 * secondSlot)));
                        }
                    }
                }
            }
        }
        else
            throw std::invalid_argument("drive support must be one or two");
        return result;
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
        : "explorations/riffle_dp_g4_g2_multi_packet_node_turnoffs.json";
    const AutonomousMap step;
    const auto onePacket = drives(1);
    const auto twoPacket = drives(2);
    std::array<std::array<std::uint64_t, 3>, 3> algebraicCounts{};
    std::array<std::array<std::uint64_t, 3>, 3> belowCounts{};
    std::vector<Witness> witnesses;

    for (unsigned initialPackets = 1; initialPackets <= 2; ++initialPackets)
    {
        const auto& initialDrives = initialPackets == 1 ? onePacket : twoPacket;
        for (const auto initialDrive : initialDrives)
        {
            auto state = initialDrive;
            std::uint64_t emittedWeight = 0;
            for (unsigned resetNode = 1; resetNode < InnerNodes; ++resetNode)
            {
                emittedWeight += std::popcount(accumulate(state));
                state = step(state);
                const auto resetPackets = nibbleSupport(state);
                if (resetPackets == 1 || resetPackets == 2)
                {
                    ++algebraicCounts[initialPackets][resetPackets];
                    if (emittedWeight <= DistanceThreshold)
                    {
                        ++belowCounts[initialPackets][resetPackets];
                        witnesses.push_back(Witness{
                            initialDrive,
                            state,
                            initialPackets,
                            resetPackets,
                            resetNode,
                            static_cast<unsigned>(emittedWeight),
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
        << "  \"schema\": \"riffle-dp-g4-g2-multi-packet-node-turnoffs-v1\",\n"
        << "  \"candidate\": \"Riffle DP g=4@g0-v1\",\n"
        << "  \"evidence_label\": \"EXACT\",\n"
        << "  \"maximum_reset_node\": " << InnerNodes - 1 << ",\n"
        << "  \"initial_drive_counts\": {\"1\":" << onePacket.size()
        << ",\"2\":" << twoPacket.size() << "},\n"
        << "  \"algebraic_counts\": {"
        << "\"1_to_1\":" << algebraicCounts[1][1] << ','
        << "\"1_to_2\":" << algebraicCounts[1][2] << ','
        << "\"2_to_1\":" << algebraicCounts[2][1] << ','
        << "\"2_to_2\":" << algebraicCounts[2][2] << "},\n"
        << "  \"below_threshold_counts\": {"
        << "\"1_to_1\":" << belowCounts[1][1] << ','
        << "\"1_to_2\":" << belowCounts[1][2] << ','
        << "\"2_to_1\":" << belowCounts[2][1] << ','
        << "\"2_to_2\":" << belowCounts[2][2] << "},\n"
        << "  \"below_threshold_witnesses\": [";
    for (std::size_t index = 0; index < witnesses.size(); ++index)
    {
        const auto& witness = witnesses[index];
        output << (index ? ",\n    {" : "\n    {")
            << "\"initial_packets\":" << witness.initialPackets
            << ",\"reset_packets\":" << witness.resetPackets
            << ",\"initial_nibbles\":";
        writeNibbles(output, witness.initialDrive);
        output << ",\"reset_nibbles\":";
        writeNibbles(output, witness.resetDrive);
        output << ",\"reset_node\":" << witness.resetNode
            << ",\"emitted_weight_before_turnoff\":" << witness.emittedWeight
            << '}';
    }
    if (!witnesses.empty())
        output << '\n';
    output << "  ]\n}\n";
    output.close();

    std::cout << "candidate=Riffle DP g=4@g0-v1\n";
    std::cout << "initial_support_1=" << onePacket.size() << '\n';
    std::cout << "initial_support_2=" << twoPacket.size() << '\n';
    for (unsigned initial = 1; initial <= 2; ++initial)
    {
        for (unsigned reset = 1; reset <= 2; ++reset)
        {
            std::cout << initial << "_to_" << reset
                << "_algebraic=" << algebraicCounts[initial][reset]
                << " below_threshold=" << belowCounts[initial][reset] << '\n';
        }
    }
    std::cout << "output=" << outputPath << '\n';
    std::cout << "status=EXACT_G2_MULTI_PACKET_NODE_TURNOFF_SEARCH\n";
}
