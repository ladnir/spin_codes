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
    using Columns = std::array<std::uint64_t, 64>;
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
        unsigned initialPackets = 0;
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
            throw std::runtime_error("six-packet search: generator inverse failed");
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

    Columns systematicRightColumns()
    {
        const auto inverse = generatorInverse();
        Columns columns{};
        for (unsigned column = 0; column < 64; ++column)
        {
            const auto message = carrylessMultiplyLow(
                std::uint64_t{1} << column,
                inverse);
            const auto word = encode(message);
            if (word.low != (std::uint64_t{1} << column))
                throw std::runtime_error("six-packet search: systematic columns failed");
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

    std::uint64_t applyColumns(const Columns& columns, std::uint64_t value)
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

    Columns forwardColumns()
    {
        const auto right = systematicRightColumns();
        Columns result{};
        for (unsigned bit = 0; bit < 64; ++bit)
            result[bit] = applyColumns(right, accumulate(std::uint64_t{1} << bit));
        return result;
    }

    Columns inverseColumns(const Columns& columns)
    {
        std::array<std::uint64_t, 64> rows{};
        std::array<std::uint64_t, 64> inverseRows{};
        for (unsigned row = 0; row < 64; ++row)
        {
            for (unsigned column = 0; column < 64; ++column)
                rows[row] |= ((columns[column] >> row) & 1) << column;
            inverseRows[row] = std::uint64_t{1} << row;
        }
        for (unsigned column = 0; column < 64; ++column)
        {
            unsigned pivot = column;
            while (pivot < 64 && ((rows[pivot] >> column) & 1) == 0)
                ++pivot;
            if (pivot == 64)
                throw std::runtime_error("six-packet search: singular autonomous map");
            std::swap(rows[column], rows[pivot]);
            std::swap(inverseRows[column], inverseRows[pivot]);
            for (unsigned row = 0; row < 64; ++row)
            {
                if (row != column && ((rows[row] >> column) & 1))
                {
                    rows[row] ^= rows[column];
                    inverseRows[row] ^= inverseRows[column];
                }
            }
        }
        Columns result{};
        for (unsigned column = 0; column < 64; ++column)
        {
            for (unsigned row = 0; row < 64; ++row)
                result[column] |= ((inverseRows[row] >> column) & 1) << row;
        }
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            if (applyColumns(columns, result[bit]) != (std::uint64_t{1} << bit))
                throw std::runtime_error("six-packet search: inverse check failed");
        }
        return result;
    }

    class ByteMap
    {
    public:
        explicit ByteMap(const Columns& columns)
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

    std::vector<std::uint64_t> drives(
        unsigned firstValue,
        unsigned secondValue,
        unsigned packetCount,
        unsigned firstCount)
    {
        if (packetCount == 0 || packetCount > 5 || firstCount > packetCount)
            throw std::invalid_argument("six-packet search: invalid drive profile");
        std::vector<std::uint64_t> result;
        const auto slotMasks = std::uint32_t{1} << 16;
        for (std::uint32_t occupied = 1; occupied < slotMasks; ++occupied)
        {
            if (static_cast<unsigned>(std::popcount(occupied)) != packetCount)
                continue;
            std::array<unsigned, 5> slots{};
            unsigned cursor = 0;
            for (unsigned slot = 0; slot < 16; ++slot)
            {
                if ((occupied >> slot) & 1)
                    slots[cursor++] = slot;
            }
            const auto assignments = std::uint32_t{1} << packetCount;
            for (std::uint32_t assignment = 0; assignment < assignments; ++assignment)
            {
                if (static_cast<unsigned>(std::popcount(assignment)) != firstCount)
                    continue;
                std::uint64_t drive = 0;
                for (unsigned index = 0; index < packetCount; ++index)
                {
                    const auto value = ((assignment >> index) & 1)
                        ? firstValue
                        : secondValue;
                    drive |= std::uint64_t{value} << (4 * slots[index]);
                }
                result.push_back(drive);
            }
        }
        return result;
    }

    bool hasCounts(
        std::uint64_t drive,
        unsigned firstValue,
        unsigned secondValue,
        unsigned firstCount,
        unsigned secondCount)
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
        return observedFirst == firstCount && observedSecond == secondCount;
    }

    std::pair<unsigned, std::uint64_t> runLap(
        const std::vector<std::uint64_t>& inputs,
        std::uint64_t initialState,
        const Columns& rightColumns)
    {
        unsigned weight = 0;
        auto state = initialState;
        for (const auto input : inputs)
        {
            const auto output = accumulate(state ^ input);
            weight += std::popcount(output);
            state = applyColumns(rightColumns, output);
        }
        return {weight, state};
    }

    Witness replay(
        unsigned firstValue,
        unsigned secondValue,
        unsigned initialPackets,
        std::uint64_t initialDrive,
        std::uint64_t resetDrive,
        unsigned resetNode,
        const Columns& rightColumns)
    {
        std::vector<std::uint64_t> inputs(resetNode + 1);
        inputs[0] = initialDrive;
        inputs[resetNode] = resetDrive;
        auto firstState = std::uint64_t{0};
        std::vector<std::uint64_t> firstOutputs;
        firstOutputs.reserve(inputs.size());
        unsigned firstWeight = 0;
        for (const auto input : inputs)
        {
            const auto output = accumulate(firstState ^ input);
            firstOutputs.push_back(output);
            firstWeight += std::popcount(output);
            firstState = applyColumns(rightColumns, output);
        }
        const auto [secondWeight, secondState] = runLap(
            firstOutputs,
            firstState,
            rightColumns);
        if (firstState != 0)
            throw std::runtime_error("six-packet search: first turnoff replay failed");
        return Witness{
            firstValue,
            secondValue,
            initialPackets,
            initialDrive,
            resetDrive,
            resetNode,
            firstWeight,
            secondState == 0 ? secondWeight : DistanceThreshold,
        };
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
        : "explorations/riffle_dp_2lap_g4_six_packet_two_node_turnoffs.json";
    const auto rightColumns = systematicRightColumns();
    const auto forwardMapColumns = forwardColumns();
    const ByteMap forward(forwardMapColumns);
    const ByteMap backward(inverseColumns(forwardMapColumns));
    const std::array<std::array<unsigned, 2>, 2> profiles{{
        {10, 13},
        {5, 13},
    }};
    std::array<std::array<std::uint64_t, 6>, 2> testedSteps{};
    std::array<std::array<std::uint64_t, 6>, 2> firstTurnoffs{};
    std::array<std::array<std::uint64_t, 6>, 2> doubleTurnoffs{};
    std::vector<Witness> witnesses;

    for (unsigned profileIndex = 0; profileIndex < profiles.size(); ++profileIndex)
    {
        const auto firstValue = profiles[profileIndex][0];
        const auto secondValue = profiles[profileIndex][1];

        // Forward search covers initial packet counts one through three.
        for (unsigned initialPackets = 1; initialPackets <= 3; ++initialPackets)
        {
            for (unsigned initialFirst = 0; initialFirst <= 3; ++initialFirst)
            {
                const auto initialSecond = initialPackets - initialFirst;
                if (initialFirst > initialPackets || initialSecond > 3)
                    continue;
                const auto initialDrives = drives(
                    firstValue, secondValue, initialPackets, initialFirst);
                for (const auto initialDrive : initialDrives)
                {
                    auto resetDrive = initialDrive;
                    for (unsigned resetNode = 1; resetNode < InnerNodes; ++resetNode)
                    {
                        resetDrive = forward(resetDrive);
                        ++testedSteps[profileIndex][initialPackets];
                        if (!hasCounts(
                                resetDrive,
                                firstValue,
                                secondValue,
                                3 - initialFirst,
                                3 - initialSecond))
                            continue;
                        ++firstTurnoffs[profileIndex][initialPackets];
                        const auto witness = replay(
                            firstValue,
                            secondValue,
                            initialPackets,
                            initialDrive,
                            resetDrive,
                            resetNode,
                            rightColumns);
                        if (witness.secondWeight < DistanceThreshold)
                        {
                            ++doubleTurnoffs[profileIndex][initialPackets];
                            witnesses.push_back(witness);
                        }
                    }
                }
            }
        }

        // Reverse search covers initial packet counts four and five cheaply.
        for (unsigned resetPackets = 1; resetPackets <= 2; ++resetPackets)
        {
            const auto initialPackets = 6 - resetPackets;
            for (unsigned resetFirst = 0; resetFirst <= 3; ++resetFirst)
            {
                const auto resetSecond = resetPackets - resetFirst;
                if (resetFirst > resetPackets || resetSecond > 3)
                    continue;
                const auto resetDrives = drives(
                    firstValue, secondValue, resetPackets, resetFirst);
                for (const auto resetDrive : resetDrives)
                {
                    auto initialDrive = resetDrive;
                    for (unsigned resetNode = 1; resetNode < InnerNodes; ++resetNode)
                    {
                        initialDrive = backward(initialDrive);
                        ++testedSteps[profileIndex][initialPackets];
                        if (!hasCounts(
                                initialDrive,
                                firstValue,
                                secondValue,
                                3 - resetFirst,
                                3 - resetSecond))
                            continue;
                        ++firstTurnoffs[profileIndex][initialPackets];
                        const auto witness = replay(
                            firstValue,
                            secondValue,
                            initialPackets,
                            initialDrive,
                            resetDrive,
                            resetNode,
                            rightColumns);
                        if (witness.secondWeight < DistanceThreshold)
                        {
                            ++doubleTurnoffs[profileIndex][initialPackets];
                            witnesses.push_back(witness);
                        }
                    }
                }
            }
        }
    }

    std::ofstream output(outputPath, std::ios::trunc);
    if (!output)
        throw std::runtime_error("six-packet search: output open failed");
    output << "{\n"
        << "  \"schema\": \"riffle-dp-2lap-g4-six-packet-two-node-turnoffs-v1\",\n"
        << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n"
        << "  \"evidence_label\": \"EXACT_BOUNDED_SEARCH\",\n"
        << "  \"coverage\": \"Exactly two active nodes containing all six packets; all splits 1+5 through 5+1; value multisets 10^3 13^3 and 5^3 13^3; every reset offset through 32771\",\n"
        << "  \"maximum_reset_node\": " << InnerNodes - 1 << ",\n"
        << "  \"distance_threshold\": " << DistanceThreshold << ",\n"
        << "  \"profiles\": [";
    for (unsigned profileIndex = 0; profileIndex < profiles.size(); ++profileIndex)
    {
        output << (profileIndex ? ",\n    {" : "\n    {")
            << "\"packet_value_counts\":{\"" << profiles[profileIndex][0]
            << "\":3,\"" << profiles[profileIndex][1] << "\":3},"
            << "\"splits\":[";
        for (unsigned initialPackets = 1; initialPackets <= 5; ++initialPackets)
        {
            output << (initialPackets == 1 ? "" : ",")
                << "{\"initial_packets\":" << initialPackets
                << ",\"reset_packets\":" << 6 - initialPackets
                << ",\"tested_orbit_steps\":" << testedSteps[profileIndex][initialPackets]
                << ",\"first_lap_turnoffs\":" << firstTurnoffs[profileIndex][initialPackets]
                << ",\"below_threshold_double_turnoffs\":"
                << doubleTurnoffs[profileIndex][initialPackets] << '}';
        }
        output << "]}";
    }
    output << "\n  ],\n"
        << "  \"below_threshold_double_turnoffs\": " << witnesses.size() << ",\n"
        << "  \"witnesses\": [";
    for (std::size_t index = 0; index < witnesses.size(); ++index)
    {
        const auto& witness = witnesses[index];
        output << (index ? ",\n    {" : "\n    {")
            << "\"packet_value_counts\":{\"" << witness.firstValue
            << "\":3,\"" << witness.secondValue << "\":3},"
            << "\"initial_packets\":" << witness.initialPackets
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
        << "  \"scope_limitation\": \"The search excludes six-packet episodes with three or more active nodes and excludes other value multisets.\"\n"
        << "}\n";

    std::uint64_t totalSteps = 0;
    std::uint64_t totalFirstTurnoffs = 0;
    for (unsigned profileIndex = 0; profileIndex < profiles.size(); ++profileIndex)
    {
        for (unsigned initialPackets = 1; initialPackets <= 5; ++initialPackets)
        {
            totalSteps += testedSteps[profileIndex][initialPackets];
            totalFirstTurnoffs += firstTurnoffs[profileIndex][initialPackets];
        }
    }
    std::cout << "candidate=Riffle DP-2Lap g=4\n";
    std::cout << "tested_orbit_steps=" << totalSteps << '\n';
    std::cout << "first_lap_turnoffs=" << totalFirstTurnoffs << '\n';
    std::cout << "below_threshold_double_turnoffs=" << witnesses.size() << '\n';
    std::cout << "output=" << outputPath << '\n';
    std::cout << "status=EXACT_BOUNDED_SIX_PACKET_TWO_NODE_SEARCH\n";
}
