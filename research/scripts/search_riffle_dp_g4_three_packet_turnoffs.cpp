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
    constexpr unsigned DriveCount = 16 * 15;
    constexpr unsigned DistanceThreshold = 188'766;

    struct Word128
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
    };

    struct Witness
    {
        unsigned first = 0;
        unsigned second = 0;
        unsigned third = 0;
        unsigned secondNode = 0;
        unsigned thirdNode = 0;
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

        std::uint64_t operator()(std::uint64_t value) const
        {
            value = accumulate(value);
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

    std::uint64_t mix(std::uint64_t value)
    {
        value += 0x9e3779b97f4a7c15ULL;
        value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
        value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }

    class MultiMap
    {
    public:
        explicit MultiMap(std::size_t maximumEntries)
        {
            std::size_t capacity = 1;
            while (capacity < maximumEntries * 2)
                capacity <<= 1;
            mKeys.assign(capacity, 0);
            mValues.resize(capacity);
            mMask = capacity - 1;
        }

        void insert(std::uint64_t key, std::uint32_t value)
        {
            if (key == 0)
                throw std::runtime_error("zero autonomous orbit state");
            auto slot = static_cast<std::size_t>(mix(key)) & mMask;
            while (mKeys[slot] != 0)
                slot = (slot + 1) & mMask;
            mKeys[slot] = key;
            mValues[slot] = value;
            ++mSize;
        }

        template<typename Callback>
        void findAll(std::uint64_t key, Callback&& callback) const
        {
            if (key == 0)
                return;
            auto slot = static_cast<std::size_t>(mix(key)) & mMask;
            while (mKeys[slot] != 0)
            {
                if (mKeys[slot] == key)
                    callback(mValues[slot]);
                slot = (slot + 1) & mMask;
            }
        }

        std::size_t size() const { return mSize; }
        std::size_t capacity() const { return mKeys.size(); }

    private:
        std::vector<std::uint64_t> mKeys;
        std::vector<std::uint32_t> mValues;
        std::size_t mMask = 0;
        std::size_t mSize = 0;
    };

    void writeDrive(std::ostream& output, unsigned index)
    {
        output << "{\"slot\":" << index / 15
            << ",\"value\":" << index % 15 + 1 << "}";
    }
}

int main(int argc, char** argv)
{
    unsigned maximumSpan = 6000;
    std::string outputPath =
        "explorations/riffle_dp_g4_g2_three_packet_turnoffs_span6000.json";
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        if (argument == "--max-span" && index + 1 < argc)
            maximumSpan = static_cast<unsigned>(std::stoul(argv[++index]));
        else if (argument == "--output" && index + 1 < argc)
            outputPath = argv[++index];
        else
            throw std::invalid_argument("unknown or incomplete argument: " + argument);
    }
    if (maximumSpan < 2 || maximumSpan > 32771)
        throw std::invalid_argument("maximum span must be in [2,32771]");
    if (maximumSpan > std::numeric_limits<std::uint32_t>::max() / DriveCount)
        throw std::invalid_argument("maximum span does not fit packed metadata");

    std::array<std::uint64_t, DriveCount> drives{};
    for (unsigned slot = 0; slot < 16; ++slot)
    {
        for (unsigned value = 1; value < 16; ++value)
            drives[slot * 15 + value - 1] = std::uint64_t{value} << (4 * slot);
    }

    const AutonomousMap step;
    const auto stride = static_cast<std::size_t>(maximumSpan) + 1;
    std::vector<std::uint64_t> orbits(DriveCount * stride);
    std::vector<std::uint32_t> prefixWeights(DriveCount * stride);
    for (unsigned drive = 0; drive < DriveCount; ++drive)
    {
        auto state = drives[drive];
        orbits[drive * stride] = state;
        for (unsigned exponent = 0; exponent < maximumSpan; ++exponent)
        {
            prefixWeights[drive * stride + exponent + 1] =
                prefixWeights[drive * stride + exponent] +
                std::popcount(accumulate(state));
            state = step(state);
            orbits[drive * stride + exponent + 1] = state;
        }
    }

    MultiMap left(static_cast<std::size_t>(DriveCount) * maximumSpan);
    std::uint64_t queryCount = 0;
    std::uint64_t algebraicSolutions = 0;
    std::vector<Witness> belowThreshold;
    for (unsigned b = maximumSpan - 1; b >= 1; --b)
    {
        const auto h = b + 1;
        for (unsigned first = 0; first < DriveCount; ++first)
        {
            left.insert(
                orbits[first * stride + h],
                static_cast<std::uint32_t>(h * DriveCount + first));
        }
        for (unsigned second = 0; second < DriveCount; ++second)
        {
            const auto secondOrbit = orbits[second * stride + b];
            for (unsigned third = 0; third < DriveCount; ++third)
            {
                ++queryCount;
                const auto target = secondOrbit ^ drives[third];
                left.findAll(target, [&](std::uint32_t packed)
                {
                    ++algebraicSolutions;
                    const auto matchedH = packed / DriveCount;
                    const auto first = packed % DriveCount;
                    const auto a = matchedH - b;
                    const auto secondDrive =
                        orbits[first * stride + a] ^ drives[second];
                    auto state = secondDrive;
                    std::uint32_t weight = prefixWeights[first * stride + a];
                    for (unsigned offset = 0; offset < b; ++offset)
                    {
                        weight += std::popcount(accumulate(state));
                        state = step(state);
                    }
                    if (state != drives[third])
                        throw std::runtime_error("turnoff equation replay failed");
                    if (weight <= DistanceThreshold)
                    {
                        belowThreshold.push_back(Witness{
                            first,
                            second,
                            third,
                            a,
                            matchedH,
                            weight,
                        });
                    }
                });
            }
        }
        if (b == 1)
            break;
    }

    std::ofstream output(outputPath, std::ios::trunc);
    if (!output)
        throw std::runtime_error("failed to open output path");
    output << "{\n"
        << "  \"schema\": \"riffle-dp-g4-g2-three-packet-turnoffs-v1\",\n"
        << "  \"candidate\": \"Riffle DP g=4@g0-v1\",\n"
        << "  \"evidence_label\": \"EXACT_BOUNDED_SEARCH\",\n"
        << "  \"maximum_third_node\": " << maximumSpan << ",\n"
        << "  \"coverage\": \"three nonzero one-nibble inputs in distinct nodes 0<a<h<=maximum_third_node\",\n"
        << "  \"distance_threshold\": " << DistanceThreshold << ",\n"
        << "  \"left_entries\": " << left.size() << ",\n"
        << "  \"hash_capacity\": " << left.capacity() << ",\n"
        << "  \"lookup_queries\": " << queryCount << ",\n"
        << "  \"algebraic_turnoff_solutions\": " << algebraicSolutions << ",\n"
        << "  \"below_threshold_count\": " << belowThreshold.size() << ",\n"
        << "  \"below_threshold_witnesses\": [";
    for (std::size_t index = 0; index < belowThreshold.size(); ++index)
    {
        const auto& witness = belowThreshold[index];
        output << (index ? ",\n    {" : "\n    {") << "\"first\":";
        writeDrive(output, witness.first);
        output << ",\"second\":";
        writeDrive(output, witness.second);
        output << ",\"third\":";
        writeDrive(output, witness.third);
        output << ",\"second_node\":" << witness.secondNode
            << ",\"third_node\":" << witness.thirdNode
            << ",\"emitted_weight_before_turnoff\":" << witness.emittedWeight
            << "}";
    }
    if (!belowThreshold.empty())
        output << '\n';
    output << "  ],\n"
        << "  \"scope_limitation\": \"The search excludes packet collisions within one node and third-node offsets above the bound.\"\n"
        << "}\n";
    output.close();

    std::cout << "candidate=Riffle DP g=4@g0-v1\n";
    std::cout << "maximum_third_node=" << maximumSpan << '\n';
    std::cout << "lookup_queries=" << queryCount << '\n';
    std::cout << "algebraic_turnoff_solutions=" << algebraicSolutions << '\n';
    std::cout << "below_threshold_count=" << belowThreshold.size() << '\n';
    std::cout << "output=" << outputPath << '\n';
    std::cout << "status=EXACT_BOUNDED_THREE_PACKET_TURNOFF_SEARCH\n";
}
