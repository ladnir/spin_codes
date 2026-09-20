#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned Nodes = 32772;
    constexpr std::uint64_t ResponseLength = 64ULL * Nodes;
    constexpr std::uint64_t BiasThreshold = 106802;
    constexpr unsigned RestartsPerSupport = 8;
    constexpr std::uint64_t RngSeed = 0x6A09E667F3BCC909ULL;
    constexpr std::array<unsigned, 7> Degrees{1, 2, 4, 9, 10, 18, 20};
    constexpr std::array<std::uint64_t, 7> Factors{
        0x3ULL, 0x7ULL, 0x13ULL, 0x373ULL, 0x519ULL, 0x7C9C3ULL, 0x1E1FFFULL};

    std::uint64_t carrylessMultiplyLow(std::uint64_t left, std::uint64_t right)
    {
        std::uint64_t result = 0;
        while (right)
        {
            const unsigned bit = std::countr_zero(right);
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
            throw std::runtime_error("high-support bias search: generator inverse failed");
        return inverse;
    }

    struct Word128
    {
        std::uint64_t low = 0;
        std::uint64_t high = 0;
    };

    Word128 encodeBch(std::uint64_t message)
    {
        Word128 result{};
        while (message)
        {
            const unsigned bit = std::countr_zero(message);
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
        const std::uint64_t inverse = generatorInverse();
        std::array<std::uint64_t, 64> columns{};
        for (unsigned column = 0; column < 64; ++column)
        {
            const std::uint64_t message = carrylessMultiplyLow(
                std::uint64_t{1} << column, inverse);
            const Word128 word = encodeBch(message);
            if (word.low != (std::uint64_t{1} << column))
                throw std::runtime_error("high-support bias search: systematic BCH failed");
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

    std::uint64_t applyPolynomial(
        const ByteMap& step,
        std::uint64_t polynomial,
        std::uint64_t value)
    {
        std::uint64_t result = 0;
        std::uint64_t current = value;
        while (polynomial)
        {
            if (polynomial & 1)
                result ^= current;
            polynomial >>= 1;
            current = step(current);
        }
        return result;
    }

    std::vector<std::uint64_t> kernelBasis(
        const ByteMap& step,
        std::uint64_t polynomial,
        unsigned expectedDimension)
    {
        std::array<std::uint64_t, 64> pivotValues{};
        std::array<std::uint64_t, 64> pivotRepresentations{};
        std::vector<std::uint64_t> basis;
        for (unsigned column = 0; column < 64; ++column)
        {
            std::uint64_t value = applyPolynomial(
                step, polynomial, std::uint64_t{1} << column);
            std::uint64_t representation = std::uint64_t{1} << column;
            while (value)
            {
                const unsigned pivot = std::bit_width(value) - 1;
                if (pivotValues[pivot])
                {
                    value ^= pivotValues[pivot];
                    representation ^= pivotRepresentations[pivot];
                }
                else
                {
                    pivotValues[pivot] = value;
                    pivotRepresentations[pivot] = representation;
                    break;
                }
            }
            if (value == 0)
                basis.push_back(representation);
        }
        if (basis.size() != expectedDimension)
            throw std::runtime_error("high-support bias search: kernel dimension mismatch");
        return basis;
    }

    std::uint64_t xorshift64star(std::uint64_t& state)
    {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        return state * 0x2545F4914F6CDD1DULL;
    }

    struct Witness
    {
        std::uint64_t coordinates = 0;
        std::uint64_t physicalState = 0;
        std::uint64_t responseWeight = 0;
        std::int64_t signedBias = 0;
        unsigned iterations = 0;
    };

    struct SupportResult
    {
        unsigned mask = 0;
        unsigned dimension = 0;
        std::uint64_t evaluations = 0;
        Witness minimumWeight{};
        Witness maximumWeight{};
    };

    std::uint64_t physicalState(
        std::uint64_t coordinates,
        const std::array<std::uint64_t, 64>& basis)
    {
        std::uint64_t state = 0;
        while (coordinates)
        {
            const unsigned bit = std::countr_zero(coordinates);
            state ^= basis[bit];
            coordinates &= coordinates - 1;
        }
        return state;
    }

    std::uint64_t responseWeight(const std::vector<std::uint64_t>& word)
    {
        std::uint64_t result = 0;
        for (const std::uint64_t node : word)
            result += std::popcount(node);
        return result;
    }

    std::uint64_t flippedWeight(
        const std::vector<std::uint64_t>& word,
        const std::vector<std::uint64_t>& row)
    {
        std::uint64_t result = 0;
        for (unsigned node = 0; node < Nodes; ++node)
            result += std::popcount(word[node] ^ row[node]);
        return result;
    }

    bool hasExactMask(
        std::uint64_t coordinates,
        unsigned supportMask,
        const std::array<std::uint64_t, 7>& segmentMasks)
    {
        for (unsigned component = 0; component < 7; ++component)
        {
            const bool present = (coordinates & segmentMasks[component]) != 0;
            if (present != (((supportMask >> component) & 1) != 0))
                return false;
        }
        return true;
    }

    Witness greedyClimb(
        std::uint64_t start,
        bool minimize,
        unsigned supportMask,
        const std::array<std::uint64_t, 7>& segmentMasks,
        const std::array<std::uint64_t, 64>& basis,
        const std::array<std::vector<std::uint64_t>, 64>& rows,
        std::uint64_t& evaluations)
    {
        if (!hasExactMask(start, supportMask, segmentMasks))
            throw std::runtime_error("high-support bias search: invalid start mask");
        std::vector<std::uint64_t> word(Nodes, 0);
        std::uint64_t active = start;
        while (active)
        {
            const unsigned bit = std::countr_zero(active);
            for (unsigned node = 0; node < Nodes; ++node)
                word[node] ^= rows[bit][node];
            active &= active - 1;
        }
        std::uint64_t weight = responseWeight(word);
        unsigned iterations = 0;
        for (;;)
        {
            std::uint64_t bestWeight = weight;
            unsigned bestBit = 64;
            for (unsigned bit = 0; bit < 64; ++bit)
            {
                const std::uint64_t candidate = start ^ (std::uint64_t{1} << bit);
                if (!hasExactMask(candidate, supportMask, segmentMasks))
                    continue;
                const std::uint64_t candidateWeight = flippedWeight(word, rows[bit]);
                ++evaluations;
                const bool improves = minimize
                    ? candidateWeight < bestWeight
                    : candidateWeight > bestWeight;
                if (improves)
                {
                    bestWeight = candidateWeight;
                    bestBit = bit;
                }
            }
            if (bestBit == 64)
                break;
            start ^= std::uint64_t{1} << bestBit;
            for (unsigned node = 0; node < Nodes; ++node)
                word[node] ^= rows[bestBit][node];
            weight = bestWeight;
            ++iterations;
        }
        const std::int64_t bias = static_cast<std::int64_t>(ResponseLength)
            - 2 * static_cast<std::int64_t>(weight);
        return Witness{start, physicalState(start, basis), weight, bias, iterations};
    }

    std::uint64_t startCoordinates(
        unsigned supportMask,
        unsigned restart,
        const std::array<unsigned, 7>& offsets,
        const std::array<std::uint64_t, 7>& segmentMasks,
        std::uint64_t& rng)
    {
        std::uint64_t result = 0;
        for (unsigned component = 0; component < 7; ++component)
        {
            if (!((supportMask >> component) & 1))
                continue;
            std::uint64_t local = 0;
            if (restart == 0)
                local = 1;
            else if (restart == 1)
                local = (std::uint64_t{1} << Degrees[component]) - 1;
            else
            {
                const std::uint64_t localMask =
                    (std::uint64_t{1} << Degrees[component]) - 1;
                do
                {
                    local = xorshift64star(rng) & localMask;
                }
                while (local == 0);
            }
            result |= local << offsets[component];
        }
        if (!hasExactMask(result, supportMask, segmentMasks))
            throw std::runtime_error("high-support bias search: start construction failed");
        return result;
    }

    void writeHex(std::ostream& output, std::uint64_t value)
    {
        output << "\"0x" << std::hex << value << std::dec << "\"";
    }

    void writeWitness(std::ostream& output, const Witness& witness, unsigned indent)
    {
        const std::string spaces(indent, ' ');
        output << "{\n"
               << spaces << "  \"component_coordinates_hex\": ";
        writeHex(output, witness.coordinates);
        output << ",\n" << spaces << "  \"physical_state_hex\": ";
        writeHex(output, witness.physicalState);
        output << ",\n"
               << spaces << "  \"response_weight\": " << witness.responseWeight << ",\n"
               << spaces << "  \"signed_response_bias\": " << witness.signedBias << ",\n"
               << spaces << "  \"absolute_response_bias\": "
               << std::uint64_t(witness.signedBias < 0 ? -witness.signedBias : witness.signedBias) << ",\n"
               << spaces << "  \"greedy_iterations\": " << witness.iterations << "\n"
               << spaces << "}";
    }
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
            throw std::runtime_error("usage: search.exe OUTPUT.json");
        const auto started = std::chrono::steady_clock::now();

        const ByteMap parityMap(systematicRightColumns());
        std::array<std::uint64_t, 64> stepColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            stepColumns[bit] = parityMap(accumulate(std::uint64_t{1} << bit));
        const ByteMap step(stepColumns);

        std::array<std::vector<std::uint64_t>, 7> componentBases;
        std::array<std::uint64_t, 64> basis{};
        std::array<unsigned, 7> offsets{};
        std::array<std::uint64_t, 7> segmentMasks{};
        unsigned basisOffset = 0;
        for (unsigned component = 0; component < 7; ++component)
        {
            componentBases[component] = kernelBasis(step, Factors[component], Degrees[component]);
            offsets[component] = basisOffset;
            segmentMasks[component] =
                ((std::uint64_t{1} << Degrees[component]) - 1) << basisOffset;
            for (const std::uint64_t state : componentBases[component])
                basis[basisOffset++] = state;
        }
        if (basisOffset != 64)
            throw std::runtime_error("high-support bias search: component dimension mismatch");

        std::array<std::vector<std::uint64_t>, 64> responseRows;
        for (unsigned bit = 0; bit < 64; ++bit)
        {
            responseRows[bit].resize(Nodes);
            std::uint64_t state = basis[bit];
            for (unsigned node = 0; node < Nodes; ++node)
            {
                responseRows[bit][node] = accumulate(state);
                state = step(state);
            }
        }

        std::uint64_t rng = RngSeed;
        std::vector<SupportResult> results;
        std::uint64_t totalEvaluations = 0;
        Witness globalBest{};
        unsigned globalBestMask = 0;
        for (unsigned supportMask = 1; supportMask < 128; ++supportMask)
        {
            unsigned dimension = 0;
            for (unsigned component = 0; component < 7; ++component)
            {
                if ((supportMask >> component) & 1)
                    dimension += Degrees[component];
            }
            if (dimension < 23)
                continue;

            SupportResult result{};
            result.mask = supportMask;
            result.dimension = dimension;
            result.minimumWeight.responseWeight = std::numeric_limits<std::uint64_t>::max();
            result.maximumWeight.responseWeight = 0;
            for (unsigned restart = 0; restart < RestartsPerSupport; ++restart)
            {
                const std::uint64_t start = startCoordinates(
                    supportMask, restart, offsets, segmentMasks, rng);
                const Witness minimum = greedyClimb(
                    start, true, supportMask, segmentMasks, basis, responseRows, result.evaluations);
                const Witness maximum = greedyClimb(
                    start, false, supportMask, segmentMasks, basis, responseRows, result.evaluations);
                if (minimum.responseWeight < result.minimumWeight.responseWeight)
                    result.minimumWeight = minimum;
                if (maximum.responseWeight > result.maximumWeight.responseWeight)
                    result.maximumWeight = maximum;
            }
            totalEvaluations += result.evaluations;
            for (const Witness* witness : {&result.minimumWeight, &result.maximumWeight})
            {
                const std::uint64_t absoluteBias = witness->signedBias < 0
                    ? std::uint64_t(-witness->signedBias)
                    : std::uint64_t(witness->signedBias);
                const std::uint64_t globalAbsolute = globalBest.signedBias < 0
                    ? std::uint64_t(-globalBest.signedBias)
                    : std::uint64_t(globalBest.signedBias);
                if (absoluteBias > globalAbsolute)
                {
                    globalBest = *witness;
                    globalBestMask = supportMask;
                }
            }
            std::cerr << "support=0x" << std::hex << supportMask << std::dec
                      << " dimension=" << dimension
                      << " min=" << result.minimumWeight.responseWeight
                      << " max=" << result.maximumWeight.responseWeight
                      << " max_abs_bias="
                      << std::max(
                             std::uint64_t(result.minimumWeight.signedBias < 0 ? -result.minimumWeight.signedBias : result.minimumWeight.signedBias),
                             std::uint64_t(result.maximumWeight.signedBias < 0 ? -result.maximumWeight.signedBias : result.maximumWeight.signedBias))
                      << '\n';
            results.push_back(result);
        }

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        const std::uint64_t globalAbsolute = globalBest.signedBias < 0
            ? std::uint64_t(-globalBest.signedBias)
            : std::uint64_t(globalBest.signedBias);

        std::ofstream output(argv[1]);
        if (!output)
            throw std::runtime_error("high-support bias search: cannot open output");
        output << "{\n"
               << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-g4-goal02-high-support-bias-search-v1\",\n"
               << "  \"candidate\": \"Riffle PacketMul-WrapMul-2Lap g=4\",\n"
               << "  \"evidence_label\": \"BOUNDED_DIAGNOSTIC_HIGH_SUPPORT_BIAS_SEARCH\",\n"
               << "  \"response_nodes\": " << Nodes << ",\n"
               << "  \"response_length\": " << ResponseLength << ",\n"
               << "  \"bias_threshold\": " << BiasThreshold << ",\n"
               << "  \"first_violating_even_bias\": " << (BiasThreshold + 2) << ",\n"
               << "  \"component_degrees\": [1, 2, 4, 9, 10, 18, 20],\n"
               << "  \"searched_support_dimension_minimum\": 23,\n"
               << "  \"searched_support_count\": " << results.size() << ",\n"
               << "  \"restart_law\": {\n"
               << "    \"restarts_per_support_per_direction\": " << RestartsPerSupport << ",\n"
               << "    \"deterministic_starts\": [\"one low coordinate per active component\", \"all coordinates in every active component\"],\n"
               << "    \"remaining_starts\": \"independent nonzero component coordinates from xorshift64star\",\n"
               << "    \"rng_seed_hex\": \"0x" << std::hex << RngSeed << std::dec << "\",\n"
               << "    \"move_rule\": \"strict best single-coordinate flip preserving the exact component-support mask\"\n"
               << "  },\n"
               << "  \"exact_bias_evaluations\": " << totalEvaluations << ",\n"
               << "  \"support_rows\": [\n";
        for (unsigned index = 0; index < results.size(); ++index)
        {
            const SupportResult& result = results[index];
            output << "    {\n"
                   << "      \"component_support_mask_hex\": \"0x" << std::hex << result.mask << std::dec << "\",\n"
                   << "      \"dimension\": " << result.dimension << ",\n"
                   << "      \"exact_bias_evaluations\": " << result.evaluations << ",\n"
                   << "      \"minimum_weight_witness\": ";
            writeWitness(output, result.minimumWeight, 6);
            output << ",\n      \"maximum_weight_witness\": ";
            writeWitness(output, result.maximumWeight, 6);
            output << "\n    }" << (index + 1 == results.size() ? "\n" : ",\n");
        }
        output << "  ],\n"
               << "  \"global_best_support_mask_hex\": \"0x" << std::hex << globalBestMask << std::dec << "\",\n"
               << "  \"global_best_witness\": ";
        writeWitness(output, globalBest, 2);
        output << ",\n"
               << "  \"maximum_found_absolute_bias\": " << globalAbsolute << ",\n"
               << "  \"violation_found\": " << (globalAbsolute > BiasThreshold ? "true" : "false") << ",\n"
               << "  \"elapsed_seconds\": " << std::setprecision(17) << elapsed << ",\n"
               << "  \"scope_limit\": \"Every reported bias is exact, but greedy restarts do not exhaust any support with dimension at least 23.\",\n"
               << "  \"result\": \"" << (globalAbsolute > BiasThreshold
                      ? "DIAGNOSTIC_BIAS_VIOLATION_FOUND"
                      : "NO_VIOLATION_FOUND_IN_BOUNDED_SEARCH") << "\"\n"
               << "}\n";
        output.close();

        std::cout << "output=" << argv[1] << '\n'
                  << "supports=" << results.size() << '\n'
                  << "evaluations=" << totalEvaluations << '\n'
                  << "maximum_absolute_bias=" << globalAbsolute << '\n'
                  << "violation=" << (globalAbsolute > BiasThreshold ? "yes" : "no") << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
