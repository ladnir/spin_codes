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
    constexpr std::uint32_t ResponseNodes = 32772;
    constexpr std::uint64_t LowWeightThreshold = 377530;
    constexpr unsigned ComponentDimensionCap = 22;
    constexpr std::array<unsigned, 7> ComponentDegrees{1, 2, 4, 9, 10, 18, 20};
    constexpr std::array<std::uint64_t, 7> ComponentFactors{
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

    std::uint64_t polynomialMultiply(std::uint64_t left, std::uint64_t right)
    {
        return carrylessMultiplyLow(left, right);
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
            throw std::runtime_error("low components: generator inverse failed");
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
                throw std::runtime_error("low components: systematic BCH failed");
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
            throw std::runtime_error("low components: kernel dimension mismatch");
        return basis;
    }

    struct CoordinateSolver
    {
        std::array<std::uint64_t, 64> pivotValues{};
        std::array<std::uint64_t, 64> pivotCoordinates{};

        explicit CoordinateSolver(const std::vector<std::uint64_t>& basis)
        {
            for (unsigned coordinate = 0; coordinate < basis.size(); ++coordinate)
            {
                std::uint64_t value = basis[coordinate];
                std::uint64_t representation = std::uint64_t{1} << coordinate;
                while (value)
                {
                    const unsigned pivot = std::bit_width(value) - 1;
                    if (pivotValues[pivot])
                    {
                        value ^= pivotValues[pivot];
                        representation ^= pivotCoordinates[pivot];
                    }
                    else
                    {
                        pivotValues[pivot] = value;
                        pivotCoordinates[pivot] = representation;
                        break;
                    }
                }
                if (value == 0)
                    throw std::runtime_error("low components: dependent combined basis");
            }
        }

        std::uint64_t coordinates(std::uint64_t value) const
        {
            std::uint64_t result = 0;
            while (value)
            {
                const unsigned pivot = std::bit_width(value) - 1;
                if (!pivotValues[pivot])
                    throw std::runtime_error("low components: state outside component basis");
                value ^= pivotValues[pivot];
                result ^= pivotCoordinates[pivot];
            }
            return result;
        }
    };

    std::uint64_t stateFromCoordinates(
        std::uint32_t coordinates,
        const std::vector<std::uint64_t>& basis)
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

    struct Row
    {
        unsigned supportMask = 0;
        unsigned dimension = 0;
        std::uint64_t exactStateCount = 0;
        std::uint64_t visitedStateCount = 0;
        std::uint64_t cycleCount = 0;
        std::uint64_t lowStateCount = 0;
        std::uint64_t minimumResponseWeight = std::numeric_limits<std::uint64_t>::max();
        std::uint32_t minimumCoordinates = 0;
        std::uint64_t minimumState = 0;
        std::uint64_t minimumCyclePeriod = 0;
        std::uint64_t maximumAbsoluteResponseBias = 0;
        std::uint32_t maximumBiasCoordinates = 0;
        std::uint64_t maximumBiasState = 0;
        std::uint64_t maximumBiasResponseWeight = 0;
    };

    bool hasExactSupport(
        std::uint32_t coordinates,
        const std::vector<std::uint32_t>& segmentMasks)
    {
        for (const std::uint32_t mask : segmentMasks)
        {
            if ((coordinates & mask) == 0)
                return false;
        }
        return true;
    }

    Row enumerateSupport(
        unsigned supportMask,
        const ByteMap& step,
        const std::array<std::vector<std::uint64_t>, 7>& componentBases)
    {
        Row row{};
        row.supportMask = supportMask;
        std::vector<std::uint64_t> basis;
        std::vector<std::uint32_t> segmentMasks;
        for (unsigned component = 0; component < ComponentDegrees.size(); ++component)
        {
            if (!((supportMask >> component) & 1))
                continue;
            const unsigned offset = static_cast<unsigned>(basis.size());
            basis.insert(
                basis.end(),
                componentBases[component].begin(),
                componentBases[component].end());
            const unsigned degree = ComponentDegrees[component];
            segmentMasks.push_back(
                ((std::uint32_t{1} << degree) - 1) << offset);
            row.exactStateCount = row.exactStateCount == 0
                ? ((std::uint64_t{1} << degree) - 1)
                : row.exactStateCount * ((std::uint64_t{1} << degree) - 1);
        }
        row.dimension = static_cast<unsigned>(basis.size());
        const std::uint32_t size = std::uint32_t{1} << row.dimension;
        const CoordinateSolver solver(basis);

        std::vector<std::uint32_t> nextColumns(row.dimension);
        for (unsigned bit = 0; bit < row.dimension; ++bit)
            nextColumns[bit] = static_cast<std::uint32_t>(solver.coordinates(step(basis[bit])));

        std::vector<std::uint32_t> next(size);
        std::vector<std::uint8_t> emittedWeights(size);
        std::vector<std::uint64_t> states(size);
        for (std::uint32_t coordinates = 1; coordinates < size; ++coordinates)
        {
            const unsigned bit = std::countr_zero(coordinates);
            const std::uint32_t rest = coordinates & (coordinates - 1);
            next[coordinates] = next[rest] ^ nextColumns[bit];
            states[coordinates] = states[rest] ^ basis[bit];
            emittedWeights[coordinates] = static_cast<std::uint8_t>(
                std::popcount(accumulate(states[coordinates])));
        }

        std::vector<std::uint8_t> visited(size);
        std::vector<std::uint32_t> cycle;
        cycle.reserve(size);
        for (std::uint32_t start = 1; start < size; ++start)
        {
            if (visited[start] || !hasExactSupport(start, segmentMasks))
                continue;
            cycle.clear();
            std::uint32_t current = start;
            std::uint64_t cycleTotal = 0;
            do
            {
                if (visited[current] || !hasExactSupport(current, segmentMasks))
                    throw std::runtime_error("low components: exact-support orbit escaped or merged");
                visited[current] = 1;
                cycle.push_back(current);
                cycleTotal += emittedWeights[current];
                current = next[current];
            }
            while (current != start);

            ++row.cycleCount;
            row.visitedStateCount += cycle.size();
            const std::uint64_t period = cycle.size();
            const std::uint64_t quotient = ResponseNodes / period;
            const std::uint64_t remainder = ResponseNodes % period;
            std::uint64_t rolling = 0;
            for (std::uint64_t index = 0; index < remainder; ++index)
                rolling += emittedWeights[cycle[index]];
            for (std::uint64_t index = 0; index < period; ++index)
            {
                const std::uint64_t responseWeight = quotient * cycleTotal + rolling;
                if (responseWeight <= LowWeightThreshold)
                    ++row.lowStateCount;
                if (responseWeight < row.minimumResponseWeight)
                {
                    row.minimumResponseWeight = responseWeight;
                    row.minimumCoordinates = cycle[index];
                    row.minimumState = states[cycle[index]];
                    row.minimumCyclePeriod = period;
                }
                const std::uint64_t responseLength = 64 * std::uint64_t{ResponseNodes};
                const std::uint64_t twiceWeight = 2 * responseWeight;
                const std::uint64_t absoluteBias = twiceWeight <= responseLength
                    ? responseLength - twiceWeight
                    : twiceWeight - responseLength;
                if (absoluteBias > row.maximumAbsoluteResponseBias)
                {
                    row.maximumAbsoluteResponseBias = absoluteBias;
                    row.maximumBiasCoordinates = cycle[index];
                    row.maximumBiasState = states[cycle[index]];
                    row.maximumBiasResponseWeight = responseWeight;
                }
                if (remainder)
                {
                    rolling -= emittedWeights[cycle[index]];
                    rolling += emittedWeights[cycle[(index + remainder) % period]];
                }
            }
        }
        if (row.visitedStateCount != row.exactStateCount)
            throw std::runtime_error("low components: exact state count mismatch");
        return row;
    }

    void writeHex(std::ostream& output, std::uint64_t value)
    {
        output << "\"0x" << std::hex << value << std::dec << "\"";
    }
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
            throw std::runtime_error("usage: certificate.exe OUTPUT.json");
        const auto started = std::chrono::steady_clock::now();
        const ByteMap parityMap(systematicRightColumns());
        std::array<std::uint64_t, 64> stepColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            stepColumns[bit] = parityMap(accumulate(std::uint64_t{1} << bit));
        const ByteMap step(stepColumns);

        std::array<std::vector<std::uint64_t>, 7> componentBases;
        for (unsigned component = 0; component < ComponentDegrees.size(); ++component)
            componentBases[component] = kernelBasis(
                step, ComponentFactors[component], ComponentDegrees[component]);
        std::vector<std::uint64_t> fullComponentBasis;
        for (const auto& componentBasis : componentBases)
            fullComponentBasis.insert(
                fullComponentBasis.end(), componentBasis.begin(), componentBasis.end());
        if (fullComponentBasis.size() != 64)
            throw std::runtime_error("low components: component degrees do not sum to 64");
        const CoordinateSolver fullComponentSolver(fullComponentBasis);
        for (const std::uint64_t basisState : fullComponentBasis)
            (void)fullComponentSolver.coordinates(step(basisState));
        const std::vector<std::uint64_t> stepColumnVector(
            stepColumns.begin(), stepColumns.end());
        const CoordinateSolver invertibilityCheck(stepColumnVector);
        (void)invertibilityCheck;

        std::vector<Row> rows;
        std::uint64_t totalExactStates = 0;
        std::uint64_t totalLowStates = 0;
        std::uint64_t globalMinimum = std::numeric_limits<std::uint64_t>::max();
        std::uint64_t globalMaximumAbsoluteBias = 0;
        for (unsigned supportMask = 1; supportMask < (1u << ComponentDegrees.size()); ++supportMask)
        {
            unsigned dimension = 0;
            for (unsigned component = 0; component < ComponentDegrees.size(); ++component)
            {
                if ((supportMask >> component) & 1)
                    dimension += ComponentDegrees[component];
            }
            if (dimension > ComponentDimensionCap)
                continue;
            Row row = enumerateSupport(supportMask, step, componentBases);
            totalExactStates += row.exactStateCount;
            totalLowStates += row.lowStateCount;
            globalMinimum = std::min(globalMinimum, row.minimumResponseWeight);
            globalMaximumAbsoluteBias = std::max(
                globalMaximumAbsoluteBias, row.maximumAbsoluteResponseBias);
            std::cerr << "support=0x" << std::hex << supportMask << std::dec
                      << " dimension=" << dimension
                      << " states=" << row.exactStateCount
                      << " cycles=" << row.cycleCount
                      << " minimum=" << row.minimumResponseWeight
                      << " max_abs_bias=" << row.maximumAbsoluteResponseBias
                      << " low=" << row.lowStateCount << '\n';
            rows.push_back(row);
        }
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();

        std::ofstream output(argv[1]);
        if (!output)
            throw std::runtime_error("low components: cannot open output receipt");
        output << "{\n"
               << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-g4-goal02-low-components-v1\",\n"
               << "  \"candidate\": \"Riffle PacketMul-WrapMul-2Lap g=4\",\n"
               << "  \"evidence_label\": \"EXACT_EXHAUSTIVE_COMPONENT_SPECTRUM_CERTIFICATE\",\n"
               << "  \"autonomous_stored_state_map\": \"T=P_after_systematic_BCH o Acc\",\n"
               << "  \"emission_map\": \"Acc\",\n"
               << "  \"response_nodes\": " << ResponseNodes << ",\n"
               << "  \"low_weight_threshold\": " << LowWeightThreshold << ",\n"
               << "  \"component_dimension_cap\": " << ComponentDimensionCap << ",\n"
               << "  \"autonomous_map_rank\": 64,\n"
               << "  \"direct_sum_component_dimension\": 64,\n"
               << "  \"component_degrees\": [1, 2, 4, 9, 10, 18, 20],\n"
               << "  \"component_factors_hex\": [\"0x3\", \"0x7\", \"0x13\", \"0x373\", \"0x519\", \"0x7c9c3\", \"0x1e1fff\"],\n"
               << "  \"support_count\": " << rows.size() << ",\n"
               << "  \"total_exact_nonzero_states\": " << totalExactStates << ",\n"
               << "  \"total_low_response_states\": " << totalLowStates << ",\n"
               << "  \"global_minimum_response_weight\": " << globalMinimum << ",\n"
               << "  \"global_maximum_absolute_response_bias\": " << globalMaximumAbsoluteBias << ",\n"
               << "  \"rows\": [\n";
        for (unsigned rowIndex = 0; rowIndex < rows.size(); ++rowIndex)
        {
            const Row& row = rows[rowIndex];
            output << "    {\n"
                   << "      \"component_support_mask_hex\": \"0x" << std::hex << row.supportMask << std::dec << "\",\n"
                   << "      \"component_degrees\": [";
            bool first = true;
            for (unsigned component = 0; component < ComponentDegrees.size(); ++component)
            {
                if ((row.supportMask >> component) & 1)
                {
                    if (!first)
                        output << ", ";
                    output << ComponentDegrees[component];
                    first = false;
                }
            }
            output << "],\n"
                   << "      \"dimension\": " << row.dimension << ",\n"
                   << "      \"exact_nonzero_state_count\": " << row.exactStateCount << ",\n"
                   << "      \"visited_state_count\": " << row.visitedStateCount << ",\n"
                   << "      \"cycle_count\": " << row.cycleCount << ",\n"
                   << "      \"low_response_state_count\": " << row.lowStateCount << ",\n"
                   << "      \"minimum_response_weight\": " << row.minimumResponseWeight << ",\n"
                   << "      \"minimum_witness_coordinates_hex\": \"0x" << std::hex << row.minimumCoordinates << std::dec << "\",\n"
                   << "      \"minimum_witness_state_hex\": ";
            writeHex(output, row.minimumState);
            output << ",\n"
                   << "      \"minimum_witness_cycle_period\": " << row.minimumCyclePeriod << "\n"
                   << "      ,\"maximum_absolute_response_bias\": " << row.maximumAbsoluteResponseBias << ",\n"
                   << "      \"maximum_bias_response_weight\": " << row.maximumBiasResponseWeight << ",\n"
                   << "      \"maximum_bias_witness_coordinates_hex\": \"0x" << std::hex << row.maximumBiasCoordinates << std::dec << "\",\n"
                   << "      \"maximum_bias_witness_state_hex\": ";
            writeHex(output, row.maximumBiasState);
            output << "\n"
                   << "    }" << (rowIndex + 1 == rows.size() ? "\n" : ",\n");
        }
        output << "  ],\n"
               << "  \"elapsed_seconds\": " << std::fixed << std::setprecision(6) << elapsed << ",\n"
               << "  \"result\": \"" << (totalLowStates == 0 ? "PASS" : "LOW_RESPONSE_STATES_FOUND") << "\"\n"
               << "}\n";
        if (!output)
            throw std::runtime_error("low components: receipt write failed");
        std::cerr << "output=" << argv[1] << '\n'
                  << "total_states=" << totalExactStates << '\n'
                  << "total_low_states=" << totalLowStates << '\n'
                  << "global_minimum=" << globalMinimum << '\n'
                  << "global_maximum_absolute_bias=" << globalMaximumAbsoluteBias << '\n'
                  << "elapsed_seconds=" << elapsed << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
