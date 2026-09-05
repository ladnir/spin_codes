#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned StateBits = 64;
    constexpr unsigned Nibbles = 16;
    constexpr unsigned RequiredBranch = 12;

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
        for (unsigned bit = 1; bit < StateBits; ++bit)
        {
            if ((carrylessMultiplyLow(Generator, inverse) >> bit) & 1)
                inverse |= std::uint64_t{1} << bit;
        }
        if (carrylessMultiplyLow(Generator, inverse) != 1)
            throw std::runtime_error("pair branch: generator inverse failed");
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
            const auto bit = std::countr_zero(message);
            result.low ^= Generator << bit;
            result.high ^=
                (bit == 0 ? 0 : Generator >> (StateBits - bit)) |
                (std::uint64_t{1} << 63);
            message &= message - 1;
        }
        return result;
    }

    std::array<std::uint64_t, StateBits> systematicStateColumns()
    {
        const auto inverse = generatorInverse();
        std::array<std::uint64_t, StateBits> columns{};
        for (unsigned column = 0; column < StateBits; ++column)
        {
            const auto message = carrylessMultiplyLow(
                std::uint64_t{1} << column,
                inverse);
            const auto word = encodeBch(message);
            if (word.low != (std::uint64_t{1} << column))
                throw std::runtime_error("pair branch: systematic BCH failed");
            columns[column] = word.high;
        }
        return columns;
    }

    std::uint64_t applyColumns(
        const std::array<std::uint64_t, StateBits>& columns,
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

    std::array<std::uint64_t, StateBits> transposeColumns(
        const std::array<std::uint64_t, StateBits>& columns)
    {
        std::array<std::uint64_t, StateBits> transpose{};
        for (unsigned output = 0; output < StateBits; ++output)
        {
            for (unsigned input = 0; input < StateBits; ++input)
                transpose[output] |= ((columns[input] >> output) & 1) << input;
        }
        return transpose;
    }

    std::array<std::uint64_t, StateBits> transposeStepColumns()
    {
        const auto stateColumns = systematicStateColumns();
        std::array<std::uint64_t, StateBits> stepColumns{};
        for (unsigned input = 0; input < StateBits; ++input)
        {
            stepColumns[input] = applyColumns(
                stateColumns,
                accumulate(std::uint64_t{1} << input));
        }
        return transposeColumns(stepColumns);
    }

    std::array<std::uint64_t, StateBits> inverseColumns(
        const std::array<std::uint64_t, StateBits>& columns)
    {
        struct AugmentedRow
        {
            std::uint64_t left = 0;
            std::uint64_t right = 0;
        };
        std::array<AugmentedRow, StateBits> rows{};
        for (unsigned output = 0; output < StateBits; ++output)
        {
            for (unsigned input = 0; input < StateBits; ++input)
                rows[output].left |= ((columns[input] >> output) & 1) << input;
            rows[output].right = std::uint64_t{1} << output;
        }
        for (unsigned column = 0; column < StateBits; ++column)
        {
            unsigned pivot = column;
            while (pivot < StateBits && !((rows[pivot].left >> column) & 1))
                ++pivot;
            if (pivot == StateBits)
                throw std::runtime_error("pair branch: singular transpose step");
            std::swap(rows[column], rows[pivot]);
            for (unsigned row = 0; row < StateBits; ++row)
            {
                if (row != column && ((rows[row].left >> column) & 1))
                {
                    rows[row].left ^= rows[column].left;
                    rows[row].right ^= rows[column].right;
                }
            }
        }
        std::array<std::uint64_t, StateBits> inverse{};
        for (unsigned output = 0; output < StateBits; ++output)
        {
            for (unsigned input = 0; input < StateBits; ++input)
            {
                inverse[output] |=
                    ((rows[input].right >> output) & 1)
                    << input;
            }
        }
        for (unsigned bit = 0; bit < StateBits; ++bit)
        {
            if (applyColumns(columns, inverse[bit]) != (std::uint64_t{1} << bit))
                throw std::runtime_error("pair branch: inverse verification failed");
        }
        return inverse;
    }

    unsigned nibbleWeight(std::uint64_t value)
    {
        unsigned result = 0;
        for (unsigned nibble = 0; nibble < Nibbles; ++nibble)
            result += ((value >> (4 * nibble)) & 0xf) != 0;
        return result;
    }

    std::uint64_t expandNibbleMask(std::uint16_t mask)
    {
        std::uint64_t result = 0;
        while (mask)
        {
            const auto nibble = std::countr_zero(mask);
            result |= std::uint64_t{0xf} << (4 * nibble);
            mask &= mask - 1;
        }
        return result;
    }

    struct Dependency
    {
        bool dependent = false;
        std::uint32_t representation = 0;
    };

    Dependency dependency(
        const std::array<std::uint64_t, 20>& columns,
        unsigned columnCount,
        std::uint64_t keepMask)
    {
        static std::array<std::uint64_t, StateBits> basis{};
        static std::array<std::uint32_t, StateBits> representations{};
        std::array<unsigned, 20> usedPivots{};
        unsigned used = 0;
        const auto clear = [&]()
        {
            for (unsigned index = 0; index < used; ++index)
            {
                basis[usedPivots[index]] = 0;
                representations[usedPivots[index]] = 0;
            }
        };
        for (unsigned column = 0; column < columnCount; ++column)
        {
            auto value = columns[column] & keepMask;
            std::uint32_t representation = std::uint32_t{1} << column;
            while (value)
            {
                const auto pivot = std::bit_width(value) - 1;
                if (basis[pivot])
                {
                    value ^= basis[pivot];
                    representation ^= representations[pivot];
                }
                else
                {
                    basis[pivot] = value;
                    representations[pivot] = representation;
                    usedPivots[used++] = pivot;
                    break;
                }
            }
            if (!value)
            {
                clear();
                return {true, representation};
            }
        }
        clear();
        return {};
    }

    std::vector<std::uint16_t> masksOfWeight(unsigned weight)
    {
        std::vector<std::uint16_t> result;
        for (unsigned mask = 0; mask < (1u << Nibbles); ++mask)
        {
            if (static_cast<unsigned>(std::popcount(mask)) == weight)
                result.push_back(static_cast<std::uint16_t>(mask));
        }
        return result;
    }

    struct DirectionResult
    {
        std::uint64_t rankChecks = 0;
        bool passed = true;
        std::uint64_t witnessInput = 0;
        std::uint64_t witnessOutput = 0;
    };

    DirectionResult certifyDirection(
        const std::array<std::uint64_t, StateBits>& columns,
        const std::string& name)
    {
        DirectionResult result{};
        for (unsigned inputWeight = 1; inputWeight <= (RequiredBranch - 1) / 2; ++inputWeight)
        {
            const unsigned allowedOutputWeight = RequiredBranch - 1 - inputWeight;
            const auto inputMasks = masksOfWeight(inputWeight);
            const auto outputMasks = masksOfWeight(allowedOutputWeight);
            const auto expected =
                static_cast<std::uint64_t>(inputMasks.size()) * outputMasks.size();
            std::uint64_t trancheChecks = 0;
            for (const auto inputMask : inputMasks)
            {
                std::array<unsigned, 20> inputBits{};
                std::array<std::uint64_t, 20> selectedColumns{};
                unsigned columnCount = 0;
                for (unsigned nibble = 0; nibble < Nibbles; ++nibble)
                {
                    if (!((inputMask >> nibble) & 1))
                        continue;
                    for (unsigned bit = 0; bit < 4; ++bit)
                    {
                        inputBits[columnCount] = 4 * nibble + bit;
                        selectedColumns[columnCount] = columns[inputBits[columnCount]];
                        ++columnCount;
                    }
                }
                for (const auto outputMask : outputMasks)
                {
                    const auto keepMask = ~expandNibbleMask(outputMask);
                    const auto found = dependency(
                        selectedColumns,
                        columnCount,
                        keepMask);
                    ++trancheChecks;
                    if (found.dependent)
                    {
                        std::uint64_t input = 0;
                        for (unsigned column = 0; column < columnCount; ++column)
                        {
                            if ((found.representation >> column) & 1)
                                input |= std::uint64_t{1} << inputBits[column];
                        }
                        const auto output = applyColumns(columns, input);
                        if (!input || nibbleWeight(input) > inputWeight ||
                            nibbleWeight(output) > allowedOutputWeight)
                            throw std::runtime_error("pair branch: invalid rank witness");
                        result.rankChecks += trancheChecks;
                        result.passed = false;
                        result.witnessInput = input;
                        result.witnessOutput = output;
                        return result;
                    }
                }
            }
            if (trancheChecks != expected)
                throw std::runtime_error("pair branch: rank-check count mismatch");
            result.rankChecks += trancheChecks;
            std::cout
                << "direction=" << name
                << " input_weight=" << inputWeight
                << " output_support_size=" << allowedOutputWeight
                << " rank_checks=" << trancheChecks
                << std::endl;
        }
        return result;
    }

    void writeResult(
        const std::string& outputPath,
        const DirectionResult& forward,
        const DirectionResult& inverse)
    {
        const bool passed = forward.passed && inverse.passed;
        std::ofstream output(outputPath);
        if (!output)
            throw std::runtime_error("pair branch: cannot open output");
        output << "{\n";
        output << "  \"schema\": \"riffle-packetmul-2lap-g4-goal02-pair-branch-primary-v1\",\n";
        output << "  \"candidate\": \"Riffle PacketMul-2Lap g=4\",\n";
        output << "  \"candidate_id\": \"riffle_packetmul_2lap_g4\",\n";
        output << "  \"evidence_label\": \""
               << (passed ? "EXACT_EXHAUSTIVE_PRIMARY" : "EXACT_COUNTEREXAMPLE")
               << "\",\n";
        output << "  \"required_branch_number\": " << RequiredBranch << ",\n";
        output << "  \"forward_rank_checks\": " << forward.rankChecks << ",\n";
        output << "  \"inverse_rank_checks\": " << inverse.rankChecks << ",\n";
        output << "  \"all_rank_checks_full\": " << (passed ? "true" : "false") << ",\n";
        output << "  \"forward_witness_input_hex\": \"0x" << std::hex << forward.witnessInput << "\",\n";
        output << "  \"forward_witness_output_hex\": \"0x" << forward.witnessOutput << "\",\n";
        output << std::dec;
        output << "  \"forward_witness_input_nibble_weight\": " << nibbleWeight(forward.witnessInput) << ",\n";
        output << "  \"forward_witness_output_nibble_weight\": " << nibbleWeight(forward.witnessOutput) << ",\n";
        output << std::hex;
        output << "  \"inverse_witness_input_hex\": \"0x" << inverse.witnessInput << "\",\n";
        output << "  \"inverse_witness_output_hex\": \"0x" << inverse.witnessOutput << "\",\n";
        output << std::dec;
        output << "  \"proof_rule\": \"A pair of total nibble weight at most 11 has one side of weight at most 5. For each exact support I of that side and each support O of size 11-|I| on the other side, the restricted columns outside O must lose rank. The forward enumeration covers the case wt(x)<=5. The inverse enumeration covers wt(Sx)<=5. Full rank in every case excludes every nonzero pair of total weight below 12.\",\n";
        output << "  \"pair_branch_lower_bound\": " << (passed ? 12 : 0) << ",\n";
        output << "  \"orbit_consequence\": \"Because 32772 is even and S is invertible, a branch lower bound of 12 implies sum_t wt_4(S^t chi) >= 6*32772 for every nonzero chi.\",\n";
        output << "  \"status\": \"" << (passed ? "EXACT_PAIR_BRANCH_CERTIFICATE" : "EXACT_LOCAL_COUNTEREXAMPLE") << "\"\n";
        output << "}\n";
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::string outputPath = argc > 1
            ? argv[1]
            : "constructions/riffle_packetmul_2lap_g4/receipts/goal02_pair_branch_primary.json";
        const auto forwardColumns = transposeStepColumns();
        const auto backwardColumns = inverseColumns(forwardColumns);
        const auto forward = certifyDirection(forwardColumns, "forward");
        DirectionResult inverse{};
        if (forward.passed)
            inverse = certifyDirection(backwardColumns, "inverse");
        writeResult(outputPath, forward, inverse);
        std::cout << "output=" << outputPath << std::endl;
        std::cout << "status="
                  << (forward.passed && inverse.passed
                      ? "EXACT_PAIR_BRANCH_CERTIFICATE"
                      : "EXACT_LOCAL_COUNTEREXAMPLE")
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
