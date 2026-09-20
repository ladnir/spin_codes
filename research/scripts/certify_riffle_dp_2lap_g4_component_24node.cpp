#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
#ifndef COMPONENT_DIMENSION
#error COMPONENT_DIMENSION must be defined at compile time
#endif
#ifndef CERTIFICATE_WINDOW_NODES
#define CERTIFICATE_WINDOW_NODES 24
#endif
#ifndef INDEPENDENT_VERIFIER
#define INDEPENDENT_VERIFIER 0
#endif
    constexpr std::array<unsigned, 7> ComponentDegrees{1, 2, 4, 9, 10, 18, 20};
    constexpr std::array<std::uint64_t, 7> ComponentFactors{
        0x3ULL, 0x7ULL, 0x13ULL, 0x373ULL, 0x519ULL, 0x7C9C3ULL, 0x1E1FFFULL};
    constexpr unsigned Dimension = COMPONENT_DIMENSION;
    constexpr unsigned WindowNodes = CERTIFICATE_WINDOW_NODES;
    static_assert(WindowNodes == 12 || WindowNodes == 24);
    constexpr unsigned Slots = 16;
    constexpr unsigned Length = WindowNodes * Slots;
    constexpr unsigned InformationSetCount = Length / Dimension;
    constexpr unsigned LeftDimension = Dimension / 2;
    constexpr unsigned RightDimension = Dimension - LeftDimension;
    constexpr std::uint64_t LeftSize = std::uint64_t{1} << LeftDimension;
    constexpr std::uint64_t RightSize = std::uint64_t{1} << RightDimension;

    struct Word384
    {
        std::uint64_t limb0 = 0;
        std::uint64_t limb1 = 0;
        std::uint64_t limb2 = 0;
        std::uint64_t limb3 = 0;
        std::uint64_t limb4 = 0;
        std::uint64_t limb5 = 0;
    };

    Word384 operator^(Word384 left, const Word384& right)
    {
        left.limb0 ^= right.limb0;
        left.limb1 ^= right.limb1;
        left.limb2 ^= right.limb2;
        left.limb3 ^= right.limb3;
        left.limb4 ^= right.limb4;
        left.limb5 ^= right.limb5;
        return left;
    }

    Word384& operator^=(Word384& left, const Word384& right)
    {
        left.limb0 ^= right.limb0;
        left.limb1 ^= right.limb1;
        left.limb2 ^= right.limb2;
        left.limb3 ^= right.limb3;
        left.limb4 ^= right.limb4;
        left.limb5 ^= right.limb5;
        return left;
    }

    unsigned weight(const Word384& word)
    {
        return
            std::popcount(word.limb0) +
            std::popcount(word.limb1) +
            std::popcount(word.limb2) +
            std::popcount(word.limb3) +
            std::popcount(word.limb4) +
            std::popcount(word.limb5);
    }

    bool coordinate(const Word384& word, unsigned index)
    {
        if (index < 64)
            return (word.limb0 >> index) & 1;
        if (index < 128)
            return (word.limb1 >> (index - 64)) & 1;
        if (index < 192)
            return (word.limb2 >> (index - 128)) & 1;
        if (index < 256)
            return (word.limb3 >> (index - 192)) & 1;
        if (index < 320)
            return (word.limb4 >> (index - 256)) & 1;
        return (word.limb5 >> (index - 320)) & 1;
    }

    void setCoordinate(Word384& word, unsigned index)
    {
        if (index < 64)
            word.limb0 |= std::uint64_t{1} << index;
        else if (index < 128)
            word.limb1 |= std::uint64_t{1} << (index - 64);
        else if (index < 192)
            word.limb2 |= std::uint64_t{1} << (index - 128);
        else if (index < 256)
            word.limb3 |= std::uint64_t{1} << (index - 192);
        else if (index < 320)
            word.limb4 |= std::uint64_t{1} << (index - 256);
        else
            word.limb5 |= std::uint64_t{1} << (index - 320);
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

    std::uint64_t polynomialMultiply(std::uint64_t left, std::uint64_t right)
    {
        return carrylessMultiplyLow(left, right);
    }

    unsigned supportDimension(unsigned supportMask)
    {
        unsigned result = 0;
        for (unsigned index = 0; index < ComponentDegrees.size(); ++index)
        {
            if ((supportMask >> index) & 1)
                result += ComponentDegrees[index];
        }
        return result;
    }

    std::uint64_t supportAnnihilator(unsigned supportMask)
    {
        std::uint64_t result = 1;
        for (unsigned index = 0; index < ComponentFactors.size(); ++index)
        {
            if ((supportMask >> index) & 1)
                result = polynomialMultiply(result, ComponentFactors[index]);
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
            throw std::runtime_error("component 24-node certificate: generator inverse failed");
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
            const auto word = encodeBch(message);
            if (word.low != (std::uint64_t{1} << column))
                throw std::runtime_error("component 24-node certificate: systematic BCH failed");
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

    std::array<std::uint64_t, 64> transposeColumns(
        const std::array<std::uint64_t, 64>& columns)
    {
        std::array<std::uint64_t, 64> transpose{};
        for (unsigned output = 0; output < 64; ++output)
        {
            for (unsigned input = 0; input < 64; ++input)
                transpose[output] |= ((columns[input] >> output) & 1) << input;
        }
        return transpose;
    }

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

    std::array<std::uint64_t, Dimension> componentBasis(
        const ByteMap& innerStep,
        std::uint64_t annihilator)
    {
        std::array<std::uint64_t, 64> stepColumns{};
        for (unsigned bit = 0; bit < 64; ++bit)
            stepColumns[bit] = innerStep(accumulate(std::uint64_t{1} << bit));
        const ByteMap transposeStep(transposeColumns(stepColumns));

        std::array<std::uint64_t, 64> pivotValues{};
        std::array<std::uint64_t, 64> pivotRepresentations{};
        std::array<std::uint64_t, Dimension> basis{};
        unsigned basisSize = 0;
        for (unsigned column = 0; column < 64; ++column)
        {
            auto value = applyPolynomial(
                transposeStep,
                annihilator,
                std::uint64_t{1} << column);
            std::uint64_t representation = std::uint64_t{1} << column;
            while (value)
            {
                const auto pivot = std::bit_width(value) - 1;
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
            {
                if (basisSize >= Dimension)
                    throw std::runtime_error("component 24-node certificate: oversized kernel");
                basis[basisSize++] = representation;
            }
        }
        if (basisSize != Dimension)
            throw std::runtime_error("component 24-node certificate: kernel dimension changed");
        return basis;
    }

    std::array<Word384, Dimension> codeRows(
        unsigned packetValue,
        const ByteMap& innerStep,
        const std::array<std::uint64_t, Dimension>& characters)
    {
        std::array<std::uint64_t, Length> observations{};
        unsigned coordinateIndex = 0;
        for (unsigned slot = 0; slot < Slots; ++slot)
        {
            std::uint64_t state = std::uint64_t{packetValue} << (4 * slot);
            for (unsigned node = 0; node < WindowNodes; ++node)
            {
                state = innerStep(accumulate(state));
                observations[coordinateIndex++] = state;
            }
        }
        std::array<Word384, Dimension> rows{};
        for (unsigned row = 0; row < Dimension; ++row)
        {
            for (unsigned coordinateIndex2 = 0; coordinateIndex2 < Length; ++coordinateIndex2)
            {
                if (std::popcount(characters[row] & observations[coordinateIndex2]) & 1)
                    setCoordinate(rows[row], coordinateIndex2);
            }
        }
        return rows;
    }

    std::array<std::uint64_t, Length> coordinateColumns(
        const std::array<Word384, Dimension>& rows)
    {
        std::array<std::uint64_t, Length> columns{};
        for (unsigned coordinateIndex = 0; coordinateIndex < Length; ++coordinateIndex)
        {
            for (unsigned row = 0; row < Dimension; ++row)
                columns[coordinateIndex] |= std::uint64_t{coordinate(rows[row], coordinateIndex)} << row;
        }
        return columns;
    }

    bool insertIndependent(std::array<std::uint64_t, Dimension>& pivots, std::uint64_t value)
    {
        while (value)
        {
            const auto pivot = std::bit_width(value) - 1;
            if (pivots[pivot])
                value ^= pivots[pivot];
            else
            {
                pivots[pivot] = value;
                return true;
            }
        }
        return false;
    }

    using InformationSet = std::array<unsigned, Dimension>;

    struct InformationPartition
    {
        std::array<InformationSet, InformationSetCount> sets{};
        std::array<unsigned, Length - InformationSetCount * Dimension> unused{};
        unsigned restarts = 0;
    };

    InformationPartition findInformationSets(
        unsigned supportMask,
        unsigned packetValue,
        const std::array<std::uint64_t, Length>& columns)
    {
        std::array<unsigned, Length> order{};
        std::iota(order.begin(), order.end(), 0);
        if constexpr (INDEPENDENT_VERIFIER)
            std::reverse(order.begin(), order.end());
        std::mt19937_64 random(
            (INDEPENDENT_VERIFIER ? 0x12A11D17BEEFULL : 0x24C0DEC0FFEEULL) ^
            (std::uint64_t{WindowNodes} << 48) ^
            (std::uint64_t{supportMask} << 16) ^ packetValue);
        for (unsigned attempt = 0; attempt < 10'000; ++attempt)
        {
            if (attempt != 0)
                std::shuffle(order.begin(), order.end(), random);
            std::array<unsigned, Length> available = order;
            unsigned availableSize = Length;
            InformationPartition result{};
            bool failed = false;
            for (unsigned setIndex = 0; setIndex < InformationSetCount; ++setIndex)
            {
                std::array<std::uint64_t, Dimension> pivots{};
                std::array<unsigned, Length> remaining{};
                unsigned selectedSize = 0;
                unsigned remainingSize = 0;
                for (unsigned index = 0; index < availableSize; ++index)
                {
                    const auto coordinateIndex = available[index];
                    if (insertIndependent(pivots, columns[coordinateIndex]))
                        result.sets[setIndex][selectedSize++] = coordinateIndex;
                    else
                        remaining[remainingSize++] = coordinateIndex;
                }
                if (selectedSize != Dimension)
                {
                    failed = true;
                    break;
                }
                available = remaining;
                availableSize = remainingSize;
            }
            if (!failed)
            {
                if (availableSize != result.unused.size())
                    throw std::runtime_error("component 24-node certificate: unused coordinate count");
                std::copy_n(available.begin(), availableSize, result.unused.begin());
                result.restarts = attempt;
                return result;
            }
        }
        throw std::runtime_error("component 24-node certificate: information sets not found");
    }

    struct SystematicRow
    {
        Word384 word{};
        std::uint64_t character = 0;
    };

    std::array<SystematicRow, Dimension> systematicRows(
        const std::array<Word384, Dimension>& rows,
        const std::array<std::uint64_t, Dimension>& characters,
        const InformationSet& informationSet)
    {
        std::array<SystematicRow, Dimension> systematic{};
        for (unsigned row = 0; row < Dimension; ++row)
            systematic[row] = SystematicRow{rows[row], characters[row]};
        for (unsigned column = 0; column < Dimension; ++column)
        {
            unsigned pivot = column;
            while (pivot < Dimension && !coordinate(systematic[pivot].word, informationSet[column]))
                ++pivot;
            if (pivot == Dimension)
                throw std::runtime_error("component 24-node certificate: singular information set");
            std::swap(systematic[column], systematic[pivot]);
            for (unsigned row = 0; row < Dimension; ++row)
            {
                if (row != column && coordinate(systematic[row].word, informationSet[column]))
                {
                    systematic[row].word ^= systematic[column].word;
                    systematic[row].character ^= systematic[column].character;
                }
            }
        }
        for (unsigned row = 0; row < Dimension; ++row)
        {
            for (unsigned column = 0; column < Dimension; ++column)
            {
                if (coordinate(systematic[row].word, informationSet[column]) != (row == column))
                    throw std::runtime_error("component 24-node certificate: systematic verification failed");
            }
        }
        return systematic;
    }

    struct TableEntry
    {
        Word384 word{};
        std::uint64_t character = 0;
    };

    std::vector<TableEntry> buildTable(
        const std::array<SystematicRow, Dimension>& rows,
        unsigned offset,
        unsigned halfDimension)
    {
        const auto halfSize = std::uint64_t{1} << halfDimension;
        std::vector<TableEntry> table(halfSize);
        for (std::uint64_t mask = 1; mask < halfSize; ++mask)
        {
            const auto bitValue = mask & (~mask + 1);
            const auto bit = std::countr_zero(bitValue);
            table[mask] = table[mask ^ bitValue];
            table[mask].word ^= rows[offset + bit].word;
            table[mask].character ^= rows[offset + bit].character;
        }
        return table;
    }

    template<unsigned Bits>
    std::array<std::vector<std::uint32_t>, Bits + 1> masksByWeight()
    {
        std::array<std::vector<std::uint32_t>, Bits + 1> result{};
        const auto size = std::uint32_t{1} << Bits;
        for (std::uint32_t mask = 0; mask < size; ++mask)
            result[std::popcount(mask)].push_back(mask);
        return result;
    }

    std::uint64_t choose(unsigned n, unsigned k)
    {
        if (k > n)
            return 0;
        k = std::min(k, n - k);
        std::uint64_t result = 1;
        for (unsigned i = 1; i <= k; ++i)
            result = result * (n - k + i) / i;
        return result;
    }

    std::uint64_t ballSize(unsigned radius)
    {
        std::uint64_t result = 0;
        for (unsigned weight2 = 0; weight2 <= radius; ++weight2)
            result += choose(Dimension, weight2);
        return result;
    }

    struct SearchResult
    {
        bool passed = true;
        unsigned minimumLow = Length;
        unsigned minimumComplement = Length;
        std::uint64_t violatingCharacter = 0;
        unsigned violatingDistance = Length;
        const char* violatingSide = "none";
        std::uint64_t candidates = 0;
    };

    SearchResult searchInformationSet(
        const std::array<SystematicRow, Dimension>& rows,
        unsigned radius,
        unsigned requiredDistance,
        const std::array<std::vector<std::uint32_t>, LeftDimension + 1>& leftMasksByWeight,
        const std::array<std::vector<std::uint32_t>, RightDimension + 1>& rightMasksByWeight)
    {
        const auto left = buildTable(rows, 0, LeftDimension);
        const auto right = buildTable(rows, LeftDimension, RightDimension);
        TableEntry allOnes{};
        for (const auto& row : rows)
        {
            allOnes.word ^= row.word;
            allOnes.character ^= row.character;
        }

        SearchResult result{};
        for (unsigned leftWeight = 0; leftWeight <= radius; ++leftWeight)
        {
            for (unsigned rightWeight = 0; rightWeight + leftWeight <= radius; ++rightWeight)
            {
                if (leftWeight > LeftDimension || rightWeight > RightDimension)
                    continue;
                const auto& leftMasks = leftMasksByWeight[leftWeight];
                const auto& rightMasks = rightMasksByWeight[rightWeight];
                for (const auto leftMask : leftMasks)
                {
                    const auto& leftEntry = left[leftMask];
                    for (const auto rightMask : rightMasks)
                    {
                        const auto& rightEntry = right[rightMask];
                        const Word384 word = leftEntry.word ^ rightEntry.word;
                        const auto character = leftEntry.character ^ rightEntry.character;
                        ++result.candidates;
                        if (character != 0)
                        {
                            const auto lowWeight = weight(word);
                            result.minimumLow = std::min(result.minimumLow, lowWeight);
                            if (lowWeight < requiredDistance)
                            {
                                result.passed = false;
                                result.violatingCharacter = character;
                                result.violatingDistance = lowWeight;
                                result.violatingSide = "weight";
                                return result;
                            }
                        }
                        const Word384 highWord = allOnes.word ^ word;
                        const auto complementWeight = Length - weight(highWord);
                        result.minimumComplement = std::min(
                            result.minimumComplement,
                            complementWeight);
                        if (complementWeight < requiredDistance)
                        {
                            result.passed = false;
                            result.violatingCharacter = allOnes.character ^ character;
                            result.violatingDistance = complementWeight;
                            result.violatingSide = "complement_weight";
                            return result;
                        }
                    }
                }
            }
        }
        if (result.candidates != ballSize(radius))
            throw std::runtime_error("component 24-node certificate: Hamming ball count mismatch");
        return result;
    }

    void inspectCandidateIndependent(
        const Word384& word,
        std::uint64_t character,
        const Word384& allOnesWord,
        std::uint64_t allOnesCharacter,
        unsigned requiredDistance,
        SearchResult& result)
    {
        ++result.candidates;
        if (character != 0)
        {
            const auto lowWeight = weight(word);
            result.minimumLow = std::min(result.minimumLow, lowWeight);
            if (lowWeight < requiredDistance)
            {
                result.passed = false;
                result.violatingCharacter = character;
                result.violatingDistance = lowWeight;
                result.violatingSide = "weight";
                return;
            }
        }
        const auto highCharacter = allOnesCharacter ^ character;
        if (highCharacter != 0)
        {
            const auto complementWeight = Length - weight(allOnesWord ^ word);
            result.minimumComplement = std::min(
                result.minimumComplement,
                complementWeight);
            if (complementWeight < requiredDistance)
            {
                result.passed = false;
                result.violatingCharacter = highCharacter;
                result.violatingDistance = complementWeight;
                result.violatingSide = "complement_weight";
            }
        }
    }

    void enumerateExactWeightIndependent(
        const std::array<SystematicRow, Dimension>& rows,
        unsigned nextRow,
        unsigned remaining,
        Word384 word,
        std::uint64_t character,
        const Word384& allOnesWord,
        std::uint64_t allOnesCharacter,
        unsigned requiredDistance,
        SearchResult& result)
    {
        if (!result.passed)
            return;
        if (remaining == 0)
        {
            inspectCandidateIndependent(
                word,
                character,
                allOnesWord,
                allOnesCharacter,
                requiredDistance,
                result);
            return;
        }
        const unsigned last = Dimension - remaining;
        for (unsigned row = nextRow; row <= last; ++row)
        {
            enumerateExactWeightIndependent(
                rows,
                row + 1,
                remaining - 1,
                word ^ rows[row].word,
                character ^ rows[row].character,
                allOnesWord,
                allOnesCharacter,
                requiredDistance,
                result);
            if (!result.passed)
                return;
        }
    }

    SearchResult searchInformationSetIndependent(
        const std::array<SystematicRow, Dimension>& rows,
        unsigned radius,
        unsigned requiredDistance)
    {
        Word384 allOnesWord{};
        std::uint64_t allOnesCharacter = 0;
        for (const auto& row : rows)
        {
            allOnesWord ^= row.word;
            allOnesCharacter ^= row.character;
        }
        SearchResult result{};
        for (unsigned exactWeight = 0; exactWeight <= radius; ++exactWeight)
        {
            enumerateExactWeightIndependent(
                rows,
                0,
                exactWeight,
                Word384{},
                0,
                allOnesWord,
                allOnesCharacter,
                requiredDistance,
                result);
            if (!result.passed)
                return result;
        }
        if (result.candidates != ballSize(radius))
            throw std::runtime_error("component certificate: independent Hamming ball count mismatch");
        return result;
    }

    void writeInformationSet(std::ostream& output, const InformationSet& set)
    {
        output << '[';
        for (unsigned index = 0; index < Dimension; ++index)
        {
            if (index)
                output << ',';
            output << set[index];
        }
        output << ']';
    }

    void writeReceipt(
        const std::string& path,
        unsigned supportMask,
        std::uint64_t annihilator,
        unsigned packetValue,
        unsigned requiredDistance,
        unsigned radius,
        const InformationPartition& partition,
        const std::array<SearchResult, InformationSetCount>& results,
        double elapsedSeconds)
    {
        bool passed = true;
        std::uint64_t totalCandidates = 0;
        unsigned minimumLow = Length;
        unsigned minimumComplement = Length;
        for (const auto& result : results)
        {
            passed = passed && result.passed;
            totalCandidates += result.candidates;
            minimumLow = std::min(minimumLow, result.minimumLow);
            minimumComplement = std::min(minimumComplement, result.minimumComplement);
        }
        std::ofstream output(path);
        if (!output)
            throw std::runtime_error("component 24-node certificate: cannot open receipt");
        output << "{\n";
        output << "  \"schema\": \"riffle-dp-2lap-g4-component-endpoint-certificate-v2\",\n";
        output << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n";
        output << "  \"evidence_label\": \""
               << (passed
                    ? (INDEPENDENT_VERIFIER
                        ? "EXACT_INDEPENDENT_EXHAUSTIVE_VERIFICATION"
                        : "EXACT_EXHAUSTIVE_CERTIFICATE")
                    : (INDEPENDENT_VERIFIER
                        ? "EXACT_INDEPENDENT_COUNTEREXAMPLE"
                        : "EXACT_COUNTEREXAMPLE"))
               << "\",\n";
        output << "  \"verification_mode\": \""
               << (INDEPENDENT_VERIFIER ? "independent_recursive" : "primary_split_table")
               << "\",\n";
        output << "  \"packet_value\": " << packetValue << ",\n";
        output << "  \"component_support_mask_hex\": \"0x" << std::hex
               << supportMask << std::dec << "\",\n";
        output << "  \"code_length\": " << Length << ",\n";
        output << "  \"window_nodes\": " << WindowNodes << ",\n";
        output << "  \"code_dimension\": " << Dimension << ",\n";
        output << "  \"component_degrees\": [";
        bool firstDegree = true;
        for (unsigned index = 0; index < ComponentDegrees.size(); ++index)
        {
            if ((supportMask >> index) & 1)
            {
                if (!firstDegree)
                    output << ',';
                output << ComponentDegrees[index];
                firstDegree = false;
            }
        }
        output << "],\n";
        output << "  \"component_annihilator_hex\": \"0x" << std::hex
               << annihilator << std::dec << "\",\n";
        output << "  \"required_two_sided_distance\": " << requiredDistance << ",\n";
        output << "  \"information_radius\": " << radius << ",\n";
        output << "  \"hamming_ball_size\": " << ballSize(radius) << ",\n";
        output << "  \"information_set_count\": " << InformationSetCount << ",\n";
        output << "  \"information_partition_restarts\": " << partition.restarts << ",\n";
        output << "  \"information_sets\": [\n";
        for (unsigned setIndex = 0; setIndex < InformationSetCount; ++setIndex)
        {
            output << "    ";
            writeInformationSet(output, partition.sets[setIndex]);
            output << (setIndex + 1 == InformationSetCount ? "\n" : ",\n");
        }
        output << "  ],\n";
        output << "  \"unused_coordinates\": [";
        for (unsigned index = 0; index < partition.unused.size(); ++index)
        {
            if (index)
                output << ',';
            output << partition.unused[index];
        }
        output << "],\n";
        output << "  \"total_information_vectors_checked\": " << totalCandidates << ",\n";
        output << "  \"minimum_enumerated_weight\": " << minimumLow << ",\n";
        output << "  \"minimum_enumerated_complement_weight\": " << minimumComplement << ",\n";
        output << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsedSeconds << ",\n";
        output << "  \"passed\": " << (passed ? "true" : "false") << ",\n";
        output << "  \"result\": \"" << (passed ? "PASS" : "COUNTEREXAMPLE") << "\",\n";
        if (!passed)
        {
            const auto failing = std::find_if(
                results.begin(),
                results.end(),
                [](const SearchResult& result) { return !result.passed; });
            output << "  \"violating_character_hex\": \"0x" << std::hex
                   << failing->violatingCharacter << std::dec << "\",\n";
            output << "  \"violating_distance\": " << failing->violatingDistance << ",\n";
            output << "  \"violating_side\": \"" << failing->violatingSide << "\",\n";
        }
        output << "  \"proof_rule\": \"The listed disjoint information sets cover "
               << InformationSetCount * Dimension
               << " coordinates. A word or complement below the required distance restricts to radius floor((d-1)/information_set_count) in one set. Each restriction determines one codeword. The search exhausts both balls around zero and one using "
               << (INDEPENDENT_VERIFIER
                    ? "direct recursive fixed-weight enumeration and an independently seeded reversed partition"
                    : "split-half lookup tables")
               << ".\",\n";
        output << "  \"scope_limitation\": \"This receipt covers one packet value and one listed irreducible-component subcode.\"\n";
        output << "}\n";
    }
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 4)
        {
            std::cerr << "usage: component_endpoint_certificate <support-mask> <packet-value> <output.json>\n";
            return 2;
        }
        const unsigned supportMask = std::stoul(argv[1], nullptr, 0);
        if (supportMask < 1 || supportMask >= (1u << ComponentDegrees.size()))
            throw std::runtime_error("component 24-node certificate: support mask outside 1..127");
        if (supportDimension(supportMask) != Dimension)
            throw std::runtime_error("component 24-node certificate: executable dimension mismatch");
        const auto annihilator = supportAnnihilator(supportMask);
        const unsigned packetValue = std::stoul(argv[2]);
        if (packetValue < 1 || packetValue > 15)
            throw std::runtime_error("component 24-node certificate: packet value outside 1..15");
        const unsigned requiredDistance = WindowNodes == 24
            ? (packetValue == 15 ? 72 : 97)
            : (packetValue == 15 ? 36 : 49);
        const unsigned radius = (requiredDistance - 1) / InformationSetCount;

        const ByteMap innerStep(systematicRightColumns());
        const auto characters = componentBasis(innerStep, annihilator);
        const auto rows = codeRows(packetValue, innerStep, characters);
        const auto partition = findInformationSets(
            supportMask,
            packetValue,
            coordinateColumns(rows));
#if !INDEPENDENT_VERIFIER
        const auto leftMasks = masksByWeight<LeftDimension>();
        const auto rightMasks = masksByWeight<RightDimension>();
#endif
        std::array<SearchResult, InformationSetCount> results{};
        const auto started = std::chrono::steady_clock::now();
        for (unsigned setIndex = 0; setIndex < InformationSetCount; ++setIndex)
        {
            const auto systematic = systematicRows(
                rows,
                characters,
                partition.sets[setIndex]);
#if INDEPENDENT_VERIFIER
            results[setIndex] = searchInformationSetIndependent(
                systematic,
                radius,
                requiredDistance);
#else
            results[setIndex] = searchInformationSet(
                systematic,
                radius,
                requiredDistance,
                leftMasks,
                rightMasks);
#endif
            std::cout
                << "value=" << packetValue
                << " information_set=" << setIndex
                << " candidates=" << results[setIndex].candidates
                << " minimum_weight=" << results[setIndex].minimumLow
                << " minimum_complement=" << results[setIndex].minimumComplement
                << " passed=" << std::boolalpha << results[setIndex].passed
                << '\n';
            if (!results[setIndex].passed)
                break;
        }
        const auto elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        writeReceipt(
            argv[3],
            supportMask,
            annihilator,
            packetValue,
            requiredDistance,
            radius,
            partition,
            results,
            elapsed);
        const bool passed = std::all_of(
            results.begin(),
            results.end(),
            [](const SearchResult& result) { return result.passed; });
        std::cout << "output=" << argv[3] << '\n';
        std::cout << "elapsed_seconds=" << std::setprecision(12) << elapsed << '\n';
        std::cout << "status=" << (passed ? "PASS" : "COUNTEREXAMPLE") << '\n';
        return passed ? 0 : 1;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 2;
    }
}
