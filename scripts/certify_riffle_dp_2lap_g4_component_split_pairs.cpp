#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

namespace
{
    struct Word384
    {
        std::array<std::uint64_t, 6> limb{};
    };

    Word384 operator^(Word384 left, const Word384& right)
    {
        left.limb[0] ^= right.limb[0];
        left.limb[1] ^= right.limb[1];
        left.limb[2] ^= right.limb[2];
        left.limb[3] ^= right.limb[3];
        left.limb[4] ^= right.limb[4];
        left.limb[5] ^= right.limb[5];
        return left;
    }

    Word384& operator^=(Word384& left, const Word384& right)
    {
        left.limb[0] ^= right.limb[0];
        left.limb[1] ^= right.limb[1];
        left.limb[2] ^= right.limb[2];
        left.limb[3] ^= right.limb[3];
        left.limb[4] ^= right.limb[4];
        left.limb[5] ^= right.limb[5];
        return left;
    }

    unsigned weight(const Word384& word)
    {
        return
            std::popcount(word.limb[0]) +
            std::popcount(word.limb[1]) +
            std::popcount(word.limb[2]) +
            std::popcount(word.limb[3]) +
            std::popcount(word.limb[4]) +
            std::popcount(word.limb[5]);
    }

    Word384 linearCombination(
        const std::vector<Word384>& basis,
        std::uint64_t information)
    {
        Word384 result;
        while (information)
        {
            const unsigned bit = std::countr_zero(information);
            result ^= basis[bit];
            information &= information - 1;
        }
        return result;
    }

    struct Input
    {
        unsigned support = 0;
        unsigned packetValue = 0;
        unsigned leftMask = 0;
        unsigned rightMask = 0;
        unsigned leftDimension = 0;
        unsigned rightDimension = 0;
        std::vector<Word384> leftBasis;
        std::vector<Word384> rightBasis;
    };

    Input readInput(const std::string& path)
    {
        std::ifstream input(path);
        if (!input)
            throw std::runtime_error("component-split pairs: cannot open input");
        std::string magic;
        input >> magic;
        if (magic != "RIFFLE_DP_2LAP_G4_COMPONENT_SPLIT_WORDS_V1")
            throw std::runtime_error("component-split pairs: input schema mismatch");
        Input result;
        input >> std::hex >> result.support >> std::dec >> result.packetValue
              >> std::hex >> result.leftMask >> result.rightMask >> std::dec
              >> result.leftDimension >> result.rightDimension;
        result.leftBasis.resize(result.leftDimension);
        result.rightBasis.resize(result.rightDimension);
        for (unsigned index = 0; index < result.leftDimension + result.rightDimension; ++index)
        {
            char side = 0;
            input >> side;
            Word384 word;
            for (auto& limb : word.limb)
                input >> std::hex >> limb;
            if (!input)
                throw std::runtime_error("component-split pairs: truncated basis");
            if (index < result.leftDimension)
            {
                if (side != 'L')
                    throw std::runtime_error("component-split pairs: left basis tag changed");
                result.leftBasis[index] = word;
            }
            else
            {
                if (side != 'R')
                    throw std::runtime_error("component-split pairs: right basis tag changed");
                result.rightBasis[index - result.leftDimension] = word;
            }
        }
        return result;
    }

    struct RightEntry
    {
        Word384 word;
        std::uint64_t information = 0;
    };

    struct Minimum
    {
        unsigned side = std::numeric_limits<unsigned>::max();
        unsigned wordWeight = 0;
        std::uint64_t leftInformation = 0;
        std::uint64_t rightInformation = 0;
    };

    bool better(const Minimum& left, const Minimum& right)
    {
        return std::tuple{
            left.side,
            left.wordWeight,
            left.leftInformation,
            left.rightInformation}
            < std::tuple{
                right.side,
                right.wordWeight,
                right.leftInformation,
                right.rightInformation};
    }
}

int main(int argc, char** argv)
{
    try
    {
        if (argc != 4)
        {
            std::cerr << "usage: component_split_pairs input.txt threads output.json\n";
            return 2;
        }
        const Input input = readInput(argv[1]);
        const unsigned threadCount = std::stoul(argv[2]);
        if (!threadCount || input.leftDimension >= 63 || input.rightDimension >= 63)
            throw std::runtime_error("component-split pairs: unsupported dimensions or threads");

        const std::uint64_t leftSize = std::uint64_t{1} << input.leftDimension;
        const std::uint64_t rightSize = std::uint64_t{1} << input.rightDimension;
        std::vector<RightEntry> right(rightSize);
        for (std::uint64_t information = 0; information < rightSize; ++information)
        {
            right[information].information = information;
            right[information].word = linearCombination(
                input.rightBasis,
                information ^ (information >> 1));
        }

        std::vector<Minimum> minima(threadCount);
        std::vector<std::thread> workers;
        const auto started = std::chrono::steady_clock::now();
        for (unsigned threadIndex = 0; threadIndex < threadCount; ++threadIndex)
        {
            const std::uint64_t begin = leftSize * threadIndex / threadCount;
            const std::uint64_t end = leftSize * (threadIndex + 1) / threadCount;
            workers.emplace_back([&, threadIndex, begin, end]() {
                Minimum local;
                for (std::uint64_t information = begin; information < end; ++information)
                {
                    const Word384 left = linearCombination(
                        input.leftBasis,
                        information ^ (information >> 1));
                    for (const RightEntry& entry : right)
                    {
                        if (!information && !entry.information)
                            continue;
                        const unsigned wordWeight = weight(left ^ entry.word);
                        const Minimum candidate{
                            std::min(wordWeight, 384U - wordWeight),
                            wordWeight,
                            information,
                            entry.information};
                        if (better(candidate, local))
                            local = candidate;
                    }
                }
                minima[threadIndex] = local;
            });
        }
        for (auto& worker : workers)
            worker.join();
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();

        Minimum minimum;
        for (const Minimum& candidate : minima)
            if (better(candidate, minimum))
                minimum = candidate;
        const std::uint64_t pairCount = leftSize * rightSize - 1;
        const Word384 witness =
            linearCombination(
                input.leftBasis,
                minimum.leftInformation ^ (minimum.leftInformation >> 1)) ^
            linearCombination(
                input.rightBasis,
                minimum.rightInformation ^ (minimum.rightInformation >> 1));
        if (weight(witness) != minimum.wordWeight)
            throw std::runtime_error("component-split pairs: witness replay failed");
        const bool passed = minimum.side >= 97;

        std::ofstream output(argv[3]);
        if (!output)
            throw std::runtime_error("component-split pairs: cannot open output");
        output << "{\n"
               << "  \"schema\": \"riffle-dp-2lap-g4-component-split-pairs-v1\",\n"
               << "  \"candidate\": \"Riffle DP-2Lap g=4\",\n"
               << "  \"evidence_label\": \"EXACT_EXHAUSTIVE_COMPONENT_SPLIT_PAIR_ENUMERATION\",\n"
               << "  \"support_mask_hex\": \"0x" << std::hex << input.support << std::dec << "\",\n"
               << "  \"packet_value\": " << input.packetValue << ",\n"
               << "  \"left_support_mask_hex\": \"0x" << std::hex << input.leftMask << std::dec << "\",\n"
               << "  \"right_support_mask_hex\": \"0x" << std::hex << input.rightMask << std::dec << "\",\n"
               << "  \"left_dimension\": " << input.leftDimension << ",\n"
               << "  \"right_dimension\": " << input.rightDimension << ",\n"
               << "  \"thread_count\": " << threadCount << ",\n"
               << "  \"pair_count\": " << pairCount << ",\n"
               << "  \"minimum_two_sided_weight\": " << minimum.side << ",\n"
               << "  \"witness_weight\": " << minimum.wordWeight << ",\n"
               << "  \"left_gray_index_hex\": \"0x" << std::hex << minimum.leftInformation << std::dec << "\",\n"
               << "  \"right_gray_index_hex\": \"0x" << std::hex << minimum.rightInformation << std::dec << "\",\n"
               << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsed << ",\n"
               << "  \"pairs_per_second\": " << std::setprecision(12) << pairCount / elapsed << ",\n"
               << "  \"result\": \"" << (passed ? "PASS" : "COUNTEREXAMPLE") << "\",\n"
               << "  \"scope_limitation\": \"This exact run covers one support, packet value, and 24-node window.\"\n"
               << "}\n";
        std::cout << "pairs=" << pairCount << '\n'
                  << "minimum_two_sided_weight=" << minimum.side << '\n'
                  << "elapsed_seconds=" << elapsed << '\n'
                  << "output=" << argv[3] << '\n'
                  << "status=" << (passed ? "PASS" : "COUNTEREXAMPLE") << '\n';
        return passed ? 0 : 1;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
