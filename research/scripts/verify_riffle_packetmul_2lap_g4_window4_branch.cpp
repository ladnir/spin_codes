#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{
    constexpr std::uint64_t Generator = 0xF4845518B9582A1FULL;
    constexpr unsigned StateBits = 64;
    constexpr unsigned Nibbles = 16;
    constexpr unsigned Window = 4;
    constexpr unsigned RequiredWeight = 24;
    constexpr std::uint64_t NibbleLowBits = 0x1111111111111111ULL;

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

    std::array<std::uint64_t, StateBits> inverseColumns(
        const std::array<std::uint64_t, StateBits>& columns)
    {
        struct Row
        {
            std::uint64_t left = 0;
            std::uint64_t right = 0;
        };
        std::array<Row, StateBits> rows{};
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
                throw std::runtime_error("window4 verifier: singular matrix");
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
        for (unsigned column = 0; column < StateBits; ++column)
        {
            for (unsigned row = 0; row < StateBits; ++row)
                inverse[column] |= ((rows[row].right >> column) & 1) << row;
        }
        for (unsigned bit = 0; bit < StateBits; ++bit)
        {
            if (applyColumns(columns, inverse[bit]) != (std::uint64_t{1} << bit))
                throw std::runtime_error("window4 verifier: inverse check failed");
        }
        return inverse;
    }

    std::array<std::uint64_t, StateBits> systematicStateColumnsIndependent()
    {
        std::array<std::uint64_t, StateBits> left{};
        std::array<std::uint64_t, StateBits> right{};
        for (unsigned messageBit = 0; messageBit < StateBits; ++messageBit)
        {
            left[messageBit] = Generator << messageBit;
            right[messageBit] =
                (messageBit == 0 ? 0 : Generator >> (StateBits - messageBit)) |
                (std::uint64_t{1} << 63);
        }
        const auto inverseLeft = inverseColumns(left);
        std::array<std::uint64_t, StateBits> state{};
        for (unsigned column = 0; column < StateBits; ++column)
        {
            if (applyColumns(left, inverseLeft[column]) != (std::uint64_t{1} << column))
                throw std::runtime_error("window4 verifier: systematic left check failed");
            state[column] = applyColumns(right, inverseLeft[column]);
        }
        return state;
    }

    std::array<std::uint64_t, StateBits> transposeColumns(
        const std::array<std::uint64_t, StateBits>& columns)
    {
        std::array<std::uint64_t, StateBits> result{};
        for (unsigned newInput = 0; newInput < StateBits; ++newInput)
        {
            for (unsigned oldInput = 0; oldInput < StateBits; ++oldInput)
                result[newInput] |= ((columns[oldInput] >> newInput) & 1) << oldInput;
        }
        return result;
    }

    std::array<std::uint64_t, StateBits> transposeStepColumnsIndependent()
    {
        const auto state = systematicStateColumnsIndependent();
        std::array<std::uint64_t, StateBits> forwardStep{};
        for (unsigned input = 0; input < StateBits; ++input)
        {
            // The accumulator sends basis bit i to bits i through 63.
            const auto accumulatedBasis = ~std::uint64_t{0} << input;
            forwardStep[input] = applyColumns(state, accumulatedBasis);
        }
        return transposeColumns(forwardStep);
    }

    unsigned nibbleWeight(std::uint64_t value)
    {
        const auto occupied = value | (value >> 1) | (value >> 2) | (value >> 3);
        return std::popcount(occupied & NibbleLowBits);
    }

    struct OrbitImages
    {
        std::array<std::uint64_t, 2 * Window - 1> states{};
    };

    OrbitImages& operator^=(OrbitImages& left, const OrbitImages& right)
    {
        for (unsigned index = 0; index < left.states.size(); ++index)
            left.states[index] ^= right.states[index];
        return left;
    }

    struct Minimum
    {
        unsigned weight = Window * Nibbles + 1;
        std::uint64_t anchor = 0;
        unsigned anchorPosition = 0;
    };

    struct Shared
    {
        std::atomic<unsigned> nextSupport{0};
        std::atomic<std::uint64_t> candidates{0};
        std::atomic<unsigned> minimumWeight{Window * Nibbles + 1};
        std::atomic<bool> counterexample{false};
        std::mutex minimumMutex;
        Minimum minimum{};
    };

    void updateMinimum(
        Shared& shared,
        unsigned weight,
        std::uint64_t anchor,
        unsigned anchorPosition)
    {
        std::lock_guard lock(shared.minimumMutex);
        if (weight < shared.minimum.weight ||
            (weight == shared.minimum.weight && anchor < shared.minimum.anchor))
        {
            shared.minimum = {weight, anchor, anchorPosition};
            shared.minimumWeight.store(weight, std::memory_order_relaxed);
        }
        if (weight < RequiredWeight)
            shared.counterexample.store(true, std::memory_order_relaxed);
    }

    void evaluate(
        const OrbitImages& orbit,
        std::uint64_t anchor,
        Shared& shared)
    {
        std::array<unsigned, 2 * Window - 1> weights{};
        for (unsigned index = 0; index < weights.size(); ++index)
            weights[index] = nibbleWeight(orbit.states[index]);
        unsigned rolling = 0;
        for (unsigned index = 0; index < Window; ++index)
            rolling += weights[index];
        unsigned best = rolling;
        unsigned bestAnchorPosition = Window - 1;
        for (unsigned start = 1; start < Window; ++start)
        {
            rolling += weights[start + Window - 1] - weights[start - 1];
            if (rolling < best)
            {
                best = rolling;
                bestAnchorPosition = Window - 1 - start;
            }
        }
        if (best <= shared.minimumWeight.load(std::memory_order_relaxed))
            updateMinimum(shared, best, anchor, bestAnchorPosition);
    }

    void worker(
        const std::vector<std::uint16_t>& supports,
        const std::array<OrbitImages, StateBits>& bitImages,
        Shared& shared)
    {
        while (!shared.counterexample.load(std::memory_order_relaxed))
        {
            const auto supportIndex = shared.nextSupport.fetch_add(1, std::memory_order_relaxed);
            if (supportIndex >= supports.size())
                break;
            const auto support = supports[supportIndex];
            std::array<unsigned, 20> bits{};
            unsigned bitCount = 0;
            for (unsigned nibble = 0; nibble < Nibbles; ++nibble)
            {
                if (!((support >> nibble) & 1))
                    continue;
                for (unsigned bit = 0; bit < 4; ++bit)
                    bits[bitCount++] = 4 * nibble + bit;
            }
            if (bitCount != 20)
                throw std::runtime_error("window4 verifier: support dimension differs");
            OrbitImages current{};
            for (std::uint64_t grayIndex = 1; grayIndex < (std::uint64_t{1} << bitCount); ++grayIndex)
            {
                const auto changed = std::countr_zero(grayIndex);
                current ^= bitImages[bits[changed]];
                evaluate(current, current.states[Window - 1], shared);
                if (shared.counterexample.load(std::memory_order_relaxed))
                    break;
            }
            if (!shared.counterexample.load(std::memory_order_relaxed))
                shared.candidates.fetch_add((std::uint64_t{1} << bitCount) - 1, std::memory_order_relaxed);
        }
    }

    std::array<std::uint64_t, Window> replayWindow(
        const Minimum& minimum,
        const std::array<std::uint64_t, StateBits>& forward,
        const std::array<std::uint64_t, StateBits>& backward)
    {
        std::array<std::uint64_t, Window> result{};
        result[minimum.anchorPosition] = minimum.anchor;
        for (unsigned index = minimum.anchorPosition + 1; index < Window; ++index)
            result[index] = applyColumns(forward, result[index - 1]);
        for (unsigned index = minimum.anchorPosition; index > 0; --index)
            result[index - 1] = applyColumns(backward, result[index]);
        return result;
    }
}

int main(int argc, char** argv)
{
    try
    {
        const std::string outputPath = argc > 1
            ? argv[1]
            : "constructions/riffle_packetmul_2lap_g4/receipts/goal02_window4_branch_independent.json";
        const auto forward = transposeStepColumnsIndependent();
        const auto backward = inverseColumns(forward);
        std::array<OrbitImages, StateBits> bitImages{};
        for (unsigned bit = 0; bit < StateBits; ++bit)
        {
            auto& image = bitImages[bit];
            image.states[Window - 1] = std::uint64_t{1} << bit;
            for (unsigned index = Window; index < image.states.size(); ++index)
                image.states[index] = applyColumns(forward, image.states[index - 1]);
            for (unsigned index = Window - 1; index > 0; --index)
                image.states[index - 1] = applyColumns(backward, image.states[index]);
        }

        std::vector<std::uint16_t> supports;
        for (unsigned mask = 0; mask < (1u << Nibbles); ++mask)
        {
            if (std::popcount(mask) == 5)
                supports.push_back(static_cast<std::uint16_t>(mask));
        }
        if (supports.size() != 4368)
            throw std::runtime_error("window4 verifier: five-nibble support count differs");

        Shared shared{};
        const unsigned threadCount = std::max(
            1u,
            std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> threads;
        for (unsigned index = 0; index < threadCount; ++index)
            threads.emplace_back(worker, std::cref(supports), std::cref(bitImages), std::ref(shared));
        for (auto& thread : threads)
            thread.join();

        const auto expected =
            static_cast<std::uint64_t>(supports.size()) * ((std::uint64_t{1} << 20) - 1);
        const auto checked = shared.candidates.load();
        const bool passed = !shared.counterexample.load();
        if (passed && checked != expected)
            throw std::runtime_error("window4 verifier: candidate count differs");
        const auto window = replayWindow(shared.minimum, forward, backward);
        std::array<unsigned, Window> weights{};
        unsigned total = 0;
        for (unsigned index = 0; index < Window; ++index)
        {
            weights[index] = nibbleWeight(window[index]);
            total += weights[index];
            if (index && applyColumns(forward, window[index - 1]) != window[index])
                throw std::runtime_error("window4 verifier: replay is not consecutive");
        }
        if (total != shared.minimum.weight)
            throw std::runtime_error("window4 verifier: replay weight differs");

        std::ofstream output(outputPath);
        if (!output)
            throw std::runtime_error("window4 verifier: cannot open output");
        output << "{\n";
        output << "  \"schema\": \"riffle-packetmul-2lap-g4-goal02-window4-branch-independent-v1\",\n";
        output << "  \"candidate\": \"Riffle PacketMul-2Lap g=4\",\n";
        output << "  \"candidate_id\": \"riffle_packetmul_2lap_g4\",\n";
        output << "  \"evidence_label\": \""
               << (passed ? "EXACT_EXHAUSTIVE_INDEPENDENT" : "EXACT_COUNTEREXAMPLE")
               << "\",\n";
        output << "  \"window_states\": 4,\n";
        output << "  \"required_window_weight\": 24,\n";
        output << "  \"enumeration\": \"Every five-nibble support and every nonzero 20-bit word on that support, using binary reflected Gray order. States of smaller support occur redundantly.\",\n";
        output << "  \"five_nibble_supports\": " << supports.size() << ",\n";
        output << "  \"words_per_support\": " << ((std::uint64_t{1} << 20) - 1) << ",\n";
        output << "  \"expected_words_checked_with_duplicates\": " << expected << ",\n";
        output << "  \"words_checked_with_duplicates\": " << checked << ",\n";
        output << "  \"thread_count\": " << threadCount << ",\n";
        output << "  \"minimum_found_weight\": " << shared.minimum.weight << ",\n";
        output << "  \"minimum_anchor_position\": " << shared.minimum.anchorPosition << ",\n";
        output << "  \"minimum_anchor_hex\": \"0x" << std::hex << shared.minimum.anchor << "\",\n";
        output << "  \"minimum_window_hex\": [";
        for (unsigned index = 0; index < Window; ++index)
            output << (index ? ", \"0x" : "\"0x") << window[index] << "\"";
        output << "],\n";
        output << std::dec;
        output << "  \"minimum_window_nibble_weights\": [";
        for (unsigned index = 0; index < Window; ++index)
            output << (index ? ", " : "") << weights[index];
        output << "],\n";
        output << "  \"coverage_rule\": \"Every state with at most five nonzero nibbles is supported inside at least one five-nibble set. Every violating four-state window has an anchor state of weight at most five. The four possible anchor positions cover every violating window.\",\n";
        output << "  \"status\": \""
               << (passed ? "EXACT_WINDOW4_BRANCH_INDEPENDENT" : "EXACT_LOCAL_COUNTEREXAMPLE")
               << "\"\n";
        output << "}\n";

        std::cout << "words_checked_with_duplicates=" << checked << std::endl;
        std::cout << "minimum_found_weight=" << shared.minimum.weight << std::endl;
        std::cout << "output=" << outputPath << std::endl;
        std::cout << "status="
                  << (passed ? "EXACT_WINDOW4_BRANCH_INDEPENDENT" : "EXACT_LOCAL_COUNTEREXAMPLE")
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
