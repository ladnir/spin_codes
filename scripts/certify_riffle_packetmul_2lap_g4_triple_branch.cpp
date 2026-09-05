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
#ifndef CERTIFICATE_WINDOW
#define CERTIFICATE_WINDOW 3
#endif
    constexpr unsigned Window = CERTIFICATE_WINDOW;
    static_assert(Window >= 2);
    constexpr unsigned RequiredWeight = 6 * Window;
    constexpr std::uint64_t NibbleLowBits = 0x1111111111111111ULL;

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
            throw std::runtime_error("triple branch: generator inverse failed");
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
                throw std::runtime_error("triple branch: systematic BCH failed");
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
                throw std::runtime_error("triple branch: singular transpose step");
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
                throw std::runtime_error("triple branch: inverse verification failed");
        }
        return inverse;
    }

    unsigned nibbleWeight(std::uint64_t value)
    {
        const auto occupied = value | (value >> 1) | (value >> 2) | (value >> 3);
        return std::popcount(occupied & NibbleLowBits);
    }

    std::uint64_t integerPower(std::uint64_t base, unsigned exponent)
    {
        std::uint64_t result = 1;
        while (exponent--)
            result *= base;
        return result;
    }

    std::uint64_t choose(unsigned n, unsigned k)
    {
        std::uint64_t result = 1;
        for (unsigned index = 1; index <= k; ++index)
            result = result * (n + 1 - index) / index;
        return result;
    }

    struct Images
    {
        std::array<std::uint64_t, 2 * Window - 1> orbit{};
    };

    Images operator^(const Images& left, const Images& right)
    {
        Images result{};
        for (unsigned index = 0; index < result.orbit.size(); ++index)
            result.orbit[index] = left.orbit[index] ^ right.orbit[index];
        return result;
    }

    Images& operator^=(Images& left, const Images& right)
    {
        for (unsigned index = 0; index < left.orbit.size(); ++index)
            left.orbit[index] ^= right.orbit[index];
        return left;
    }

    struct Task
    {
        std::uint16_t support = 0;
        unsigned weight = 0;
    };

    struct Minimum
    {
        unsigned weight = Window * Nibbles + 1;
        std::uint64_t anchor = 0;
        unsigned phase = 0;
    };

    struct Shared
    {
        std::atomic<std::size_t> nextTask{0};
        std::atomic<std::uint64_t> candidates{0};
        std::atomic<bool> counterexample{false};
        std::atomic<unsigned> minimumWeight{Window * Nibbles + 1};
        std::mutex minimumMutex;
        Minimum minimum{};
    };

    void updateMinimum(Shared& shared, unsigned weight, std::uint64_t anchor, unsigned phase)
    {
        std::lock_guard lock(shared.minimumMutex);
        if (weight < shared.minimum.weight ||
            (weight == shared.minimum.weight && anchor < shared.minimum.anchor))
        {
            shared.minimum = {weight, anchor, phase};
            shared.minimumWeight.store(weight, std::memory_order_relaxed);
        }
        if (weight < RequiredWeight)
            shared.counterexample.store(true, std::memory_order_relaxed);
    }

    void processTask(
        const Task& task,
        const std::array<std::array<Images, 16>, Nibbles>& images,
        Shared& shared)
    {
        std::array<unsigned, 5> positions{};
        unsigned positionCount = 0;
        for (unsigned nibble = 0; nibble < Nibbles; ++nibble)
        {
            if ((task.support >> nibble) & 1)
                positions[positionCount++] = nibble;
        }
        if (positionCount != task.weight)
            throw std::runtime_error("triple branch: task support weight mismatch");

        std::array<unsigned, 5> values{};
        Images current{};
        for (unsigned digit = 0; digit < positionCount; ++digit)
        {
            values[digit] = 1;
            current ^= images[positions[digit]][1];
        }
        const auto taskCandidates = integerPower(15, positionCount);
        for (std::uint64_t candidate = 0; candidate < taskCandidates; ++candidate)
        {
            std::array<unsigned, 2 * Window - 1> orbitWeights{};
            for (unsigned index = 0; index < orbitWeights.size(); ++index)
                orbitWeights[index] = nibbleWeight(current.orbit[index]);
            unsigned rollingWeight = 0;
            for (unsigned index = 0; index < Window; ++index)
                rollingWeight += orbitWeights[index];
            unsigned bestWeight = rollingWeight;
            unsigned bestPhase = Window - 1;
            for (unsigned start = 1; start < Window; ++start)
            {
                rollingWeight += orbitWeights[start + Window - 1];
                rollingWeight -= orbitWeights[start - 1];
                const unsigned anchorPhase = Window - 1 - start;
                if (rollingWeight < bestWeight)
                {
                    bestWeight = rollingWeight;
                    bestPhase = anchorPhase;
                }
            }
            if (bestWeight <= shared.minimumWeight.load(std::memory_order_relaxed))
                updateMinimum(shared, bestWeight, current.orbit[Window - 1], bestPhase);
            if (bestWeight < RequiredWeight)
                break;
            if (candidate + 1 == taskCandidates)
                break;

            for (unsigned digit = 0; digit < positionCount; ++digit)
            {
                const auto oldValue = values[digit];
                const auto newValue = oldValue == 15 ? 1 : oldValue + 1;
                current ^= images[positions[digit]][oldValue];
                current ^= images[positions[digit]][newValue];
                values[digit] = newValue;
                if (oldValue != 15)
                    break;
            }
        }
        if (!shared.counterexample.load(std::memory_order_relaxed))
            shared.candidates.fetch_add(taskCandidates, std::memory_order_relaxed);
    }

    void worker(
        const std::vector<Task>& tasks,
        const std::array<std::array<Images, 16>, Nibbles>& images,
        Shared& shared)
    {
        while (!shared.counterexample.load(std::memory_order_relaxed))
        {
            const auto index = shared.nextTask.fetch_add(1, std::memory_order_relaxed);
            if (index >= tasks.size())
                break;
            processTask(tasks[index], images, shared);
        }
    }

    std::array<std::uint64_t, Window> windowFor(
        const Minimum& minimum,
        const std::array<std::uint64_t, StateBits>& forward,
        const std::array<std::uint64_t, StateBits>& backward)
    {
        const auto x = minimum.anchor;
        std::array<std::uint64_t, Window> result{};
        result[minimum.phase] = x;
        for (unsigned index = minimum.phase + 1; index < Window; ++index)
            result[index] = applyColumns(forward, result[index - 1]);
        for (unsigned index = minimum.phase; index > 0; --index)
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
            : "constructions/riffle_packetmul_2lap_g4/receipts/goal02_triple_branch_primary.json";
        const auto forward = transposeStepColumns();
        const auto backward = inverseColumns(forward);
        std::array<std::array<Images, 16>, Nibbles> images{};
        for (unsigned nibble = 0; nibble < Nibbles; ++nibble)
        {
            for (unsigned value = 1; value < 16; ++value)
            {
                const auto state = std::uint64_t{value} << (4 * nibble);
                auto& image = images[nibble][value];
                image.orbit[Window - 1] = state;
                for (unsigned index = Window; index < image.orbit.size(); ++index)
                    image.orbit[index] = applyColumns(forward, image.orbit[index - 1]);
                for (unsigned index = Window - 1; index > 0; --index)
                    image.orbit[index - 1] = applyColumns(backward, image.orbit[index]);
            }
        }

        std::vector<Task> tasks;
        for (unsigned weight = 5; weight >= 1; --weight)
        {
            for (unsigned mask = 1; mask < (1u << Nibbles); ++mask)
            {
                if (static_cast<unsigned>(std::popcount(mask)) == weight)
                    tasks.push_back({static_cast<std::uint16_t>(mask), weight});
            }
        }
        const auto expectedCandidates = [&]()
        {
            std::uint64_t result = 0;
            for (unsigned weight = 1; weight <= 5; ++weight)
                result += choose(Nibbles, weight) * integerPower(15, weight);
            return result;
        }();

        Shared shared{};
        const unsigned threadCount = std::max(
            1u,
            std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> threads;
        for (unsigned index = 0; index < threadCount; ++index)
            threads.emplace_back(worker, std::cref(tasks), std::cref(images), std::ref(shared));
        for (auto& thread : threads)
            thread.join();

        const auto checked = shared.candidates.load();
        const bool passed = !shared.counterexample.load();
        if (passed && checked != expectedCandidates)
            throw std::runtime_error("triple branch: candidate count mismatch");
        const auto minimumWindow = windowFor(shared.minimum, forward, backward);
        std::array<unsigned, Window> weights{};
        unsigned replayWeight = 0;
        for (unsigned index = 0; index < Window; ++index)
        {
            weights[index] = nibbleWeight(minimumWindow[index]);
            replayWeight += weights[index];
        }
        if (replayWeight != shared.minimum.weight)
            throw std::runtime_error("triple branch: minimum replay mismatch");
        for (unsigned index = 1; index < Window; ++index)
        {
            if (applyColumns(forward, minimumWindow[index - 1]) != minimumWindow[index])
                throw std::runtime_error("triple branch: witness is not consecutive");
        }

        std::ofstream output(outputPath);
        if (!output)
            throw std::runtime_error("triple branch: cannot open output");
        output << "{\n";
        output << "  \"schema\": \"riffle-packetmul-2lap-g4-goal02-window-branch-primary-v1\",\n";
        output << "  \"candidate\": \"Riffle PacketMul-2Lap g=4\",\n";
        output << "  \"candidate_id\": \"riffle_packetmul_2lap_g4\",\n";
        output << "  \"evidence_label\": \""
               << (passed ? "EXACT_EXHAUSTIVE_PRIMARY" : "EXACT_COUNTEREXAMPLE")
               << "\",\n";
        output << "  \"window_states\": " << Window << ",\n";
        output << "  \"required_window_weight\": " << RequiredWeight << ",\n";
        output << "  \"low_anchor_maximum_weight\": 5,\n";
        output << "  \"thread_count\": " << threadCount << ",\n";
        output << "  \"expected_low_anchor_states\": " << expectedCandidates << ",\n";
        output << "  \"low_anchor_states_checked\": " << checked << ",\n";
        output << "  \"minimum_found_weight\": " << shared.minimum.weight << ",\n";
        output << "  \"minimum_anchor_phase\": " << shared.minimum.phase << ",\n";
        output << "  \"minimum_anchor_hex\": \"0x" << std::hex << shared.minimum.anchor << "\",\n";
        output << "  \"minimum_window_hex\": [";
        for (unsigned index = 0; index < Window; ++index)
        {
            if (index)
                output << ", ";
            output << "\"0x" << minimumWindow[index] << "\"";
        }
        output << "],\n";
        output << std::dec;
        output << "  \"minimum_window_nibble_weights\": [";
        for (unsigned index = 0; index < Window; ++index)
        {
            if (index)
                output << ", ";
            output << weights[index];
        }
        output << "],\n";
        output << "  \"proof_rule\": \"If a consecutive w-state window has total nibble weight at most 6w-1, one of its states has weight at most 5. Re-anchor at that state. The exhaustive search enumerates every nonzero 64-bit state of nibble weight at most 5 exactly once and evaluates all w anchor phases. Absence of a weight-below-6w window proves the local bound.\",\n";
        output << "  \"orbit_consequence\": \"When w divides 32772, invertibility and a window weight lower bound of 6w imply sum_t wt_4(S^t chi) >= 6*32772 for every nonzero chi.\",\n";
        output << "  \"status\": \""
               << (passed ? "EXACT_WINDOW_BRANCH_CERTIFICATE" : "EXACT_LOCAL_COUNTEREXAMPLE")
               << "\"\n";
        output << "}\n";

        std::cout << "low_anchor_states_checked=" << checked << std::endl;
        std::cout << "minimum_found_weight=" << shared.minimum.weight << std::endl;
        std::cout << "minimum_window_weights=";
        for (unsigned index = 0; index < Window; ++index)
            std::cout << (index ? "," : "") << weights[index];
        std::cout << std::endl;
        std::cout << "output=" << outputPath << std::endl;
        std::cout << "status="
                  << (passed ? "EXACT_WINDOW_BRANCH_CERTIFICATE" : "EXACT_LOCAL_COUNTEREXAMPLE")
                  << std::endl;
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
