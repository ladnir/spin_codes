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
#include <vector>

namespace {

constexpr unsigned Width = 16;
constexpr unsigned Nodes = 24;
constexpr unsigned StateBits = 32;
constexpr unsigned Words = 6;
constexpr unsigned AnchorNode = 1;
constexpr unsigned RestrictionNode = 2;
constexpr unsigned GroupCap = 5;
constexpr std::uint32_t KnownWitness = 0xcd7e70d6U;

constexpr std::array<std::uint16_t, Width> ParityColumns = {
    0x8faf, 0x9f5e, 0xbebc, 0xfd78,
    0x755f, 0xe511, 0x458d, 0x84b5,
    0x896a, 0x92d4, 0xa5a8, 0xcb50,
    0x190f, 0x321e, 0x643c, 0xc7d7,
};

struct Observation {
    std::array<std::uint64_t, Words> word{};
};

struct Candidate {
    Observation observation{};
    std::uint32_t state = 0;
    std::uint16_t selector = 0;
};

struct Minimum {
    unsigned weight = std::numeric_limits<unsigned>::max();
    std::uint64_t count = 0;
    std::uint32_t state = 0;
    std::uint16_t left = 0;
    std::uint16_t right = 0;
};

std::uint16_t accumulate(std::uint16_t value) {
    value ^= static_cast<std::uint16_t>(value << 1);
    value ^= static_cast<std::uint16_t>(value << 2);
    value ^= static_cast<std::uint16_t>(value << 4);
    value ^= static_cast<std::uint16_t>(value << 8);
    return value;
}

std::uint16_t parityMap(std::uint16_t value) {
    std::uint16_t result = 0;
    while (value != 0) {
        const unsigned bit = std::countr_zero(value);
        result ^= ParityColumns[bit];
        value &= static_cast<std::uint16_t>(value - 1);
    }
    return result;
}

Observation observe(std::uint32_t state) {
    std::uint16_t a = static_cast<std::uint16_t>(state);
    std::uint16_t b = static_cast<std::uint16_t>(state >> Width);
    Observation result;
    for (unsigned node = 0; node < Nodes; ++node) {
        const std::uint16_t accA = accumulate(a);
        const std::uint16_t accB = accumulate(b);
        const std::uint16_t acc2A = accumulate(accA);
        const std::uint16_t output = acc2A ^ accB;
        result.word[node / 4] |=
            std::uint64_t{output} << (Width * (node % 4));
        a = parityMap(accA);
        b = parityMap(acc2A) ^ parityMap(accB);
    }
    return result;
}

Observation xorObservation(Observation left, const Observation& right) {
    for (unsigned word = 0; word < Words; ++word) {
        left.word[word] ^= right.word[word];
    }
    return left;
}

unsigned weight(const Observation& observation) {
    return std::popcount(observation.word[0])
        + std::popcount(observation.word[1])
        + std::popcount(observation.word[2])
        + std::popcount(observation.word[3])
        + std::popcount(observation.word[4])
        + std::popcount(observation.word[5]);
}

std::uint16_t packet(const Observation& observation, unsigned node) {
    return static_cast<std::uint16_t>(
        observation.word[node / 4] >> (Width * (node % 4))
    );
}

std::uint32_t chartImage(const Observation& observation) {
    return std::uint32_t{packet(observation, AnchorNode)}
        | (std::uint32_t{packet(observation, RestrictionNode)} << Width);
}

std::array<std::uint32_t, StateBits> invertChart(
    const std::array<Observation, StateBits>& stateGenerators
) {
    std::array<std::uint32_t, StateBits> pivotValues{};
    std::array<std::uint32_t, StateBits> pivotStates{};
    for (unsigned bit = 0; bit < StateBits; ++bit) {
        std::uint32_t value = chartImage(stateGenerators[bit]);
        std::uint32_t state = std::uint32_t{1} << bit;
        while (value != 0) {
            const unsigned pivot = std::bit_width(value) - 1;
            if (pivotValues[pivot] != 0) {
                value ^= pivotValues[pivot];
                state ^= pivotStates[pivot];
            } else {
                pivotValues[pivot] = value;
                pivotStates[pivot] = state;
                break;
            }
        }
        if (value == 0) {
            throw std::runtime_error("chart map is singular");
        }
    }

    std::array<std::uint32_t, StateBits> result{};
    for (unsigned coordinate = 0; coordinate < StateBits; ++coordinate) {
        std::uint32_t value = std::uint32_t{1} << coordinate;
        std::uint32_t state = 0;
        while (value != 0) {
            const unsigned pivot = std::bit_width(value) - 1;
            value ^= pivotValues[pivot];
            state ^= pivotStates[pivot];
        }
        result[coordinate] = state;
        if (chartImage(observe(state)) != (std::uint32_t{1} << coordinate)) {
            throw std::runtime_error("chart inverse replay failed");
        }
    }
    return result;
}

std::vector<Candidate> buildCandidates(
    const std::array<std::uint32_t, StateBits>& chartBasis,
    unsigned offset
) {
    std::array<Observation, Width> generators{};
    for (unsigned bit = 0; bit < Width; ++bit) {
        generators[bit] = observe(chartBasis[offset + bit]);
    }
    std::vector<Candidate> result;
    result.reserve(6885);
    Observation current{};
    std::uint32_t state = 0;
    std::uint16_t previousGray = 0;
    for (std::uint32_t index = 0; index < (std::uint32_t{1} << Width); ++index) {
        const std::uint16_t gray = static_cast<std::uint16_t>(index ^ (index >> 1));
        if (index != 0) {
            const std::uint16_t difference = gray ^ previousGray;
            const unsigned bit = std::countr_zero(difference);
            current = xorObservation(current, generators[bit]);
            state ^= chartBasis[offset + bit];
        }
        if (std::popcount(gray) <= GroupCap) {
            result.push_back(Candidate{current, state, gray});
        }
        previousGray = gray;
    }
    if (result.size() != 6885) {
        throw std::runtime_error("candidate-list size mismatch");
    }
    return result;
}

void update(
    Minimum& minimum,
    unsigned candidateWeight,
    const Candidate& left,
    const Candidate& right
) {
    const std::uint32_t state = left.state ^ right.state;
    if (state == 0) {
        return;
    }
    if (candidateWeight < minimum.weight) {
        minimum.weight = candidateWeight;
        minimum.count = 1;
        minimum.state = state;
        minimum.left = left.selector;
        minimum.right = right.selector;
    } else if (candidateWeight == minimum.weight) {
        ++minimum.count;
        if (state < minimum.state) {
            minimum.state = state;
            minimum.left = left.selector;
            minimum.right = right.selector;
        }
    }
}

void merge(Minimum& target, const Minimum& source) {
    if (source.count == 0) {
        return;
    }
    if (source.weight < target.weight) {
        target = source;
    } else if (source.weight == target.weight) {
        target.count += source.count;
        if (source.state < target.state) {
            target.state = source.state;
            target.left = source.left;
            target.right = source.right;
        }
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3) {
            throw std::runtime_error(
                "usage: analyze_riffle_packetmul_wrapmul_2lap_g4_goal09_chart_exact_width16.exe OUTPUT.json [THREADS]"
            );
        }
        unsigned threadCount = std::thread::hardware_concurrency();
        if (argc == 3) {
            threadCount = static_cast<unsigned>(std::stoul(argv[2]));
        }
        threadCount = std::max(1u, threadCount);

        std::array<Observation, StateBits> stateGenerators{};
        for (unsigned bit = 0; bit < StateBits; ++bit) {
            stateGenerators[bit] = observe(std::uint32_t{1} << bit);
        }
        const auto chartBasis = invertChart(stateGenerators);
        const auto left = buildCandidates(chartBasis, 0);
        const auto right = buildCandidates(chartBasis, Width);
        const Observation knownObservation = observe(KnownWitness);
        if (weight(knownObservation) != 121
            || std::popcount(packet(knownObservation, AnchorNode)) != 4
            || std::popcount(packet(knownObservation, RestrictionNode)) != 4) {
            throw std::runtime_error("known minimum witness replay changed");
        }

        std::atomic<std::size_t> nextLeft{0};
        std::vector<Minimum> localMinima(threadCount);
        std::vector<std::thread> workers;
        workers.reserve(threadCount);
        const auto started = std::chrono::steady_clock::now();
        for (unsigned thread = 0; thread < threadCount; ++thread) {
            workers.emplace_back([&, thread] {
                Minimum& local = localMinima[thread];
                while (true) {
                    const std::size_t leftIndex = nextLeft.fetch_add(1);
                    if (leftIndex >= left.size()) {
                        break;
                    }
                    const Candidate& leftCandidate = left[leftIndex];
                    for (const Candidate& rightCandidate : right) {
                        const Observation joined = xorObservation(
                            leftCandidate.observation,
                            rightCandidate.observation
                        );
                        update(
                            local,
                            weight(joined),
                            leftCandidate,
                            rightCandidate
                        );
                    }
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
        const auto stopped = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(stopped - started).count();
        Minimum minimum;
        for (const Minimum& local : localMinima) {
            merge(minimum, local);
        }
        const Observation witnessObservation = observe(minimum.state);
        if (weight(witnessObservation) != minimum.weight
            || chartImage(witnessObservation)
                != (std::uint32_t{minimum.left}
                    | (std::uint32_t{minimum.right} << Width))) {
            throw std::runtime_error("minimum witness replay failed");
        }

        std::ofstream out(argv[1], std::ios::binary);
        if (!out) {
            throw std::runtime_error("cannot open output file");
        }
        out << "{\n"
            << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-g4-goal09-chart-exact-width16-v1\",\n"
            << "  \"result\": \"EXHAUSTED\",\n"
            << "  \"chart\": {\"anchor_node\":1,\"restriction_node\":2,"
               "\"rank\":32,\"group_weight_cap\":5},\n"
            << "  \"list_size_each\": " << left.size() << ",\n"
            << "  \"cartesian_pairs\": " << (left.size() * right.size()) << ",\n"
            << "  \"nonzero_pairs\": " << (left.size() * right.size() - 1) << ",\n"
            << "  \"minimum_weight\": " << minimum.weight << ",\n"
            << "  \"minimum_count\": " << minimum.count << ",\n"
            << "  \"witness_hex\": \"0x" << std::hex << std::setw(8)
            << std::setfill('0') << minimum.state << std::dec << std::setfill(' ')
            << "\",\n"
            << "  \"witness_anchor_weight\": "
            << std::popcount(packet(witnessObservation, AnchorNode)) << ",\n"
            << "  \"witness_restriction_weight\": "
            << std::popcount(packet(witnessObservation, RestrictionNode)) << ",\n"
            << "  \"known_global_minimum_recovered\": "
            << (minimum.weight == 121 ? "true" : "false") << ",\n"
            << "  \"threads\": " << threadCount << ",\n"
            << "  \"elapsed_seconds\": " << std::fixed << std::setprecision(6)
            << seconds << "\n"
            << "}\n";
        std::cout << "minimum_weight=" << minimum.weight << '\n'
                  << "minimum_count=" << minimum.count << '\n'
                  << "witness_hex=0x" << std::hex << std::setw(8)
                  << std::setfill('0') << minimum.state << std::dec << '\n'
                  << "pairs=" << (left.size() * right.size()) << '\n'
                  << "elapsed_seconds=" << std::fixed << std::setprecision(6)
                  << seconds << '\n'
                  << "status=EXACT_SYSTEMATIC_CHART_ENUMERATION\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error=" << error.what() << '\n';
        return 1;
    }
}
