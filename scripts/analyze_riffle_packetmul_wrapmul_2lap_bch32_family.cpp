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
constexpr unsigned StateBits = 2 * Width;
constexpr unsigned Words = Nodes * Width / 64;
constexpr unsigned LowBits = 20;
constexpr unsigned HighBits = StateBits - LowBits;
constexpr std::uint32_t LowCount = std::uint32_t{1} << LowBits;
constexpr std::uint32_t TaskCount = std::uint32_t{1} << HighBits;
constexpr unsigned MaxWeight = Nodes * Width;
constexpr unsigned Inf = std::numeric_limits<unsigned>::max();

// Canonical left-systematic parity columns of the extended primitive
// BCH [32,16,8] code. The underlying primitive polynomial is 0x25 and the
// cyclic BCH generator polynomial is 0x8faf.
constexpr std::array<std::uint16_t, Width> ParityColumns = {
    0x8faf, 0x9f5e, 0xbebc, 0xfd78,
    0x755f, 0xe511, 0x458d, 0x84b5,
    0x896a, 0x92d4, 0xa5a8, 0xcb50,
    0x190f, 0x321e, 0x643c, 0xc7d7,
};

struct Observation {
    std::array<std::uint64_t, Words> word{};
};

struct Minimum {
    unsigned weight = Inf;
    std::uint64_t count = 0;
    std::uint32_t witness = 0;
};

struct LocalResult {
    std::array<Minimum, Nodes> prefix{};
    std::array<std::array<Minimum, Width + 1>, Nodes> anchor{};
    std::array<std::uint64_t, MaxWeight + 1> spectrum{};
    std::uint64_t states = 0;
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
        const std::uint16_t q = acc2A ^ accB;
        result.word[node / 4] |=
            std::uint64_t{q} << (Width * (node % 4));

        const std::uint16_t nextA = parityMap(accA);
        const std::uint16_t nextB = parityMap(acc2A) ^ parityMap(accB);
        a = nextA;
        b = nextB;
    }
    return result;
}

Observation xorObservation(Observation left, const Observation& right) {
    for (unsigned i = 0; i < Words; ++i) {
        left.word[i] ^= right.word[i];
    }
    return left;
}

void updateMinimum(Minimum& minimum, unsigned weight, std::uint32_t state) {
    if (weight < minimum.weight) {
        minimum.weight = weight;
        minimum.count = 1;
        minimum.witness = state;
    } else if (weight == minimum.weight) {
        ++minimum.count;
        minimum.witness = std::min(minimum.witness, state);
    }
}

void processObservation(
    const Observation& observation,
    std::uint32_t state,
    LocalResult& result
) {
    std::array<unsigned char, Nodes> weights{};
    unsigned total = 0;
    for (unsigned node = 0; node < Nodes; ++node) {
        const std::uint16_t block = static_cast<std::uint16_t>(
            observation.word[node / 4] >> (Width * (node % 4))
        );
        const unsigned weight = std::popcount(block);
        weights[node] = static_cast<unsigned char>(weight);
        total += weight;
        updateMinimum(result.prefix[node], total, state);
    }
    ++result.spectrum[total];
    for (unsigned node = 0; node < Nodes; ++node) {
        updateMinimum(result.anchor[node][weights[node]], total, state);
    }
    ++result.states;
}

void mergeMinimum(Minimum& target, const Minimum& source) {
    if (source.count == 0) {
        return;
    }
    if (source.weight < target.weight) {
        target = source;
    } else if (source.weight == target.weight) {
        target.count += source.count;
        target.witness = std::min(target.witness, source.witness);
    }
}

std::array<unsigned, Nodes> replayWeights(std::uint32_t state) {
    const Observation observation = observe(state);
    std::array<unsigned, Nodes> result{};
    for (unsigned node = 0; node < Nodes; ++node) {
        const std::uint16_t block = static_cast<std::uint16_t>(
            observation.word[node / 4] >> (Width * (node % 4))
        );
        result[node] = std::popcount(block);
    }
    return result;
}

void writeMinimum(std::ostream& out, const Minimum& minimum) {
    if (minimum.count == 0) {
        out << "null";
        return;
    }
    out << "{\"weight\":" << minimum.weight
        << ",\"count\":" << minimum.count
        << ",\"witness_hex\":\"0x"
        << std::hex << std::setw(8) << std::setfill('0') << minimum.witness
        << std::dec << std::setfill(' ') << "\"}";
}

void writeWeights(std::ostream& out, const std::array<unsigned, Nodes>& weights) {
    out << '[';
    for (unsigned i = 0; i < Nodes; ++i) {
        if (i != 0) {
            out << ',';
        }
        out << weights[i];
    }
    out << ']';
}

} // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 2 || argc > 3) {
            throw std::runtime_error(
                "usage: analyze_riffle_packetmul_wrapmul_2lap_bch32_family.exe OUTPUT.json [THREADS]"
            );
        }
        unsigned threadCount = std::thread::hardware_concurrency();
        if (argc == 3) {
            threadCount = static_cast<unsigned>(std::stoul(argv[2]));
        }
        threadCount = std::max(1u, threadCount);

        std::array<Observation, StateBits> generators{};
        for (unsigned bit = 0; bit < StateBits; ++bit) {
            generators[bit] = observe(std::uint32_t{1} << bit);
        }

        std::atomic<std::uint32_t> nextTask{0};
        std::vector<LocalResult> locals(threadCount);
        std::vector<std::thread> workers;
        workers.reserve(threadCount);

        const auto start = std::chrono::steady_clock::now();
        for (unsigned thread = 0; thread < threadCount; ++thread) {
            workers.emplace_back([&, thread] {
                LocalResult& local = locals[thread];
                while (true) {
                    const std::uint32_t high = nextTask.fetch_add(1);
                    if (high >= TaskCount) {
                        break;
                    }
                    Observation current{};
                    std::uint32_t h = high;
                    while (h != 0) {
                        const unsigned bit = std::countr_zero(h);
                        current = xorObservation(current, generators[LowBits + bit]);
                        h &= h - 1;
                    }

                    const std::uint32_t highState = high << LowBits;
                    if (highState != 0) {
                        processObservation(current, highState, local);
                    }
                    std::uint32_t previousGray = 0;
                    for (std::uint32_t index = 1; index < LowCount; ++index) {
                        const std::uint32_t gray = index ^ (index >> 1);
                        const std::uint32_t difference = gray ^ previousGray;
                        current = xorObservation(
                            current,
                            generators[std::countr_zero(difference)]
                        );
                        processObservation(current, highState | gray, local);
                        previousGray = gray;
                    }
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
        const auto stop = std::chrono::steady_clock::now();
        const double seconds = std::chrono::duration<double>(stop - start).count();

        LocalResult total;
        for (const LocalResult& local : locals) {
            total.states += local.states;
            for (unsigned weight = 0; weight <= MaxWeight; ++weight) {
                total.spectrum[weight] += local.spectrum[weight];
            }
            for (unsigned node = 0; node < Nodes; ++node) {
                mergeMinimum(total.prefix[node], local.prefix[node]);
                for (unsigned anchor = 0; anchor <= Width; ++anchor) {
                    mergeMinimum(total.anchor[node][anchor], local.anchor[node][anchor]);
                }
            }
        }

        if (total.states != std::numeric_limits<std::uint32_t>::max()) {
            throw std::runtime_error("exact state count mismatch");
        }

        std::ofstream out(argv[1], std::ios::binary);
        if (!out) {
            throw std::runtime_error("cannot open output file");
        }
        out << "{\n"
            << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-bch32-family-v1\",\n"
            << "  \"result\": \"EXHAUSTED\",\n"
            << "  \"instance\": {\"m\":5,\"primitive_polynomial_hex\":\"0x25\","
               "\"designed_distance\":7,\"cyclic_parameters\":\"[31,16,7]\","
               "\"extended_parameters\":\"[32,16,8]\","
               "\"generator_polynomial_hex\":\"0x8faf\","
               "\"state_bits\":32,\"output_bits_per_node\":16,\"nodes\":24,"
               "\"parity_columns_hex\":[";
        for (unsigned bit = 0; bit < Width; ++bit) {
            if (bit != 0) {
                out << ',';
            }
            out << "\"0x" << std::hex << ParityColumns[bit] << std::dec << "\"";
        }
        out << "]},\n"
            << "  \"enumeration\": {\"nonzero_states\":" << total.states
            << ",\"threads\":" << threadCount
            << ",\"seconds\":" << std::fixed << std::setprecision(6) << seconds
            << "},\n";

        out << "  \"prefix_minima\": [";
        for (unsigned node = 0; node < Nodes; ++node) {
            if (node != 0) {
                out << ',';
            }
            writeMinimum(out, total.prefix[node]);
        }
        out << "],\n";

        out << "  \"weight24_spectrum\": {";
        bool first = true;
        for (unsigned weight = 0; weight <= MaxWeight; ++weight) {
            if (total.spectrum[weight] == 0) {
                continue;
            }
            if (!first) {
                out << ',';
            }
            first = false;
            out << '\"' << weight << "\":" << total.spectrum[weight];
        }
        out << "},\n";

        out << "  \"anchor_minima\": [\n";
        for (unsigned node = 0; node < Nodes; ++node) {
            out << "    [";
            for (unsigned anchor = 0; anchor <= Width; ++anchor) {
                if (anchor != 0) {
                    out << ',';
                }
                writeMinimum(out, total.anchor[node][anchor]);
            }
            out << ']' << (node + 1 == Nodes ? "\n" : ",\n");
        }
        out << "  ],\n";

        const Minimum& global = total.prefix[Nodes - 1];
        out << "  \"global_witness_node_weights\": ";
        writeWeights(out, replayWeights(global.witness));
        out << ",\n";

        out << "  \"selected_anchor_witness_weights\": {\n";
        bool firstSelected = true;
        for (unsigned node : {0u, 12u, 23u}) {
            for (unsigned anchor : {0u, 1u}) {
                const Minimum& minimum = total.anchor[node][anchor];
                if (minimum.count == 0) {
                    continue;
                }
                if (!firstSelected) {
                    out << ",\n";
                }
                firstSelected = false;
                out << "    \"node" << node << "_weight" << anchor << "\": ";
                writeWeights(out, replayWeights(minimum.witness));
            }
        }
        out << "\n  }\n}\n";
        out.close();

        std::cout << "EXHAUSTED states=" << total.states
                  << " d24=" << global.weight
                  << " count=" << global.count
                  << " witness=0x" << std::hex << global.witness << std::dec
                  << " seconds=" << seconds << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
