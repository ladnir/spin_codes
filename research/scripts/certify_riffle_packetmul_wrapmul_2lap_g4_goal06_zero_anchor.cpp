#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr unsigned kMaximumSets = 22;
constexpr unsigned kDimension = 64;
constexpr unsigned kWords = 23;
constexpr unsigned kBound = 143;

struct Row {
    std::array<std::uint64_t, kWords> word{};
    std::array<std::uint64_t, 2> state{};
};

struct Input {
    std::vector<std::array<Row, kDimension>> rows;
    unsigned base_weight = 0;
    unsigned base_weight_set_count = 0;
    std::uint64_t expected_candidates = 0;
};

template <class T>
void read_exact(std::ifstream& input, T* destination, std::size_t count = 1) {
    input.read(reinterpret_cast<char*>(destination),
               static_cast<std::streamsize>(sizeof(T) * count));
    if (!input) {
        throw std::runtime_error("goal06 zero anchor: truncated input");
    }
}

Input load_input(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("goal06 zero anchor: cannot open input");
    }
    std::array<char, 8> magic{};
    read_exact(input, magic.data(), magic.size());
    if (std::string(magic.data(), magic.size()) != std::string("RZ24Z02\0", 8)) {
        throw std::runtime_error("goal06 zero anchor: bad input magic");
    }
    std::uint32_t sets = 0;
    std::uint32_t dimension = 0;
    std::uint32_t words = 0;
    std::uint32_t bound = 0;
    std::uint32_t base_weight = 0;
    std::uint32_t base_weight_set_count = 0;
    std::uint64_t expected_candidates = 0;
    read_exact(input, &sets);
    read_exact(input, &dimension);
    read_exact(input, &words);
    read_exact(input, &bound);
    read_exact(input, &base_weight);
    read_exact(input, &base_weight_set_count);
    read_exact(input, &expected_candidates);
    if (sets == 0 || sets > kMaximumSets || dimension != kDimension ||
        words != kWords || bound != kBound || base_weight > 7 ||
        base_weight_set_count == 0 || base_weight_set_count > sets) {
        throw std::runtime_error("goal06 zero anchor: input dimensions changed");
    }
    Input result;
    result.rows.resize(sets);
    result.base_weight = base_weight;
    result.base_weight_set_count = base_weight_set_count;
    result.expected_candidates = expected_candidates;
    for (auto& set_rows : result.rows) {
        for (auto& row : set_rows) {
            read_exact(input, row.word.data(), row.word.size());
            read_exact(input, row.state.data(), row.state.size());
        }
    }
    char extra = 0;
    if (input.read(&extra, 1)) {
        throw std::runtime_error("goal06 zero anchor: trailing input data");
    }
    return result;
}

struct Search {
    const Input& input;
    const std::array<Row, kDimension>* active_rows = nullptr;
    std::uint64_t candidates = 0;
    bool counterexample = false;
    std::uint16_t counterexample_weight = 0;
    std::array<std::uint64_t, 2> counterexample_state{};

    void accept(std::uint64_t selector, std::uint16_t weight) {
        std::array<std::uint64_t, 2> state{};
        std::uint64_t remaining = selector;
        while (remaining != 0) {
            const unsigned bit = std::countr_zero(remaining);
            state[0] ^= (*active_rows)[bit].state[0];
            state[1] ^= (*active_rows)[bit].state[1];
            remaining &= remaining - 1;
        }
        if ((state[0] | state[1]) == 0) {
            return;
        }
        counterexample = true;
        counterexample_weight = weight;
        counterexample_state = state;
    }

    void evaluate_zero() {
        ++candidates;
    }

    void evaluate_last(
        std::uint64_t selector,
        unsigned bit,
        const std::array<std::uint64_t, kWords>& current) {
        ++candidates;
        std::uint16_t weight = 0;
        const Row& row = (*active_rows)[bit];
        for (unsigned word = 0; word < kWords; ++word) {
            weight = static_cast<std::uint16_t>(
                weight + std::popcount(current[word] ^ row.word[word]));
            if (weight > kBound) {
                return;
            }
        }
        accept(selector | (std::uint64_t{1} << bit), weight);
    }

    template <unsigned Remaining>
    void enumerate_exact(
        unsigned start,
        std::uint64_t selector,
        std::array<std::uint64_t, kWords>& current) {
        if constexpr (Remaining == 0) {
            evaluate_zero();
        } else if constexpr (Remaining == 1) {
            for (unsigned bit = start; bit < kDimension && !counterexample; ++bit) {
                evaluate_last(selector, bit, current);
            }
        } else {
            const unsigned last = kDimension - Remaining;
            for (unsigned bit = start; bit <= last && !counterexample; ++bit) {
                for (unsigned word = 0; word < kWords; ++word) {
                    current[word] ^= (*active_rows)[bit].word[word];
                }
                enumerate_exact<Remaining - 1>(
                    bit + 1,
                    selector | (std::uint64_t{1} << bit),
                    current);
                for (unsigned word = 0; word < kWords; ++word) {
                    current[word] ^= (*active_rows)[bit].word[word];
                }
            }
        }
    }

    void enumerate_weight(unsigned target) {
        std::array<std::uint64_t, kWords> current{};
        switch (target) {
            case 0: enumerate_exact<0>(0, 0, current); break;
            case 1: enumerate_exact<1>(0, 0, current); break;
            case 2: enumerate_exact<2>(0, 0, current); break;
            case 3: enumerate_exact<3>(0, 0, current); break;
            case 4: enumerate_exact<4>(0, 0, current); break;
            case 5: enumerate_exact<5>(0, 0, current); break;
            case 6: enumerate_exact<6>(0, 0, current); break;
            case 7: enumerate_exact<7>(0, 0, current); break;
            default:
                throw std::runtime_error("goal06 zero anchor: target above six");
        }
    }

    void run() {
        for (unsigned set = 0; set < input.rows.size() && !counterexample; ++set) {
            active_rows = &input.rows[set];
            for (unsigned target = 0;
                 target < input.base_weight && !counterexample;
                 ++target) {
                enumerate_weight(target);
            }
            std::cerr << "small_set=" << (set + 1)
                      << " candidates=" << candidates << '\n';
        }
        for (unsigned set = 0;
             set < input.base_weight_set_count && !counterexample;
             ++set) {
            active_rows = &input.rows[set];
            enumerate_weight(input.base_weight);
            std::cerr << "base_weight_set=" << (set + 1)
                      << " candidates=" << candidates << '\n';
        }
    }
};

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 3) {
            std::cerr << "usage: certificate INPUT.bin OUTPUT.json\n";
            return 2;
        }
        const auto started = std::chrono::steady_clock::now();
        const Input input = load_input(argv[1]);
        Search search{input};
        search.run();
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        if (!search.counterexample && search.candidates != input.expected_candidates) {
            throw std::runtime_error("goal06 zero anchor: candidate count changed");
        }

        std::ofstream output(argv[2]);
        if (!output) {
            throw std::runtime_error("goal06 zero anchor: cannot open output");
        }
        output << "{\n"
               << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-g4-goal06-zero-anchor-raw-v1\",\n"
               << "  \"result\": \""
               << (search.counterexample ? "COUNTEREXAMPLE" : "EXHAUSTED") << "\",\n"
               << "  \"candidates\": " << search.candidates << ",\n"
               << "  \"expected_candidates\": " << input.expected_candidates << ",\n"
               << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsed;
        if (search.counterexample) {
            output << ",\n  \"counterexample\": {\n"
                   << "    \"weight\": " << search.counterexample_weight << ",\n"
                   << "    \"initial_a_hex\": \"0x" << std::hex
                   << std::setw(16) << std::setfill('0') << search.counterexample_state[0]
                   << "\",\n"
                   << "    \"initial_b_hex\": \"0x" << std::setw(16)
                   << search.counterexample_state[1] << "\"\n"
                   << "  }";
        }
        output << "\n}\n";
        std::cout << "result="
                  << (search.counterexample ? "COUNTEREXAMPLE" : "EXHAUSTED") << '\n'
                  << "candidates=" << std::dec << search.candidates << '\n'
                  << "elapsed_seconds=" << elapsed << '\n'
                  << "output=" << argv[2] << '\n';
        return search.counterexample ? 1 : 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 2;
    }
}
