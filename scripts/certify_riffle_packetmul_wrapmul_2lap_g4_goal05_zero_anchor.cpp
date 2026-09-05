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

constexpr std::uint32_t kSets = 14;
constexpr std::uint32_t kDimension = 40;
constexpr std::uint32_t kWords = 9;
constexpr std::uint32_t kThreshold = 105;

struct Row {
    std::array<std::uint64_t, kWords> word{};
    std::array<std::uint64_t, 2> state{};
};

struct Coset {
    std::uint16_t even_weight = 0;
    std::array<std::uint64_t, kWords> odd_word{};
    std::array<std::uint64_t, 2> state{};
};

struct Input {
    std::array<std::array<std::uint16_t, kDimension>, kSets> information_sets{};
    std::array<std::array<Row, kDimension>, kSets> rows{};
    std::vector<Coset> cosets;
};

template <class T>
void read_exact(std::ifstream& input, T* destination, std::size_t count = 1) {
    input.read(reinterpret_cast<char*>(destination),
               static_cast<std::streamsize>(sizeof(T) * count));
    if (!input) {
        throw std::runtime_error("zero-anchor certificate: truncated input");
    }
}

Input load_input(const std::string& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("zero-anchor certificate: cannot open input");
    }
    std::array<char, 8> magic{};
    read_exact(input, magic.data(), magic.size());
    if (std::string(magic.data(), magic.size()) != std::string("RZOC01\0\0", 8)) {
        throw std::runtime_error("zero-anchor certificate: bad input magic");
    }
    std::uint32_t sets = 0;
    std::uint32_t dimension = 0;
    std::uint32_t words = 0;
    std::uint32_t cosets = 0;
    std::uint32_t threshold = 0;
    read_exact(input, &sets);
    read_exact(input, &dimension);
    read_exact(input, &words);
    read_exact(input, &cosets);
    read_exact(input, &threshold);
    if (sets != kSets || dimension != kDimension || words != kWords ||
        threshold != kThreshold || cosets != 1025) {
        throw std::runtime_error("zero-anchor certificate: input dimensions changed");
    }

    Input result;
    for (auto& information_set : result.information_sets) {
        read_exact(input, information_set.data(), information_set.size());
    }
    for (auto& set_rows : result.rows) {
        for (auto& row : set_rows) {
            read_exact(input, row.word.data(), row.word.size());
            read_exact(input, row.state.data(), row.state.size());
        }
    }
    result.cosets.resize(cosets);
    for (auto& coset : result.cosets) {
        read_exact(input, &coset.even_weight);
        std::array<char, 6> padding{};
        read_exact(input, padding.data(), padding.size());
        read_exact(input, coset.odd_word.data(), coset.odd_word.size());
        read_exact(input, coset.state.data(), coset.state.size());
    }
    char extra = 0;
    if (input.read(&extra, 1)) {
        throw std::runtime_error("zero-anchor certificate: trailing input data");
    }
    return result;
}

struct Search {
    const Input& input;
    std::uint64_t candidates = 0;
    bool counterexample = false;
    std::uint16_t counterexample_even_weight = 0;
    std::uint16_t counterexample_odd_weight = 0;
    std::array<std::uint64_t, 2> counterexample_state{};

    const std::array<Row, kDimension>* active_rows = nullptr;
    std::array<std::uint64_t, kWords> base_word{};
    std::array<std::uint64_t, 2> base_state{};
    std::uint16_t active_even_weight = 0;
    std::uint16_t active_odd_limit = 0;

    static std::uint16_t weight(
        const std::array<std::uint64_t, kWords>& word,
        std::uint16_t limit) {
        std::uint16_t result = 0;
        for (std::uint32_t output_word = 0; output_word < kWords; ++output_word) {
            result = static_cast<std::uint16_t>(
                result + std::popcount(word[output_word]));
            if (result > limit) {
                break;
            }
        }
        return result;
    }

    void accept_if_nonzero(std::uint64_t selector, std::uint16_t odd_weight) {
        const auto& rows = *active_rows;
        std::array<std::uint64_t, 2> state = base_state;
        std::uint64_t remaining = selector;
        while (remaining != 0) {
            const unsigned bit = std::countr_zero(remaining);
            state[0] ^= rows[bit].state[0];
            state[1] ^= rows[bit].state[1];
            remaining &= remaining - 1;
        }
        if ((state[0] | state[1]) == 0) {
            return;
        }
        counterexample = true;
        counterexample_even_weight = active_even_weight;
        counterexample_odd_weight = odd_weight;
        counterexample_state = state;
    }

    void evaluate_zero() {
        ++candidates;
        const std::uint16_t odd_weight = weight(base_word, active_odd_limit);
        if (odd_weight <= active_odd_limit) {
            accept_if_nonzero(0, odd_weight);
        }
    }

    void evaluate_last(
        std::uint64_t selector,
        unsigned bit,
        const std::array<std::uint64_t, kWords>& current) {
        ++candidates;
        const auto& row = (*active_rows)[bit];
        std::uint16_t odd_weight = 0;
        for (unsigned output_word = 0; output_word < kWords; ++output_word) {
            odd_weight = static_cast<std::uint16_t>(
                odd_weight + std::popcount(current[output_word] ^ row.word[output_word]));
            if (odd_weight > active_odd_limit) {
                return;
            }
        }
        accept_if_nonzero(selector | (std::uint64_t{1} << bit), odd_weight);
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
                    bit + 1, selector | (std::uint64_t{1} << bit), current);
                for (unsigned word = 0; word < kWords; ++word) {
                    current[word] ^= (*active_rows)[bit].word[word];
                }
            }
        }
    }

    void enumerate_weight(unsigned target) {
        auto current = base_word;
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
                throw std::runtime_error("zero-anchor certificate: target above seven");
        }
    }

    void prepare_set(const Coset& coset, unsigned set) {
        active_rows = &input.rows[set];
        base_word = coset.odd_word;
        base_state = coset.state;
        std::uint64_t selected = 0;
        for (unsigned bit = 0; bit < kDimension; ++bit) {
            const unsigned coordinate = input.information_sets[set][bit];
            const unsigned word = coordinate / 64;
            const unsigned offset = coordinate % 64;
            selected |= ((coset.odd_word[word] >> offset) & 1) << bit;
        }
        std::uint64_t remaining = selected;
        while (remaining != 0) {
            const unsigned bit = std::countr_zero(remaining);
            for (unsigned word = 0; word < kWords; ++word) {
                base_word[word] ^= (*active_rows)[bit].word[word];
            }
            base_state[0] ^= (*active_rows)[bit].state[0];
            base_state[1] ^= (*active_rows)[bit].state[1];
            remaining &= remaining - 1;
        }
        for (unsigned bit = 0; bit < kDimension; ++bit) {
            const unsigned coordinate = input.information_sets[set][bit];
            if ((base_word[coordinate / 64] >> (coordinate % 64)) & 1) {
                throw std::runtime_error("zero-anchor certificate: coset zeroing failed");
            }
        }
    }

    void run() {
        for (std::size_t coset_index = 0;
             coset_index < input.cosets.size() && !counterexample;
             ++coset_index) {
            const Coset& coset = input.cosets[coset_index];
            active_even_weight = coset.even_weight;
            active_odd_limit = static_cast<std::uint16_t>(kThreshold - coset.even_weight);
            const unsigned base = active_odd_limit / kSets;
            const unsigned remainder = active_odd_limit % kSets;

            for (unsigned set = 0; set < kSets && !counterexample; ++set) {
                prepare_set(coset, set);
                for (unsigned target = 0; target < base && !counterexample; ++target) {
                    enumerate_weight(target);
                }
            }
            for (unsigned set = 0; set <= remainder && !counterexample; ++set) {
                prepare_set(coset, set);
                enumerate_weight(base);
            }
            if ((coset_index + 1) % 100 == 0) {
                std::cerr << "cosets=" << (coset_index + 1)
                          << " candidates=" << candidates << '\n';
            }
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

        std::ofstream output(argv[2]);
        if (!output) {
            throw std::runtime_error("zero-anchor certificate: cannot open output");
        }
        output << "{\n"
               << "  \"schema\": \"riffle-packetmul-wrapmul-2lap-g4-goal05-zero-anchor-certificate-raw-v1\",\n"
               << "  \"result\": \""
               << (search.counterexample ? "COUNTEREXAMPLE" : "EXHAUSTED") << "\",\n"
               << "  \"cosets\": " << input.cosets.size() << ",\n"
               << "  \"candidates\": " << search.candidates << ",\n"
               << "  \"elapsed_seconds\": " << std::setprecision(12) << elapsed;
        if (search.counterexample) {
            output << ",\n  \"counterexample\": {\n"
                   << "    \"even_weight\": " << search.counterexample_even_weight << ",\n"
                   << "    \"odd_weight\": " << search.counterexample_odd_weight << ",\n"
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
