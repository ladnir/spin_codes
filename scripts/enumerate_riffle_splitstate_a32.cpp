#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <omp.h>

namespace {

struct Word128 {
    std::uint64_t low{};
    std::uint64_t high{};
};

Word128& operator^=(Word128& left, const Word128 right) {
    left.low ^= right.low;
    left.high ^= right.high;
    return left;
}

Word128 parseWord(const std::string& value) {
    if (value.size() != 32) {
        throw std::runtime_error("generator word must have 32 hex digits");
    }
    return Word128{
        std::stoull(value.substr(16), nullptr, 16),
        std::stoull(value.substr(0, 16), nullptr, 16),
    };
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: enumerate_a32 generators.txt histogram.csv\n";
        return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) {
        throw std::runtime_error("failed to open generator input");
    }
    std::array<Word128, 32> generators{};
    std::string token;
    for (auto& generator : generators) {
        if (!(input >> token)) {
            throw std::runtime_error("generator input has fewer than 32 words");
        }
        generator = parseWord(token);
    }
    if (input >> token) {
        throw std::runtime_error("generator input has more than 32 words");
    }

    constexpr unsigned lowBits = 20;
    constexpr std::uint32_t prefixCount = std::uint32_t{1} << (32 - lowBits);
    constexpr std::uint32_t suffixCount = std::uint32_t{1} << lowBits;
    std::array<unsigned long long, 129> histogram{};

#pragma omp parallel
    {
        std::array<unsigned long long, 129> local{};
#pragma omp for schedule(static)
        for (std::uint32_t prefix = 0; prefix < prefixCount; ++prefix) {
            Word128 current{};
            for (unsigned bit = 0; bit < 32 - lowBits; ++bit) {
                if ((prefix >> bit) & 1U) {
                    current ^= generators[lowBits + bit];
                }
            }
            auto observe = [&]() {
                const unsigned weight = std::popcount(current.low)
                    + std::popcount(current.high);
                ++local[weight];
            };
            observe();
            for (std::uint32_t index = 1; index < suffixCount; ++index) {
                current ^= generators[std::countr_zero(index)];
                observe();
            }
        }
#pragma omp critical
        {
            for (unsigned weight = 0; weight <= 128; ++weight) {
                histogram[weight] += local[weight];
            }
        }
    }

    unsigned long long mass = 0;
    for (const auto count : histogram) {
        mass += count;
    }
    if (mass != (std::uint64_t{1} << 32)) {
        throw std::runtime_error("A histogram mass mismatch");
    }
    std::ofstream output(argv[2]);
    output << "weight,count\n";
    for (unsigned weight = 0; weight <= 128; ++weight) {
        if (histogram[weight]) {
            output << weight << ',' << histogram[weight] << '\n';
        }
    }
    std::cout << "threads," << omp_get_max_threads() << '\n';
    std::cout << "mass," << mass << '\n';
    std::cout << "wrote," << argv[2] << '\n';
    return 0;
}

