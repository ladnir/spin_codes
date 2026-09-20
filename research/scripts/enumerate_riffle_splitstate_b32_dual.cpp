#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

namespace {

struct Word128 {
    std::uint64_t low{};
    std::uint64_t high{};
};

Word128 operator^(Word128 left, const Word128 right) {
    left.low ^= right.low;
    left.high ^= right.high;
    return left;
}

Word128& operator^=(Word128& left, const Word128 right) {
    left.low ^= right.low;
    left.high ^= right.high;
    return left;
}

std::uint32_t parseHex32(const std::string& value) {
    return static_cast<std::uint32_t>(std::stoul(value, nullptr, 16));
}

} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: enumerate_b32_dual columns.txt histogram.csv\n";
        return 2;
    }

    std::ifstream input(argv[1]);
    if (!input) {
        throw std::runtime_error("failed to open column input");
    }
    std::vector<std::uint32_t> columns;
    std::string token;
    while (input >> token) {
        columns.push_back(parseHex32(token));
    }
    if (columns.size() != 128) {
        throw std::runtime_error("column input must contain 128 words");
    }

    std::array<Word128, 32> rows{};
    for (std::size_t column = 0; column < columns.size(); ++column) {
        for (unsigned row = 0; row < 32; ++row) {
            if ((columns[column] >> row) & 1U) {
                if (column < 64) {
                    rows[row].low |= std::uint64_t{1} << column;
                } else {
                    rows[row].high |= std::uint64_t{1} << (column - 64);
                }
            }
        }
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
                    current ^= rows[lowBits + bit];
                }
            }
            auto observe = [&]() {
                const unsigned weight = std::popcount(current.low)
                    + std::popcount(current.high);
                ++local[weight];
            };
            observe();
            for (std::uint32_t index = 1; index < suffixCount; ++index) {
                current ^= rows[std::countr_zero(index)];
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
        throw std::runtime_error("dual histogram mass mismatch");
    }

    std::ofstream output(argv[2]);
    if (!output) {
        throw std::runtime_error("failed to open histogram output");
    }
    output << "weight,count\n";
    for (unsigned weight = 0; weight <= 128; ++weight) {
        if (histogram[weight] != 0) {
            output << weight << ',' << histogram[weight] << '\n';
        }
    }
    std::cout << "threads," << omp_get_max_threads() << '\n';
    std::cout << "mass," << mass << '\n';
    std::cout << "wrote," << argv[2] << '\n';
    return 0;
}

