#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
    constexpr std::size_t dimension = 32;
    constexpr std::size_t length = 64;
    constexpr std::size_t lowDimension = 20;
    constexpr std::size_t highDimension = dimension - lowDimension;

    std::array<std::uint64_t, dimension> readRows(const std::string& path)
    {
        std::ifstream input(path);
        if (!input)
            throw std::runtime_error("could not open generator-row file");
        std::array<std::uint64_t, dimension> rows{};
        for (auto& row : rows)
        {
            input >> std::hex >> row;
            if (!input)
                throw std::runtime_error("generator-row file is truncated");
        }
        std::string extra;
        if (input >> extra)
            throw std::runtime_error("generator-row file has extra rows");
        return rows;
    }
}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
        std::cerr << "usage: enumerate_binary64_code_spectrum ROWS OUT.csv\n";
        return 2;
    }
    const auto rows = readRows(argv[1]);
    int threadCount = 1;
#ifdef _OPENMP
    threadCount = omp_get_max_threads();
#endif
    std::vector<std::array<std::uint64_t, length + 1>> local(
        static_cast<std::size_t>(threadCount));

    constexpr std::uint64_t highCount = std::uint64_t{1} << highDimension;
    constexpr std::uint64_t lowCount = std::uint64_t{1} << lowDimension;
#pragma omp parallel for schedule(static)
    for (std::int64_t high = 0; high < static_cast<std::int64_t>(highCount); ++high)
    {
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        auto& histogram = local[static_cast<std::size_t>(thread)];
        std::uint64_t highWord = 0;
        for (std::size_t bit = 0; bit < highDimension; ++bit)
            if ((static_cast<std::uint64_t>(high) >> bit) & 1)
                highWord ^= rows[lowDimension + bit];

        std::uint64_t word = highWord;
        ++histogram[std::popcount(word)];
        for (std::uint64_t index = 1; index < lowCount; ++index)
        {
            word ^= rows[std::countr_zero(index)];
            ++histogram[std::popcount(word)];
        }
    }

    std::array<std::uint64_t, length + 1> histogram{};
    for (const auto& one : local)
        for (std::size_t weight = 0; weight <= length; ++weight)
            histogram[weight] += one[weight];
    std::uint64_t mass = 0;
    for (auto count : histogram)
        mass += count;
    if (mass != (std::uint64_t{1} << dimension))
        throw std::runtime_error("enumerated spectrum has the wrong mass");

    std::ofstream output(argv[2]);
    if (!output)
        throw std::runtime_error("could not open output file");
    output << "weight,count\n";
    for (std::size_t weight = 0; weight <= length; ++weight)
        if (histogram[weight])
            output << weight << ',' << histogram[weight] << '\n';
    std::cout << "threads," << threadCount << "\n";
    std::cout << "mass," << mass << "\n";
    std::cout << "minimum_weight,";
    for (std::size_t weight = 1; weight <= length; ++weight)
        if (histogram[weight])
        {
            std::cout << weight << "\n";
            break;
        }
    std::cout << "output," << argv[2] << "\n";
}
