// Generate the complete EBCH [128,64] weight-24 message shell from affine orbits.

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr unsigned kWeight = 24;
constexpr std::size_t kExpectedCount = 6855968;
constexpr std::uint64_t kDecoderInverse = 0x92a451cd307d0e63ULL;

struct Codeword {
    std::uint64_t low;
    std::uint64_t high;
    bool operator==(const Codeword&) const = default;
};

struct CodewordHash {
    std::size_t operator()(const Codeword& value) const noexcept {
        std::uint64_t x = value.low ^ std::rotl(value.high, 27);
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return static_cast<std::size_t>(x);
    }
};

constexpr std::array<Codeword, 64> kGeneratorRows{{
    {0xf4845518b9582a1fULL,0x8000000000000000ULL},{0xe908aa3172b0543eULL,0x8000000000000001ULL},
    {0xd2115462e560a87cULL,0x8000000000000003ULL},{0xa422a8c5cac150f8ULL,0x8000000000000007ULL},
    {0x4845518b9582a1f0ULL,0x800000000000000fULL},{0x908aa3172b0543e0ULL,0x800000000000001eULL},
    {0x2115462e560a87c0ULL,0x800000000000003dULL},{0x422a8c5cac150f80ULL,0x800000000000007aULL},
    {0x845518b9582a1f00ULL,0x80000000000000f4ULL},{0x08aa3172b0543e00ULL,0x80000000000001e9ULL},
    {0x115462e560a87c00ULL,0x80000000000003d2ULL},{0x22a8c5cac150f800ULL,0x80000000000007a4ULL},
    {0x45518b9582a1f000ULL,0x8000000000000f48ULL},{0x8aa3172b0543e000ULL,0x8000000000001e90ULL},
    {0x15462e560a87c000ULL,0x8000000000003d21ULL},{0x2a8c5cac150f8000ULL,0x8000000000007a42ULL},
    {0x5518b9582a1f0000ULL,0x800000000000f484ULL},{0xaa3172b0543e0000ULL,0x800000000001e908ULL},
    {0x5462e560a87c0000ULL,0x800000000003d211ULL},{0xa8c5cac150f80000ULL,0x800000000007a422ULL},
    {0x518b9582a1f00000ULL,0x80000000000f4845ULL},{0xa3172b0543e00000ULL,0x80000000001e908aULL},
    {0x462e560a87c00000ULL,0x80000000003d2115ULL},{0x8c5cac150f800000ULL,0x80000000007a422aULL},
    {0x18b9582a1f000000ULL,0x8000000000f48455ULL},{0x3172b0543e000000ULL,0x8000000001e908aaULL},
    {0x62e560a87c000000ULL,0x8000000003d21154ULL},{0xc5cac150f8000000ULL,0x8000000007a422a8ULL},
    {0x8b9582a1f0000000ULL,0x800000000f484551ULL},{0x172b0543e0000000ULL,0x800000001e908aa3ULL},
    {0x2e560a87c0000000ULL,0x800000003d211546ULL},{0x5cac150f80000000ULL,0x800000007a422a8cULL},
    {0xb9582a1f00000000ULL,0x80000000f4845518ULL},{0x72b0543e00000000ULL,0x80000001e908aa31ULL},
    {0xe560a87c00000000ULL,0x80000003d2115462ULL},{0xcac150f800000000ULL,0x80000007a422a8c5ULL},
    {0x9582a1f000000000ULL,0x8000000f4845518bULL},{0x2b0543e000000000ULL,0x8000001e908aa317ULL},
    {0x560a87c000000000ULL,0x8000003d2115462eULL},{0xac150f8000000000ULL,0x8000007a422a8c5cULL},
    {0x582a1f0000000000ULL,0x800000f4845518b9ULL},{0xb0543e0000000000ULL,0x800001e908aa3172ULL},
    {0x60a87c0000000000ULL,0x800003d2115462e5ULL},{0xc150f80000000000ULL,0x800007a422a8c5caULL},
    {0x82a1f00000000000ULL,0x80000f4845518b95ULL},{0x0543e00000000000ULL,0x80001e908aa3172bULL},
    {0x0a87c00000000000ULL,0x80003d2115462e56ULL},{0x150f800000000000ULL,0x80007a422a8c5cacULL},
    {0x2a1f000000000000ULL,0x8000f4845518b958ULL},{0x543e000000000000ULL,0x8001e908aa3172b0ULL},
    {0xa87c000000000000ULL,0x8003d2115462e560ULL},{0x50f8000000000000ULL,0x8007a422a8c5cac1ULL},
    {0xa1f0000000000000ULL,0x800f4845518b9582ULL},{0x43e0000000000000ULL,0x801e908aa3172b05ULL},
    {0x87c0000000000000ULL,0x803d2115462e560aULL},{0x0f80000000000000ULL,0x807a422a8c5cac15ULL},
    {0x1f00000000000000ULL,0x80f4845518b9582aULL},{0x3e00000000000000ULL,0x81e908aa3172b054ULL},
    {0x7c00000000000000ULL,0x83d2115462e560a8ULL},{0xf800000000000000ULL,0x87a422a8c5cac150ULL},
    {0xf000000000000000ULL,0x8f4845518b9582a1ULL},{0xe000000000000000ULL,0x9e908aa3172b0543ULL},
    {0xc000000000000000ULL,0xbd2115462e560a87ULL},{0x8000000000000000ULL,0xfa422a8c5cac150fULL},
}};

Codeword encode(std::uint64_t message) {
    Codeword result{0, 0};
    while (message != 0) {
        const unsigned bit = std::countr_zero(message);
        result.low ^= kGeneratorRows[bit].low;
        result.high ^= kGeneratorRows[bit].high;
        message &= message - 1;
    }
    return result;
}

unsigned codeword_weight(const Codeword value) {
    return std::popcount(value.low) + std::popcount(value.high);
}

std::uint8_t gf_multiply(std::uint8_t left, std::uint8_t right) {
    std::uint8_t result = 0;
    for (unsigned bit = 0; bit < 7; ++bit) {
        if ((right & 1U) != 0) {
            result ^= left;
        }
        right >>= 1;
        const bool top = (left & 0x40U) != 0;
        left = static_cast<std::uint8_t>((left << 1) & 0x7fU);
        if (top) {
            left ^= 0x03U;  // x^7+x+1
        }
    }
    return result;
}

std::array<std::uint8_t, 128> coordinate_to_point() {
    std::array<std::uint8_t, 128> result{};
    std::uint8_t point = 1;
    for (unsigned coordinate = 0; coordinate < 127; ++coordinate) {
        result[coordinate] = point;
        point = gf_multiply(point, 2);
    }
    result[127] = 0;
    return result;
}

std::array<std::uint8_t, 128> point_to_coordinate(
    const std::array<std::uint8_t, 128>& coordinates) {
    std::array<std::uint8_t, 128> result{};
    for (unsigned coordinate = 0; coordinate < 128; ++coordinate) {
        result[coordinates[coordinate]] = static_cast<std::uint8_t>(coordinate);
    }
    return result;
}

std::vector<std::uint8_t> support_points(
    Codeword codeword, const std::array<std::uint8_t, 128>& coordinates) {
    std::vector<std::uint8_t> points;
    points.reserve(kWeight);
    while (codeword.low != 0) {
        const unsigned coordinate = std::countr_zero(codeword.low);
        points.push_back(coordinates[coordinate]);
        codeword.low &= codeword.low - 1;
    }
    while (codeword.high != 0) {
        const unsigned bit = std::countr_zero(codeword.high);
        points.push_back(coordinates[64 + bit]);
        codeword.high &= codeword.high - 1;
    }
    if (points.size() != kWeight) {
        throw std::runtime_error("affine representative has the wrong weight");
    }
    return points;
}

std::vector<Codeword> affine_orbit(
    Codeword representative,
    const std::array<std::uint8_t, 128>& coordinates,
    const std::array<std::uint8_t, 128>& inverse_coordinates,
    const std::array<std::array<std::uint8_t, 128>, 128>& products) {
    const auto points = support_points(representative, coordinates);
    std::unordered_set<Codeword, CodewordHash> unique;
    unique.reserve(17000);
    for (unsigned multiplier = 1; multiplier < 128; ++multiplier) {
        for (unsigned translation = 0; translation < 128; ++translation) {
            Codeword image{0, 0};
            for (const auto point : points) {
                const unsigned coordinate =
                    inverse_coordinates[products[multiplier][point] ^ translation];
                if (coordinate < 64) {
                    image.low |= 1ULL << coordinate;
                } else {
                    image.high |= 1ULL << (coordinate - 64);
                }
            }
            unique.insert(image);
        }
    }
    return {unique.begin(), unique.end()};
}

std::uint64_t decode_low_half(std::uint64_t codeword_low) {
    std::uint64_t result = 0;
    std::uint64_t multiplier = kDecoderInverse;
    while (multiplier != 0) {
        const unsigned bit = std::countr_zero(multiplier);
        result ^= codeword_low << bit;
        multiplier &= multiplier - 1;
    }
    return result;
}

std::unordered_set<Codeword, CodewordHash> read_seed_codewords(
    const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw std::runtime_error("failed to open authenticated low-support records");
    }
    std::unordered_set<Codeword, CodewordHash> result;
    while (true) {
        std::uint8_t support = 0;
        std::uint64_t message = 0;
        input.read(reinterpret_cast<char*>(&support), 1);
        if (input.eof()) {
            break;
        }
        input.read(reinterpret_cast<char*>(&message), sizeof(message));
        if (!input) {
            throw std::runtime_error("authenticated record file has a partial record");
        }
        const Codeword codeword = encode(message);
        if (codeword_weight(codeword) == kWeight) {
            result.insert(codeword);
        }
    }
    return result;
}

}  // namespace

int main(int argc, char** argv) try {
    if (argc != 3) {
        std::cerr << "usage: enumerator authenticated_records.bin output_messages.bin\n";
        return 2;
    }
    auto uncovered = read_seed_codewords(argv[1]);
    const std::size_t authenticated_seeds = uncovered.size();
    const auto coordinates = coordinate_to_point();
    const auto inverse_coordinates = point_to_coordinate(coordinates);
    std::array<std::array<std::uint8_t, 128>, 128> products{};
    for (unsigned multiplier = 1; multiplier < 128; ++multiplier) {
        for (unsigned point = 0; point < 128; ++point) {
            products[multiplier][point] = gf_multiply(
                static_cast<std::uint8_t>(multiplier), static_cast<std::uint8_t>(point));
        }
    }

    std::vector<Codeword> complete;
    complete.reserve(kExpectedCount);
    std::vector<std::size_t> orbit_sizes;
    while (!uncovered.empty()) {
        const Codeword representative = *uncovered.begin();
        auto orbit = affine_orbit(
            representative, coordinates, inverse_coordinates, products);
        for (const auto codeword : orbit) {
            uncovered.erase(codeword);
        }
        orbit_sizes.push_back(orbit.size());
        complete.insert(complete.end(), orbit.begin(), orbit.end());
        if (orbit_sizes.size() % 64 == 0) {
            std::cerr << "orbits=" << orbit_sizes.size()
                      << " generated=" << complete.size()
                      << " uncovered_seeds=" << uncovered.size() << '\n';
        }
    }
    if (complete.size() != kExpectedCount) {
        throw std::runtime_error(
            "authenticated affine orbits do not cover committed A_24: generated "
            + std::to_string(complete.size()));
    }

    std::vector<std::uint64_t> messages;
    messages.reserve(complete.size());
    for (const auto codeword : complete) {
        const std::uint64_t message = decode_low_half(codeword.low);
        if (encode(message) != codeword) {
            throw std::runtime_error("affine image failed BCH re-encoding");
        }
        messages.push_back(message);
    }
    std::sort(messages.begin(), messages.end());
    if (std::adjacent_find(messages.begin(), messages.end()) != messages.end()) {
        throw std::runtime_error("generated weight-24 message set contains duplicates");
    }
    std::ofstream output(argv[2], std::ios::binary);
    output.write(
        reinterpret_cast<const char*>(messages.data()),
        static_cast<std::streamsize>(messages.size() * sizeof(std::uint64_t)));
    if (!output) {
        throw std::runtime_error("failed to write complete weight-24 message set");
    }

    std::sort(orbit_sizes.begin(), orbit_sizes.end());
    std::cout << "{\n"
              << "  \"schema\": \"ebch128-weight24-affine-enumeration-v1\",\n"
              << "  \"evidence_label\": \"EXACT_AFFINE_ORBIT_ENUMERATION_WITH_REENCODING\",\n"
              << "  \"authenticated_weight24_seeds\": " << authenticated_seeds << ",\n"
              << "  \"affine_orbits\": " << orbit_sizes.size() << ",\n"
              << "  \"generated_messages\": " << messages.size() << ",\n"
              << "  \"matches_committed_A_24\": true,\n"
              << "  \"all_images_reencoded\": true,\n"
              << "  \"orbit_size_histogram\": {";
    bool first = true;
    for (std::size_t begin = 0; begin < orbit_sizes.size();) {
        const auto size = orbit_sizes[begin];
        const auto end = std::upper_bound(orbit_sizes.begin() + begin, orbit_sizes.end(), size);
        if (!first) {
            std::cout << ',';
        }
        first = false;
        std::cout << "\n    \"" << size << "\": " << (end - (orbit_sizes.begin() + begin));
        begin = static_cast<std::size_t>(end - orbit_sizes.begin());
    }
    std::cout << "\n  },\n  \"output\": \"" << argv[2] << "\"\n}\n";
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
