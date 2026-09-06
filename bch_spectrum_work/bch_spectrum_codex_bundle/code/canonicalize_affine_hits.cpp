#include <array>
#include <bit>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

struct Word256 {
    std::array<std::uint64_t, 4> x{};
    bool operator==(const Word256& other) const noexcept { return x == other.x; }
};

struct WordHash {
    std::size_t operator()(const Word256& word) const noexcept {
        std::uint64_t h = 0x9e3779b97f4a7c15ULL;
        for (std::uint64_t value : word.x) {
            value ^= value >> 30;
            value *= 0xbf58476d1ce4e5b9ULL;
            value ^= value >> 27;
            value *= 0x94d049bb133111ebULL;
            value ^= value >> 31;
            h ^= value + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};

static bool numeric_less(const Word256& a, const Word256& b) noexcept {
    for (int limb = 3; limb >= 0; --limb) {
        if (a.x[limb] != b.x[limb]) return a.x[limb] < b.x[limb];
    }
    return false;
}

struct Field {
    unsigned q = 0;
    unsigned length = 0;
    unsigned m = 0;
    unsigned modulus = 0;
    std::array<std::array<std::uint8_t, 256>, 256> multiply{};
    std::array<std::uint8_t, 256> powers{};
    std::array<std::uint8_t, 256> logs{};
    std::array<std::uint8_t, 256> inverse{};

    explicit Field(unsigned total_length) {
        q = total_length;
        length = q - 1;
        if (q == 8) { m = 3; modulus = 0xB; }
        else if (q == 32) { m = 5; modulus = 0x25; }
        else if (q == 128) { m = 7; modulus = 0x83; }
        else if (q == 256) { m = 8; modulus = 0x14D; }
        else throw std::runtime_error("unsupported field size");

        for (unsigned a = 0; a < q; ++a) {
            for (unsigned b = 0; b < q; ++b) {
                unsigned left = a, right = b, result = 0;
                while (right) {
                    if (right & 1) result ^= left;
                    right >>= 1;
                    left <<= 1;
                    if (left & q) left ^= modulus;
                }
                multiply[a][b] = static_cast<std::uint8_t>(result & length);
            }
        }
        unsigned value = 1;
        for (unsigned exponent = 0; exponent < length; ++exponent) {
            powers[exponent] = static_cast<std::uint8_t>(value);
            logs[value] = static_cast<std::uint8_t>(exponent);
            value = multiply[value][2];
        }
        if (value != 1) throw std::runtime_error("field generator is not primitive");
        for (unsigned value2 = 1; value2 < q; ++value2) {
            for (unsigned candidate = 1; candidate < q; ++candidate) {
                if (multiply[value2][candidate] == 1) {
                    inverse[value2] = static_cast<std::uint8_t>(candidate);
                    break;
                }
            }
            if (!inverse[value2]) throw std::runtime_error("missing field inverse");
        }
    }
};

struct CanonicalResult {
    Word256 word{};
    unsigned stabilizer = 0;
};

static CanonicalResult canonicalize(const Word256& word, const Field& field) {
    std::array<std::uint8_t, 256> support{};
    unsigned support_size = 0;
    for (unsigned limb = 0; limb != 4; ++limb) {
        std::uint64_t bits = word.x[limb];
        while (bits) {
            const unsigned source = 64 * limb + std::countr_zero(bits);
            if (source >= field.q) throw std::runtime_error("nonzero padding bit");
            support[support_size++] = source == field.length
                ? 0 : field.powers[source];
            bits &= bits - 1;
        }
    }

    CanonicalResult result{};
    bool initialized = false;
    for (unsigned ix = 0; ix < support_size; ++ix) {
        const unsigned x = support[ix];
        for (unsigned iy = 0; iy < support_size; ++iy) {
            if (ix == iy) continue;
            const unsigned y = support[iy];
            const unsigned a = field.inverse[x ^ y];
            const unsigned b = field.multiply[a][x];
            Word256 image{};
            for (unsigned iz = 0; iz < support_size; ++iz) {
                const unsigned target = field.multiply[a][support[iz]] ^ b;
                const unsigned coordinate = target == 0 ? field.length : field.logs[target];
                image.x[coordinate >> 6] |= std::uint64_t{1} << (coordinate & 63);
            }
            if (!initialized || numeric_less(image, result.word)) {
                result.word = image;
                result.stabilizer = 1;
                initialized = true;
            } else if (image == result.word) {
                ++result.stabilizer;
            }
        }
    }
    if (!initialized || !result.stabilizer) throw std::runtime_error("empty support");
    return result;
}

struct OrbitRecord {
    std::uint64_t hits = 0;
    unsigned stabilizer = 0;
};

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

int main(int argc, char** argv) try {
    if (argc != 4) {
        std::cerr << "usage: canonicalize_affine_hits HITS.tsv ORBITS.tsv TOTAL_LENGTH\n";
        return 2;
    }
    const Field field(static_cast<unsigned>(std::stoul(argv[3])));
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open hit file");
    std::string line;
    std::getline(input, line);
    if (line != "weight\tlimb0\tlimb1\tlimb2\tlimb3")
        throw std::runtime_error("unexpected hit-file header");

    std::unordered_map<Word256, OrbitRecord, WordHash> orbits;
    std::uint64_t input_words = 0;
    unsigned shell_weight = 0;
    while (std::getline(input, line)) {
        if (line.empty()) continue;
        const auto fields = split_tabs(line);
        if (fields.size() != 5) throw std::runtime_error("bad hit row");
        const unsigned weight = static_cast<unsigned>(std::stoul(fields[0]));
        if (shell_weight && shell_weight != weight)
            throw std::runtime_error("canonicalizer expects one weight shell");
        shell_weight = weight;
        Word256 word{};
        for (unsigned limb = 0; limb != 4; ++limb)
            word.x[limb] = std::stoull(fields[limb + 1], nullptr, 16);
        const CanonicalResult canonical = canonicalize(word, field);
        OrbitRecord& record = orbits[canonical.word];
        if (record.hits && record.stabilizer != canonical.stabilizer)
            throw std::runtime_error("inconsistent stabilizer");
        ++record.hits;
        record.stabilizer = canonical.stabilizer;
        ++input_words;
    }

    std::ofstream output(argv[2]);
    if (!output) throw std::runtime_error("cannot open orbit file");
    output << "weight\tlimb0\tlimb1\tlimb2\tlimb3\torbit_size\tstabilizer_size\tsearch_hits\n";
    output << std::hex << std::setfill('0');
    std::uint64_t lower_bound = 0;
    for (const auto& [word, record] : orbits) {
        const unsigned orbit_size = field.q * field.length / record.stabilizer;
        lower_bound += orbit_size;
        output << std::dec << shell_weight << std::hex;
        for (std::uint64_t limb : word.x) output << '\t' << std::setw(16) << limb;
        output << std::dec << '\t' << orbit_size << '\t' << record.stabilizer
               << '\t' << record.hits << '\n';
    }
    std::cout << "weight=" << shell_weight << " input_words=" << input_words
              << " affine_orbits=" << orbits.size()
              << " rigorous_lower_bound=" << lower_bound << '\n';
    std::cout << "wrote=" << argv[2] << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
