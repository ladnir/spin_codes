#include <array>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
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

static inline Word256 xor2(const Word256& a, const Word256& b) noexcept {
    return {{{a.x[0] ^ b.x[0], a.x[1] ^ b.x[1],
              a.x[2] ^ b.x[2], a.x[3] ^ b.x[3]}}};
}

static inline unsigned weight(const Word256& a) noexcept {
    return std::popcount(a.x[0]) + std::popcount(a.x[1]) +
           std::popcount(a.x[2]) + std::popcount(a.x[3]);
}

static Word256 permute_to_standard(
    const Word256& word, unsigned exponent, unsigned total_length) {
    Word256 result{};
    for (unsigned limb = 0; limb != 4; ++limb) {
        std::uint64_t bits = word.x[limb];
        while (bits) {
            const unsigned offset = std::countr_zero(bits);
            const unsigned source = 64 * limb + offset;
            if (source >= total_length) throw std::runtime_error("nonzero padding bit");
            const unsigned punctured_length = total_length - 1;
            const unsigned target = source == punctured_length
                ? punctured_length : (exponent * source) % punctured_length;
            result.x[target >> 6] |= std::uint64_t{1} << (target & 63);
            bits &= bits - 1;
        }
    }
    return result;
}

template <class T>
static void read_exact(std::ifstream& stream, T& value) {
    stream.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!stream) throw std::runtime_error("truncated basis file");
}

struct Basis {
    std::uint32_t exponent = 0;
    std::vector<Word256> rows;
};

static std::vector<Basis> read_bases(
    const std::string& path, std::uint32_t& total_length, std::uint32_t& dimension) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot open basis file: " + path);
    std::array<char, 8> magic{};
    stream.read(magic.data(), magic.size());
    if (std::string(magic.data(), magic.size()) != "BCHBAS1\n")
        throw std::runtime_error("bad basis-file magic");
    std::uint32_t variants = 0;
    read_exact(stream, total_length);
    read_exact(stream, dimension);
    read_exact(stream, variants);
    if (total_length < 8 || total_length > 256 || dimension == 0 ||
        dimension >= total_length || variants == 0)
        throw std::runtime_error("unexpected basis-file dimensions");

    std::vector<Basis> bases(variants);
    for (Basis& basis : bases) {
        read_exact(stream, basis.exponent);
        basis.rows.resize(dimension);
        for (Word256& row : basis.rows) read_exact(stream, row);
    }
    return bases;
}

struct Search {
    unsigned minimum_weight;
    unsigned maximum_weight;
    unsigned total_length;
    std::array<std::uint64_t, 257> raw_hits{};
    std::unordered_set<Word256, WordHash> unique_hits;

    inline void inspect(const Word256& candidate, unsigned exponent) {
        const unsigned w = weight(candidate);
        if (w < minimum_weight || w > maximum_weight || (w & 1)) return;
        ++raw_hits[w];
        unique_hits.insert(permute_to_standard(candidate, exponent, total_length));
    }
};

static std::uint64_t choose(std::uint64_t n, unsigned k) {
    std::uint64_t result = 1;
    for (unsigned i = 1; i <= k; ++i) result = result * (n + 1 - i) / i;
    return result;
}

static void scan_basis(const Basis& basis, unsigned max_rows, Search& search) {
    const auto& r = basis.rows;
    const unsigned n = static_cast<unsigned>(r.size());
    if (max_rows >= 1) {
        for (unsigned i = 0; i < n; ++i) search.inspect(r[i], basis.exponent);
    }
    if (max_rows >= 2) {
        for (unsigned i = 0; i + 1 < n; ++i)
            for (unsigned j = i + 1; j < n; ++j)
                search.inspect(xor2(r[i], r[j]), basis.exponent);
    }
    if (max_rows >= 3) {
        for (unsigned i = 0; i + 2 < n; ++i)
            for (unsigned j = i + 1; j + 1 < n; ++j) {
                const Word256 pair = xor2(r[i], r[j]);
                for (unsigned k = j + 1; k < n; ++k)
                    search.inspect(xor2(pair, r[k]), basis.exponent);
            }
    }
    if (max_rows >= 4) {
        for (unsigned i = 0; i + 3 < n; ++i)
            for (unsigned j = i + 1; j + 2 < n; ++j) {
                const Word256 pair = xor2(r[i], r[j]);
                for (unsigned k = j + 1; k + 1 < n; ++k) {
                    const Word256 triple = xor2(pair, r[k]);
                    for (unsigned l = k + 1; l < n; ++l)
                        search.inspect(xor2(triple, r[l]), basis.exponent);
                }
            }
    }
    if (max_rows >= 5) {
        for (unsigned i = 0; i + 4 < n; ++i)
            for (unsigned j = i + 1; j + 3 < n; ++j) {
                const Word256 pair = xor2(r[i], r[j]);
                for (unsigned k = j + 1; k + 2 < n; ++k) {
                    const Word256 triple = xor2(pair, r[k]);
                    for (unsigned l = k + 1; l + 1 < n; ++l) {
                        const Word256 quadruple = xor2(triple, r[l]);
                        for (unsigned m = l + 1; m < n; ++m)
                            search.inspect(xor2(quadruple, r[m]), basis.exponent);
                    }
                }
            }
    }
}

int main(int argc, char** argv) try {
    if (argc != 6) {
        std::cerr << "usage: enumerate_wambach_rows BASES.bin HITS.tsv MAX_ROWS MIN_WEIGHT MAX_WEIGHT\n";
        return 2;
    }
    const std::string bases_path = argv[1];
    const std::string hits_path = argv[2];
    const unsigned max_rows = static_cast<unsigned>(std::stoul(argv[3]));
    const unsigned minimum_weight = static_cast<unsigned>(std::stoul(argv[4]));
    const unsigned maximum_weight = static_cast<unsigned>(std::stoul(argv[5]));
    if (max_rows < 1 || max_rows > 5 || minimum_weight > maximum_weight || maximum_weight > 256)
        throw std::runtime_error("invalid numeric argument");

    std::uint32_t total_length = 0, dimension = 0;
    const std::vector<Basis> bases = read_bases(bases_path, total_length, dimension);
    Search search{minimum_weight, maximum_weight, total_length};
    std::uint64_t combinations_per_basis = 0;
    for (unsigned rows = 1; rows <= max_rows; ++rows)
        combinations_per_basis += choose(dimension, rows);

    const auto start = std::chrono::steady_clock::now();
    for (std::size_t index = 0; index < bases.size(); ++index) {
        scan_basis(bases[index], max_rows, search);
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        std::cerr << "basis " << (index + 1) << '/' << bases.size()
                  << " exponent=" << bases[index].exponent
                  << " elapsed=" << std::fixed << std::setprecision(3) << seconds << "s\n";
    }

    std::ofstream output(hits_path);
    if (!output) throw std::runtime_error("cannot open hit file: " + hits_path);
    output << "weight\tlimb0\tlimb1\tlimb2\tlimb3\n" << std::hex << std::setfill('0');
    std::array<std::uint64_t, 257> unique_by_weight{};
    for (const Word256& word : search.unique_hits) {
        const unsigned w = weight(word);
        ++unique_by_weight[w];
        output << std::dec << w << std::hex;
        for (std::uint64_t limb : word.x) output << '\t' << std::setw(16) << limb;
        output << '\n';
    }

    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start).count();
    std::cout << "length=" << total_length
              << " bases=" << bases.size()
              << " dimension=" << dimension
              << " combinations=" << combinations_per_basis * bases.size()
              << " seconds=" << std::fixed << std::setprecision(6) << seconds << '\n';
    for (unsigned w = minimum_weight; w <= maximum_weight; ++w) {
        if (search.raw_hits[w] || unique_by_weight[w])
            std::cout << "weight=" << w << " raw_hits=" << search.raw_hits[w]
                      << " distinct_words=" << unique_by_weight[w] << '\n';
    }
    std::cout << "wrote=" << hits_path << '\n';
    return 0;
} catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
}
