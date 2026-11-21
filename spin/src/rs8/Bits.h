#pragma once
#include "Rs8.h"
#include <cstdint>
#include <vector>

// Bit-vector transposed encoding: the Turbo transpose 2K -> K applied to one
// bit per position. It equals every bit lane of the 128-bit record transpose
// (records are 128 independent copies of the same binary map). Bit i is bit
// i%64 of word i/64. Phases, each the bit form of the record oracle:
//   inner  : reverse epochs, y = x ^ C s, s' = F x ^ R_e s (22-bit state, in
//            the fast-circuit basis of plan.updates, so banked plans work)
//   route  : each group's 64 nibble slots read their routed packet of y
//   symbols: scalar * w in the GF(2^32) tower, nine GF(2^8) products
//   rs     : the fixed RS[8,4] map, 256 -> 128 bits per group
namespace spin::detail::rs8::bits {
// Plan-independent tables, built once.
struct Constants {
    std::uint64_t ct[3][256];      // C' s' by state byte (C' = C * inverseBasis)
    std::uint32_t ft[8][256];      // F' x by input byte (F' = basis * F)
    std::uint64_t gc[3][8];        // GFNI 8x8 blocks of C': input byte k, output byte o
    std::uint64_t rs[32][256][2];  // RS map by mixed byte
    std::uint8_t log[256], exp[512];
};
const Constants& constants();

// Plan-dependent tables. Built from routeBytes, updates and field, which every
// setup path (and banked refresh) keeps current. Rebuilding with the same
// geometry reuses the allocations.
struct Tables {
    bool fast = false;
    std::vector<std::uint32_t> src;        // per slot g*64+c: y byte offset << 3 | 4*(nibble&1)
    std::vector<std::uint32_t> f1, f2, f3; // fast: tower coefficients per symbol, lane layout
    std::vector<std::uint32_t> cols;       // portable: columns of R'_e, epochs x 22
    void build(const Plan&, bool fast);
    std::size_t bytes() const noexcept;
};

// Per-call scratch (one per concurrent call).
struct Scratch {
    std::vector<std::uint64_t> y;
    std::vector<std::uint32_t> w, mixed;
    explicit Scratch(const Plan& p)
        : y(p.n / 64), w(8 * p.groups + 16), mixed(8 * p.groups + 16) {}
    std::size_t bytes() const noexcept { return y.size() * 8 + (w.size() + mixed.size()) * 4; }
};

void transposePortable(const std::uint64_t* in, std::uint64_t* out, const Plan&, const Tables&, Scratch&);
#if SPIN_AVX512
// Requires cpu_packet512(). Tables must be built with fast = true.
void transposeFast(const std::uint64_t* in, std::uint64_t* out, const Plan&, const Tables&, Scratch&);
#endif
}
