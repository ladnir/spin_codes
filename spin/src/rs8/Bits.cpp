#include "Bits.h"
#include "../certified22/productionMap.h"
#include <array>
#include <bit>
#include <cstring>

namespace spin::detail::rs8::bits {
namespace {
namespace map = certified22::map;

// Row-mask matrix times a 22-bit vector: bit i = parity(m[i] & v).
std::uint32_t apply(const std::array<std::uint32_t, 22>& m, std::uint32_t v) {
    std::uint32_t r = 0;
    for (unsigned i = 0; i < 22; ++i) r |= std::uint32_t(std::popcount(m[i] & v) & 1) << i;
    return r;
}
// GF(16) RS rows exactly as SetupScalar.cpp builds them.
std::array<std::uint32_t, 16> rsRows() {
    auto mul = [](unsigned a, unsigned b) {
        unsigned x = 0;
        for (unsigned i = 0; i < 4; ++i) { if (b & 1U) x ^= a; b >>= 1; a <<= 1; if (a & 16U) a ^= 0x13U; }
        return x;
    };
    auto inv = [&](unsigned a) { unsigned x = 1; for (unsigned i = 0; i < 14; ++i) x = mul(x, a); return x; };
    std::array<std::uint32_t, 16> r{};
    for (unsigned mb = 0; mb < 16; ++mb)
        for (unsigned sym = 0; sym < 8; ++sym) {
            unsigned num = 1, den = 1;
            for (unsigned o = 0; o < 4; ++o)
                if (o != mb / 4) { num = mul(num, sym ^ o); den = mul(den, (mb / 4) ^ o); }
            r[mb] |= std::uint32_t(mul(mul(num, inv(den)), 1U << (mb % 4))) << (4 * sym);
        }
    return r;
}
Constants buildConstants() {
    Constants c{};
    // C' row j (a mask over the fast-basis state s' = basis s): C[j] * inverseBasis.
    std::array<std::uint32_t, 64> cRows{};
    for (unsigned j = 0; j < 64; ++j)
        for (auto m = map::oldColumns[j]; m; m &= m - 1) cRows[j] ^= map::inverseBasis[std::countr_zero(m)];
    for (unsigned t = 0; t < 3; ++t)
        for (unsigned b = 0; b < 256; ++b)
            for (unsigned j = 0; j < 64; ++j)
                if (std::popcount(cRows[j] & (b << (8 * t))) & 1) c.ct[t][b] ^= std::uint64_t(1) << j;
    for (unsigned t = 0; t < 8; ++t)
        for (unsigned b = 0; b < 256; ++b) {
            std::uint32_t f = 0;
            for (unsigned i = 0; i < 8; ++i)
                if ((b >> i) & 1U) f ^= map::oldColumns[8 * t + i];
            c.ft[t][b] = apply(map::basis, f);
        }
    // GFNI: result bit i of output byte o is parity(block byte (7-i) & input byte k).
    for (unsigned k = 0; k < 3; ++k)
        for (unsigned o = 0; o < 8; ++o) {
            std::uint64_t q = 0;
            for (unsigned i = 0; i < 8; ++i) q |= std::uint64_t((cRows[8 * o + i] >> (8 * k)) & 255U) << (8 * (7 - i));
            c.gc[k][o] = q;
        }
    // Mixed bit 32*sym + 4*row + t feeds output bit 16*row + m when RS row m has bit 4*sym + t.
    const auto rows = rsRows();
    for (unsigned byte = 0; byte < 32; ++byte)
        for (unsigned v = 0; v < 256; ++v) {
            std::uint64_t out[2] = {0, 0};
            for (unsigned b = 0; b < 8; ++b) {
                if (!((v >> b) & 1U)) continue;
                const unsigned bit = 8 * byte + b, sym = bit / 32, row = (bit % 32) / 4, t = bit % 4;
                for (unsigned m = 0; m < 16; ++m)
                    if ((rows[m] >> (4 * sym + t)) & 1U) out[(16 * row + m) / 64] ^= std::uint64_t(1) << ((16 * row + m) % 64);
            }
            c.rs[byte][v][0] = out[0]; c.rs[byte][v][1] = out[1];
        }
    // GF(2^8) mod 0x11b, generator 3.
    unsigned x = 1;
    for (unsigned i = 0; i < 255; ++i) {
        c.exp[i] = c.exp[i + 255] = std::uint8_t(x);
        c.log[x] = std::uint8_t(i);
        x ^= x << 1; if (x & 256U) x ^= 0x11bU;
    }
    return c;
}
}

const Constants& constants() {
    static const Constants c = buildConstants();
    return c;
}

void Tables::build(const Plan& p, bool useFast) {
    fast = useFast;
    const auto& c = *p.certified22Plan;
    const std::size_t slots = 64 * p.groups;
    src.assign(slots, 0);
    for (std::size_t i = 0; i < p.n / 4; ++i) {
        const std::uint32_t off = c.routeBytes[i] / 16, g = off / std::uint32_t(groupStride),
                            slot = (off % std::uint32_t(groupStride)) / 4;
        // Packet i is nibble i%16 of y word i/16: byte (i/16)*8 + (i%16)/2 = i/2.
        src[64 * std::size_t(g) + slot] = std::uint32_t((i >> 1) << 3) | std::uint32_t((i & 1) << 2);
    }
    const std::size_t symbols = 8 * p.groups;
    if (fast) {
        f1.assign(symbols + 16, 0); f2.assign(symbols + 16, 0); f3.assign(symbols + 16, 0);
        for (std::size_t s = 0; s < symbols; ++s) {
            const auto& f = p.field[s];
            f1[s] = f[0] | std::uint32_t(f[1]) << 8 | std::uint32_t(f[2]) << 16 | std::uint32_t(f[3]) << 24;
            f2[s] = f[4] | std::uint32_t(f[5]) << 8 | std::uint32_t(f[6]) << 16 | std::uint32_t(f[7]) << 24;
            f3[s] = f[8];
        }
        cols.clear();
    } else {
        // R'_e rows from the packed blocks: block (o, k) byte j holds row 8o+j's input byte k.
        const auto epochs = p.n / 64;
        cols.assign(22 * epochs, 0);
        for (std::size_t e = 0; e < epochs; ++e)
            for (unsigned o = 0; o < 3; ++o)
                for (unsigned k = 0; k < 3; ++k) {
                    const auto q = c.updates[9 * e + 3 * o + k];
                    for (unsigned j = 0; j < 8 && 8 * o + j < 22; ++j) {
                        const auto byte = unsigned(q >> (8 * j)) & 255U;
                        for (unsigned b = 0; b < 8; ++b)
                            if ((byte >> b) & 1U) cols[22 * e + 8 * k + b] |= 1U << (8 * o + j);
                    }
                }
        f1.clear(); f2.clear(); f3.clear();
    }
}
std::size_t Tables::bytes() const noexcept {
    return (src.capacity() + f1.capacity() + f2.capacity() + f3.capacity() + cols.capacity()) * 4;
}

void transposePortable(const std::uint64_t* x, std::uint64_t* out, const Plan& p, const Tables& tb, Scratch& s) {
    const auto& k = constants();
    const std::size_t epochs = p.n / 64;
    std::uint64_t* y = s.y.data();
    std::uint32_t st = 0;
    for (std::size_t step = 0; step < epochs; ++step) {
        const auto e = epochs - 1 - step;
        const auto xe = x[e];
        y[e] = xe ^ k.ct[0][st & 255U] ^ k.ct[1][(st >> 8) & 255U] ^ k.ct[2][st >> 16];
        if (!e) break;
        std::uint32_t fb = 0;
        if (xe)
            for (unsigned t = 0; t < 8; ++t) fb ^= k.ft[t][(xe >> (8 * t)) & 255U];
        const auto* col = tb.cols.data() + 22 * e;
        std::uint32_t a[22];
        for (unsigned c = 0; c < 22; ++c) a[c] = col[c] & (0U - ((st >> c) & 1U));
        std::uint32_t r = fb;
        for (unsigned c = 0; c < 22; c += 2) r ^= a[c] ^ a[c + 1];
        st = r;
    }
    const auto* yb = reinterpret_cast<const std::uint8_t*>(y);
    for (std::size_t g = 0; g < p.groups; ++g) {
        std::uint64_t acc[4];
        for (unsigned q = 0; q < 4; ++q) {
            std::uint64_t a = 0;
            for (unsigned t = 0; t < 16; ++t) {
                const auto i = tb.src[64 * g + 16 * q + t];
                a |= std::uint64_t((yb[i >> 3] >> (i & 4U)) & 15U) << (4 * t);
            }
            acc[q] = a;
        }
        std::memcpy(s.w.data() + 8 * g, acc, 32);
    }
    // Multiplication by the symbol scalar in the tower basis: with the nine
    // coefficients f, products p_k0 = x0 f[3k], p_k1 = x1 f[3k+1],
    // p_k2 = (x0^x1) f[3k+2] for the 16-bit halves X0, X1, X0^X1.
    auto mul = [&](std::uint8_t a, std::uint8_t b) -> std::uint8_t {
        return (a && b) ? k.exp[k.log[a] + k.log[b]] : 0;
    };
    for (std::size_t sym = 0; sym < 8 * p.groups; ++sym) {
        const auto wv = s.w[sym];
        const auto& f = p.field[sym];
        const std::uint8_t b0 = std::uint8_t(wv), b1 = std::uint8_t(wv >> 8), b2 = std::uint8_t(wv >> 16), b3 = std::uint8_t(wv >> 24);
        const std::uint8_t X[3][2] = {{b0, b1}, {b2, b3}, {std::uint8_t(b0 ^ b2), std::uint8_t(b1 ^ b3)}};
        std::uint8_t lo[3], hi[3];
        for (unsigned t = 0; t < 3; ++t) {
            const auto p0 = mul(X[t][0], f[3 * t]), p1 = mul(X[t][1], f[3 * t + 1]),
                       p2 = mul(std::uint8_t(X[t][0] ^ X[t][1]), f[3 * t + 2]);
            lo[t] = std::uint8_t(p0 ^ p1); hi[t] = std::uint8_t(p0 ^ p2);
        }
        s.mixed[sym] = std::uint32_t(lo[0] ^ lo[1]) | std::uint32_t(hi[0] ^ hi[1]) << 8 |
                       std::uint32_t(lo[0] ^ lo[2]) << 16 | std::uint32_t(hi[0] ^ hi[2]) << 24;
    }
    for (std::size_t g = 0; g < p.groups; ++g) {
        const auto* b = reinterpret_cast<const std::uint8_t*>(s.mixed.data() + 8 * g);
        std::uint64_t o0 = 0, o1 = 0;
        for (unsigned j = 0; j < 32; ++j) { o0 ^= k.rs[j][b[j]][0]; o1 ^= k.rs[j][b[j]][1]; }
        out[2 * g] = o0; out[2 * g + 1] = o1;
    }
}
}
