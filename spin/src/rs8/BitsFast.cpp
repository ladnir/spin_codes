// AVX-512 + GFNI bit-vector transpose (see Bits.h). Same values as
// transposePortable; spin_bits_test checks both against the record transpose.
#include "Bits.h"
#include <immintrin.h>
#include <cstring>

namespace spin::detail::rs8::bits {
namespace {
// Inner. The 22-bit state lives in bytes 0..2 of an XMM. Per epoch, R'_e s and
// C' s are three GFNI 8x8-block products each (one per state byte), reduced
// with a ternary XOR; vpmovqb takes output byte o from qword lane o.
// plan.updates packs block (o, k) at qword 3o+k with row 8o+j in byte j; GFNI
// wants row j in byte 7-j, hence the per-qword byte reversal.
void inner(const Plan& p, const std::uint64_t* x, std::uint64_t* y) {
    const auto& k = constants();
    const auto* up = p.certified22Plan->updates.data();
    const std::size_t epochs = p.n / 64;
    alignas(64) std::uint8_t revBytes[64];
    for (unsigned i = 0; i < 64; ++i) revBytes[i] = std::uint8_t((i & 8U) | (7U - (i & 7U))); // lane-local
    const __m512i rev = _mm512_load_si512(revBytes);
    const __m512i i0 = _mm512_set_epi64(0, 0, 0, 0, 0, 6, 3, 0), i1 = _mm512_set_epi64(0, 0, 0, 0, 0, 7, 4, 1),
                  i2 = _mm512_set_epi64(0, 0, 0, 0, 0, 8, 5, 2);
    const __m512i gc0 = _mm512_loadu_si512(k.gc[0]), gc1 = _mm512_loadu_si512(k.gc[1]), gc2 = _mm512_loadu_si512(k.gc[2]);
    __m128i sv = _mm_setzero_si128();
    for (std::size_t step = 0; step < epochs; ++step) {
        const auto e = epochs - 1 - step;
        const auto xe = x[e];
        const __m128i s1 = _mm_srli_si128(sv, 1), s2 = _mm_srli_si128(sv, 2);
        const __m512i z0 = _mm512_broadcastb_epi8(sv), z1 = _mm512_broadcastb_epi8(s1), z2 = _mm512_broadcastb_epi8(s2);
        const __m512i cy = _mm512_ternarylogic_epi64(_mm512_gf2p8affine_epi64_epi8(z0, gc0, 0),
            _mm512_gf2p8affine_epi64_epi8(z1, gc1, 0), _mm512_gf2p8affine_epi64_epi8(z2, gc2, 0), 0x96);
        y[e] = xe ^ std::uint64_t(_mm_cvtsi128_si64(_mm512_cvtepi64_epi8(cy)));
        if (!e) break;
        std::uint32_t fb = 0;
        if (xe)
            fb = k.ft[0][xe & 255U] ^ k.ft[1][(xe >> 8) & 255U] ^ k.ft[2][(xe >> 16) & 255U] ^ k.ft[3][(xe >> 24) & 255U]
               ^ k.ft[4][(xe >> 32) & 255U] ^ k.ft[5][(xe >> 40) & 255U] ^ k.ft[6][(xe >> 48) & 255U] ^ k.ft[7][xe >> 56];
        const __m512i u = _mm512_shuffle_epi8(_mm512_loadu_si512(up + 9 * e), rev);
        const __m512i u8 = _mm512_shuffle_epi8(_mm512_set1_epi64((long long)up[9 * e + 8]), rev);
        const __m256i m0 = _mm512_castsi512_si256(_mm512_maskz_permutex2var_epi64(7, u, i0, u8));
        const __m256i m1 = _mm512_castsi512_si256(_mm512_maskz_permutex2var_epi64(7, u, i1, u8));
        const __m256i m2 = _mm512_castsi512_si256(_mm512_maskz_permutex2var_epi64(7, u, i2, u8));
        const __m256i b0 = _mm256_broadcastb_epi8(sv), b1 = _mm256_broadcastb_epi8(s1), b2 = _mm256_broadcastb_epi8(s2);
        const __m256i r = _mm256_ternarylogic_epi64(_mm256_gf2p8affine_epi64_epi8(b0, m0, 0),
            _mm256_gf2p8affine_epi64_epi8(b1, m1, 0), _mm256_gf2p8affine_epi64_epi8(b2, m2, 0), 0x96);
        sv = _mm_xor_si128(_mm256_cvtepi64_epi8(r), _mm_cvtsi32_si128(int(fb)));
    }
}

void route(const Plan& p, const Tables& tb, const std::uint64_t* y, std::uint32_t* w) {
    const auto* yb = reinterpret_cast<const std::uint8_t*>(y);
    const auto* sb = tb.src.data();
    for (std::size_t g = 0; g < p.groups; ++g) {
        std::uint64_t acc[4];
        for (unsigned q = 0; q < 4; ++q) {
            std::uint64_t a0 = 0, a1 = 0;
            for (unsigned t = 0; t < 16; t += 2) {
                const auto j0 = sb[64 * g + 16 * q + t], j1 = sb[64 * g + 16 * q + t + 1];
                a0 |= std::uint64_t((yb[j0 >> 3] >> (j0 & 4U)) & 15U) << (4 * t);
                a1 |= std::uint64_t((yb[j1 >> 3] >> (j1 & 4U)) & 15U) << (4 * t + 4);
            }
            acc[q] = a0 | a1;
        }
        std::memcpy(w + 8 * g, acc, 32);
    }
}

// Tower multiplication, 16 symbols per ZMM. Per symbol (bytes b0..b3) with
// coefficient lanes F1 = f0..f3, F2 = f4..f7, F3 = f8:
//   X1 = (b0, b1, b0^b1, b2), X2 = (b3, b2^b3, b0^b2, b1^b3), X3 = b0^b1^b2^b3
//   P1 = X1 F1 = (p00,p01,p02,p10), P2 = X2 F2 = (p11,p12,p20,p21), P3 = p22
//   out = (p00^p01^p10^p11, p00^p02^p10^p12, p00^p01^p20^p21, p00^p02^p20^p22)
__m512i dwordPattern(int a, int b, int c, int d) {
    alignas(64) std::uint8_t m[64];
    for (int l = 0; l < 16; ++l) {
        const int v[4] = {a, b, c, d};
        for (int i = 0; i < 4; ++i) m[4 * l + i] = std::uint8_t(v[i] < 0 ? 0x80 : v[i] + 4 * (l % 4));
    }
    return _mm512_load_si512(m);
}
void symbols(const Plan& p, const Tables& tb, const std::uint32_t* w, std::uint32_t* mixed) {
    const __m512i s1a = dwordPattern(0, 1, 0, 2), s1b = dwordPattern(-1, -1, 1, -1);
    const __m512i s2a = dwordPattern(3, 2, 0, 1), s2b = dwordPattern(-1, 3, 2, 3);
    const __m512i s3a = dwordPattern(2, -1, -1, -1), s3b = dwordPattern(3, -1, -1, -1);
    const __m512i v1 = dwordPattern(0, 0, 0, 0), v2 = dwordPattern(3, 3, -1, -1), v3 = dwordPattern(-1, -1, 2, 2);
    const __m512i v4 = dwordPattern(1, 2, 1, 2), v5 = dwordPattern(0, 1, 3, -1), v6 = dwordPattern(-1, -1, -1, 0);
    const std::size_t count = 8 * p.groups;
    for (std::size_t i = 0; i < count; i += 16) {
        const __m512i W = _mm512_loadu_si512(w + i);
        const __m512i X1 = _mm512_xor_si512(_mm512_shuffle_epi8(W, s1a), _mm512_shuffle_epi8(W, s1b));
        const __m512i X2 = _mm512_xor_si512(_mm512_shuffle_epi8(W, s2a), _mm512_shuffle_epi8(W, s2b));
        const __m512i X3 = _mm512_xor_si512(_mm512_shuffle_epi8(X2, s3a), _mm512_shuffle_epi8(X2, s3b));
        const __m512i P1 = _mm512_gf2p8mul_epi8(X1, _mm512_loadu_si512(tb.f1.data() + i));
        const __m512i P2 = _mm512_gf2p8mul_epi8(X2, _mm512_loadu_si512(tb.f2.data() + i));
        const __m512i P3 = _mm512_gf2p8mul_epi8(X3, _mm512_loadu_si512(tb.f3.data() + i));
        __m512i r = _mm512_ternarylogic_epi64(_mm512_shuffle_epi8(P1, v1), _mm512_shuffle_epi8(P1, v2),
                                              _mm512_shuffle_epi8(P2, v3), 0x96);
        r = _mm512_ternarylogic_epi64(r, _mm512_shuffle_epi8(P1, v4), _mm512_shuffle_epi8(P2, v5), 0x96);
        _mm512_storeu_si512(mixed + i, _mm512_xor_si512(r, _mm512_shuffle_epi8(P3, v6)));
    }
}

void rs(const Plan& p, const std::uint32_t* mixed, std::uint64_t* out) {
    const auto& k = constants();
    for (std::size_t g = 0; g < p.groups; ++g) {
        const auto* b = reinterpret_cast<const std::uint8_t*>(mixed + 8 * g);
        __m128i a[4] = {_mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128(), _mm_setzero_si128()};
        for (unsigned j = 0; j < 32; ++j)
            a[j & 3] = _mm_xor_si128(a[j & 3], _mm_loadu_si128(reinterpret_cast<const __m128i*>(k.rs[j][b[j]])));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(out + 2 * g),
                         _mm_xor_si128(_mm_xor_si128(a[0], a[1]), _mm_xor_si128(a[2], a[3])));
    }
}
}

void transposeFast(const std::uint64_t* x, std::uint64_t* out, const Plan& p, const Tables& tb, Scratch& s) {
    inner(p, x, s.y.data());
    route(p, tb, s.y.data(), s.w.data());
    symbols(p, tb, s.w.data(), s.mixed.data());
    rs(p, s.mixed.data(), out);
}
}
