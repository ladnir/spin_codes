#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <immintrin.h>

namespace spin::research::char2::medium32 {

// These are two representations of GF(2^32), not four independent byte fields.
// Polynomial representation: GF(2)[X]/(X^32+X^7+X^3+X^2+1).
inline constexpr std::uint64_t polynomial = 0x10000008dULL;

constexpr std::uint64_t polynomialRemainder(std::uint64_t a, std::uint64_t b) noexcept {
    unsigned degreeB = 0;
    for (auto t = b; t >>= 1;) ++degreeB;
    for (unsigned i = 64; i-- > degreeB;)
        if ((a >> i) & 1U) a ^= b << (i - degreeB);
    return a;
}

constexpr std::uint64_t polynomialGcd(std::uint64_t a, std::uint64_t b) noexcept {
    while (b) { const auto r = polynomialRemainder(a, b); a = b; b = r; }
    return a;
}

// Independent carryless schoolbook product and polynomial long division.
constexpr std::uint32_t multiplyOracle(std::uint32_t a, std::uint32_t b) noexcept {
    std::uint64_t product = 0;
    for (unsigned i = 0; i != 32; ++i)
        if ((b >> i) & 1U) product ^= std::uint64_t(a) << i;
    return std::uint32_t(polynomialRemainder(product, polynomial));
}

constexpr std::uint32_t multiplyScalar(std::uint32_t a, std::uint32_t b) noexcept {
    std::uint32_t result = 0;
    for (unsigned i = 0; i != 32; ++i) {
        result ^= a & (0U - (b & 1U));
        a = (a << 1) ^ (0x8dU & (0U - (a >> 31)));
        b >>= 1;
    }
    return result;
}

constexpr std::uint32_t frobenius(unsigned iterations) noexcept {
    std::uint32_t x = 2;
    for (unsigned i = 0; i != iterations; ++i) x = multiplyOracle(x, x);
    return x;
}

// Rabin's criterion needs only prime divisor 2 of degree 32.
static_assert(frobenius(16) == 0x9127c27aU);
static_assert(polynomialGcd(frobenius(16) ^ 2U, polynomial) == 1);
static_assert(frobenius(32) == 2);

constexpr std::uint8_t byteMultiply(std::uint8_t a, std::uint8_t b) noexcept {
    std::uint8_t result = 0;
    for (unsigned i = 0; i != 8; ++i) {
        result ^= std::uint8_t(a & (0U - (b & 1U)));
        a = std::uint8_t((unsigned(a) << 1) ^ (0x1bU & (0U - (a >> 7))));
        b >>= 1;
    }
    return result;
}

constexpr std::uint8_t byteOracle(std::uint8_t a, std::uint8_t b) noexcept {
    std::uint16_t product = 0;
    for (unsigned i = 0; i != 8; ++i)
        if ((b >> i) & 1U) product ^= std::uint16_t(a) << i;
    return std::uint8_t(polynomialRemainder(product, 0x11b));
}

constexpr std::uint8_t bytePower(std::uint8_t a, unsigned n) noexcept {
    std::uint8_t result = 1;
    while (n) { if (n & 1U) result = byteMultiply(result, a); a = byteMultiply(a, a); n >>= 1; }
    return result;
}

constexpr std::uint8_t byteTrace(std::uint8_t a) noexcept {
    auto result = a;
    for (unsigned i = 1; i != 8; ++i) { a = byteMultiply(a, a); result ^= a; }
    return result;
}

// Tower representation, byte order [1,u,v,uv]:
// F = AES GF256, u^2+u+delta=0, v^2+v+delta*u=0, delta=0x20.
// Tr_F/GF2(delta)=1 proves the first quadratic irreducible. In F(u),
// Tr_F(u)/F(delta*u)=delta, so trace transitivity proves the second.
// AES-byte multiplication applied to all four bytes is EXACTLY multiplication
// by the embedded F scalar. Arbitrary 32-bit tower scalars mix all four bytes.
inline constexpr std::uint8_t delta = 0x20;
inline constexpr std::uint8_t delta2 = byteMultiply(delta, delta);
inline constexpr std::uint8_t delta3 = byteMultiply(delta2, delta);
inline constexpr std::uint8_t deltaInverse = bytePower(delta, 254);
static_assert(byteTrace(delta) == 1 && delta2 == 0x6c);
static_assert(byteMultiply(delta, deltaInverse) == 1);

constexpr std::uint16_t quadraticMultiply(std::uint16_t a, std::uint16_t b) noexcept {
    const auto a0 = std::uint8_t(a), a1 = std::uint8_t(a >> 8);
    const auto b0 = std::uint8_t(b), b1 = std::uint8_t(b >> 8);
    const auto p = byteMultiply(a0, b0), q = byteMultiply(a1, b1);
    const auto r = byteMultiply(a0 ^ a1, b0 ^ b1);
    return std::uint16_t((p ^ byteMultiply(delta, q)) | (unsigned(r ^ p) << 8));
}

constexpr std::uint32_t towerMultiplyScalar(std::uint32_t a, std::uint32_t b) noexcept {
    const auto a0 = std::uint16_t(a), a1 = std::uint16_t(a >> 16);
    const auto b0 = std::uint16_t(b), b1 = std::uint16_t(b >> 16);
    const auto p = quadraticMultiply(a0, b0), q = quadraticMultiply(a1, b1);
    const auto r = quadraticMultiply(a0 ^ a1, b0 ^ b1);
    return std::uint32_t(p ^ quadraticMultiply(std::uint16_t(delta) << 8, q)) |
           (std::uint32_t(r ^ p) << 16);
}

// Independent oracle: convert to powers of v, schoolbook multiply, reduce
// modulo v^4+(1+delta)*v^2+delta*v+delta^3, convert back. This quartic follows
// by eliminating u=(v^2+v)/delta from the first quadratic.
constexpr std::uint32_t towerMultiplyOracle(std::uint32_t a, std::uint32_t b) noexcept {
    const auto convert = [](std::uint32_t x) {
        const auto x0 = std::uint8_t(x), x1 = std::uint8_t(x >> 8);
        const auto x2 = std::uint8_t(x >> 16), x3 = std::uint8_t(x >> 24);
        return std::array<std::uint8_t, 4>{x0, std::uint8_t(x2 ^ byteOracle(x1, deltaInverse)),
            byteOracle(x1 ^ x3, deltaInverse), byteOracle(x3, deltaInverse)};
    };
    const auto av = convert(a), bv = convert(b);
    std::array<std::uint8_t, 7> product{};
    for (unsigned i = 0; i != 4; ++i)
        for (unsigned j = 0; j != 4; ++j) product[i + j] ^= byteOracle(av[i], bv[j]);
    for (unsigned i = 7; i-- > 4;) {
        const auto top = product[i];
        product[i - 4] ^= byteOracle(top, delta3);
        product[i - 3] ^= byteOracle(top, delta);
        product[i - 2] ^= byteOracle(top, 1 ^ delta);
    }
    return std::uint32_t(product[0]) |
           (std::uint32_t(byteOracle(product[2] ^ product[3], delta)) << 8) |
           (std::uint32_t(product[1] ^ product[2] ^ product[3]) << 16) |
           (std::uint32_t(byteOracle(product[3], delta)) << 24);
}

struct alignas(16) Record128 {
    std::array<std::uint32_t, 4> lane{};
    bool operator==(const Record128& other) const noexcept { return lane == other.lane; }
};
static_assert(sizeof(Record128) == 16);

// One coefficient per 128-bit record, replicated over its four GF32 lanes.
// The 16-byte coefficient load is exact: no padding or overread is required.
inline __m512i broadcastFour32(const std::uint32_t* coefficients) noexcept {
    const auto packed = _mm512_broadcast_i32x4(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(coefficients)));
    return _mm512_permutexvar_epi32(_mm512_setr_epi32(
        0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3), packed);
}

inline __m512i broadcastFour8(const std::uint8_t* coefficients) noexcept {
    std::uint32_t packed;
    std::memcpy(&packed, coefficients, sizeof(packed));
    return _mm512_shuffle_epi8(_mm512_set1_epi32(std::int32_t(packed)),
        _mm512_setr_epi64(0, 0, 0x0101010101010101LL, 0x0101010101010101LL,
            0x0202020202020202LL, 0x0202020202020202LL,
            0x0303030303030303LL, 0x0303030303030303LL));
}

inline __m512i xor3(__m512i a, __m512i b, __m512i c) noexcept {
    return _mm512_ternarylogic_epi32(a, b, c, 0x96);
}

// For P=L+X^32 H, first fold H*(1+X^2+X^3+X^7). Its overflow has degree <=6;
// the second fold therefore has degree <=13 and cannot overflow again.
inline __m512i reduce512(__m512i low, __m512i high) noexcept {
    const auto overflow = xor3(_mm512_srli_epi32(high, 30),
        _mm512_srli_epi32(high, 29), _mm512_srli_epi32(high, 25));
    const auto folded = xor3(high, _mm512_slli_epi32(high, 2), _mm512_slli_epi32(high, 3));
    const auto second = xor3(overflow, _mm512_slli_epi32(overflow, 2), _mm512_slli_epi32(overflow, 3));
    return xor3(xor3(low, folded, _mm512_slli_epi32(high, 7)),
        second, _mm512_slli_epi32(overflow, 7));
}

// Four VPCLMUL instructions calculate sixteen independent 32x32 products.
// Zeroing the neighboring dword is essential: CLMUL of packed 64-bit pairs
// without that separation would entangle different GF32 lanes.
inline __m512i multiply512(__m512i a, __m512i b) noexcept {
    const auto mask = _mm512_set1_epi64(0xffffffffULL);
    const auto ae = _mm512_and_si512(a, mask), be = _mm512_and_si512(b, mask);
    const auto ao = _mm512_srli_epi64(a, 32), bo = _mm512_srli_epi64(b, 32);
    const auto p0 = _mm512_clmulepi64_epi128(ae, be, 0x00);
    const auto p1 = _mm512_clmulepi64_epi128(ao, bo, 0x00);
    const auto p2 = _mm512_clmulepi64_epi128(ae, be, 0x11);
    const auto p3 = _mm512_clmulepi64_epi128(ao, bo, 0x11);
    const auto p01 = _mm512_unpacklo_epi32(p0, p1), p23 = _mm512_unpacklo_epi32(p2, p3);
    return reduce512(_mm512_unpacklo_epi64(p01, p23), _mm512_unpackhi_epi64(p01, p23));
}

template<unsigned A, unsigned B, unsigned C, unsigned D>
inline __m512i selectBytes(__m512i x) noexcept {
    constexpr auto indices = [](unsigned base) {
        const auto index = [base](unsigned i) { return i == 128 ? 128U : base + i; };
        return std::int32_t(index(A) | (index(B) << 8) | (index(C) << 16) | (index(D) << 24));
    };
    return _mm512_shuffle_epi8(x, _mm512_setr_epi32(indices(0), indices(4), indices(8), indices(12),
        indices(0), indices(4), indices(8), indices(12), indices(0), indices(4), indices(8), indices(12),
        indices(0), indices(4), indices(8), indices(12)));
}

inline __m512i towerCrossInput(__m512i x) noexcept {
    const auto halves = _mm512_xor_si512(x, selectBytes<2, 3, 0, 1>(x));
    const auto total = _mm512_xor_si512(halves, selectBytes<1, 0, 3, 2>(halves));
    return _mm512_mask_blend_epi16(0xaaaaaaaaU, halves, total);
}

// Packed Karatsuba: p=[t0,t1,t2,t3], q=[t4,t4,t5,t5], r=[t6,t7,t8,t8].
// These nine variable byte products take three GFNI instructions. Reduction
// takes two additional GFNI multiplications by fixed delta and delta^2.
inline __m512i towerMultiply512(__m512i a, __m512i b) noexcept {
    const auto p = _mm512_gf2p8mul_epi8(a, b);
    const auto q = _mm512_gf2p8mul_epi8(_mm512_xor_si512(a, selectBytes<1, 0, 3, 2>(a)),
        _mm512_xor_si512(b, selectBytes<1, 0, 3, 2>(b)));
    const auto r = _mm512_gf2p8mul_epi8(towerCrossInput(a), towerCrossInput(b));
    const auto base = _mm512_xor_si512(selectBytes<0, 0, 0, 0>(p),
        xor3(selectBytes<128, 0, 128, 0>(q), selectBytes<128, 128, 0, 0>(r),
             selectBytes<128, 128, 128, 2>(r)));
    const auto linear = xor3(selectBytes<1, 128, 1, 128>(p),
        selectBytes<128, 2, 128, 128>(q), selectBytes<128, 128, 1, 128>(r));
    const auto quadratic = _mm512_xor_si512(selectBytes<2, 3, 128, 128>(p),
        selectBytes<2, 128, 128, 128>(q));
    return xor3(base, _mm512_gf2p8mul_epi8(linear, _mm512_set1_epi8(delta)),
        _mm512_gf2p8mul_epi8(quadratic, _mm512_set1_epi8(delta2)));
}

struct SmallMatrices { std::uint64_t low, high; };
constexpr SmallMatrices smallMatrices(std::uint8_t c) noexcept {
    SmallMatrices result{};
    for (unsigned output = 0; output != 8; ++output)
        for (unsigned input = 0; input != 8; ++input) {
            const auto product = unsigned(c) << input;
            result.low |= std::uint64_t((product >> output) & 1U) << (8 * (7 - output) + input);
            result.high |= std::uint64_t((product >> (output + 8)) & 1U) << (8 * (7 - output) + input);
        }
    return result;
}
inline constexpr auto smallMatrixTable = [] {
    std::array<SmallMatrices, 256> result{};
    for (unsigned i = 0; i != 256; ++i) result[i] = smallMatrices(std::uint8_t(i));
    return result;
}();

// Restricted polynomial-byte coefficients, NOT the embedded GF256 subfield
// and NOT the uniform nonzero GF32 ensemble. Included only as a cost control.
inline __m512i smallMultiply512(__m512i x, const std::uint8_t* coefficients) noexcept {
    const auto a = smallMatrixTable[coefficients[0]], b = smallMatrixTable[coefficients[1]];
    const auto c = smallMatrixTable[coefficients[2]], d = smallMatrixTable[coefficients[3]];
    const auto ml = _mm512_setr_epi64(a.low, a.low, b.low, b.low, c.low, c.low, d.low, d.low);
    const auto mh = _mm512_setr_epi64(a.high, a.high, b.high, b.high, c.high, c.high, d.high, d.high);
    const auto low = _mm512_gf2p8affine_epi64_epi8(x, ml, 0);
    const auto high = _mm512_gf2p8affine_epi64_epi8(x, mh, 0);
    const auto overflow = _mm512_srli_epi32(high, 24);
    return xor3(_mm512_xor_si512(low, _mm512_slli_epi32(high, 8)), overflow,
        xor3(_mm512_slli_epi32(overflow, 2), _mm512_slli_epi32(overflow, 3),
             _mm512_slli_epi32(overflow, 7)));
}

enum class Kernel { Gf8, Poly32, Tower32, Poly32Small };

#if defined(_MSC_VER)
#define SPIN_MEDIUM_NOINLINE __declspec(noinline)
#else
#define SPIN_MEDIUM_NOINLINE __attribute__((noinline))
#endif

// Exact in-place input/output alias is supported. Coefficients cannot alias
// output. All setup preparation and allocation are outside these hot loops.
template<Kernel K>
SPIN_MEDIUM_NOINLINE inline void diagonal(const Record128* input, const std::uint32_t* coefficients32,
    const std::uint8_t* coefficients8, Record128* output, std::size_t count) noexcept {
    std::size_t i = 0;
    for (; i + 4 <= count; i += 4) {
        const auto x = _mm512_loadu_si512(input + i);
        __m512i y;
        if constexpr (K == Kernel::Gf8) y = _mm512_gf2p8mul_epi8(x, broadcastFour8(coefficients8 + i));
        else if constexpr (K == Kernel::Poly32Small) y = smallMultiply512(x, coefficients8 + i);
        else if constexpr (K == Kernel::Poly32) y = multiply512(x, broadcastFour32(coefficients32 + i));
        else y = towerMultiply512(x, broadcastFour32(coefficients32 + i));
        _mm512_storeu_si512(output + i, y);
    }
    for (; i != count; ++i)
        for (unsigned j = 0; j != 4; ++j) {
            const auto x = input[i].lane[j];
            if constexpr (K == Kernel::Poly32) output[i].lane[j] = multiplyScalar(x, coefficients32[i]);
            else if constexpr (K == Kernel::Tower32) output[i].lane[j] = towerMultiplyScalar(x, coefficients32[i]);
            else if constexpr (K == Kernel::Poly32Small) output[i].lane[j] = multiplyScalar(x, coefficients8[i]);
            else {
                const auto c = coefficients8[i];
                output[i].lane[j] = std::uint32_t(byteMultiply(std::uint8_t(x), c)) |
                    (std::uint32_t(byteMultiply(std::uint8_t(x >> 8), c)) << 8) |
                    (std::uint32_t(byteMultiply(std::uint8_t(x >> 16), c)) << 16) |
                    (std::uint32_t(byteMultiply(std::uint8_t(x >> 24), c)) << 24);
            }
        }
}

#undef SPIN_MEDIUM_NOINLINE
} // namespace spin::research::char2::medium32
