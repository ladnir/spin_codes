#pragma once
#include "Cauchy256Inner.h"

namespace spin::research::char2 {

// Isolated A/B candidate: the SAME certified map and AES-byte record basis.
// This trades the generated coefficient-basis circuit for more packed GFNI
// products and fewer packing/XOR operations. No speed advantage is assumed.
// Same alias, diagonal, persistent-state, and CPU-feature contracts as the
// generated inner. Input diagonal remains external to both packed hooks.

template<std::uint8_t A, std::uint8_t B, std::uint8_t C, std::uint8_t D>
inline __m512i c256DenseCoefficients() noexcept {
    constexpr auto repeat = [](std::uint8_t x) {
        return static_cast<long long>(0x0101010101010101ULL*x);
    };
    return _mm512_setr_epi64(repeat(A), repeat(A), repeat(B), repeat(B),
                            repeat(C), repeat(C), repeat(D), repeat(D));
}

template<unsigned Lane>
inline __m512i c256DenseBroadcast(__m512i records) noexcept {
    static_assert(Lane < 4);
    return _mm512_shuffle_i32x4(records, records, Lane*0x55);
}

template<unsigned Phase, unsigned Batch, unsigned State, std::uint8_t Scale>
inline __m512i c256DenseEmissionCoefficient() noexcept {
    static_assert(Phase < 4 && Batch < 4 && State < 8);
    return c256DenseCoefficients<sfMultiply(Scale, c256Coefficient(Phase, 4*Batch, State)),
        sfMultiply(Scale, c256Coefficient(Phase, 4*Batch+1, State)),
        sfMultiply(Scale, c256Coefficient(Phase, 4*Batch+2, State)),
        sfMultiply(Scale, c256Coefficient(Phase, 4*Batch+3, State))>();
}

template<unsigned Phase, unsigned Batch, unsigned Input, std::uint8_t Scale>
inline __m512i c256DenseInjectionCoefficient() noexcept {
    static_assert(Phase < 4 && Batch < 2 && Input < 16);
    return c256DenseCoefficients<sfMultiply(Scale, c256Coefficient(Phase, Input, 4*Batch)),
        sfMultiply(Scale, c256Coefficient(Phase, Input, 4*Batch+1)),
        sfMultiply(Scale, c256Coefficient(Phase, Input, 4*Batch+2)),
        sfMultiply(Scale, c256Coefficient(Phase, Input, 4*Batch+3))>();
}

template<unsigned Phase, unsigned Batch, unsigned First, std::uint8_t Scale>
inline __m512i c256DenseEmissionPair(__m512i a, __m512i b, __m512i previous) noexcept {
    const auto pa = _mm512_gf2p8mul_epi8(a, c256DenseEmissionCoefficient<Phase, Batch, First, Scale>());
    const auto pb = _mm512_gf2p8mul_epi8(b, c256DenseEmissionCoefficient<Phase, Batch, First+1, Scale>());
    if constexpr (First == 0) return _mm512_xor_si512(pa, pb);
    else return _mm512_ternarylogic_epi32(previous, pa, pb, 0x96);
}

template<unsigned Phase, unsigned First, std::uint8_t Scale>
inline void c256DenseEmissionAccumulate(__m512i records, __m512i& r0, __m512i& r1,
                                        __m512i& r2, __m512i& r3) noexcept {
    static_assert(First < 8 && First%2 == 0);
    const auto a = c256DenseBroadcast<First%4>(records);
    const auto b = c256DenseBroadcast<(First+1)%4>(records);
    r0 = c256DenseEmissionPair<Phase, 0, First, Scale>(a, b, r0);
    r1 = c256DenseEmissionPair<Phase, 1, First, Scale>(a, b, r1);
    r2 = c256DenseEmissionPair<Phase, 2, First, Scale>(a, b, r2);
    r3 = c256DenseEmissionPair<Phase, 3, First, Scale>(a, b, r3);
}

template<unsigned Phase, std::uint8_t Scale = 1>
inline void c256EmitDensePacked(const __m512i* input2, __m512i* output4) noexcept {
    const auto h0 = input2[0], h1 = input2[1];
    auto r0 = _mm512_setzero_si512(), r1 = r0, r2 = r0, r3 = r0;
    c256DenseEmissionAccumulate<Phase, 0, Scale>(h0, r0, r1, r2, r3);
    c256DenseEmissionAccumulate<Phase, 2, Scale>(h0, r0, r1, r2, r3);
    c256DenseEmissionAccumulate<Phase, 4, Scale>(h1, r0, r1, r2, r3);
    c256DenseEmissionAccumulate<Phase, 6, Scale>(h1, r0, r1, r2, r3);
    output4[0] = r0;
    output4[1] = r1;
    output4[2] = r2;
    output4[3] = r3;
}

template<unsigned Phase, unsigned Batch, unsigned First, std::uint8_t Scale>
inline __m512i c256DenseInjectionPair(__m512i a, __m512i b, __m512i previous) noexcept {
    const auto pa = _mm512_gf2p8mul_epi8(a, c256DenseInjectionCoefficient<Phase, Batch, First, Scale>());
    const auto pb = _mm512_gf2p8mul_epi8(b, c256DenseInjectionCoefficient<Phase, Batch, First+1, Scale>());
    if constexpr (First == 0) return _mm512_xor_si512(pa, pb);
    else return _mm512_ternarylogic_epi32(previous, pa, pb, 0x96);
}

template<unsigned Phase, unsigned First, std::uint8_t Scale>
inline void c256DenseInjectionAccumulate(__m512i records, __m512i& r0, __m512i& r1) noexcept {
    static_assert(First < 16 && First%2 == 0);
    const auto a = c256DenseBroadcast<First%4>(records);
    const auto b = c256DenseBroadcast<(First+1)%4>(records);
    r0 = c256DenseInjectionPair<Phase, 0, First, Scale>(a, b, r0);
    r1 = c256DenseInjectionPair<Phase, 1, First, Scale>(a, b, r1);
}

template<unsigned Phase, std::uint8_t Scale = 1>
inline void c256InjectDensePacked(const __m512i* input4, __m512i* output2) noexcept {
    const auto v0 = input4[0], v1 = input4[1], v2 = input4[2], v3 = input4[3];
    auto r0 = _mm512_setzero_si512(), r1 = r0;
    c256DenseInjectionAccumulate<Phase, 0, Scale>(v0, r0, r1);
    c256DenseInjectionAccumulate<Phase, 2, Scale>(v0, r0, r1);
    c256DenseInjectionAccumulate<Phase, 4, Scale>(v1, r0, r1);
    c256DenseInjectionAccumulate<Phase, 6, Scale>(v1, r0, r1);
    c256DenseInjectionAccumulate<Phase, 8, Scale>(v2, r0, r1);
    c256DenseInjectionAccumulate<Phase, 10, Scale>(v2, r0, r1);
    c256DenseInjectionAccumulate<Phase, 12, Scale>(v3, r0, r1);
    c256DenseInjectionAccumulate<Phase, 14, Scale>(v3, r0, r1);
    output2[0] = r0;
    output2[1] = r1;
}

// 64 map GFNI512 + 2 diagonal products, independent of phase. Eta is folded
// into fixed map coefficients at compile time, not a separate multiplication.
// Map data movement: 24 reused record broadcasts; 32 XOR/XOR3 operations.
// Recurrence: u=v+A*h, previous h=D_s*(h+2*A^T*v).
template<unsigned Phase>
inline void c256InnerTransposeDensePacked(const __m512i* input4, __m512i* output4,
                                         __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {input4[0], input4[1], input4[2], input4[3]};
    const __m512i h[2] = {state2[0], state2[1]};
    __m512i moments[2], emitted[4];
    c256InjectDensePacked<Phase, 2>(input, moments);
    c256EmitDensePacked<Phase>(h, emitted);
    output4[0] = _mm512_xor_si512(input[0], emitted[0]);
    output4[1] = _mm512_xor_si512(input[1], emitted[1]);
    output4[2] = _mm512_xor_si512(input[2], emitted[2]);
    output4[3] = _mm512_xor_si512(input[3], emitted[3]);
    state2[0] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h[0], moments[0]), diagonal2[0]);
    state2[1] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h[1], moments[1]), diagonal2[1]);
}

template<unsigned Phase>
inline void c256InnerForwardDensePacked(const __m512i* input4, __m512i* output4,
                                       __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {input4[0], input4[1], input4[2], input4[3]};
    const __m512i r[2] = {_mm512_gf2p8mul_epi8(state2[0], diagonal2[0]),
                          _mm512_gf2p8mul_epi8(state2[1], diagonal2[1])};
    __m512i moments[2], emitted[4];
    c256InjectDensePacked<Phase>(input, moments);
    c256EmitDensePacked<Phase, 2>(r, emitted);
    output4[0] = _mm512_xor_si512(input[0], emitted[0]);
    output4[1] = _mm512_xor_si512(input[1], emitted[1]);
    output4[2] = _mm512_xor_si512(input[2], emitted[2]);
    output4[3] = _mm512_xor_si512(input[3], emitted[3]);
    state2[0] = _mm512_xor_si512(r[0], moments[0]);
    state2[1] = _mm512_xor_si512(r[1], moments[1]);
}

} // namespace spin::research::char2
