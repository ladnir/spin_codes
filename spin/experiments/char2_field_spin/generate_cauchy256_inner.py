"""Emit Cauchy256Inner.h to stdout; --check only checks scalar circuits.

The generated hot paths have no runtime loops, coefficient lookup, or dense
GFNI fallback. Basis products are packed four records per GFNI512, then
immediately accumulated to keep their live ranges short. No file is written.
"""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
from random import Random
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3] /
                      'research/workstreams/char2_field_spin/theory'))
import state_capacity_reference as ref


def circuit(phase):
    gamma = [ref.zeta(ref.coefficients(ref.LABELS[phase][first] ^ ref.BETA, True))
             for first in (0, 8)]
    bases, masks = [], []
    for source in range(8):
        coefficients = [gamma[block][subset] for block in range(2)
                        for subset in range(1, 8) if not source & subset]
        frequency = Counter(coefficients)
        basis, span = [1], {0: 0, 1: 1}
        # Frequent actual coefficients shorten the boolean reconstruction.
        for value in sorted(frequency, key=lambda x: (-frequency[x], x)):
            if value not in span:
                bit = 1 << len(basis)
                span.update({x ^ value: mask ^ bit for x, mask in list(span.items())})
                basis.append(value)
        bases.append(basis)
        masks.append(span)
    # matrix[output][source] is a bitmask of source coefficient-basis products.
    matrix = []
    for block in range(2):
        for destination in range(8):
            matrix.append([masks[source][gamma[block][destination ^ source]]
                           if source & destination == source else 0
                           for source in range(8)])
    products = [(source, bit, value) for source, basis in enumerate(bases)
                for bit, value in enumerate(basis) if bit]
    return bases, matrix, products


def zeta(values, transpose=False):
    result = list(values)
    for bit in ((4, 2, 1) if transpose else (1, 2, 4)):
        for mask in range(8):
            if not mask & bit:
                destination, source = (mask | bit, mask) if transpose else (mask, mask | bit)
                result[destination] ^= result[source]
    return result


def evaluate(phase, values, transpose=False):
    bases, matrix, _ = circuit(phase)
    if not transpose:
        transformed = zeta(values)
        result = [ref.xor_sum(ref.mul(coefficient, transformed[source])
                             for source, basis in enumerate(bases)
                             for bit, coefficient in enumerate(basis)
                             if row[source] >> bit & 1) for row in matrix]
        return tuple(zeta(result[:8])+zeta(result[8:]))
    transformed = zeta(values[:8], True)+zeta(values[8:], True)
    result = [ref.xor_sum(ref.mul(coefficient,
                    ref.xor_sum(transformed[j] for j, row in enumerate(matrix)
                                if row[source] >> bit & 1))
                         for bit, coefficient in enumerate(basis))
              for source, basis in enumerate(bases)]
    return tuple(zeta(result, True))


def check():
    checked = 0
    rng = Random(0xC256)
    for phase in range(4):
        a, c = ref.matrices(phase, symmetric=True)
        for transpose, dimension, matrix in ((False, 8, a), (True, 16, c)):
            vectors = [tuple(int(i == j) for i in range(dimension)) for j in range(dimension)]
            vectors += [tuple(rng.randrange(256) for _ in range(dimension)) for _ in range(64)]
            for vector in vectors:
                assert evaluate(phase, vector, transpose) == ref.matvec(matrix, vector)
                checked += 1
        for _ in range(64):
            x, state = [tuple(rng.randrange(256) for _ in range(n)) for n in (16, 8)]
            d = tuple(rng.randrange(1, 256) for _ in range(8))
            r = tuple(ref.mul(di, si) for di, si in zip(d, state))
            y = tuple(xi ^ ref.mul(2, ai) for xi, ai in zip(x, evaluate(phase, r)))
            following = tuple(ri ^ ci for ri, ci in zip(r, evaluate(phase, x, True)))
            assert (y, following) == ref.step(x, state, d, phase, symmetric=True)
            u = tuple(xi ^ ai for xi, ai in zip(x, evaluate(phase, state)))
            previous = tuple(ref.mul(di, si ^ ref.mul(2, ci))
                             for di, si, ci in zip(d, state, evaluate(phase, x, True)))
            assert (u, previous) == ref.adjoint(x, state, d, phase, symmetric=True)
            checked += 2
    return checked


def xor_expression(terms):
    if not terms:
        return '_mm_setzero_si128()'
    if len(terms) == 1:
        return terms[0]
    if len(terms) == 2:
        return f'_mm_xor_si128({terms[0]}, {terms[1]})'
    groups = [terms[i::3] for i in range(3)]
    return '_mm_ternarylogic_epi32('+', '.join(xor_expression(g) for g in groups)+', 0x96)'


def zeta_code(names, transpose=False):
    lines = []
    for bit in ((4, 2, 1) if transpose else (1, 2, 4)):
        for mask in range(8):
            if not mask & bit:
                destination, source = (mask | bit, mask) if transpose else (mask, mask | bit)
                lines.append(f'{names[destination]} = _mm_xor_si128({names[destination]}, {names[source]});')
    return lines


def map_body(phase, transpose):
    _, matrix, products = circuit(phase)
    size, outputs = (16, 8) if transpose else (8, 16)
    lines = [f'// {len(products)} shared scalar products; {(len(products)+3)//4} GFNI512 instructions.']
    for i in range(size):
        expression = f'_mm512_extracti32x4_epi32(input[{i//4}], {i%4})' if i%4 else f'_mm512_castsi512_si128(input[{i//4}])'
        lines.append(f'auto a{i} = {expression};')
    for first in range(0, size, 8):
        lines += zeta_code([f'a{i}' for i in range(first, first+8)], transpose)
    for destination in range(outputs):
        terms = ([f'a{j}' for j, row in enumerate(matrix) if row[destination] & 1]
                 if transpose else [f'a{j}' for j, mask in enumerate(matrix[destination]) if mask & 1])
        lines.append(f'auto r{destination} = {xor_expression(terms)};')
    for first in range(0, len(products), 4):
        group = products[first:first+4]
        values = []
        for source, bit, _ in group:
            terms = ([f'a{j}' for j, row in enumerate(matrix) if row[source] >> bit & 1]
                     if transpose else [f'a{source}'])
            values.append(xor_expression(terms))
        values += ['_mm_setzero_si128()']*(4-len(group))
        constants = [str(value) for _, _, value in group]+['0']*(4-len(group))
        lines.append('{')
        lines.append('    const auto product = c256MultiplyFour<'+', '.join(constants)+'>('+', '.join(values)+');')
        for lane in range(len(group)):
            expression = f'_mm512_extracti32x4_epi32(product, {lane})' if lane else '_mm512_castsi512_si128(product)'
            lines.append(f'    const auto p{lane} = {expression};')
        for destination in range(outputs):
            terms = [f'p{lane}' for lane, (source, bit, _) in enumerate(group)
                     if (destination == source if transpose else matrix[destination][source] >> bit & 1)]
            if terms:
                lines.append(f'    r{destination} = {xor_expression([f"r{destination}"]+terms)};')
        lines.append('}')
    for first in range(0, outputs, 8):
        lines += zeta_code([f'r{i}' for i in range(first, first+8)], transpose)
    for first in range(0, outputs, 4):
        lines.append(f'output[{first//4}] = c256PackFour('+', '.join(f'r{i}' for i in range(first, first+4))+');')
    return lines


PREFIX = '''// Generated by generate_cauchy256_inner.py; do not hand-edit the map circuits.
#pragma once
#include "Subfield.h"

namespace spin::research::char2 {

// Certified symmetric state-eight geometry in the AES GF256 basis. A record
// is sixteen independent bytes for these subfield scalar operations.
// Requires AVX512F/BW/VL and GFNI. Input/output may alias; state and diagonal
// storage must be disjoint from input/output. State is persistent across regions.
// Each diagonal coefficient must be nonzero; the hot hooks do not validate it.
// The maps use two eight-point XOR circulants. Their nilpotent-basis products
// share GF2 coefficient bases, and injection is the exact transposed circuit.
// Map product batches are consumed immediately, not stored in a large array.
inline constexpr unsigned c256InnerWidth = 16, c256InnerState = 8;
inline constexpr std::uint8_t c256Eta = 2;
'''


SUFFIX = '''
// State h is the incoming ADJOINT of the following packet state.
// output = v + A*h; previous h = D_s*(h + eta*A^T*v).
// The independent INPUT diagonal belongs outside this hook, after output.
template<unsigned Phase>
inline void c256InnerTransposePacked(const __m512i* input4, __m512i* output4,
                                    __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {input4[0], input4[1], input4[2], input4[3]};
    const __m512i h[2] = {state2[0], state2[1]};
    __m512i moments[2], emitted[4];
    c256InjectPacked<Phase>(input, moments);
    c256EmitPacked<Phase>(h, emitted);
    output4[0] = _mm512_xor_si512(input[0], emitted[0]);
    output4[1] = _mm512_xor_si512(input[1], emitted[1]);
    output4[2] = _mm512_xor_si512(input[2], emitted[2]);
    output4[3] = _mm512_xor_si512(input[3], emitted[3]);
    state2[0] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h[0], sfConstant4<2>(moments[0])), diagonal2[0]);
    state2[1] = _mm512_gf2p8mul_epi8(_mm512_xor_si512(h[1], sfConstant4<2>(moments[1])), diagonal2[1]);
}

// Input has already received its independent input diagonal.
// r=D_s*state; output=input+eta*A*r; following state=r+A^T*input.
template<unsigned Phase>
inline void c256InnerForwardPacked(const __m512i* input4, __m512i* output4,
                                  __m512i* state2, const __m512i* diagonal2) noexcept {
    static_assert(Phase < 4);
    const __m512i input[4] = {input4[0], input4[1], input4[2], input4[3]};
    const __m512i r[2] = {_mm512_gf2p8mul_epi8(state2[0], diagonal2[0]),
                          _mm512_gf2p8mul_epi8(state2[1], diagonal2[1])};
    const __m512i scaled[2] = {sfConstant4<2>(r[0]), sfConstant4<2>(r[1])};
    __m512i moments[2], emitted[4];
    c256InjectPacked<Phase>(input, moments);
    c256EmitPacked<Phase>(scaled, emitted);
    output4[0] = _mm512_xor_si512(input[0], emitted[0]);
    output4[1] = _mm512_xor_si512(input[1], emitted[1]);
    output4[2] = _mm512_xor_si512(input[2], emitted[2]);
    output4[3] = _mm512_xor_si512(input[3], emitted[3]);
    state2[0] = _mm512_xor_si512(r[0], moments[0]);
    state2[1] = _mm512_xor_si512(r[1], moments[1]);
}

template<unsigned Phase, bool Transpose>
inline void c256InnerBlocks(const Block* input16, Block* output16, Block* state8,
                            const std::uint8_t* diagonal8) noexcept {
    const __m512i input[4] = {_mm512_loadu_si512(input16), _mm512_loadu_si512(input16+4),
                              _mm512_loadu_si512(input16+8), _mm512_loadu_si512(input16+12)};
    __m512i output[4], state[2] = {_mm512_loadu_si512(state8), _mm512_loadu_si512(state8+4)};
    const __m512i diagonal[2] = {sfBroadcastFour(diagonal8), sfBroadcastFour(diagonal8+4)};
    if constexpr (Transpose) c256InnerTransposePacked<Phase>(input, output, state, diagonal);
    else c256InnerForwardPacked<Phase>(input, output, state, diagonal);
    _mm512_storeu_si512(output16, output[0]);
    _mm512_storeu_si512(output16+4, output[1]);
    _mm512_storeu_si512(output16+8, output[2]);
    _mm512_storeu_si512(output16+12, output[3]);
    _mm512_storeu_si512(state8, state[0]);
    _mm512_storeu_si512(state8+4, state[1]);
}

inline void c256InnerTranspose(const Block* input16, Block* output16, Block* state8,
                               const std::uint8_t* diagonal8, unsigned phase) {
    switch (phase) {
    case 0: c256InnerBlocks<0, true>(input16, output16, state8, diagonal8); break;
    case 1: c256InnerBlocks<1, true>(input16, output16, state8, diagonal8); break;
    case 2: c256InnerBlocks<2, true>(input16, output16, state8, diagonal8); break;
    case 3: c256InnerBlocks<3, true>(input16, output16, state8, diagonal8); break;
    default: throw std::invalid_argument("Cauchy256 packet phase");
    }
}

inline void c256InnerForward(const Block* input16, Block* output16, Block* state8,
                             const std::uint8_t* diagonal8, unsigned phase) {
    switch (phase) {
    case 0: c256InnerBlocks<0, false>(input16, output16, state8, diagonal8); break;
    case 1: c256InnerBlocks<1, false>(input16, output16, state8, diagonal8); break;
    case 2: c256InnerBlocks<2, false>(input16, output16, state8, diagonal8); break;
    case 3: c256InnerBlocks<3, false>(input16, output16, state8, diagonal8); break;
    default: throw std::invalid_argument("Cauchy256 packet phase");
    }
}

constexpr std::uint8_t c256Coefficient(unsigned phase, unsigned output, unsigned state) noexcept {
    const auto x = std::uint8_t(c256Labels[phase][output] ^ c256Nodes[state]);
    const auto x2 = sfMultiply(x, x), x3 = sfMultiply(x2, x), x4 = sfMultiply(x2, x2);
    return std::uint8_t(sfMultiply(x4, x3) ^ x3 ^ x ^ 1U);
}

// Slow one-byte reference hooks only. These are not fallbacks for a hot path.
// They permit input/output alias, but require distinct state storage.
inline void c256InnerForwardScalar(const std::uint8_t* input16, std::uint8_t* output16,
                                   std::uint8_t* state8, const std::uint8_t* diagonal8,
                                   unsigned phase) {
    if (phase >= 4) throw std::invalid_argument("Cauchy256 scalar phase");
    std::uint8_t x[16], r[8], following[8];
    std::memcpy(x, input16, 16);
    for (unsigned j=0; j!=8; ++j) following[j] = r[j] = sfMultiply(state8[j], diagonal8[j]);
    for (unsigned i=0; i!=16; ++i) {
        auto y = x[i];
        for (unsigned j=0; j!=8; ++j) {
            const auto a = c256Coefficient(phase, i, j);
            y ^= sfMultiply(2, sfMultiply(a, r[j]));
            following[j] ^= sfMultiply(a, x[i]);
        }
        output16[i] = y;
    }
    std::memcpy(state8, following, 8);
}

inline void c256InnerTransposeScalar(const std::uint8_t* input16, std::uint8_t* output16,
                                     std::uint8_t* state8, const std::uint8_t* diagonal8,
                                     unsigned phase) {
    if (phase >= 4) throw std::invalid_argument("Cauchy256 scalar phase");
    std::uint8_t x[16], h[8], previous[8];
    std::memcpy(x, input16, 16);
    std::memcpy(h, state8, 8);
    std::memcpy(previous, h, 8);
    for (unsigned i=0; i!=16; ++i) {
        auto u = x[i];
        for (unsigned j=0; j!=8; ++j) {
            const auto a = c256Coefficient(phase, i, j);
            u ^= sfMultiply(a, h[j]);
            previous[j] ^= sfMultiply(2, sfMultiply(a, x[i]));
        }
        output16[i] = u;
    }
    for (unsigned j=0; j!=8; ++j) state8[j] = sfMultiply(diagonal8[j], previous[j]);
}

} // namespace spin::research::char2
'''


def generate():
    lines = [PREFIX]
    lines.append('inline constexpr std::uint8_t c256Labels[4][16] = {')
    lines += ['    {'+', '.join(map(str, row))+'},' for row in ref.LABELS]
    lines.append('};')
    lines.append('inline constexpr std::uint8_t c256Nodes[8] = {'+', '.join(map(str, ref.NODES))+'};')
    lines += ['''
inline __m512i c256PackFour(__m128i a, __m128i b, __m128i c, __m128i d) noexcept {
    auto result = _mm512_castsi128_si512(a);
    result = _mm512_inserti32x4(result, b, 1);
    result = _mm512_inserti32x4(result, c, 2);
    return _mm512_inserti32x4(result, d, 3);
}

template<std::uint8_t A, std::uint8_t B, std::uint8_t C, std::uint8_t D>
inline __m512i c256MultiplyFour(__m128i a, __m128i b, __m128i c, __m128i d) noexcept {
    constexpr auto repeat = [](std::uint8_t x) { return static_cast<long long>(0x0101010101010101ULL*x); };
    const auto coefficients = _mm512_setr_epi64(repeat(A), repeat(A), repeat(B), repeat(B),
                                               repeat(C), repeat(C), repeat(D), repeat(D));
    return _mm512_gf2p8mul_epi8(c256PackFour(a, b, c, d), coefficients);
}
''']
    for transpose, name in ((False, 'c256EmitPacked'), (True, 'c256InjectPacked')):
        lines += ['template<unsigned Phase>', f'inline void {name}(const __m512i* input, __m512i* output) noexcept {{',
                  '    static_assert(Phase < 4);']
        for phase in range(4):
            lines.append(f'    {"if" if phase == 0 else "else if"} constexpr (Phase == {phase}) {{')
            lines += ['        '+line for line in map_body(phase, transpose)]
            lines.append('    }')
        lines.append('}\n')
    lines.append(SUFFIX)
    return '\n'.join(lines)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    checks = check()
    if args.check:
        print(dict(scalar_circuit_checks=checks,
                   scalar_products=[len(circuit(p)[2]) for p in range(4)],
                   packed_products=[(len(circuit(p)[2])+3)//4 for p in range(4)]))
    else:
        print(generate())
