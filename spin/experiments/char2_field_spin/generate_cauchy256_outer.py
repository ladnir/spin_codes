"""Generate/check the fixed AES-field GRS[256,176] outer.

The generated native circuit is not a dense matrix loop. It implements
three 64-point XOR-circulants and eleven 16-point XOR-circulants, sharing
binary coefficient bases before bounded sixteen-output stages. Forward
uses the transposed coefficient-basis circuit, not another dense map.
--patch emits an apply_patch fragment; this script does not write files.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
import random


def multiply(a, b):
    result = 0
    for _ in range(8):
        if b & 1:
            result ^= a
        a = ((a << 1) ^ (0x11b if a & 128 else 0)) & 255
        b >>= 1
    return result


MUL = tuple(bytes(multiply(a, b) for b in range(256)) for a in range(256))
INV = (0,) + tuple(next(b for b in range(1, 256) if MUL[a][b] == 1)
                   for a in range(1, 256))


def interpolation(dimension, targets):
    weights = []
    for source in range(dimension):
        denominator = 1
        for other in range(dimension):
            if source != other:
                denominator = MUL[denominator][source ^ other]
        weights.append(INV[denominator])
    result = []
    for target in targets:
        numerator = 1
        for source in range(dimension):
            numerator = MUL[numerator][target ^ source]
        result.append([MUL[MUL[numerator][INV[target ^ source]]][weights[source]]
                       for source in range(dimension)])
    return result


def zeta(values, transpose=False):
    values = list(values)
    for shift in range((len(values)-1).bit_length()):
        bit = 1 << shift
        for mask in range(len(values)):
            if not mask & bit:
                destination, source = (mask | bit, mask) if transpose else (mask, mask | bit)
                values[destination] ^= values[source]
    return values


def basis_for(coefficients):
    basis, span = [1], {0: 0, 1: 1}
    for coefficient in sorted(set(coefficients)-{0}):
        if coefficient not in span:
            bit = 1 << len(basis)
            span.update({a ^ coefficient: mask | bit for a, mask in tuple(span.items())})
            basis.append(coefficient)
    return basis, span


@dataclass
class Map:
    dimension: int
    size: int

    def __post_init__(self):
        self.matrix = interpolation(self.dimension, range(self.dimension, self.dimension+self.size))
        self.blocks = self.dimension//self.size
        self.kernels = [self.matrix[0][j*self.size:(j+1)*self.size]
                        for j in range(self.blocks)]
        assert all(self.matrix[row][column] ==
                   self.kernels[column//self.size][row ^ (column % self.size)]
                   for row in range(self.size) for column in range(self.dimension))
        self.gamma = [zeta(kernel) for kernel in self.kernels]
        self.bases, self.spans, self.locations = [], [], []
        self.products = 0
        for source in range(self.size):
            needed = [g[subset] for g in self.gamma for subset in range(self.size)
                      if not source & subset]
            basis, span = basis_for(needed)
            self.bases.append(basis)
            self.spans.append(span)
            self.locations.append([source] + list(range(self.size+self.products,
                                                        self.size+self.products+len(basis)-1)))
            self.products += len(basis)-1
        self.scratch = self.size+self.products
        self.tiles = []
        for first in range(0, self.dimension, 16):
            sources = []
            for source in range(self.size):
                by_coefficient = {}
                for destination in range(first, first+16):
                    mask = destination % self.size
                    if source & mask == source:
                        coefficient = self.gamma[destination//self.size][mask ^ source]
                        if coefficient:
                            by_coefficient.setdefault(coefficient, []).append(destination-first)
                if by_coefficient:
                    sources.append((source, sorted(by_coefficient.items())))
            self.tiles.append((first, sources))

    def wide(self, values):
        """The same prepared-basis, tiled schedule emitted for native code."""
        prepared = zeta(values) + [0]*self.products
        for source, (basis, locations) in enumerate(zip(self.bases, self.locations)):
            for coefficient, location in zip(basis[1:], locations[1:]):
                prepared[location] = MUL[coefficient][prepared[source]]
        out = [0]*self.dimension
        for first, sources in self.tiles:
            accumulators = [0]*16
            for source, groups in sources:
                for coefficient, destinations in groups:
                    mask, value = self.spans[source][coefficient], 0
                    for bit, location in enumerate(self.locations[source]):
                        if mask >> bit & 1:
                            value ^= prepared[location]
                    for destination in destinations:
                        accumulators[destination] ^= value
            out[first:first+16] = accumulators
        return sum((zeta(out[first:first+self.size])
                    for first in range(0, self.dimension, self.size)), [])

    def narrow(self, values):
        """Field transpose of wide: reverse each coefficient-basis fan-out."""
        transformed = sum((zeta(values[first:first+self.size], True)
                           for first in range(0, self.dimension, self.size)), [])
        out = []
        for source, basis in enumerate(self.bases):
            value = 0
            for bit, coefficient in enumerate(basis):
                term = 0
                for destination in range(self.dimension):
                    mask = destination % self.size
                    if source & mask == source:
                        c = self.gamma[destination//self.size][mask ^ source]
                        if self.spans[source][c] >> bit & 1:
                            term ^= transformed[destination]
                value ^= MUL[coefficient][term]
            out.append(value)
        return zeta(out, True)


def load(pointer, index):
    return f'_mm512_loadu_si512({pointer}+{index})'


def store(pointer, index, value):
    return f'_mm512_storeu_si512({pointer}+{index}, {value});'


def emit_sum(lines, name, terms, counter, indent='    '):
    assert terms
    lines.append(f'{indent}auto {name} = {terms[0]};')
    for start in range(1, len(terms), 2):
        if start+1 < len(terms):
            lines.append(f'{indent}{name} = _mm512_ternarylogic_epi64({name}, {terms[start]}, {terms[start+1]}, 0x96);')
        else:
            lines.append(f'{indent}{name} = _mm512_xor_si512({name}, {terms[start]});')
        counter['boolean'] += 1


def emit_zeta_helpers(lines):
    for transpose in (False, True):
        suffix = 'T' if transpose else 'N'
        lines += [f'static C256_NOINLINE void zeta16{suffix}(const __m512i* input, __m512i* output) {{']
        lines += [f'    auto x{i} = {load("input", i)};' for i in range(16)]
        for bit in (1, 2, 4, 8):
            for mask in range(16):
                if not mask & bit:
                    destination, source = (mask | bit, mask) if transpose else (mask, mask | bit)
                    lines.append(f'    x{destination} = _mm512_xor_si512(x{destination}, x{source});')
        lines += ['    '+store('output', i, f'x{i}') for i in range(16)]
        lines += ['}', f'static C256_NOINLINE void zetaHigh{suffix}(__m512i* lower, __m512i* upper) {{']
        destination = 'upper' if transpose else 'lower'
        for i in range(16):
            lines.append('    '+store(destination, i, f'_mm512_xor_si512({load("lower", i)}, {load("upper", i)})'))
        lines += ['}', f'static C256_NOINLINE void zeta64{suffix}(const __m512i* input, __m512i* output) {{']
        lines += [f'    zeta16{suffix}(input+{i}, output+{i});' for i in range(0, 64, 16)]
        lines += [f'    zetaHigh{suffix}(output+{a}, output+{b});'
                  for a, b in ((0, 16), (32, 48), (0, 32), (16, 48))]
        lines.append('}')


def emit_map(lines, name, mapping):
    wide, narrow = Counter(), Counter()
    size = mapping.size
    lines.append(f'static C256_NOINLINE void {name}Prepare(__m512i* products) {{')
    for source, basis in enumerate(mapping.bases):
        if len(basis) == 1:
            continue
        lines += ['    {', f'        const auto x = {load("products", source)};']
        for coefficient, location in zip(basis[1:], mapping.locations[source][1:]):
            lines.append('        '+store('products', location,
                f'_mm512_gf2p8mul_epi8(x, _mm512_set1_epi8(static_cast<char>({coefficient})))'))
        lines.append('    }')
    lines.append('}')
    for tile_number, (first, sources) in enumerate(mapping.tiles):
        lines += [f'static C256_NOINLINE void {name}Tile{tile_number}(const __m512i* products, __m512i* output) {{',
                  '    __m512i '+', '.join(f'y{i}' for i in range(16))+';']
        initialized = set()
        for source, groups in sources:
            masks = [mapping.spans[source][coefficient] for coefficient, _ in groups]
            used = {i for mask in masks for i in range(len(mapping.bases[source])) if mask >> i & 1}
            lines.append('    {')
            for bit in sorted(used):
                lines.append(f'        const auto b{bit} = {load("products", mapping.locations[source][bit])};')
                wide['prepared_loads'] += 1
            for coefficient, destinations in groups:
                mask = mapping.spans[source][coefficient]
                terms = [f'b{i}' for i in sorted(used) if mask >> i & 1]
                lines.append('        {')
                emit_sum(lines, 'term', terms, wide, '            ')
                for destination in destinations:
                    if destination in initialized:
                        lines.append(f'            y{destination} = _mm512_xor_si512(y{destination}, term);')
                        wide['boolean'] += 1
                    else:
                        initialized.add(destination)
                        lines.append(f'            y{destination} = term;')
                lines.append('        }')
            lines.append('    }')
        assert len(initialized) == 16
        lines += ['    '+store('output', i, f'y{i}') for i in range(16)]
        lines.append('}')
    lines += [f'static C256_NOINLINE void {name}Wide(const __m512i* input, __m512i* output) {{',
              f'    alignas(64) __m512i products[{mapping.scratch}];',
              f'    zeta{size}N(input, products);', f'    {name}Prepare(products);']
    lines += [f'    {name}Tile{j}(products, output+{first});'
              for j, (first, _) in enumerate(mapping.tiles)]
    lines += [f'    zeta{size}N(output+{first}, output+{first});'
              for first in range(0, mapping.dimension, size)]
    lines.append('}')
    for source, basis in enumerate(mapping.bases):
        lines.append(f'static C256_NOINLINE __m512i {name}Narrow{source}(const __m512i* input) {{')
        for bit, coefficient in enumerate(basis):
            terms = []
            for destination in range(mapping.dimension):
                mask = destination % size
                if source & mask == source:
                    c = mapping.gamma[destination//size][mask ^ source]
                    if mapping.spans[source][c] >> bit & 1:
                        terms.append(load('input', destination))
            narrow['prepared_loads'] += len(terms)
            assert terms
            lines.append('    {') if bit else None
            emit_sum(lines, 'term' if bit else 'result', terms, narrow, '        ' if bit else '    ')
            if bit:
                lines.append(f'        result = _mm512_xor_si512(result, _mm512_gf2p8mul_epi8(term, _mm512_set1_epi8(static_cast<char>({coefficient}))));')
                narrow['boolean'] += 1
                lines.append('    }')
        lines += ['    return result;', '}']
    lines += [f'static C256_NOINLINE void {name}Narrow(const __m512i* input, __m512i* output) {{',
              f'    alignas(64) __m512i transformed[{mapping.dimension}];']
    lines += [f'    zeta{size}T(input+{first}, transformed+{first});'
              for first in range(0, mapping.dimension, size)]
    lines += ['    '+store('output', source, f'{name}Narrow{source}(transformed)')
              for source in range(size)]
    lines += [f'    zeta{size}T(output, output);', '}']
    transforms = (mapping.blocks+1)*(size//2)*(size.bit_length()-1)
    for count in (wide, narrow):
        count['boolean'] += transforms
        count['products'] = mapping.products
    return dict(wide=dict(wide), narrow=dict(narrow))


def generate():
    base, shorten = Map(192, 64), Map(176, 16)
    lines = ['// Generated by generate_cauchy256_outer.py; do not edit manually.',
             '// Fixed systematic AES GF256 GRS[256,176], source points0..175.',
             '// Every vector holds four128-bit records; GF256 acts on each byte.',
             '// Caller input/output ranges must not overlap; unaligned pointers are accepted.',
             '// No caller scratch. Explicitly allocated nested scratch:42KiB transpose,12KiB forward, plus compiler frames.',
             '// All scratch is in ordinary noinline functions, never coroutines.',
             '#pragma once', '#include <immintrin.h>',
             '#if defined(_MSC_VER)', '#define C256_NOINLINE __declspec(noinline)',
             '#else', '#define C256_NOINLINE __attribute__((noinline))', '#endif',
             'namespace spin::research::char2 {', 'namespace c256_outer_detail {']
    emit_zeta_helpers(lines)
    cost = {'base': emit_map(lines, 'base', base), 'shorten': emit_map(lines, 'shorten', shorten)}
    lines += ['} // namespace c256_outer_detail',
              'static C256_NOINLINE void c256OuterTranspose(const __m512i* input256, __m512i* output176) {',
              '    alignas(64) __m512i base[192], missing[16];',
              '    c256_outer_detail::baseWide(input256+192, base);']
    lines += ['    '+store('missing', i, f'_mm512_xor_si512({load("input256", 176+i)}, {load("base", 176+i)})')
              for i in range(16)]
    lines.append('    c256_outer_detail::shortenWide(missing, output176);')
    lines += ['    '+store('output176', i, f'_mm512_ternarylogic_epi64({load("output176", i)}, {load("input256", i)}, {load("base", i)}, 0x96)')
              for i in range(176)]
    lines += ['}', 'static C256_NOINLINE void c256OuterForward(const __m512i* input176, __m512i* output256) {']
    lines += ['    '+store('output256', i, load('input176', i)) for i in range(176)]
    lines += ['    c256_outer_detail::shortenNarrow(input176, output256+176);',
              '    c256_outer_detail::baseNarrow(output256, output256+192);', '}',
              '} // namespace spin::research::char2', '#undef C256_NOINLINE']
    # Marked fragments permit apply_patch generation without shell file writes.
    parts = [lines[first:first+400] + [f'// c256-generated-part {first//400}']
             for first in range(0, len(lines), 400)]
    header = '\n'.join(sum(parts, []))+'\n'
    return header, parts, cost, base, shorten


def check(header, base, shorten):
    from functools import reduce
    xor = lambda values: reduce(int.__xor__, values, 0)
    matrix = interpolation(176, range(176, 256))
    def forward(x):
        intermediate = x+shorten.narrow(x)
        return intermediate+base.narrow(intermediate)
    def transpose(x):
        intermediate = base.wide(x[192:])
        missing = [a ^ b for a, b in zip(x[176:192], intermediate[176:])]
        correction = shorten.wide(missing)
        return [a ^ b ^ c for a, b, c in zip(x, intermediate, correction)]
    checked = 0
    for source in range(176):
        x = [int(i == source) for i in range(176)]
        assert forward(x) == x+[row[source] for row in matrix]
        checked += 1
    for source in range(256):
        x = [int(i == source) for i in range(256)]
        expected = [int(i == source) for i in range(176)] if source < 176 else matrix[source-176]
        assert transpose(x) == expected
        checked += 1
    rng = random.Random(0xC2560176)
    for _ in range(16):
        x, y = [rng.randrange(256) for _ in range(176)], [rng.randrange(256) for _ in range(256)]
        encoded = forward(x)
        expected = x+[xor(MUL[a][b] for a, b in zip(row, x)) for row in matrix]
        assert encoded == expected
        expected_t = [y[j] ^ xor(MUL[matrix[i][j]][y[176+i]] for i in range(80)) for j in range(176)]
        actual_t = transpose(y)
        assert actual_t == expected_t
        assert xor(MUL[a][b] for a, b in zip(encoded, y)) == xor(MUL[a][b] for a, b in zip(x, actual_t))
        checked += 3
    saved = Path(__file__).with_name('Cauchy256Outer.h')
    assert saved.read_text(encoding='utf-8') == header
    assert base.products == 400 and shorten.products == 103
    return dict(exact_scalar_checks=checked, products_each_direction=503,
                generated_header_matches=True, native_execution=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--emit', action='store_true')
    parser.add_argument('--patch', type=int)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    header, parts, cost, base, shorten = generate()
    if args.emit:
        print(header, end='')
    elif args.patch is not None:
        index = args.patch
        target = Path(__file__).with_name('Cauchy256Outer.h').as_posix()
        print('*** Begin Patch')
        print(f'*** {"Add" if index == 0 else "Update"} File: {target}')
        if index:
            print('@@')
            print(f' // c256-generated-part {index-1}')
        print('\n'.join('+'+line for line in parts[index]))
        print('*** End Patch')
    else:
        result = dict(cost=cost, parts=len(parts), header_bytes=len(header))
        if args.check:
            result.update(check(header, base, shorten))
        print(result)
