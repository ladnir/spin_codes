"""Generate formally checked BCH transpose circuits using fixed polynomial products.

Raw evaluates the q/p polynomial basis (same code, different message basis).
Exact reconstructs the production systematic basis using a truncated reciprocal
of q and five rank corrections. Neither implementation changes production.
The polynomial products use fixed-coefficient Karatsuba, then the complete
forward XOR graph is transposed. All final 256-bit forms are checked exactly.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re


P = 0x1C3B42FE115AA90747020D79AC738ADB
Q = 0x1F3126FF25AF19628F815685CB83DA770F


def raw_rows():
    extend = lambda value: value | ((value.bit_count() & 1) << 255)
    return [extend(Q << i) for i in range(123)] + [extend(P << i) for i in range(5)]


def reduced_basis(rows):
    rows = list(rows)
    rank = 0
    for bit in range(256):
        pivot = next((i for i in range(rank, len(rows)) if rows[i] >> bit & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(len(rows)):
            if i != rank and rows[i] >> bit & 1:
                rows[i] ^= rows[rank]
        rank += 1
        if rank == len(rows):
            break
    return rows[:rank]


class Graph:
    def __init__(self, inputs):
        self.inputs = inputs
        self.gates = []
        self.forms = [1 << i for i in range(inputs)]
        self.cache = {}
        self.by_form = {value: i for i, value in enumerate(self.forms)}

    def xor(self, a, b):
        if a == -1:
            return b
        if b == -1:
            return a
        if a == b:
            return -1
        form = self.forms[a] ^ self.forms[b]
        if form in self.by_form:
            return self.by_form[form]
        key = tuple(sorted((a, b)))
        if key not in self.cache:
            self.cache[key] = len(self.forms)
            self.gates.append(key)
            self.forms.append(form)
            self.by_form[form] = self.cache[key]
        return self.cache[key]

    def total(self, terms):
        terms = list(terms)
        if not terms:
            return -1
        while len(terms) > 1:
            terms = [self.xor(terms[i], terms[i + 1]) if i + 1 < len(terms)
                     else terms[i] for i in range(0, len(terms), 2)]
        return terms[0]

    def form(self, signal):
        return 0 if signal == -1 else self.forms[signal]


def multiply(graph, values, coefficient, leaf):
    """Power-of-two fixed-coefficient Karatsuba; zero and monomials are free."""
    n = len(values)
    assert n and not n & (n - 1) and coefficient < 1 << n
    out = [-1] * (2 * n - 1)
    if not coefficient or all(x == -1 for x in values):
        return out
    if coefficient.bit_count() == 1:
        shift = coefficient.bit_length() - 1
        out[shift:shift + n] = values
        return out
    if n <= leaf:
        for i in range(2 * n - 1):
            out[i] = graph.total(values[j] for j in range(max(0, i - n + 1), min(n, i + 1))
                                 if coefficient >> (i - j) & 1)
        return out
    half = n // 2
    low = coefficient & ((1 << half) - 1)
    high = coefficient >> half
    a = multiply(graph, values[:half], low, leaf)
    b = multiply(graph, values[half:], high, leaf)
    mixed = [graph.xor(values[i], values[half + i]) for i in range(half)]
    c = multiply(graph, mixed, low ^ high, leaf)
    for i in range(n - 1):
        out[i] = graph.xor(out[i], a[i])
        out[half + i] = graph.xor(out[half + i], graph.total((a[i], b[i], c[i])))
        out[n + i] = graph.xor(out[n + i], b[i])
    return out


def product(graph, values, coefficient, leaf):
    n = 1 << (len(values) - 1).bit_length()
    values = values + [-1] * (n - len(values))
    result = [-1] * (len(values) + coefficient.bit_length() - 1)
    for shift in range(0, coefficient.bit_length(), n):
        part = multiply(graph, values, coefficient >> shift & ((1 << n) - 1), leaf)
        for j, signal in enumerate(part):
            if shift + j < len(result):
                result[shift + j] = graph.xor(result[shift + j], signal)
    return result


def transpose(forward, outputs):
    graph = Graph(len(outputs))
    uses = [[] for _ in forward.forms]
    for i, value in enumerate(outputs):
        if value != -1:
            uses[value].append(i)
    for value in range(len(forward.forms) - 1, forward.inputs - 1, -1):
        signal = graph.total(uses[value])
        if signal != -1:
            for parent in forward.gates[value - forward.inputs]:
                uses[parent].append(signal)
    return graph, [graph.total(uses[i]) for i in range(forward.inputs)]


def reciprocal(poly, n):
    assert poly & 1
    inverse = 1
    for i in range(1, n):
        parity = sum(((poly >> j) & 1) * ((inverse >> (i - j)) & 1)
                     for j in range(1, i + 1)) & 1
        inverse |= parity << i
    return inverse


def build(rows, leaf, exact):
    graph = Graph(128)
    message = list(range(128))
    values = message[:123]
    if exact:
        values = product(graph, values, reciprocal(Q, 123), leaf)[:123]
    word = product(graph, values, Q, leaf)[:255]
    word += [-1] * (255 - len(word))
    word.append(graph.total(values) if Q.bit_count() & 1 else -1)
    if exact:
        # The last five reduced generators vanish on coordinates 0..122.
        for i, row in enumerate(rows[123:]):
            pivot = (row & -row).bit_length() - 1
            correction = graph.xor(message[123 + i], word[pivot])
            for j in range(256):
                if row >> j & 1:
                    word[j] = graph.xor(word[j], correction)
    else:
        for i in range(5):
            row = P << i
            row |= (row.bit_count() & 1) << 255
            for j in range(256):
                if row >> j & 1:
                    word[j] = graph.xor(word[j], message[123 + i])
    transposed, output = transpose(graph, word)
    target = rows if exact else raw_rows()
    assert [transposed.form(v) for v in output] == target
    return transposed, output


def emit(graph, outputs, name):
    lines = [f'SPIN_NOINLINE void bchTranspose4Algebraic{name}(const block* __restrict a,block* __restrict x) {{']
    computed = set()
    def visit(value):
        if value in computed:
            return
        if value < graph.inputs:
            lines.append(f'const auto v{value}=_mm512_loadu_si512(a+4*{value});')
        else:
            left, right = graph.gates[value - graph.inputs]
            visit(left)
            visit(right)
            lines.append(f'const auto v{value}=_mm512_xor_si512(v{left},v{right});')
        computed.add(value)
    for i, value in enumerate(outputs):
        assert value != -1
        visit(value)
        for lane in range(4):
            lines.append(f'x[{128 * lane + i}]=block(_mm512_extracti32x4_epi32(v{value},{lane}));')
    lines.append('}')
    return lines, sum(v >= graph.inputs for v in computed)


def emit_raw_reference():
    # Deliberately simple, untimed correctness reference for the alternative basis.
    lines = ['static constexpr std::uint64_t BchAlgebraicRawRows[128][4]={']
    for row in raw_rows():
        lines.append('{' + ','.join(f'0x{(row >> (64 * j)) & ((1 << 64) - 1):016x}ULL'
                                    for j in range(4)) + '},')
    lines += ['};',
              'SPIN_NOINLINE void bchTranspose4AlgebraicRawReference(const block* a,block* x) {',
              'for(unsigned row=0;row<128;++row) {',
              'auto value=_mm512_setzero_si512();',
              'for(unsigned column=0;column<256;++column)',
              'if((BchAlgebraicRawRows[row][column/64]>>(column%64))&1)',
              'value=_mm512_xor_si512(value,_mm512_loadu_si512(a+4*column));']
    for lane in range(4):
        lines.append(f'x[{128 * lane}+row]=block(_mm512_extracti32x4_epi32(value,{lane}));')
    lines += ['}', '}']
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--header', type=Path)
    parser.add_argument('--leaf', type=int, default=32, choices=(1, 2, 4, 8, 16, 32))
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.header is None:
        args.header = Path(__file__).resolve().parents[3] / 'spin/src/kernels/generated/BchCircuit.h'
    words = [int(x, 16) for x in re.findall(r'0x([0-9a-f]+)ULL', args.header.read_text())]
    assert len(words) == 512
    rows = [sum(words[4 * i + j] << (64 * j) for j in range(4)) for i in range(128)]
    production_basis = reduced_basis(rows)
    polynomial_basis = reduced_basis(raw_rows())
    assert len(production_basis) == len(polynomial_basis) == 128
    assert production_basis == polynomial_basis
    assert production_basis == rows  # The structured conversion uses this RREF.
    code = ['// Generated by bch_algebraic.py; exact formal checks passed.', '#include "Spin.h"',
            'namespace spin::detail::kernel {']
    counts = {}
    for name, exact in (('Raw', False), ('Exact', True)):
        graph, outputs = build(rows, args.leaf, exact)
        lines, count = emit(graph, outputs, name)
        code.extend(lines)
        counts[name] = dict(vector_xors=count, verified_forms=128, production_basis=exact)
    code.extend(emit_raw_reference())
    code.append('}')
    print(json.dumps(dict(leaf=args.leaf, raw_production_span_equal=True, rank=128,
                         implementations=counts)))
    if args.output:
        args.output.write_text('\n'.join(code) + '\n', encoding='utf-8')


if __name__ == '__main__':
    main()
