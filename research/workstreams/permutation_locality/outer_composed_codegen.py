"""Exact packed-byte polynomial BCH prototype, not a new code construction.

The retained GL32 preparation remains separate: fully composing its four row
lanes into every BCH block would multiply the dense GFNI work. Instead this
prototype factors the fixed BCH matrix through polynomial convolution. Its
Karatsuba leaves are 8-by-8 GFNI maps, not coordinate-wise XOR circuits.

Usage: python outer_composed_codegen.py BchCircuit.h --analyze-only
       python outer_composed_codegen.py BchCircuit.h > OuterPolynomial.cpp
       python outer_composed_codegen.py BchCircuit.h --basis raw > OuterRaw.cpp

Generation checks all 128 output linear forms against all 256 input columns.
The C++ replacement exports the ordinary coefficient entry points; only Compact
uses the prototype. Compile like the selected tile5/layout2 outer object.
With --basis raw, the retained exports remain unchanged and the DIFFERENT
message-coordinate map is exported as bchPackedCoeffRawPolynomial. Its BCH
generator has the same rank128 row span as the original, verified exactly;
its transpose outputs are not interchangeable with the original outputs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path
import re
import subprocess
import sys

import outer_layout_codegen as layout
import packed_bch_tune_codegen as tune
from packed_coeff_codegen import definition


Q = 0x1F3126FF25AF19628F815685CB83DA770F
P = 0x1C3B42FE115AA90747020D79AC738ADB
IDENTITY = 0x0102040810204080


def row(matrix, bit):
    return (matrix >> (8 * (7 - bit))) & 255


def matrix_from_rows(rows):
    return sum(value << (8 * (7 - bit)) for bit, value in enumerate(rows))


def transpose_matrix(matrix):
    return matrix_from_rows([
        sum(((row(matrix, j) >> bit) & 1) << j for j in range(8))
        for bit in range(8)])


def compose(left, right):
    result = []
    for bit in range(8):
        mask, value = row(left, bit), 0
        for j in range(8):
            if mask >> j & 1:
                value ^= row(right, j)
        result.append(value)
    return matrix_from_rows(result)


class Graph:
    """Byte-valued linear graph with exact binary forms and affine folding."""
    def __init__(self, inputs):
        self.inputs = inputs
        self.nodes = [('input', i) for i in range(inputs)]
        self.forms = [tuple(1 << (8 * i + j) for j in range(8))
                      for i in range(inputs)]
        self.by_form = {form: i for i, form in enumerate(self.forms)}

    def form(self, signal):
        return (0,) * 8 if signal == -1 else self.forms[signal]

    def add(self, form, operation):
        if not any(form):
            return -1
        if form in self.by_form:
            return self.by_form[form]
        result = len(self.nodes)
        self.nodes.append(operation)
        self.forms.append(form)
        self.by_form[form] = result
        return result

    def affine(self, signal, matrix):
        if signal == -1 or matrix == 0:
            return -1
        if matrix == IDENTITY:
            return signal
        previous = self.nodes[signal]
        if previous[0] == 'gfni':
            return self.affine(previous[1], compose(matrix, previous[2]))
        source = self.forms[signal]
        result = []
        for bit in range(8):
            value = 0
            for j in range(8):
                if row(matrix, bit) >> j & 1:
                    value ^= source[j]
            result.append(value)
        return self.add(tuple(result), ('gfni', signal, matrix))

    def xor(self, left, right):
        if left == -1:
            return right
        if right == -1:
            return left
        if left == right:
            return -1
        return self.add(tuple(a ^ b for a, b in zip(self.forms[left], self.forms[right])),
                        ('xor', min(left, right), max(left, right)))

    def total(self, signals):
        signals = list(signals)
        if not signals:
            return -1
        while len(signals) > 1:
            signals = [self.xor(signals[i], signals[i + 1])
                       if i + 1 < len(signals) else signals[i]
                       for i in range(0, len(signals), 2)]
        return signals[0]

    def reachable(self, outputs):
        seen = set()
        def visit(signal):
            if signal == -1 or signal in seen:
                return
            seen.add(signal)
            node = self.nodes[signal]
            if node[0] == 'gfni':
                visit(node[1])
            elif node[0] == 'xor':
                visit(node[1]); visit(node[2])
        for output in outputs:
            visit(output)
        return sorted(seen)


def multiply(graph, values, coefficient, leaf):
    """Fixed polynomial product in radix x^8, returning twice as many bytes."""
    n = len(values)
    assert n and not n & (n - 1) and coefficient < 1 << (8 * n)
    if not coefficient or all(value == -1 for value in values):
        return [-1] * (2 * n)
    if coefficient == 1:
        return values + [-1] * n
    if n <= leaf:
        result = []
        for byte in range(2 * n):
            terms = []
            for source, value in enumerate(values):
                matrix = matrix_from_rows([
                    sum(((coefficient >> (8 * (byte - source) + out - inp)) & 1) << inp
                        for inp in range(8) if 0 <= 8 * (byte - source) + out - inp < 8 * n)
                    for out in range(8)])
                terms.append(graph.affine(value, matrix))
            result.append(graph.total(terms))
        return result
    half = n // 2
    low = coefficient & ((1 << (8 * half)) - 1)
    high = coefficient >> (8 * half)
    a = multiply(graph, values[:half], low, leaf)
    b = multiply(graph, values[half:], high, leaf)
    c = multiply(graph, [graph.xor(values[i], values[half + i])
                         for i in range(half)], low ^ high, leaf)
    result = [-1] * (2 * n)
    for i in range(n):
        result[i] = graph.xor(result[i], a[i])
        result[half + i] = graph.xor(result[half + i], graph.total((a[i], b[i], c[i])))
        result[n + i] = graph.xor(result[n + i], b[i])
    return result


def product(graph, values, coefficient, leaf):
    n = 1 << (len(values) - 1).bit_length()
    values = values + [-1] * (n - len(values))
    count = (8 * n + coefficient.bit_length() - 1 + 7) // 8
    result = [-1] * count
    for offset in range(0, coefficient.bit_length(), 8 * n):
        part = multiply(graph, values, (coefficient >> offset) & ((1 << (8 * n)) - 1), leaf)
        for j, signal in enumerate(part):
            if offset // 8 + j < count:
                result[offset // 8 + j] = graph.xor(result[offset // 8 + j], signal)
    return result


def reciprocal(poly, count):
    assert poly & 1
    result = 1
    for bit in range(1, count):
        parity = sum(((poly >> j) & 1) * ((result >> (bit - j)) & 1)
                     for j in range(1, bit + 1)) & 1
        result |= parity << bit
    return result


def mask_bits(graph, signal, mask):
    return graph.affine(signal, matrix_from_rows([
        (1 << bit) if mask >> bit & 1 else 0 for bit in range(8)]))


def raw_rows():
    def extend(value):
        return value | ((value.bit_count() & 1) << 255)
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


def build_forward(rows, basis, leaf):
    graph = Graph(16)
    values = list(range(15)) + [mask_bits(graph, 15, 7)]
    if basis == 'exact':
        values = product(graph, values, reciprocal(Q, 123), leaf)[:16]
        values[15] = mask_bits(graph, values[15], 7)
    word = product(graph, values, Q, leaf)[:32]
    word += [-1] * (32 - len(word))
    word[31] = mask_bits(graph, word[31], 127)
    if Q.bit_count() & 1:
        parity = graph.affine(graph.total(values), matrix_from_rows([0] * 7 + [255]))
        word[31] = graph.xor(word[31], parity)
    # Gather the five systematic corrections into bits0..4 of one byte.
    # Their pivots include coordinate129, not coordinate127.
    if basis == 'exact':
        corrections = []
        for j, form in enumerate(rows[123:]):
            pivot = (form & -form).bit_length() - 1
            message = graph.affine(15, matrix_from_rows([
                (1 << (3 + j)) if bit == j else 0 for bit in range(8)]))
            previous = graph.affine(word[pivot // 8], matrix_from_rows([
                (1 << (pivot % 8)) if bit == j else 0 for bit in range(8)]))
            corrections.append(graph.xor(message, previous))
        correction = graph.total(corrections)
    else:
        correction = graph.affine(15, matrix_from_rows([
            (1 << (3 + bit)) if bit < 5 else 0 for bit in range(8)]))
    for byte in range(32):
        matrix = matrix_from_rows([
            sum(((rows[123 + j] >> (8 * byte + bit)) & 1) << j for j in range(5))
            for bit in range(8)])
        word[byte] = graph.xor(word[byte], graph.affine(correction, matrix))
    expected = [sum(((rows[j] >> bit) & 1) << j for j in range(128))
                for bit in range(256)]
    actual = [form for value in word for form in graph.form(value)]
    assert actual == expected, 'forward polynomial construction differs from BCH'
    return graph, word


def transposed(forward, outputs):
    graph = Graph(len(outputs))
    uses = [[] for _ in forward.nodes]
    for i, signal in enumerate(outputs):
        if signal != -1:
            uses[signal].append(i)
    for signal in range(len(forward.nodes) - 1, forward.inputs - 1, -1):
        value = graph.total(uses[signal])
        node = forward.nodes[signal]
        if value == -1:
            continue
        if node[0] == 'xor':
            uses[node[1]].append(value); uses[node[2]].append(value)
        else:
            assert node[0] == 'gfni'
            uses[node[1]].append(graph.affine(value, transpose_matrix(node[2])))
    return graph, [graph.total(uses[i]) for i in range(forward.inputs)]


def build(header, basis, leaf):
    words = [int(value, 16) for value in re.findall(r'0x([0-9a-f]+)ULL', header.read_text())]
    assert len(words) == 512
    rows = [sum(words[4 * i + j] << (64 * j) for j in range(4)) for i in range(128)]
    original = rows
    if basis == 'raw':
        rows = raw_rows()
        canonical = reduced_basis(original)
        assert len(canonical) == 128 and canonical == reduced_basis(rows), 'raw code image differs'
    forward, outputs = build_forward(rows, basis, leaf)
    graph, outputs = transposed(forward, outputs)
    assert [form for value in outputs for form in graph.form(value)] == rows
    return graph, outputs, rows


def emit_graph(graph, outputs):
    print('static SPIN_NOINLINE void polynomialPlane(const __m512i* __restrict src,'
          '__m512i* __restrict packed,unsigned plane) {')
    uses = Counter(outputs)
    for signal in graph.reachable(outputs):
        node = graph.nodes[signal]
        if node[0] == 'gfni':
            uses[node[1]] += 1
        elif node[0] == 'xor':
            uses[node[1]] += 1; uses[node[2]] += 1
    emitted = set()
    def visit(signal):
        if signal in emitted:
            return
        node = graph.nodes[signal]
        if node[0] == 'input':
            expression = f'_mm512_load_si512(src+8*{node[1]}+plane)'
        elif node[0] == 'xor':
            parents = list(node[1:])
            for position, parent in enumerate(parents):
                if uses[parent] == 1 and graph.nodes[parent][0] == 'xor':
                    parents[position:position + 1] = graph.nodes[parent][1:]
                    break
            for parent in parents:
                visit(parent)
            if len(parents) == 3:
                expression = '_mm512_ternarylogic_epi64(' + ','.join(f'v{i}' for i in parents) + ',0x96)'
            else:
                expression = f'_mm512_xor_si512(v{parents[0]},v{parents[1]})'
        else:
            visit(node[1])
            expression = (f'_mm512_gf2p8affine_epi64_epi8(v{node[1]},'
                          f'_mm512_set1_epi64(0x{node[2]:016x}ULL),0)')
        print(f'const auto v{signal}={expression};')
        emitted.add(signal)
    for i, signal in enumerate(outputs):
        visit(signal)
        print(f'_mm512_store_si512(packed+8*{i}+plane,v{signal});')
    print('}')


def verified_graph_source(graph, outputs):
    """Replay emitted C++ wiring on 256-bit linear forms, including ternaries."""
    captured = StringIO()
    with redirect_stdout(captured):
        emit_graph(graph, outputs)
    source = captured.getvalue()
    values, actual = {}, {}
    for line in source.splitlines():
        declaration = re.fullmatch(r'const auto v(\d+)=(.*);', line)
        if declaration:
            signal, expression = int(declaration[1]), declaration[2]
            loaded = re.fullmatch(r'_mm512_load_si512\(src\+8\*(\d+)\+plane\)', expression)
            affine = re.fullmatch(r'_mm512_gf2p8affine_epi64_epi8\(v(\d+),'
                                  r'_mm512_set1_epi64\(0x([0-9a-f]+)ULL\),0\)', expression)
            logical = re.fullmatch(r'_mm512_(xor_si512|ternarylogic_epi64)\((.*)\)', expression)
            if loaded:
                value = tuple(1 << (8 * int(loaded[1]) + bit) for bit in range(8))
            elif affine:
                parent, matrix = values[int(affine[1])], int(affine[2], 16)
                value = []
                for bit in range(8):
                    result = 0
                    for j in range(8):
                        if row(matrix, bit) >> j & 1:
                            result ^= parent[j]
                    value.append(result)
                value = tuple(value)
            else:
                assert logical, 'unrecognized emitted graph instruction'
                args = logical[2].split(',')
                if logical[1] == 'ternarylogic_epi64':
                    assert args.pop() == '0x96' and len(args) == 3
                else:
                    assert len(args) == 2
                parents = [values[int(arg[1:])] for arg in args]
                value = list(parents[0])
                for parent in parents[1:]:
                    value = [a ^ b for a, b in zip(value, parent)]
                value = tuple(value)
            values[signal] = value
        stored = re.fullmatch(r'_mm512_store_si512\(packed\+8\*(\d+)\+plane,v(\d+)\);', line)
        if stored:
            output = int(stored[1])
            assert output not in actual
            actual[output] = values[int(stored[2])]
    assert sorted(actual) == list(range(16))
    assert [actual[i] for i in range(16)] == [graph.form(signal) for signal in outputs]
    return source


def emit_raw_reference(rows):
    print('static constexpr std::uint64_t rawPolynomialRows[128][4]={')
    for form in rows:
        print('{' + ','.join(f'0x{(form >> (64*j)) & ((1 << 64)-1):016x}ULL'
                            for j in range(4)) + '},')
    print('};')
    print('SPIN_NOINLINE void bchPackedCoeffRawPolynomialScalar(const block* a,'
          'block* out,const std::uint64_t* coeff) {')
    print('alignas(64) block mixed[1024];')
    print('for(unsigned group=0;group<32;++group)for(unsigned r=0;r<32;++r) {')
    print('const unsigned lane=r/8,j=r%8;std::uint32_t mask=0;')
    print('for(unsigned d=0;d<4;++d)mask|=std::uint32_t('
          '(coeff[16*group+4*d+lane]>>(8*(7-j)))&255)<<(8*((lane+d)%4));')
    print('auto value=_mm_setzero_si128();for(unsigned c=0;c<32;++c)if(mask>>c&1)'
          'value=_mm_xor_si128(value,a[4*(8*group+c%8)+c/8].mData);')
    print('mixed[4*(8*group+j)+lane]=block(value);}\n')
    print('for(unsigned lane=0;lane<4;++lane)for(unsigned r=0;r<128;++r) {')
    print('auto value=_mm_setzero_si128();for(unsigned c=0;c<256;++c)'
          'if(rawPolynomialRows[r][c/64]>>(c%64)&1)'
          'value=_mm_xor_si128(value,mixed[4*c+lane].mData);')
    print('out[128*lane+r]=block(value);}\n}')


def emit(header, graph, outputs, basis, rows):
    reference = subprocess.run([sys.executable, str(Path(layout.__file__).resolve()),
                               str(header), '--tile-mode', '5', '--mode', '2'],
                              check=True, capture_output=True, text=True)
    if reference.stderr:
        print(reference.stderr, file=sys.stderr, end='')
    old = definition(reference.stdout, 'SPIN_NOINLINE void bchPackedCoeffCompact(')
    replacement = StringIO()
    with redirect_stdout(replacement):
        print('namespace coeff_polynomial_detail {')
        print(verified_graph_source(graph, outputs), end='')
        print('}')
        name = 'bchPackedCoeffCompact' if basis == 'exact' else 'bchPackedCoeffRawPolynomial'
        print(f'SPIN_NOINLINE void {name}(const block* __restrict a,'
              'block* __restrict out,const std::uint64_t* __restrict coeff) {')
        print('using namespace coeff_probe_detail;')
        layout.emit_preparation(2, packed_layout=True)
        print('alignas(64) __m512i packed[128];')
        print('for(unsigned plane=0;plane<8;++plane)'
              'coeff_polynomial_detail::polynomialPlane(src,packed,plane);')
        print('for(unsigned output=0;output<16;output+=2) {')
        for p in range(2):
            for j in range(8):
                print(f'auto y{8*p+j}=_mm512_load_si512(packed+8*(output+{p})+{j});')
        tune.emit_output(2, packed_layout=True)
        print('}\n}')
        if basis == 'raw':
            emit_raw_reference(rows)
    if basis == 'exact':
        print(reference.stdout.replace(old, replacement.getvalue().rstrip()), end='')
    else:
        print(reference.stdout, end='')
        print('namespace spin::detail::kernel {')
        print(replacement.getvalue(), end='')
        print('}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('header', type=Path)
    parser.add_argument('--basis', choices=('exact', 'raw'), default='exact')
    parser.add_argument('--leaf-bytes', type=int, choices=(1, 2, 4, 8, 16), default=2)
    parser.add_argument('--analyze-only', action='store_true')
    args = parser.parse_args()
    graph, outputs, rows = build(args.header, args.basis, args.leaf_bytes)
    counts = Counter(graph.nodes[i][0] for i in graph.reachable(outputs))
    report = dict(verified_output_forms=128, verified_input_coordinates=256,
                  per_plane=dict(counts), per_four_row_tile={key: value * 8 for key, value in counts.items()},
                  retained_dense_gfni_per_tile=2048, packed_output_scratch_bytes=8192,
                  basis=args.basis, leaf_bytes=args.leaf_bytes, code_image_rank=128,
                  scope=('exact fixed BCH map' if args.basis == 'exact' else
                         'different message-coordinate map; identical BCH code image'),
                  gl32_setup_unchanged=True, measured_speedup=False)
    print(json.dumps(report, sort_keys=True), file=sys.stderr)
    if not args.analyze_only:
        emit(args.header, graph, outputs, args.basis, rows)


if __name__ == '__main__':
    main()
