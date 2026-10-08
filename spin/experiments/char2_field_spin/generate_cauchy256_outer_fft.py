"""Generate and validate the additive-FFT AES-field GRS[256,176] outer.

Only emits text/apply_patch fragments; it never writes files. The independent
Lagrange oracle is imported from the original outer generator. No native build
or benchmark is performed by this script.
"""

from __future__ import annotations

import argparse
from collections import Counter
from functools import cache
from pathlib import Path
import random

from generate_cauchy256_outer import INV, MUL, interpolation


@cache
def normalized(level, offset):
    value = scale = 1
    for j in range(1 << level):
        value = MUL[value][offset ^ j]
        scale = MUL[scale][(1 << level) ^ j]
    return MUL[value][INV[scale]]


@cache
def fft_ops(n, offset, inverse=False, transpose=False):
    """Each tuple means x[d] ^= coefficient*x[s], over AES GF256."""
    if n == 1:
        return ()
    half = n // 2
    coefficient = normalized(half.bit_length()-1, offset)
    children = (fft_ops(half, offset, inverse) +
                tuple((d+half, s+half, c)
                      for d, s, c in fft_ops(half, offset ^ half, inverse)))
    pairs = []
    for i in range(half):
        if inverse:
            pairs.append((i+half, i, 1))
        if coefficient:
            pairs.append((i, i+half, coefficient))
        if not inverse:
            pairs.append((i+half, i, 1))
    result = children + tuple(pairs) if inverse else tuple(pairs) + children
    if transpose:
        result = tuple((s, d, c) for d, s, c in reversed(result))
    return result


def evaluate(ops, values):
    result = list(values)
    for destination, source, coefficient in ops:
        result[destination] ^= MUL[coefficient][result[source]]
    return result


def translate(values):
    n = len(values)
    return evaluate(fft_ops(n, n), evaluate(fft_ops(n, 0, True), values))


def parity_transpose(values):
    coefficients = evaluate(fft_ops(16, 0, True), values)
    return sum((evaluate(fft_ops(16, offset), coefficients)
                for offset in (48, 32, 16)), [])


def parity(values):
    transformed = [evaluate(fft_ops(16, offset, transpose=True),
                            values[first:first+16])
                   for first, offset in zip((0, 16, 32), (48, 32, 16))]
    combined = [transformed[0][i] ^ transformed[1][i] ^ transformed[2][i]
                for i in range(16)]
    return evaluate(fft_ops(16, 0, True, True), combined)


def forward(values):
    p = translate(values[:128])
    residual = [a ^ b for a, b in zip(values[128:], p[:48])]
    r = residual + parity(residual)
    return (list(values) + [a ^ b for a, b in zip(p[48:64], r[48:64])] +
            [a ^ b for a, b in zip(p[64:], translate(r))])


def transpose(values):
    q = translate(values[192:])
    q[48:] = [a ^ b for a, b in zip(q[48:], values[176:192])]
    u = [a ^ b for a, b in zip(parity_transpose(q[48:]), q[:48])]
    tail = [a ^ b for a, b in zip(values[128:176], u)]
    translated = translate(u + list(values[176:]))
    return [a ^ b for a, b in zip(values[:128], translated)] + tail


def counts(ops):
    return Counter(gfni=sum(c != 1 for _, _, c in ops), xor=len(ops))


def circuit_counts():
    translation = counts(fft_ops(128, 0, True) + fft_ops(128, 128) +
                         fft_ops(64, 0, True) + fft_ops(64, 64))
    narrow = counts(fft_ops(16, 0, True) + fft_ops(16, 48) +
                    fft_ops(16, 32) + fft_ops(16, 16))
    return {"forward": dict(translation+narrow+Counter(xor=32+48+80)),
            "transpose": dict(translation+narrow+Counter(xor=16+48+48+128))}


def load(pointer, index):
    return f"_mm512_loadu_si512({pointer}+{index})"


def store(pointer, index, value):
    return f"_mm512_storeu_si512({pointer}+{index}, {value});"


def xor(a, b):
    return f"_mm512_xor_si512({a}, {b})"


def operation(destination, source, coefficient):
    rhs = source
    if coefficient != 1:
        rhs = (f"_mm512_gf2p8mul_epi8({source}, "
               f"_mm512_set1_epi8(static_cast<char>({coefficient})))")
    return f"{destination} = {xor(destination, rhs)};"


class Generator:
    def __init__(self):
        self.lines = []
        self.emitted = set()

    def emit(self, line=""):
        self.lines.append(line)

    def leaf(self, offset, inverse, transposed=False):
        name = f"leaf16_{offset}_{int(inverse)}_{int(transposed)}"
        if name in self.emitted:
            return name
        self.emitted.add(name)
        self.emit(f"static C256_FFT_NOINLINE void {name}(__m512i* p) {{")
        for i in range(16):
            self.emit(f"    auto x{i} = {load('p', i)};")
        for d, s, c in fft_ops(16, offset, inverse, transposed):
            self.emit("    " + operation(f"x{d}", f"x{s}", c))
        for i in range(16):
            self.emit("    " + store("p", i, f"x{i}"))
        self.emit("}")
        return name

    def stage(self, half, coefficient, inverse):
        name = f"stage{half}_{coefficient}_{int(inverse)}"
        if name in self.emitted:
            return name
        self.emitted.add(name)
        self.emit(f"static C256_FFT_NOINLINE void {name}(__m512i* p) {{")
        for first in range(0, half, 8):
            self.emit("    {")
            for i in range(8):
                self.emit(f"        auto a{i} = {load('p', first+i)};")
                self.emit(f"        auto b{i} = {load('p', first+half+i)};")
            for i in range(8):
                if inverse:
                    self.emit("        " + operation(f"b{i}", f"a{i}", 1))
                if coefficient:
                    self.emit("        " + operation(f"a{i}", f"b{i}", coefficient))
                if not inverse:
                    self.emit("        " + operation(f"b{i}", f"a{i}", 1))
            for i in range(8):
                self.emit("        " + store("p", first+i, f"a{i}"))
                self.emit("        " + store("p", first+half+i, f"b{i}"))
            self.emit("    }")
        self.emit("}")
        return name

    def schedule(self, n, offset, inverse, start=0):
        if n == 16:
            return [f"    {self.leaf(offset, inverse)}(p+{start});"]
        half = n // 2
        c = normalized(half.bit_length()-1, offset)
        stage = [f"    {self.stage(half, c, inverse)}(p+{start});"]
        children = (self.schedule(half, offset, inverse, start) +
                    self.schedule(half, offset ^ half, inverse, start+half))
        return children+stage if inverse else stage+children

    def translation(self, n):
        schedule = self.schedule(n, 0, True) + self.schedule(n, n, False)
        self.emit(f"static C256_FFT_NOINLINE void translate{n}(__m512i* p) {{")
        self.lines.extend(schedule)
        self.emit("}")

    def copy(self, source, destination, count, source_first=0, destination_first=0):
        for i in range(count):
            self.emit("    " + store(destination, destination_first+i,
                                      load(source, source_first+i)))

    def add(self, a, b, destination, count, afirst=0, bfirst=0, dfirst=0):
        for i in range(count):
            self.emit("    " + store(destination, dfirst+i,
                                      xor(load(a, afirst+i), load(b, bfirst+i))))

    def parity_maps(self):
        inverse = self.leaf(0, True)
        normal = [self.leaf(offset, False) for offset in (48, 32, 16)]
        adjoints = [self.leaf(offset, False, True) for offset in (48, 32, 16)]
        adj_inverse = self.leaf(0, True, True)
        self.emit("static C256_FFT_NOINLINE void parityTranspose(const __m512i* input, __m512i* output) {")
        self.copy("input", "output", 16)
        self.emit(f"    {inverse}(output);")
        self.copy("output", "output", 16, destination_first=16)
        self.copy("output", "output", 16, destination_first=32)
        for first, name in zip((0, 16, 32), normal):
            self.emit(f"    {name}(output+{first});")
        self.emit("}")
        self.emit("static C256_FFT_NOINLINE void parityForward(const __m512i* input, __m512i* output) {")
        self.emit("    alignas(64) __m512i temporary[16];")
        self.copy("input", "output", 16)
        self.emit(f"    {adjoints[0]}(output);")
        for first, name in zip((16, 32), adjoints[1:]):
            self.copy("input", "temporary", 16, source_first=first)
            self.emit(f"    {name}(temporary);")
            self.add("output", "temporary", "output", 16)
        self.emit(f"    {adj_inverse}(output);")
        self.emit("}")

    def wrappers(self):
        self.emit("} // namespace c256_outer_fft_detail")
        self.emit("static C256_FFT_NOINLINE void c256OuterForwardFft(const __m512i* input176, __m512i* output256) {")
        self.emit("    alignas(64) __m512i p[128], r[64];")
        self.copy("input176", "p", 128)
        self.emit("    c256_outer_fft_detail::translate128(p);")
        self.add("input176", "p", "r", 48, afirst=128)
        self.emit("    c256_outer_fft_detail::parityForward(r, r+48);")
        self.copy("input176", "output256", 176)
        self.add("p", "r", "output256", 16, afirst=48, bfirst=48, dfirst=176)
        self.emit("    c256_outer_fft_detail::translate64(r);")
        self.add("p", "r", "output256", 64, afirst=64, dfirst=192)
        self.emit("}")
        self.emit("static C256_FFT_NOINLINE void c256OuterTransposeFft(const __m512i* input256, __m512i* output176) {")
        self.emit("    alignas(64) __m512i p[128], q[64];")
        self.copy("input256", "q", 64, source_first=192)
        self.emit("    c256_outer_fft_detail::translate64(q);")
        self.add("q", "input256", "q", 16, afirst=48, bfirst=176, dfirst=48)
        self.emit("    c256_outer_fft_detail::parityTranspose(q+48, p);")
        self.add("p", "q", "p", 48)
        self.add("input256", "p", "output176", 48, afirst=128, dfirst=128)
        self.copy("input256", "p", 80, source_first=176, destination_first=48)
        self.emit("    c256_outer_fft_detail::translate128(p);")
        self.add("input256", "p", "output176", 128)
        self.emit("}")


def generate():
    generator = Generator()
    generator.lines.extend([
        "// Generated by generate_cauchy256_outer_fft.py; do not edit manually.",
        "// Exact fixed systematic AES GF256 GRS[256,176], source points 0..175.",
        "// Each vector holds four independent 128-bit records; coefficients act per byte.",
        "// Distinct nonoverlapping public input/output ranges; unaligned pointers accepted.",
        "// Per call: 1203 GFNI; forward 2771 XOR, transpose 2851 XOR.",
        "// Explicit nested scratch: forward 13 KiB, transpose 12 KiB, plus compiler frames.",
        "// Ordinary noinline functions only: do not inline these arrays into a coroutine.",
        "// Sixteen-vector leaves / eight-pair stages bound the live SIMD working set.",
        "#pragma once", "#include <immintrin.h>",
        "#if defined(_MSC_VER)",
        "#define C256_FFT_NOINLINE __declspec(noinline)", "#else",
        "#define C256_FFT_NOINLINE __attribute__((noinline))", "#endif",
        "namespace spin::research::char2 {", "namespace c256_outer_fft_detail {",
    ])
    generator.translation(64)
    generator.translation(128)
    generator.parity_maps()
    generator.wrappers()
    generator.lines.extend(["} // namespace spin::research::char2", "#undef C256_FFT_NOINLINE"])
    chunks = [generator.lines[first:first+400]
              for first in range(0, len(generator.lines), 400)]
    parts = [chunk + [f"// c256-fft-generated-part {i}"] for i, chunk in enumerate(chunks)]
    return "\n".join(sum(parts, [])) + "\n", parts


def dot(a, b):
    result = 0
    for x, y in zip(a, b):
        result ^= MUL[x][y]
    return result


def check(header):
    total = 0
    matrix = interpolation(176, range(176, 256))
    for source in range(176):
        values = [0]*176
        values[source] = 1
        assert forward(values) == values+[row[source] for row in matrix]
        total += 1
    for source in range(256):
        values = [0]*256
        values[source] = 1
        expected = ([int(j == source) for j in range(176)] if source < 176
                    else matrix[source-176])
        assert transpose(values) == expected
        total += 1
    for n in (64, 128):
        translation_matrix = interpolation(n, range(n, 2*n))
        assert all(translation_matrix[i][j] == translation_matrix[j][i]
                   for i in range(n) for j in range(n))
        for source in range(n):
            values = [0]*n
            values[source] = 1
            actual = translate(values)
            assert actual == [row[source] for row in translation_matrix]
            assert translate(actual) == values
            total += 1
    rng = random.Random(256176)
    for _ in range(16):
        a = [rng.randrange(256) for _ in range(176)]
        b = [rng.randrange(256) for _ in range(256)]
        encoded = forward(a)
        dual = transpose(b)
        assert encoded == a+[dot(row, a) for row in matrix]
        assert dual == [b[j] ^ dot([row[j] for row in matrix], b[176:])
                        for j in range(176)]
        assert dot(encoded, b) == dot(a, dual)
        total += 3
    assert circuit_counts() == {"forward": {"gfni": 1203, "xor": 2771},
                                "transpose": {"gfni": 1203, "xor": 2851}}
    path = Path(__file__).with_name("Cauchy256OuterFft.h")
    if path.exists():
        assert path.read_text() == header, "generated header is incomplete or stale"
    print(f"PASS {total} exact basis/random/adjoint checks; counts={circuit_counts()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--emit", action="store_true")
    parser.add_argument("--patch", type=int)
    args = parser.parse_args()
    header, parts = generate()
    if args.patch is not None:
        number = args.patch
        target = Path(__file__).with_name("Cauchy256OuterFft.h").as_posix()
        print("*** Begin Patch")
        print(f"*** {'Add' if number == 0 else 'Update'} File: {target}")
        if number:
            print("@@")
            print(f" // c256-fft-generated-part {number-1}")
        print("\n".join("+"+line for line in parts[number]))
        print("*** End Patch")
    elif args.emit:
        print(header, end="")
    elif args.check:
        check(header)
    else:
        print(f"parts={len(parts)} lines={len(header.splitlines())} {circuit_counts()}")
