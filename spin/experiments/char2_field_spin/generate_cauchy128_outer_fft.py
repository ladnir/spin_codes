"""Generate/check the exact AES-field systematic GRS[128,87] outer.

The prefix decomposition is 87=64+23, 23=16+7. Completing seven evaluations
to eight needs only their XOR. The header has full forward/transpose hooks
and a parity-only transpose hook for a four-row-interleaved frontend.
This script emits text/apply_patch fragments; it performs no native work.
"""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from pathlib import Path
import random

from generate_cauchy256_outer import basis_for, interpolation
from generate_cauchy256_outer_fft import (Generator, counts, dot, evaluate,
                                         fft_ops, load, operation, store,
                                         translate, xor)


def parity23(values):
    p = translate(values[:16])
    residual = [a ^ b for a, b in zip(values[16:], p[:7])]
    last = 0
    for value in residual:
        last ^= value
    residual.append(last)
    return [p[7] ^ last]+[a ^ b for a, b in zip(p[8:], translate(residual))]


def parity23_transpose(values):
    q = translate(values[1:])
    q[7] ^= values[0]
    u = [value ^ q[7] for value in q[:7]]
    return translate(u+list(values))+u


def forward(values):
    p = translate(values[:64])
    residual = [a ^ b for a, b in zip(values[64:], p[:23])]
    residual += parity23(residual)
    return (list(values)+[a ^ b for a, b in zip(p[23:32], residual[23:])] +
            [a ^ b for a, b in zip(p[32:], translate(residual))])


def parity_transpose(values):
    q = translate(values[9:])
    q[23:] = [a ^ b for a, b in zip(q[23:], values[:9])]
    u = [a ^ b for a, b in zip(parity23_transpose(q[23:]), q[:23])]
    return translate(u+list(values))+u


def transpose(values):
    return [a ^ b for a, b in zip(values[:87], parity_transpose(values[87:]))]


def circuit_counts():
    translations = sum((counts(fft_ops(n, 0, True)+fft_ops(n, n))
                        for n in (8, 16, 32, 64)), Counter())
    return {"forward": dict(translations+Counter(xor=7+6+9+23+41)),
            "transpose": dict(translations+Counter(xor=1+7+9+23+87)),
            "parity_transpose": dict(translations+Counter(xor=1+7+9+23))}


class C128Generator(Generator):
    def translation(self, n):
        if n != 8:
            return super().translation(n)
        self.emit("static C256_FFT_NOINLINE void translate8(__m512i* p) {")
        for i in range(8):
            self.emit(f"    auto x{i} = {load('p', i)};")
        for d, s, c in fft_ops(8, 0, True)+fft_ops(8, 8):
            self.emit("    "+operation(f"x{d}", f"x{s}", c))
        for i in range(8):
            self.emit("    "+store("p", i, f"x{i}"))
        self.emit("}")

    def small_maps(self):
        self.emit("static C256_FFT_NOINLINE void parity23Forward(const __m512i* input, __m512i* output) {")
        self.emit("    alignas(64) __m512i p[16], r[8];")
        self.copy("input", "p", 16)
        self.emit("    translate16(p);")
        self.add("input", "p", "r", 7, afirst=16)
        # Balanced six-XOR tree; operation accounting uses binary XORs.
        self.emit(f"    const auto x01 = {xor(load('r', 0), load('r', 1))};")
        self.emit(f"    const auto x23 = {xor(load('r', 2), load('r', 3))};")
        self.emit(f"    const auto x45 = {xor(load('r', 4), load('r', 5))};")
        self.emit(f"    const auto x0123 = {xor('x01', 'x23')};")
        self.emit(f"    const auto x456 = {xor('x45', load('r', 6))};")
        self.emit("    "+store("r", 7, xor("x0123", "x456")))
        self.add("p", "r", "output", 1, afirst=7, bfirst=7)
        self.emit("    translate8(r);")
        self.add("p", "r", "output", 8, afirst=8, dfirst=1)
        self.emit("}")
        self.emit("static C256_FFT_NOINLINE void parity23Transpose(const __m512i* input, __m512i* output) {")
        self.emit("    alignas(64) __m512i q[8];")
        self.copy("input", "q", 8, source_first=1)
        self.emit("    translate8(q);")
        self.add("q", "input", "q", 1, afirst=7, dfirst=7)
        for i in range(7):
            self.emit(f"    const auto u{i} = {xor(load('q', i), load('q', 7))};")
            self.emit("    "+store("output", i, f"u{i}"))
            self.emit("    "+store("output", i+16, f"u{i}"))
        self.copy("input", "output", 9, destination_first=7)
        self.emit("    translate16(output);")
        self.emit("}")

    def wrappers128(self):
        self.emit("} // namespace c128_outer_fft_detail")
        self.emit("static C256_FFT_NOINLINE void c128OuterForwardFft(const __m512i* input87, __m512i* output128) {")
        self.emit("    alignas(64) __m512i p[64], r[32];")
        self.copy("input87", "p", 64)
        self.emit("    c128_outer_fft_detail::translate64(p);")
        self.add("input87", "p", "r", 23, afirst=64)
        self.emit("    c128_outer_fft_detail::parity23Forward(r, r+23);")
        self.copy("input87", "output128", 87)
        self.add("p", "r", "output128", 9, afirst=23, bfirst=23, dfirst=87)
        self.emit("    c128_outer_fft_detail::translate32(r);")
        self.add("p", "r", "output128", 32, afirst=32, dfirst=96)
        self.emit("}")
        self.emit("// Input41 are parity points87..127. Systematic identity is NOT included.")
        self.emit("static C256_FFT_NOINLINE void c128OuterParityTransposeFft(const __m512i* input41, __m512i* output87) {")
        self.emit("    alignas(64) __m512i q[32];")
        self.copy("input41", "q", 32, source_first=9)
        self.emit("    c128_outer_fft_detail::translate32(q);")
        self.add("q", "input41", "q", 9, afirst=23, dfirst=23)
        self.emit("    c128_outer_fft_detail::parity23Transpose(q+23, output87);")
        self.add("output87", "q", "output87", 23)
        self.copy("output87", "output87", 23, destination_first=64)
        self.copy("input41", "output87", 41, destination_first=23)
        self.emit("    c128_outer_fft_detail::translate64(output87);")
        self.emit("}")
        self.emit("static C256_FFT_NOINLINE void c128OuterTransposeFft(const __m512i* input128, __m512i* output87) {")
        self.emit("    c128OuterParityTransposeFft(input128+87, output87);")
        self.add("input128", "output87", "output87", 87)
        self.emit("}")


def generate():
    g = C128Generator()
    g.lines.extend([
        "// Generated by generate_cauchy128_outer_fft.py; do not edit manually.",
        "// Fixed systematic AES GF256 GRS[128,87], source points0..86.",
        "// Each vector contains four independent 128-bit records, one per outer row.",
        "// GF256 constants act independently on every byte; no large-field products.",
        "// Public input/output ranges must be distinct/nonoverlapping; unaligned accepted.",
        "// Per call:516 GFNI. XOR counts:forward1234, transpose1275, parity-only1188.",
        "// Explicit nested scratch:forward7.5KiB; transpose/parity2.5KiB, plus compiler frames.",
        "// Parity hook uses its output head as workspace. No caller scratch or heap.",
        "// Intended frontend: add87 systematic planes, then store87 aligned vectors in",
        "// FOUR-ROW-INTERLEAVED source order. The hook itself does not reorder records.",
        "// All arrays live in ordinary noinline helpers, never coroutine frames.",
        "#pragma once", "#include <immintrin.h>", "#if defined(_MSC_VER)",
        "#define C256_FFT_NOINLINE __declspec(noinline)", "#else",
        "#define C256_FFT_NOINLINE __attribute__((noinline))", "#endif",
        "namespace spin::research::char2 {", "namespace c128_outer_fft_detail {",
    ])
    for n in (8, 16, 32, 64):
        g.translation(n)
    g.small_maps()
    g.wrappers128()
    g.lines.extend(["} // namespace spin::research::char2", "#undef C256_FFT_NOINLINE"])
    lines = [line.replace("C256_FFT_NOINLINE", "C128_FFT_NOINLINE") for line in g.lines]
    chunks = [lines[first:first+400] for first in range(0, len(lines), 400)]
    parts = [chunk+[f"// c128-fft-generated-part {i}"] for i, chunk in enumerate(chunks)]
    return "\n".join(sum(parts, []))+"\n", parts


def cost_model():
    correction = interpolation(87, [87])[0]
    basis, span = basis_for(correction)
    assert len(basis)-1 == 5
    reconstruction = sum(span[c].bit_count()-1 for c in set(correction))
    assert reconstruction == 13
    output = {}
    for length, dimension, products, boolean in ((256, 176, 1203, 2851),
                                                (128, 87, 516, 1275)):
        rows = 64*((2**20+64*dimension-1)//(64*dimension))
        outer = Fraction(products, 4*dimension)
        inner = Fraction(70*length, 16*dimension)
        output[f"{length}/{dimension}"] = dict(
            rows=rows, source_pad=rows*dimension, encoded=rows*length,
            outer_products_per_source=str(outer), dense_inner_products_per_source=str(inner),
            total_products_per_source=str(outer+inner),
            total_boolean_per_source=str(Fraction(boolean, 4*dimension)+
                                          Fraction(38*length, 16*dimension)))
    output["via88"] = dict(shortening_basis=basis, shortening_coefficients=[correction[8*j] for j in range(11)],
                            fft_products=491+5, parity_binary_xor=1115+88+reconstruction,
                            coefficient_basis_products=232+5,
                            arithmetic_only=True)
    return output


def check(header):
    checks = 0
    matrix = interpolation(87, range(87, 128))
    small = interpolation(23, range(23, 32))
    for source in range(23):
        x = [int(i == source) for i in range(23)]
        assert parity23(x) == [row[source] for row in small]
        checks += 1
    for source in range(9):
        x = [int(i == source) for i in range(9)]
        assert parity23_transpose(x) == small[source]
        checks += 1
    for source in range(87):
        x = [int(i == source) for i in range(87)]
        assert forward(x) == x+[row[source] for row in matrix]
        checks += 1
    for source in range(128):
        x = [int(i == source) for i in range(128)]
        expected = ([int(i == source) for i in range(87)] if source < 87 else matrix[source-87])
        assert transpose(x) == expected
        checks += 1
    for source in range(41):
        x = [int(i == source) for i in range(41)]
        assert parity_transpose(x) == matrix[source]
        checks += 1
    # Complementary interpolation supplies a second exact description of P^T.
    dual = interpolation(41, range(41, 128))
    assert all(matrix[j][i] == dual[86-i][40-j] for i in range(87) for j in range(41))
    checks += 87*41
    rng = random.Random(12887)
    for _ in range(16):
        x = [rng.randrange(256) for _ in range(87)]
        y = [rng.randrange(256) for _ in range(128)]
        fx, ty = forward(x), transpose(y)
        assert fx == x+[dot(row, x) for row in matrix]
        assert ty == [y[i] ^ dot([row[i] for row in matrix], y[87:]) for i in range(87)]
        assert dot(fx, y) == dot(x, ty)
        assert parity_transpose(y[87:]) == [a ^ b for a, b in zip(ty, y[:87])]
        checks += 4
    expected = {"forward": {"gfni": 516, "xor": 1234},
                "transpose": {"gfni": 516, "xor": 1275},
                "parity_transpose": {"gfni": 516, "xor": 1188}}
    assert circuit_counts() == expected
    path = Path(__file__).with_name("Cauchy128OuterFft.h")
    if path.exists():
        assert path.read_text() == header, "header is incomplete or stale"
    print(f"PASS {checks} exact small/full/parity/dual/basis/random/adjoint checks; {circuit_counts()}")
    print(cost_model())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--emit", action="store_true")
    parser.add_argument("--patch", type=int)
    args = parser.parse_args()
    header, parts = generate()
    if args.patch is not None:
        i = args.patch
        target = Path(__file__).with_name("Cauchy128OuterFft.h").as_posix()
        print("*** Begin Patch")
        print(f"*** {'Add' if i == 0 else 'Update'} File: {target}")
        if i:
            print("@@")
            print(f" // c128-fft-generated-part {i-1}")
        print("\n".join("+"+line for line in parts[i]))
        print("*** End Patch")
    elif args.emit:
        print(header, end="")
    elif args.check:
        check(header)
    else:
        print(f"parts={len(parts)} lines={len(header.splitlines())} {circuit_counts()}")
