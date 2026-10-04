"""Authenticate the generated forward RS circuit against the frozen transpose.

This evaluates every column of the 64-by-64 binary byte-plane map. It does not
reuse the generator's regular expressions, transpose helper, or circuit lists.
Run with a repository checkout; the installed library does not depend on this.
"""
from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "spin/experiments/k16_codesign_100us/kernel/NativeOuterFinish.h"
FORWARD = ROOT / "spin/src/paired15/Paired15ForwardOuter.cpp"


def arguments(text):
    parts, start, depth = [], 0, 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch == "," and depth == 0:
            parts.append(text[start:i])
            start = i + 1
    return parts + [text[start:]]


def expression(text, values):
    text = text.strip()
    if text in values:
        return values[text]
    if re.fullmatch(r"(?:0x[0-9a-f]+|\d+)(?:ULL)?", text):
        return int(text.removesuffix("ULL"), 0)
    name, args = text.split("(", 1)
    args = [expression(x, values) for x in arguments(args[:-1])]
    if name == "_mm512_set1_epi64":
        return args[0]
    if name == "_mm512_xor_si512":
        return args[0] ^ args[1]
    if name == "_mm512_gf2p8affine_epi64_epi8":
        byte, matrix, constant = args
        return constant ^ sum(
            (((matrix >> (8 * (7 - i))) & byte).bit_count() & 1) << i
            for i in range(8)
        )
    raise AssertionError(f"unrecognized operation: {name}")


def execute(statements, values):
    values = values.copy()
    for statement in statements.split(";"):
        statement = statement.strip()
        if not statement:
            continue
        statement = statement.removeprefix("const auto ").removeprefix("auto ")
        left, right = statement.split("=", 1)
        if "packed[" in right or "_mm512_load_si512" in right:
            continue  # Inputs were supplied directly, without their SIMD loads.
        if left.startswith("packed["):
            continue  # Return the register tuple before stores.
        values[left] = expression(right, values)
    return values


class ForwardRsAdjoint(unittest.TestCase):
    def test_all_binary_columns(self):
        source = FROZEN.read_text()
        constants = {
            f"m{index}": int(value, 16)
            for index, value in re.findall(
                r"const auto m(\d+)=_mm512_set1_epi64\(0x([0-9a-f]+)ULL\);", source
            )
        }
        reverse = source[source.index("auto x0="):source.index("return {{")]
        source = FORWARD.read_text()
        forward = source[source.index("auto y0=packed["):source.index("\n}\nstatic SPIN_FORCEINLINE void parityPlanes")]
        reverse_columns, forward_columns = [], []
        for bit in range(64):
            rv = execute(reverse, constants | {f"x{i}": (1 << bit) >> (8 * i) & 255 for i in range(8)})
            fv = execute(forward, {f"y{i}": (1 << bit) >> (8 * i) & 255 for i in range(8)})
            reverse_columns.append(sum(rv[f"y{i}"] << (8 * i) for i in range(8)))
            forward_columns.append(sum(fv[f"x{i}"] << (8 * i) for i in range(8)))
        for bit, actual in enumerate(forward_columns):
            expected = sum(((column >> bit) & 1) << j for j, column in enumerate(reverse_columns))
            self.assertEqual(actual, expected, f"forward parity input bit {bit}")


if __name__ == "__main__":
    unittest.main()
