"""Exact scalar and source-preservation controls; not native execution/timing."""

from pathlib import Path
import random
import re
import unittest


WORKSTREAM = Path(__file__).resolve().parents[1]
REPO = WORKSTREAM.parents[2]
NATIVE = REPO / "spin/experiments/char2_field_spin"


def multiply(a, b):
    out = 0
    for _ in range(8):
        if b & 1:
            out ^= a
        a = ((a << 1) ^ (0x11B if a & 128 else 0)) & 255
        b >>= 1
    return out


def old_keys(t):
    return (
        t[3], t[3] ^ t[2], t[3] ^ t[4],
        t[5] ^ t[3], t[5] ^ t[3] ^ t[4] ^ t[2],
        t[5] ^ t[3] ^ t[6] ^ t[4],
        t[1] ^ t[3], t[1] ^ t[3] ^ t[0] ^ t[2],
        t[1] ^ t[3] ^ t[2] ^ t[4],
    )


def compact_keys(t):
    keys = old_keys(t)
    return tuple(keys[i] for i in (0, 1, 2, 3, 5, 6, 7))


def reconstruct(keys):
    a0, a1, a2, v0, v2, w0, w1 = keys
    k = a1 ^ a2
    return a0, a1, a2, v0, v0 ^ k, v2, w0, w1, w0 ^ k


def dense(t, x):
    out = []
    for i in range(4):
        y = 0
        for j in range(4):
            y ^= multiply(t[3 + i - j], x[j])
        out.append(y)
    return tuple(out)


def compact(t, x):
    a0, a1, a2, v0, v2, w0, w1 = compact_keys(t)
    k = a1 ^ a2
    low, high = x[0] ^ x[2], x[1] ^ x[3]
    shared_u = multiply(low ^ high, a0)
    u0 = shared_u ^ multiply(high, a1)
    u1 = shared_u ^ multiply(low, a2)
    shared_v = multiply(x[0] ^ x[1], v0)
    out_v0 = shared_v ^ multiply(x[1], v0 ^ k)
    out_v1 = shared_v ^ multiply(x[0], v2)
    shared_w = multiply(x[2] ^ x[3], w0)
    out_w0 = shared_w ^ multiply(x[3], w1)
    out_w1 = shared_w ^ multiply(x[2], w0 ^ k)
    return u0 ^ out_w0, u1 ^ out_w1, u0 ^ out_v0, u1 ^ out_v1


def project_compact_off(source):
    """Remove only the added Compact7 preprocessor branches, not a C parser.

    Existing preprocessor conditions remain literal. This checks that every
    old branch survives unchanged, including setup, checks and dispatch.
    """
    lines = source.splitlines(keepends=True)
    flag = "HDGROUP_TOEPLITZ_COMPACT7"
    opening = re.compile(r"#(?:if|ifdef|ifndef)\b")
    boundary = re.compile(r"#(?:elif|else|endif)\b")

    def body(index):
        out = []
        while index < len(lines):
            text = lines[index].strip()
            if boundary.match(text):
                break
            if opening.match(text):
                rendered, index = conditional(index)
                out.append(rendered)
            else:
                out.append(lines[index])
                index += 1
        return "".join(out), index

    def conditional(index):
        branches = []
        while True:
            header = lines[index]
            content, index = body(index + 1)
            branches.append((header, content))
            if lines[index].strip() == "#endif":
                end = lines[index]
                index += 1
                break
        first = branches[0][0].strip()
        if first == f"#ifndef {flag}":
            assert len(branches) == 1
            assert branches[0][1].strip() == f"#define {flag} 0"
            return "", index
        if first == f"#if {flag}":
            assert len(branches) <= 2
            if len(branches) == 2:
                assert branches[1][0].strip() == "#else"
                return branches[1][1], index
            return "", index
        kept = [(h, b) for h, b in branches if h.strip() != f"#elif {flag}"]
        assert all(flag not in h for h, _ in kept)
        return "".join(h + b for h, b in kept) + end, index

    result, stop = body(0)
    assert stop == len(lines)
    return result


class Compact7Tests(unittest.TestCase):
    def test_reconstructed_keys_all_coefficient_bits(self):
        for bit in range(56):
            t = [0] * 7
            t[bit // 8] = 1 << (bit % 8)
            self.assertEqual(reconstruct(compact_keys(t)), old_keys(t))

    def test_bilinear_basis_matches_dense(self):
        # Both circuits are GF2-bilinear. These 56*32 comparisons prove the
        # scalar identity for every seven-byte Toeplitz key and input column.
        for coefficient in range(56):
            t = [0] * 7
            t[coefficient // 8] = 1 << (coefficient % 8)
            for direction in range(32):
                x = [0] * 4
                x[direction // 8] = 1 << (direction % 8)
                self.assertEqual(compact(t, x), dense(t, x))

    def test_random_dense_and_singular_matrices(self):
        rng = random.Random(9128)
        for _ in range(128):
            t = [rng.randrange(256) for _ in range(7)]
            x = [rng.randrange(256) for _ in range(4)]
            self.assertEqual(reconstruct(compact_keys(t)), old_keys(t))
            self.assertEqual(compact(t, x), dense(t, x))
        self.assertEqual(compact([0] * 7, [1, 2, 3, 4]), (0, 0, 0, 0))

    def test_expansion_lane_order_and_xor(self):
        # vpshufb indexes each 128-bit lane of the broadcast dword.
        selectors = (0,) * 16 + (1,) * 16 + (2,) * 16 + (3,) * 16

        def expand(keys):
            packed = bytes(keys) * 16
            return bytes(packed[16 * (i // 16) + s] for i, s in enumerate(selectors))

        for bit in range(32):
            keys = [0] * 4
            keys[bit // 8] = 1 << (bit % 8)
            self.assertEqual(expand(keys), b"".join(bytes([x]) * 16 for x in keys))
        rng = random.Random(213)
        for _ in range(64):
            a, b = ([rng.randrange(256) for _ in range(4)] for _ in range(2))
            self.assertEqual(expand([x ^ y for x, y in zip(a, b)]),
                             bytes(x ^ y for x, y in zip(expand(a), expand(b))))

    def test_packed_stride_covers_exact_table(self):
        rows = 12
        occupied = set()
        for group in range(0, rows, 4):
            for column in range(64):
                for key in range(7):
                    for lane in range(4):
                        index = group * 448 + 4 * (column * 7 + key) + lane
                        self.assertNotIn(index, occupied)
                        occupied.add(index)
                        self.assertEqual(index, group * 448 + 28 * column + 4 * key + lane)
        self.assertEqual(occupied, set(range(rows * 448)))
        self.assertEqual((576 - 448) * 8192, 1 << 20)

    def test_native_static_operation_counts(self):
        header = (NATIVE / "HdGroupedOuter.h").read_text()
        kernel = header.split("HDGROUP_INLINE void compactToeplitz4(", 1)[1].split("\n}\n", 1)[0]
        self.assertEqual(kernel.count("compactToeplitzKey("), 7)
        self.assertEqual(kernel.count("_mm512_gf2p8mul_epi8("), 9)
        self.assertEqual(kernel.count("_mm512_xor_si512("), 18)
        self.assertNotIn("for(", kernel)
        self.assertNotIn("std::function", kernel)

    def test_flag_off_preserves_old_source_branches(self):
        for name in ("HdGroupedOuter.h", "hd_grouped_forward.cpp"):
            old = (WORKSTREAM / "compact7_baseline" / name).read_text()
            new = project_compact_off((NATIVE / name).read_text())
            # Added conditional blocks can leave blank separators only.
            normalize = lambda x: "\n".join(line for line in x.splitlines() if line.strip())
            self.assertEqual(normalize(new), normalize(old), name)


if __name__ == "__main__":
    unittest.main()
