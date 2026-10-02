"""Independent algebra and bit-layout review of the native field16 outer.

This models the declared intrinsics and setup formulas. It does not execute
the compiled SIMD kernel; the C++ scalar/native checks cover that boundary.
"""
import random
import unittest


def mul8(a, b):
    value = 0
    for bit in range(8):
        if (b >> bit) & 1:
            value ^= a << bit
    for bit in range(14, 7, -1):
        if (value >> bit) & 1:
            value ^= 0x11b << (bit-8)
    return value


def mul16(c, x):
    c0, c1, x0, x1 = c & 255, c >> 8, x & 255, x >> 8
    lo = mul8(c0, x0) ^ mul8(0x20, mul8(c1, x1))
    hi = mul8(c0, x1) ^ mul8(c1, x0) ^ mul8(c1, x1)
    return lo | (hi << 8)


def circuit(c, x):
    c0, c1, x0, x1 = c & 255, c >> 8, x & 255, x >> 8
    p0 = mul8(c0, x0)
    p1 = mul8(mul8(0x20, c1), x1)
    p2 = mul8(c0 ^ c1, x0 ^ x1)
    return (p0 ^ p1) | ((p0 ^ p2) << 8)


def transpose_coordinate(c):
    return 4*(c % 4)+c//4


def permute_word(value):
    return sum(((value >> c) & 1) << transpose_coordinate(c) for c in range(16))


def apply(rows, value):
    return sum(((row & value).bit_count() & 1) << i for i, row in enumerate(rows))


def multiplication_rows(c):
    columns = [mul16(c, 1 << bit) for bit in range(16)]
    return [sum(((column >> row) & 1) << bit for bit, column in enumerate(columns))
            for row in range(16)]


def pack_native(raw):
    """Literal packIndex and GFNI semantics for the native input packing."""
    output = []
    for pair in (0, 1):
        source = raw[128*pair:128*(pair+1)]
        for half in (0, 8):
            before = [source[16*(7-coordinate)+half+payload_byte]
                      for payload_byte in range(8) for coordinate in range(8)]
            # GFNI's matrix argument is before; matrix byte7 supplies outputbit0.
            after = []
            for base in range(0, 64, 8):
                for byte in range(8):
                    operand = 1 << (7-byte)
                    after.append(sum(((before[base+7-bit] & operand).bit_count() & 1) << bit
                                     for bit in range(8)))
            output.append(after)
    return output


class NativeFieldTests(unittest.TestCase):
    def test_tower_irreducibility_trace(self):
        trace, power = 0, 0x20
        for _ in range(8):
            trace ^= power
            power = mul8(power, power)
        self.assertEqual(trace, 1)

    def test_full_native_pack_bit_basis(self):
        for bit in range(2048):
            raw = [0]*256
            raw[bit//8] = 1 << (bit % 8)
            packed = pack_native(raw)
            coordinate, payload = divmod(bit, 128)
            expected = [[0]*64 for _ in range(4)]
            expected[2*(coordinate//8)+payload//64][8*((payload//8) % 8)+7-payload % 8] = 1 << (coordinate % 8)
            self.assertEqual(packed, expected)

    def test_compacted_gl_column_permutation_on_all_matrix_entries(self):
        # Every coefficient of an arbitrary16x16 binary reverse symbol map.
        for row in range(16):
            for column in range(16):
                old = [0]*16
                old[row] = 1 << column
                native = [permute_word(value) for value in old]
                for physical in range(16):
                    self.assertEqual(apply(native, 1 << physical),
                                     apply(old, permute_word(1 << physical)))

    def test_scalar_native_and_adjoint_basis_maps(self):
        rng = random.Random(1729)
        scalars = [1, 2, 0x100, 0x20, 0xffff, *[rng.randrange(1, 65536) for _ in range(32)]]
        for scalar in scalars:
            rows = multiplication_rows(scalar)
            reverse = [permute_word(row) for row in rows]  # M_c Q in scalarcoordinates.
            native = [permute_word(row) for row in reverse]  # M_c after nativecompaction.
            self.assertEqual(native, rows)
            for physical in range(16):
                x = 1 << physical
                self.assertEqual(apply(reverse, permute_word(x)), circuit(scalar, x))
                self.assertEqual(circuit(scalar, x), mul16(scalar, x))
            for x in (1, 7, 0x8001, 0xffff):
                adjoint = sum(((mul16(scalar, 1 << bit) & x).bit_count() & 1) << bit
                              for bit in range(16))
                forward = permute_word(adjoint)
                for y in (1, 0x61a2, 0xffff):
                    self.assertEqual((forward & y).bit_count() & 1,
                                     (x & apply(reverse, y)).bit_count() & 1)

    def test_nonzero_scalar_orbits_for_fixed_inputs(self):
        # The proof covers every nonzero input. These exhaustive scalarorbits
        # independently check four inputs, including mixed high/low fieldbytes.
        for x in (1, 2, 0x8001, 0xffff):
            basis = []
            for bit in range(16):
                c = 1 << bit
                adjoint = sum(((mul16(c, 1 << i) & x).bit_count() & 1) << i
                              for i in range(16))
                basis.append(permute_word(adjoint))
            orbit = [0]
            for value in basis:
                orbit.extend(previous ^ value for previous in orbit[:])
            self.assertEqual(len(set(orbit)), 65536)
            self.assertNotIn(0, orbit[1:])


if __name__ == '__main__':
    unittest.main()
