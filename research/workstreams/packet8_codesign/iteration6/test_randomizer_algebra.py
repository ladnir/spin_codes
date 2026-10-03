"""Independent scalar checks for the isolated native randomizer experiments.

These are algebra tests, not performance measurements or certificates.
"""
import random
import unittest


def poly_multiply(a, b, degree, polynomial):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & (1 << degree):
            a ^= polynomial
    return out


def aes(a, b):
    return poly_multiply(a, b, 8, 0x11B)


def f16(a, b):
    return poly_multiply(a, b, 4, 0x13)


def quadratic(a, b, width, mul, constant):
    mask = (1 << width) - 1
    a0, a1 = a & mask, a >> width
    b0, b1 = b & mask, b >> width
    low = mul(a0, b0) ^ mul(constant, mul(a1, b1))
    high = mul(a0, b1) ^ mul(a1, b0) ^ mul(a1, b1)
    return low | (high << width)


def wide16(a, b):
    return quadratic(a, b, 8, aes, 0x20)


def wide32(a, b):
    return quadratic(a, b, 16, wide16, 0x2000)


def tower8(a, b):
    return quadratic(a, b, 4, f16, 8)


def coefficients(scalar):
    parts = (scalar & 65535, wide16(0x2000, scalar >> 16),
             (scalar & 65535) ^ (scalar >> 16))
    return [v for p in parts for v in (p & 255, aes(0x20, p >> 8), (p & 255) ^ (p >> 8))]


def flat9(word, c, product=aes):
    x0, x1, x2, x3 = [(word >> (8*i)) & 255 for i in range(4)]
    p0, p1, p2 = product(c[0], x0), product(c[1], x1), product(c[2], x0 ^ x1)
    q0, q1, q2 = product(c[3], x2), product(c[4], x3), product(c[5], x2 ^ x3)
    r0, r1, r2 = product(c[6], x0 ^ x2), product(c[7], x1 ^ x3), product(c[8], x0 ^ x1 ^ x2 ^ x3)
    y = (p0 ^ q0 ^ p1 ^ q1, p0 ^ q0 ^ p2 ^ q2,
         p0 ^ r0 ^ p1 ^ r1, p0 ^ r0 ^ p2 ^ r2)
    return sum(value << (8*i) for i, value in enumerate(y))


def affine_matrix(scalar, mul):
    return sum(sum(((mul(scalar, 1 << c) >> r) & 1) << c for c in range(8))
               << (8*(7-r)) for r in range(8))


def gfni_emulate(matrix, byte):
    return sum((((matrix >> (8*(7-r))) & byte).bit_count() & 1) << r for r in range(8))


H = ((2, 3, 1, 1), (1, 2, 3, 1), (1, 1, 2, 3), (3, 1, 1, 2))


def structured(word, scales):
    x = [(word >> (8*i)) & 255 for i in range(4)]
    x[1:] = [tower8(scales[i], x[i+1]) for i in range(3)]
    out = []
    for row in range(4):
        value = 0
        for col in range(4):
            value ^= tower8(H[row][col], x[col])
        out.append(tower8(scales[3+row], value))
    return sum(value << (8*i) for i, value in enumerate(out))


def rows_from_images(images):
    return [sum(((image >> row) & 1) << col for col, image in enumerate(images)) for row in range(32)]


def apply_rows(rows, word):
    return sum(((row & word).bit_count() & 1) << i for i, row in enumerate(rows))


def optimized_structured(word, scales, bitwise):
    x = [(word >> (8*i)) & 255 for i in range(4)]
    for i in range(1, 4):
        x[i] = gfni_emulate(affine_matrix(scales[i-1], tower8), x[i])
    total = x[0] ^ x[1] ^ x[2] ^ x[3]
    def xtime(value):
        if not bitwise:
            return gfni_emulate(affine_matrix(2, tower8), value)
        carry = (value >> 3) & 0x11
        return ((value << 1) & 0xEE) ^ carry ^ (carry << 1)
    differences = [xtime(x[i] ^ x[i+1]) for i in range(3)]
    differences.append(differences[0] ^ differences[1] ^ differences[2])
    y = [x[i] ^ total ^ differences[i] for i in range(4)]
    return sum(gfni_emulate(affine_matrix(scales[3+i], tower8), y[i]) << (8*i) for i in range(4))


class RandomizerAlgebraTests(unittest.TestCase):
    def test_flat9_and_direct16_match_literal_field(self):
        randomizer = random.Random(0x329AFF1)
        scalars = [1, 2, 0x100, 0x10000, 0xFFFFFFFF] + [randomizer.randrange(1, 1 << 32) for _ in range(20)]
        words = [1 << i for i in range(32)] + [randomizer.getrandbits(32) for _ in range(24)]
        for scalar in scalars:
            c = coefficients(scalar)
            m = [affine_matrix(value, aes) for value in c]
            rows = rows_from_images([wide32(scalar, 1 << i) for i in range(32)])
            direct = [[sum(((rows[8*out+r] >> (8*inp)) & 255) << (8*(7-r)) for r in range(8))
                       for inp in range(4)] for out in range(4)]
            for word in words:
                expected = wide32(scalar, word)
                self.assertEqual(flat9(word, c), expected)
                self.assertEqual(flat9(word, m, gfni_emulate), expected)
                byte_values = []
                for out in range(4):
                    value = 0
                    for inp in range(4):
                        value ^= gfni_emulate(direct[out][inp], (word >> (8*inp)) & 255)
                    byte_values.append(value)
                self.assertEqual(sum(v << (8*i) for i, v in enumerate(byte_values)), expected)

    def test_tower_byte_is_field_and_affine_is_ordinary(self):
        for scalar in range(1, 256):
            self.assertEqual(len({tower8(scalar, value) for value in range(1, 256)}), 255)
            matrix = affine_matrix(scalar, tower8)
            for bit in range(8):
                self.assertEqual(gfni_emulate(matrix, 1 << bit), tower8(scalar, 1 << bit))

    def test_two_nibble_xtime(self):
        for value in range(256):
            carry = (value >> 3) & 0x11
            self.assertEqual(((value << 1) & 0xEE) ^ carry ^ (carry << 1), tower8(2, value))

    def test_structured_kernels_and_binary_adjoints(self):
        randomizer = random.Random(0x42D1A6)
        for _ in range(16):
            scales = [randomizer.randrange(1, 256) for _ in range(7)]
            columns = [structured(1 << i, scales) for i in range(32)]
            rows = rows_from_images(columns)
            for word in [1 << i for i in range(32)] + [randomizer.getrandbits(32) for _ in range(16)]:
                expected = structured(word, scales)
                self.assertEqual(apply_rows(rows, word), expected)
                self.assertEqual(optimized_structured(word, scales, False), expected)
                self.assertEqual(optimized_structured(word, scales, True), expected)
                dual = randomizer.getrandbits(32)
                transpose = apply_rows(columns, dual)
                self.assertEqual((expected & dual).bit_count() & 1, (word & transpose).bit_count() & 1)


if __name__ == "__main__":
    unittest.main()
