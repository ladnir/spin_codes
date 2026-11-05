"""Tiny multiply-by-Y oracle/models; native scalar/SIMD checks are separate."""

from random import Random
import unittest

from field128_oracle import aes_multiply
import sparse_byte_oracle as oracle


def xtime(value):
    return ((value << 1) ^ (0x1b if value >> 7 else 0)) & 255


def multiply_byte_by32(value):
    x2 = xtime(value)
    x4 = xtime(x2)
    x8 = xtime(x4)
    x16 = xtime(x8)
    x32 = xtime(x16)
    return x2 ^ x16 ^ x32


def multiply_y(value):
    top = value >> 120
    return ((value << 8) & oracle.MASK) ^ multiply_byte_by32(top) ^ (top << 24) ^ (top << 40)


def simd_lane_model(value):
    """Literal lane shifts, byte-shuffle mask, and GFNI-constant model."""
    raw = value.to_bytes(16, 'little')
    mask = (0x80800f800f808080).to_bytes(8, 'little') + bytes([0x80]*8)
    taps = bytes(0 if index & 128 else raw[index & 15] for index in mask)
    shifted = bytes([0]) + raw[:15]
    constant = bytes([aes_multiply(raw[15], 0x32)]) + bytes(15)
    return int.from_bytes(bytes(a ^ b ^ c for a, b, c in zip(shifted, constant, taps)), 'little')


class MultiplyYTests(unittest.TestCase):
    def test_all_constants_and_binary_basis(self):
        for byte in range(256):
            self.assertEqual(multiply_byte_by32(byte), aes_multiply(byte, 0x32))
        for bit in range(128):
            value = 1 << bit
            expected = oracle.multiply(value, 0x100)
            self.assertEqual(multiply_y(value), expected)
            self.assertEqual(simd_lane_model(value), expected)
        self.assertEqual(multiply_y(1 << 120), 0x10001000032)

    def test_random_and_embedded_byte_linearity(self):
        rng = Random(0x12848)
        for unused in range(256):
            a, b, byte = rng.getrandbits(128), rng.getrandbits(128), rng.randrange(256)
            expected = oracle.multiply(a, 0x100)
            self.assertEqual(multiply_y(a), expected)
            self.assertEqual(simd_lane_model(a), expected)
            self.assertEqual(multiply_y(a ^ b), multiply_y(a) ^ multiply_y(b))
            self.assertEqual(multiply_y(oracle.bytewise(a, byte)), oracle.bytewise(multiply_y(a), byte))
        self.assertEqual(multiply_y(0), 0)
        self.assertEqual(multiply_y(1), 0x100)

    def test_fullfield_linearity_and_degree16_precondition(self):
        self.assertTrue(oracle.irreducible())
        rng = Random(0x32)
        for unused in range(16):
            a, b = rng.getrandbits(128), rng.getrandbits(128)
            self.assertEqual(multiply_y(oracle.multiply(a, b)), oracle.multiply(multiply_y(a), b))


if __name__ == '__main__':
    unittest.main()
