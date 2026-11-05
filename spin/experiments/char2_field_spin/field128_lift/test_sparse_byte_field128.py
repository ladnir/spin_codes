"""Independent field and arithmetic-model tests; no timings or encoder."""

import random
import unittest

import sparse_byte_oracle as f
from field128_oracle import aes_multiply


class SparseByteField128Tests(unittest.TestCase):
    def test_irreducibility(self):
        self.assertEqual(len(f.MODULUS), 17)
        self.assertTrue(f.irreducible())

    def test_defining_relation(self):
        self.assertEqual(f.multiply(1 << 120, 1 << 8),
                         0x32 ^ (1 << 24) ^ (1 << 40))

    def test_all_binary_products_match_spaced_model(self):
        for a in range(128):
            for b in range(128):
                self.assertEqual(f.multiply_spaced(1 << a, 1 << b),
                                 f.multiply(1 << a, 1 << b))

    def test_all_embedded_byte_products(self):
        for a in range(256):
            for b in range(256):
                self.assertEqual(f.multiply(a, b), aes_multiply(a, b))

    def test_fixed_byte_gfni_is_exact_field_scaling(self):
        for bit in range(128):
            for coefficient_bit in range(8):
                value, coefficient = 1 << bit, 1 << coefficient_bit
                self.assertEqual(f.multiply(value, coefficient), f.bytewise(value, coefficient))
        rng = random.Random(0x3216)
        for coefficient in range(256):
            value = rng.getrandbits(128)
            self.assertEqual(f.multiply_spaced(value, coefficient), f.bytewise(value, coefficient))

    def test_dense_products_and_full_field_linearity(self):
        rng = random.Random(0xF12832)
        for _ in range(64):
            a, b, c, scalar = (rng.getrandbits(128) for _ in range(4))
            self.assertEqual(f.multiply_spaced(a, b), f.multiply(a, b))
            self.assertEqual(f.multiply(a, b ^ c), f.multiply(a, b) ^ f.multiply(a, c))
            self.assertEqual(f.multiply(f.multiply(a, b), scalar),
                             f.multiply(a, f.multiply(b, scalar)))
            # A fixed AES coefficient commutes with arbitrary full-field scaling.
            self.assertEqual(f.bytewise(f.multiply(a, scalar), 0x53),
                             f.multiply(f.bytewise(a, 0x53), scalar))
        for value in (0, 1, f.MASK, 1 << 127):
            self.assertEqual(f.multiply_spaced(value, 0), 0)
            self.assertEqual(f.multiply_spaced(value, 1), value)


if __name__ == "__main__":
    unittest.main()
