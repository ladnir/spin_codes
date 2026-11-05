"""Portable arithmetic/embedding checks; no encoding or performance claims."""

import random
import unittest
from pathlib import Path

import field128_oracle as f


def bytewise_multiply(value, coefficient):
    return sum(f.aes_multiply(value >> (8 * byte) & 255, coefficient) << (8 * byte)
               for byte in range(16))


class Field128LiftTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root, cls.generator, cls.exponent = f.derive_aes_root()
        cls.embedded = f.embedding(cls.root)
        cls.columns = f.byte_basis_columns(cls.embedded)
        cls.inverse = f.inverse_columns(cls.columns)

    def test_native_polynomial_irreducible(self):
        value, middle = 2, None
        for step in range(1, 129):
            value = f.multiply_oracle(value, value)
            if step == 64:
                middle = value
        self.assertEqual(value, 2)
        # Degree 128 has only the prime divisor 2: this is Rabin's test.
        a, b = f.MODULUS, middle ^ 2
        while b:
            a, b = b, f.polynomial_remainder(a, b)
        self.assertEqual(a, 1)

    def test_deterministic_embedding_pin(self):
        self.assertEqual(self.generator, 0xE6114072B8CA57AFD9DB18ED46787786)
        self.assertEqual(self.exponent, 95)
        self.assertEqual(self.root, 0x3D5BD35C94646A247573DA4A5F7710ED)
        self.assertEqual(f.aes_polynomial(self.root), 0)
        self.assertEqual(f.power(self.generator, 255), 1)
        self.assertEqual(len(set(self.embedded)), 256)
        header = Path(__file__).with_name("Field128.h").read_text()
        self.assertIn("aesRootPower=95", header)
        self.assertIn("0x7573da4a5f7710edULL,0x3d5bd35c94646a24ULL", header)

    def test_every_embedded_byte_product_and_sum(self):
        for a in range(256):
            for b in range(256):
                self.assertEqual(self.embedded[a ^ b], self.embedded[a] ^ self.embedded[b])
                self.assertEqual(f.multiply_oracle(self.embedded[a], self.embedded[b]),
                                 self.embedded[f.aes_multiply(a, b)])

    def test_independent_native_arithmetic(self):
        for left in range(128):
            for right in range(128):
                self.assertEqual(f.multiply(1 << left, 1 << right),
                                 f.multiply_oracle(1 << left, 1 << right))
        rng = random.Random(0x128113)
        for _ in range(128):
            a, b, c = (rng.getrandbits(128) for _ in range(3))
            self.assertEqual(f.multiply(a, b), f.multiply_oracle(a, b))
            self.assertEqual(f.multiply_oracle(a, b ^ c),
                             f.multiply_oracle(a, b) ^ f.multiply_oracle(a, c))
        for value in (0, 1, f.MASK, 1 << 127):
            self.assertEqual(f.multiply(value, 0), 0)
            self.assertEqual(f.multiply(value, 1), value)

    def test_full_binary_basis_roundtrips(self):
        self.assertEqual(len(self.inverse), 128)
        for bit in range(128):
            value = 1 << bit
            self.assertEqual(f.transform(f.transform(value, self.columns), self.inverse), value)
            self.assertEqual(f.transform(f.transform(value, self.inverse), self.columns), value)
        rng = random.Random(0xBA515)
        for _ in range(64):
            value = rng.getrandbits(128)
            self.assertEqual(f.transform(f.transform(value, self.columns), self.inverse), value)

    def test_bytewise_gfni_compatibility_after_conversion(self):
        # Bilinearity reduces arbitrary byte coefficients and records to
        # 8*128 basis pairs. All 256 coefficients also get a dense record.
        for coefficient_bit in range(8):
            coefficient = 1 << coefficient_bit
            for bit in range(128):
                coordinates = 1 << bit
                expected = f.transform(bytewise_multiply(coordinates, coefficient), self.columns)
                native = f.transform(coordinates, self.columns)
                self.assertEqual(f.multiply_oracle(native, self.embedded[coefficient]), expected)
        rng = random.Random(0xC0FFEE)
        for coefficient in range(256):
            coordinates = rng.getrandbits(128)
            native = f.transform(coordinates, self.columns)
            expected = bytewise_multiply(coordinates, coefficient)
            actual = f.transform(f.multiply_oracle(native, self.embedded[coefficient]), self.inverse)
            self.assertEqual(actual, expected)

    def test_unembedded_gfni_is_an_explicit_negative_control(self):
        # GFNI byte scaling by 2 on raw native records is not multiplication
        # by the embedded AES 2 and is not even GF128-linear.
        self.assertNotEqual(bytewise_multiply(1, 2), self.embedded[2])
        scalar, value = 0x80, 1
        left = bytewise_multiply(f.multiply_oracle(scalar, value), 2)
        right = f.multiply_oracle(scalar, bytewise_multiply(value, 2))
        self.assertEqual((left, right), (0x1B, 0x100))
        self.assertNotEqual(left, right)
        self.assertNotEqual(f.aes_polynomial(2), 0)


if __name__ == "__main__":
    unittest.main()
