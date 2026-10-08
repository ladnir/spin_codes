import random
import unittest

import subfield_basis as b


class SubfieldBasisTests(unittest.TestCase):
    def test_declared_moduli_are_irreducible(self):
        for bits, modulus in b.MODULI.items():
            self.assertEqual(len(modulus) - 1, bits // 8)
            self.assertTrue(b.irreducible(modulus), bits)

    def test_known_reducible_polynomials(self):
        # Every GF2 degree-four irreducible splits over GF256.
        self.assertFalse(b.irreducible((1, 1, 0, 0, 1)))
        for degree in (4, 8, 16):
            self.assertFalse(b.irreducible((1,) + (0,) * (degree - 1) + (1,)))

    def test_gf256_table_and_inverses(self):
        self.assertEqual(b.gf_mul(0x57, 0x13), 0xfe)
        for a in range(1, 256):
            self.assertEqual(b.gf_mul(a, b.gf_pow(a, 254)), 1)

    def test_small_coefficient_is_bytewise(self):
        rng = random.Random(20261003)
        for bits, modulus in b.MODULI.items():
            degree = bits // 8
            for _ in range(24):
                x = tuple(rng.randrange(256) for _ in range(degree))
                for coefficient in (0, 1, 2, 0x63, 0xa5, 255):
                    actual = b.multiply(x, (coefficient,), modulus)
                    expected = b.trim(b.gf_mul(value, coefficient) for value in x)
                    self.assertEqual(actual, expected)

    def test_general_scalar_commutes_with_subfield_matrix(self):
        rng = random.Random(20261004)
        for bits, modulus in b.MODULI.items():
            degree = bits // 8
            for _ in range(8):
                x, y, scalar = [tuple(rng.randrange(256) for _ in range(degree)) for _ in range(3)]
                c, d = rng.randrange(256), rng.randrange(256)
                def linear(a, value):
                    return b.add(b.multiply(a, (c,), modulus), b.multiply(value, (d,), modulus))
                self.assertEqual(b.multiply(linear(x, y), scalar, modulus),
                                 linear(b.multiply(x, scalar, modulus),
                                        b.multiply(y, scalar, modulus)))

    def test_frobenius_field_order_on_random_elements(self):
        rng = random.Random(20261005)
        for bits, modulus in b.MODULI.items():
            degree = bits // 8
            for _ in range(4):
                a = b.trim(rng.randrange(256) for _ in range(degree))
                value = a
                for _ in range(degree):
                    value = b.frobenius(value, modulus)
                self.assertEqual(value, a)


if __name__ == "__main__":
    unittest.main()
