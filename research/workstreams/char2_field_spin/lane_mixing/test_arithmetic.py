"""Exact checks for field definitions, byte scalar sets, and shear atoms."""

from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
import random
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import arithmetic as a


class ArithmeticTests(unittest.TestCase):
    def test_moduli_irreducible(self):
        for m, modulus in a.MODULI.items():
            self.assertEqual(modulus.bit_length() - 1, m)
            self.assertTrue(a.irreducible(modulus), m)
        self.assertFalse(a.irreducible((1 << 8) | 1))
        self.assertFalse(a.irreducible(a.carryless_product(0x13, 0x19)))

    def test_small_coefficients_are_distinct_nonzero(self):
        for m, modulus in a.MODULI.items():
            self.assertEqual(len({a.polynomial_remainder(c, modulus)
                                  for c in range(1, 256)}), 255)
            for value in (1, 2, (1 << m) - 1, 1 << (m - 1)):
                products = {a.multiply(value, c, modulus) for c in range(1, 256)}
                self.assertEqual(len(products), 255)
                self.assertNotIn(0, products)

    def test_byte_set_is_not_subfield_above_degree_eight(self):
        for m, modulus in a.MODULI.items():
            if m == 8:
                for x in range(256):
                    self.assertEqual(a.frobenius(x, 8, modulus), x)
            else:
                self.assertNotEqual(a.frobenius(2, 8, modulus), 2)
                # Both 0x80 and 2 are byte polynomials; their product is not.
                self.assertEqual(a.multiply(0x80, 2, modulus), 0x100)

    def test_horner_matches_literal_product(self):
        rng = random.Random(0xC2F1E1D)
        for m, modulus in a.MODULI.items():
            values = [0, 1, 2, 1 << (m - 1), (1 << m) - 1]
            values += [rng.getrandbits(m) for _ in range(16)]
            for value, coefficient in product(values, range(256)):
                self.assertEqual(a.multiply_byte_horner(value, coefficient, modulus),
                                 a.multiply(value, coefficient, modulus))

    def test_onefold_matches_literal_product(self):
        rng = random.Random(0xC1A0)
        for m, modulus in a.MODULI.items():
            if m == 8:
                with self.assertRaises(ValueError):
                    a.multiply_byte_onefold(128, 255, modulus)
                continue
            values = [0, 1, 2, 1 << (m - 1), (1 << m) - 1]
            values += [rng.getrandbits(m) for _ in range(16)]
            for value, coefficient in product(values, range(256)):
                self.assertEqual(a.multiply_byte_onefold(value, coefficient, modulus),
                                 a.multiply(value, coefficient, modulus))

    def test_shears_preserve_nonzero_and_are_field_linear(self):
        modulus = 0x13
        for coeffs in ((1, 1, 1), (2, 3, 7), (15, 9, 4)):
            matrix = a.shear_matrix(*coeffs, modulus)
            determinant = a.multiply(matrix[0], matrix[3], modulus) ^ a.multiply(matrix[1], matrix[2], modulus)
            self.assertEqual(determinant, 1)
            images = {a.three_shears(v, *coeffs, modulus) for v in product(range(16), repeat=2)}
            self.assertEqual(len(images), 256)
            for scalar in range(16):
                state = (5, 11)
                scaled = tuple(a.multiply(scalar, x, modulus) for x in state)
                self.assertEqual(a.three_shears(scaled, *coeffs, modulus),
                                 tuple(a.multiply(scalar, x, modulus)
                                       for x in a.three_shears(state, *coeffs, modulus)))

    def test_three_shear_vector_atoms(self):
        # Scaling a nonzero state and target by a common nonzero field scalar
        # preserves each coefficient fiber. These are all projective states.
        for modulus, coefficients in ((0x7, (1, 2, 3)),
                                      (0x13, (1, 2, 4, 7)),
                                      (0x13, tuple(range(1, 16)))):
            field_size = 1 << (modulus.bit_length() - 1)
            q = len(coefficients)
            for state in [(1, 0)] + [(x, 1) for x in range(field_size)]:
                counts = Counter(a.three_shears(state, *abc, modulus)
                                 for abc in product(coefficients, repeat=3))
                self.assertNotIn((0, 0), counts)
                self.assertLessEqual(max(counts.values()), q)
                for (x, y), count in counts.items():
                    if y == 0:
                        self.assertLessEqual(count, q)

    def test_rank_two_matrix_atoms(self):
        for modulus, coefficients in ((0x7, (1, 2, 3)),
                                      (0x13, tuple(range(1, 16)))):
            matrices = [a.shear_matrix(*abc, modulus)
                        for abc in product(coefficients, repeat=3)]
            self.assertEqual(len(set(matrices)), len(coefficients) ** 3)

    def test_two_shears_do_not_have_rank_two_vector_atoms(self):
        modulus, coefficients = 0x13, tuple(range(1, 16))
        q = len(coefficients)
        counts = Counter(a.two_shears((1, 1), x, y, modulus)
                         for x, y in product(coefficients, repeat=2))
        self.assertEqual(Fraction(counts[(0, 1)], q*q), Fraction(1, q))
        self.assertGreater(Fraction(counts[(0, 1)], q*q), Fraction(1, q*q))

    def test_independent_pivot_bound_is_not_marginal_multiplication(self):
        modulus, coefficients = 0x13, (1, 2, 4, 7)
        q = len(coefficients)
        # Three independent coefficients times rank-two forms e1,e2,e1+e2.
        independent = Counter((x ^ z, y ^ z)
                              for x, y, z in product(coefficients, repeat=3))
        self.assertLessEqual(max(independent.values()), q)
        # A shared coefficient gives both coordinates uniform marginals but
        # only one constraint: the joint target atom is 1/q, not 1/q**2.
        dependent = Counter((x, x) for x in coefficients)
        self.assertEqual(Fraction(dependent[(1, 1)], q), Fraction(1, q))
        self.assertGreater(Fraction(dependent[(1, 1)], q), Fraction(1, q*q))

    def test_binary_linear_map_need_not_preserve_field_correlation(self):
        modulus = a.MODULI[8]
        # Frobenius is an invertible GF2-linear map but not GF256-linear.
        square = lambda x: a.multiply(x, x, modulus)
        self.assertEqual(square(5 ^ 11), square(5) ^ square(11))
        self.assertNotEqual(square(a.multiply(2, 1, modulus)),
                            a.multiply(2, square(1), modulus))
        # Independent AES-field byte products are not multiplication in the
        # selected degree-32 polynomial basis.
        def bytewise_double(x):
            return sum(a.multiply((x >> (8*j)) & 255, 2, modulus) << (8*j)
                       for j in range(4))
        p32 = a.MODULI[32]
        self.assertNotEqual(bytewise_double(a.multiply(2, 128, p32)),
                            a.multiply(2, bytewise_double(128), p32))

    def test_field_transpose_uses_field_pairing(self):
        modulus = a.MODULI[8]
        coefficient, x, y = 2, 128, 1
        self.assertEqual(a.multiply(a.multiply(coefficient, x, modulus), y, modulus),
                         a.multiply(x, a.multiply(coefficient, y, modulus), modulus))
        binary_dot = lambda left, right: (left & right).bit_count() & 1
        self.assertNotEqual(binary_dot(a.multiply(coefficient, x, modulus), y),
                            binary_dot(x, a.multiply(coefficient, y, modulus)))


if __name__ == "__main__":
    unittest.main()
