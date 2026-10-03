"""Exact tiny-family checks for thinning, rebasing, and the G4 macro."""
from fractions import Fraction
from itertools import product
from math import comb
import unittest

import numpy as np
import macro_outward as m


def fixture():
    family = np.zeros((3, 4, 4))
    family[0, 0, 0] = 1
    family[1, 0, 2] = .75
    family[2, 0, 0] = .125
    family[2, 0, 3] = .5
    for j in range(3):
        for source in range(1, 4):
            family[j, source, 1] = Fraction(source+j+1, 16)
            if j:
                family[j, source, 0] = Fraction(source+j+1, 256)
    return family


def exact_macro(family, group_steps=4):
    family = m.exact_array(family)
    W, size = len(family)-1, family.shape[1]
    result = m.zero_array((W*group_steps+1, size, size))
    for occupancies in product(range(W+1), repeat=group_steps):
        matrix = m.zero_array((size, size))
        for i in range(size):
            matrix[i, i] = Fraction(1)
        numerator = 1
        for j in occupancies:
            matrix = matrix@family[j]
            numerator *= comb(W, j)
        total = sum(occupancies)
        result[total] += matrix*Fraction(numerator, comb(W*group_steps, total))
    return result


class OutwardMacro(unittest.TestCase):
    def test_exact_thinning(self):
        old = m.exact_array(fixture())
        actual = m.thin_exact(old, packet_bits=2)
        self.assertTrue(np.array_equal(actual[0], old[0]))
        self.assertTrue(np.array_equal(actual[1], old[0]/4+3*old[1]/4))
        self.assertTrue(np.array_equal(actual[2], old[0]/16+6*old[1]/16+9*old[2]/16))

    def test_exact_rebase_and_chronological_products(self):
        selected = m.upper_array(m.thin_exact(fixture(), packet_bits=2))
        new, change = m.rebase_exact(selected)
        old = m.exact_array(selected)
        self.assertTrue(all(sum(row) == 1 for row in change))
        for sequence in product(range(3), repeat=3):
            left, right = new[sequence[0]].copy(), old[sequence[0]].copy()
            for j in sequence[1:]:
                left = left@new[j]
                right = right@old[j]
            self.assertTrue(np.array_equal(left@change, change@right))

    def test_normalization_and_alpha_one_exact_macro(self):
        physical, z = fixture(), Fraction(3, 4)
        upper, a, metadata = m.prepare_physical_upper(physical, z, packet_bits=2)
        self.assertEqual(a, Fraction(49, 64))
        self.assertTrue(metadata['exact_intertwining_before_final_rounding'])
        exact = exact_macro(upper)
        bound, _ = m.macro_operators(upper, 1)
        for index in np.ndindex(exact.shape):
            self.assertGreaterEqual(Fraction(float(bound[index])), exact[index])
            self.assertEqual(bound[index] == 0, exact[index] == 0)
        # Independent ordinary conditional placement at alpha1.
        family = m.exact_array(upper)
        polynomial = [m.zero_array((4, 4)) for _ in range(9)]
        for i in range(4):
            polynomial[0][i, i] = Fraction(1)
        for step in range(4):
            following = [m.zero_array((4, 4)) for _ in range(9)]
            for total in range(2*step+1):
                for j in range(3):
                    following[total+j] += (polynomial[total]@family[j])*comb(2, j)
            polynomial = following
        for total in range(9):
            self.assertTrue(np.array_equal(exact[total], polynomial[total]/comb(8, total)))

    def test_fractional_one_step_integer_inequality(self):
        potential, _, _ = m.prepare_physical_upper(fixture(), Fraction(3, 4), packet_bits=2)
        alpha = Fraction(2, 5)
        bound, _ = m.macro_operators(potential, alpha, group_steps=1, tangent_bins=16)
        for index in np.ndindex(potential.shape):
            self.assertGreaterEqual(Fraction(float(bound[index]))**5, Fraction(float(potential[index]))**2)
            self.assertEqual(bound[index] == 0, potential[index] == 0)

    def test_invalid_structure_and_missing_receipt_rejected(self):
        broken = fixture()
        broken[1, 2, 3] = 1
        with self.assertRaises(ValueError):
            m.prepare_physical_upper(broken, Fraction(3, 4), packet_bits=2)
        with self.assertRaises(FileNotFoundError):
            m.load_local(m.HERE/'definitely-no-local-receipt.json')


if __name__ == '__main__':
    unittest.main()
