"""Tiny exact placement checks for the outward one-group contribution."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

import numpy as np
import q1_outward as q1


def exact_product(a, b):
    return [[sum(x*y for x, y in zip(row, column)) for column in zip(*b)] for row in a]


def exact_supports(local, regions, epochs, z, cutoff):
    """Enumerate region subsets and the active physical step in each region."""
    size = len(local[0])
    matrices = [[[F(float(x)) for x in row] for row in matrix] for matrix in local[:2]]
    result = []
    for support in range(regions+1):
        total = F(0)
        for active_regions in combinations(range(regions), support):
            for locations in product(range(epochs), repeat=support):
                active = {r*epochs+e for r, e in zip(active_regions, locations)}
                row = [[F(int(i == 0)) for i in range(size)]]
                for step in range(regions*epochs):
                    row = exact_product(row, matrices[int(step in active)])
                total += sum(row[0])
        result.append(min(F(1), total/(comb(regions, support)*epochs**support)/z**cutoff))
    return tuple(result)


class Q1Outward(unittest.TestCase):
    def test_tiny_noncommuting_placement_enumeration(self):
        local = np.array([[[.75, .0625], [.125, .5]],
                          [[.25, .125], [.1875, .375]]])
        for regions, epochs, cutoff in ((3, 2, 0), (2, 3, 1), (1, 1, 0)):
            expected = exact_supports(local, regions, epochs, F(3, 4), cutoff)
            upper = q1.q1_support_upper(local, F(3, 4), regions=regions, epochs=epochs, cutoff=cutoff)
            for exact, bound in zip(expected, upper):
                self.assertGreaterEqual(bound, exact)
                self.assertLessEqual(bound, exact*F(100000000001, 100000000000))
                self.assertEqual(bound.denominator & (bound.denominator-1), 0)

    def test_extreme_small_moment_uses_scaled_endpoints(self):
        local = np.array([[[2.**-100]], [[2.**-100]]])
        upper = q1.q1_support_upper(local, F(1, 2), regions=8, epochs=8, cutoff=0)
        exact = F(1, 1 << 6400)
        for bound in upper:
            self.assertGreaterEqual(bound, exact)
            self.assertLess(bound, exact*F(1000000001, 1000000000))

    def test_zero_support_and_probability_cap(self):
        local = np.array([[[1.]], [[0.]]])
        self.assertEqual(q1.q1_support_upper(local, F(1, 2), regions=3, epochs=2, cutoff=4),
                         (F(1), F(0), F(0), F(0)))

    def test_exact_shell_combination_and_witness_minimum(self):
        supports = (F(1), F(1, 8), F(1, 32))
        counts = (F(0), F(7, 3), F(5, 3))
        exact = 3*sum(x*y for x, y in zip(supports, counts))
        bound = q1.combine_supports(supports, counts=counts, groups=3)
        self.assertGreaterEqual(bound, exact)
        self.assertLess(bound-exact, F(1, 1 << 70))
        self.assertEqual(q1.min_supports(supports, (F(1), F(1, 16), F(1, 16))),
                         (F(1), F(1, 16), F(1, 32)))

    def test_wider_exact_shell_mass(self):
        counts = q1.wider_shell_counts()
        self.assertEqual(len(counts), 65)
        self.assertEqual(counts[0], 0)
        self.assertEqual(sum(counts), (1 << 256)-1)

    def test_endpoint_record_preserves_integer_mantissa(self):
        value = F((1 << 79)+1593, 1 << 211)
        record = q1.endpoint_record(value)
        self.assertNotIn('mantissa_hex', record)
        self.assertEqual(F(int(record['numerator_hex'], 16))*F(2)**record['exponent'], value)

    def test_reject_invalid_inputs(self):
        local = np.ones((2, 1, 1))
        for z in (0, 1, -1):
            with self.assertRaises(ValueError):
                q1.q1_support_upper(local, F(z))
        for kwargs in ({'regions':0}, {'epochs':0}, {'cutoff':-1}):
            with self.assertRaises(ValueError):
                q1.q1_support_upper(local, F(1, 2), **kwargs)


if __name__ == '__main__':
    unittest.main()
