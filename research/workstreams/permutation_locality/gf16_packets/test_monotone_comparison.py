import unittest
from fractions import Fraction as Q
from itertools import product

from monotone_comparison import thinned_shells


class MonotoneComparisonTests(unittest.TestCase):
    def test_binary_linear_maps_and_support_couplings(self):
        # Exhaust every binary map from two two-bit packets to two output
        # bits. This includes zero maps and maps with packet cancellation.
        for columns in product(range(4), repeat=4):
            for z in (Q(1, 4), Q(3, 4)):
                original, comparison = [], []
                for mask in range(4):
                    a = b = Q(0)
                    for left, right in product(range(4), repeat=2):
                        values = [left, right]
                        actual = relaxed = Q(1)
                        output = 0
                        for j, value in enumerate(values):
                            if mask >> j & 1:
                                actual *= Q(1, 3) if value else 0
                                relaxed *= Q(1, 6) if value else Q(1, 2)
                            else:
                                actual *= not value
                                relaxed *= not value
                            for bit in range(2):
                                if value >> bit & 1:
                                    output ^= columns[2*j+bit]
                        factor = z**output.bit_count()
                        a += actual*factor
                        b += relaxed*factor
                    original.append(a)
                    comparison.append(b)
                    self.assertLessEqual(a, b)
                for smaller in range(4):
                    for larger in range(4):
                        if smaller & larger == smaller:
                            self.assertGreaterEqual(comparison[smaller], comparison[larger])
                # Actual counts [0,2,7], upper CDF [0,5,9].
                actual_moment = original[1]+original[2]+7*original[3]
                bound = Q(5, 2)*(comparison[1]+comparison[2])+4*comparison[3]
                self.assertLessEqual(actual_moment, bound)

    def test_thinning_and_zero_label_mass(self):
        for q in (4, 16):
            shells = thinned_shells([0, 5, 9], q)
            keep, drop = Q(q-2, q), Q(2, q)
            self.assertEqual(shells, [5*drop+4*drop*drop, 5*keep+8*keep*drop, 4*keep*keep])
            self.assertEqual(sum(shells), 9)
            self.assertGreater(shells[0], 0)
            self.assertIsInstance(shells[0], Q)

    def test_invalid_cdf(self):
        for cdf, q in (([1, 2], 16), ([0, 3, 2], 16), ([0, -1, 3], 16),
                       ([0, 1], 3), ([0, 1], 2), ([0, 1], True), ([0, Q(1, 2)], 16)):
            with self.assertRaises(ValueError):
                thinned_shells(cdf, q)


if __name__ == '__main__':
    unittest.main()
