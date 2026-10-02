"""Coverage and exact outer-conversion checks for the narrow union wrapper."""
from fractions import Fraction as Q
import unittest

import packet_rs_whole as screen


class WholeUnionTests(unittest.TestCase):
    def test_complete_disjoint_coverage(self):
        screen.check_coverage([(1, 1), (2, 2), (3, 32), (33, 512)])
        for ranges in ([(1, 32), (32, 512)], [(1, 31), (33, 512)], [(0, 512)], [(1, 513)]):
            with self.assertRaises(ValueError):
                screen.check_coverage(ranges)

    def test_exact_beta_conversion(self):
        beta8, _ = screen.counts_for('rs8')
        beta16, _ = screen.counts_for('rs16')
        self.assertEqual(beta16 / beta8, Q(65537, 65535)**4)

    def test_dyadic_endpoints(self):
        self.assertEqual(screen.dyadic([3, -2]), Q(3, 4))
        for value in ([0, 0], [-1, 2], [True, 0], [1, 1000001], ['1', 0]):
            with self.assertRaises(ValueError):
                screen.dyadic(value)


if __name__ == '__main__':
    unittest.main()
