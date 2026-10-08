"""Exact checks for the bounded non-self-orthogonal Cauchy8/6 screen."""
from fractions import Fraction as F
from itertools import combinations
import math
import unittest

import numpy as np

import gf8_state6_packet_screen as target
from gf8_packet_maps import determinant, rank


class StateSixPacketTests(unittest.TestCase):
    def test_extended_invertibility_and_return_minima(self):
        for asymmetric in (False, True):
            result = target.validate_maps(asymmetric)
            self.assertEqual(result['return_rank'], 8)
            self.assertEqual(result['state_schur_rank'], 6)
            self.assertFalse(result['ca_is_zero'])

    def test_every_cauchy_minor(self):
        for asymmetric in (False, True):
            a, c = target.maps(asymmetric)
            count = 0
            for size in range(1, 7):
                for rows in combinations(range(8), size):
                    for columns in combinations(range(6), size):
                        self.assertNotEqual(determinant([[a[i][j] for j in columns] for i in rows]), 0)
                        self.assertNotEqual(determinant([[c[j][i] for i in rows] for j in columns]), 0)
                        count += 2
            self.assertEqual(count, 6004)

    def test_packet_aggregation(self):
        x, y, tau = F(2, 5), F(3, 4), F(11)
        for minima in (target.MINIMA, target.ASYMMETRIC_MINIMA):
            direct = np.asarray([[sum(math.comb(8, w)*x**w*target.local_matrix(w, y, tau, minima)[i][j]
                                      for w in range(9)) for j in range(2)] for i in range(2)], dtype=float)
            np.testing.assert_allclose(target.packet(math.log(x), math.log(y), math.log(tau), minima),
                                       direct, rtol=2e-13, atol=1e-15)


if __name__ == '__main__':
    unittest.main()
