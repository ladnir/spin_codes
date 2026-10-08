"""Consistency tests for isolated floating puncture-frontier screens."""
import math
import unittest

import numpy as np

from gf8_puncture_frontier import dense_point, middle_bound, period_matrix, profile, sparse_bound
from gf8_punctured_screen import middle_bound as old_middle, sparse_bound as old_sparse


class FrontierTests(unittest.TestCase):
    def test_homogeneous_control(self):
        for outputs, dimension in ((6, 80), (7, 84)):
            g = profile(dimension, (outputs,))
            for density in (.03, .3, .98):
                point = (-.3, -.8, 2., .4)
                actual = middle_bound(point, density, dimension, (outputs,))
                control = old_middle(point, round(density*g['input_n']), 8, outputs, 8, 128, dimension)
                self.assertAlmostEqual(actual, control, places=7)
            for active in (1, 8, 32):
                point = (6.5, 2.)
                actual = sparse_bound(point, active, dimension, (outputs,))
                control = old_sparse(active, point, 8, outputs, 8, 128, dimension)
                self.assertAlmostEqual(actual, control, places=7)

    def test_repeated_period(self):
        point = (.7, -1.1, 2.7, -.2)
        for value in (.03, .4, .98):
            self.assertAlmostEqual(middle_bound(point, value, 80, (6,)),
                                   middle_bound(point, value, 80, (6, 6, 6, 6)), places=7)
        for active in (1, 8):
            self.assertAlmostEqual(sparse_bound((6.5, 2.), active, 80, (6,)),
                                   sparse_bound((6.5, 2.), active, 80, (6, 6)), places=6)

    def test_ordered_product(self):
        point = (.3, -.4, 1.2, .7)
        left, ls = period_matrix(point, (6,))
        right, rs = period_matrix(point, (7,))
        mixed, ms = period_matrix(point, (6, 7))
        np.testing.assert_allclose(mixed*math.exp(ms), left@right*math.exp(ls+rs), rtol=1e-14)

    def test_dense_endpoint_identity(self):
        for dimension, schedule in ((80, (6,)), (81, (6,)), (84, (6, 6, 6, 7))):
            g = profile(dimension, schedule)
            p = g['output_n']/g['input_n']
            expected = (p*(8*(1-g['nominal_rate'])-1.27)-math.log2(256/255))
            actual = dense_point(1., dimension, schedule, asymptotic=True)['margin_bits']
            self.assertAlmostEqual(expected, actual, places=12)

    def test_six_output_geometry_and_join(self):
        g = profile(80, (6,))
        self.assertEqual((g['rows'], g['input_n'], g['output_n'], g['threshold']),
                         (13120, 1679360, 1259520, 125951))
        # Sparse A<=1028, followed by middle W>=ceil(.03 Nin), overlaps.
        self.assertGreaterEqual(49*1029, math.ceil(.03*g['input_n']))


if __name__ == '__main__':
    unittest.main()
