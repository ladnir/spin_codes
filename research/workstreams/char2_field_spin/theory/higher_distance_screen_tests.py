import math
import unittest
from fractions import Fraction as F

from cauchy_firstmoment_screen import geometry
from superregular_return_screen import bound as reference_bound
from higher_distance_screen import bound, output_threshold, density_bound
from state8_density_screen import bound as density_reference


class HigherDistanceTests(unittest.TestCase):
    def test_integer_threshold(self):
        self.assertEqual(output_threshold(100, F(3, 10)), 29)
        self.assertEqual(output_threshold(101, F(3, 10)), 30)
        self.assertEqual(output_threshold(1, F(1, 10)), 0)

    def test_reproduces_reference(self):
        p = (4., math.log(32))
        self.assertEqual(bound(2, p, 16, 128, 64, F(1, 10)),
                         reference_bound(2, p, 16, 128, 64, 8, True, 16))

    def test_threshold_change_has_correct_sign(self):
        p = (4., math.log(32))
        low = bound(2, p, 16, 128, 64, F(1, 10))
        high = bound(2, p, 16, 128, 64, F(3, 10))
        self.assertGreater(high, low)
        g = geometry(16, 128, 64)
        logy = -math.exp(p[0])*2/g['n']
        delta = output_threshold(g['n'], F(3, 10))-output_threshold(g['n'], F(1, 10))
        self.assertAlmostEqual(high-low, -delta*logy/math.log(2), places=10)

    def test_invalid_threshold(self):
        for n, d in ((0, F(1, 2)), (2, F(0)), (2, F(1))):
            with self.assertRaises(ValueError):
                output_threshold(n, d)

    def test_density_reference_and_monotonicity(self):
        p = (-2., -.2, math.log(16), 0.)
        reference = density_reference(p, 20000, 16, 128, 64, 16, True)
        self.assertEqual(density_bound(p, 20000, 16, 128, 64, F(1, 10)), reference)
        self.assertGreater(density_bound(p, 20000, 16, 128, 64, F(3, 10)), reference)


if __name__ == '__main__':
    unittest.main()
