"""Tiny scalar/normalization tests only; no bank, FFT, files or full objective."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
import numpy as np
import hd_field128_marker_saddles_v40 as pilot


class Tests(unittest.TestCase):
    def test_markers(self):
        original = [['1', '2'], ['3', '4']]
        extended, labels, centers = pilot.marker_extension(original, F(1, 256))
        self.assertEqual(extended[:2], original)
        self.assertEqual(original, [['1', '2'], ['3', '4']])
        self.assertIn(['22/1013', '1'], extended)
        self.assertEqual(centers['quiet_weighted'], [str(F(float(F(22, 1013))*2.**.5)), '2'])
        self.assertEqual(len(labels), 18)
        self.assertEqual(len({tuple(F(x) for x in pair) for pair in extended}), len(extended))

    def test_quiet_live_cap_formulas(self):
        local = np.zeros((65, 2, 2))
        local[:, 0, 0] = -128*np.arange(65)+.125
        local[:, 1, 1] = -64.
        caps = pilot.diagonal_caps(local, .5, .25)
        self.assertEqual(caps['quiet_residual_max'], .125)
        self.assertEqual(caps['quiet_regional'], 128*.25*.125-.25*4096*128)
        self.assertEqual(caps['live_regional_from_local'], -2048.)
        self.assertAlmostEqual(caps['live_symbolic_uniform'], -64.)
        triangle = np.array([[0., 7.], [-np.inf, 3.]])
        capped = pilot.cap_j0(triangle, caps)
        self.assertEqual(capped[0, 1], 7.)
        self.assertEqual(triangle[0, 0], 0.)
        self.assertEqual(capped[1, 1], caps['live_regional'])
        with self.assertRaises(ValueError):
            pilot.cap_j0(np.zeros((2, 2)), caps)

    def test_route_histogram_cancellation(self):
        # Four1-slot packets: fixed e=1,s=2,f=1, total occupancy4.
        alpha, quiet, live = .37, [1., .2, .03], [.3, .4, .2]
        q = 8
        c = max(log2(quiet[j])+j*log2(q) for j in range(3))
        d = max(live)
        terms = [word for word in product(range(3), repeat=4)
                 if word.count(0) == 1 and word.count(1) == 2 and word.count(2) == 1]
        h = comb(4, 2)*comb(2, 1)
        self.assertEqual(len(terms), h)
        qsum = sum(np.prod([quiet[j]**alpha for j in word]) for word in terms)/h
        vsum = sum(np.prod([live[j]**alpha for j in word]) for word in terms)/h
        self.assertLessEqual(qsum, 2.**(4*alpha*c-4*alpha*log2(q))*(1+1e-13))
        self.assertLessEqual(vsum, d**(4*alpha)*(1+1e-13))

if __name__ == '__main__':
    unittest.main()
