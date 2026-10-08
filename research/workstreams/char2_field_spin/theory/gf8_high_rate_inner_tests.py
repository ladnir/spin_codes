"""Local algebra and proof-transfer regressions; no whole-code claim."""
from fractions import Fraction as F
import math
import unittest

from flint import ctx, fmpq
import numpy as np

import gf8_high_rate_inner_screen as target
from gf8_packet_width_screen import local_matrix as coarse_local
from gf8_rank_history_screen import cone_preserved


class HighRateInnerTests(unittest.TestCase):
    def test_continuation_cone_exact(self):
        for y in (F(1, 4), F(3, 4), F(999, 1000), F(1)):
            for tau in (F(1), F(8), F(255)):
                for weight in range(9):
                    self.assertTrue(cone_preserved(np.asarray(target.local_matrix(weight, y, tau), dtype=object)))

    def test_cap_transfer_dominated_by_two_state(self):
        for y in (F(1, 4), F(3, 4), F(999, 1000)):
            for tau in (F(1), F(8), F(255)):
                for weight in range(9):
                    refined = target.local_matrix(weight, y, tau)
                    coarse = coarse_local(weight, y, tau, 8, 4, target.MINIMA)
                    self.assertEqual(refined[0][0], coarse[0][0])
                    self.assertEqual(sum(refined[0][1:]), coarse[0][1])
                    for cap in range(1, 5):
                        self.assertLessEqual(refined[cap][0], coarse[1][0])
                        self.assertLessEqual(sum(refined[cap][1:]), coarse[1][1])

    def test_aggregate_matches_positive_sum(self):
        for x, y, tau in ((.1, .9, 12.), (.4, .5, 255.), (1., .1, 1.)):
            direct = sum(math.comb(8, w)*x**w*np.asarray(target.local_matrix(w, y, tau))
                         for w in range(9))
            np.testing.assert_allclose(target.packet(math.log(x), math.log(y), math.log(tau)),
                                       direct, rtol=2e-13, atol=1e-15)

    def test_exact_polynomial_one_packet(self):
        old = ctx.prec
        try:
            ctx.prec = 256
            local = [target.local_matrix(w, F(9, 10), F(17, 2)) for w in range(9)]
            polynomial = target.polynomial_power(local, 1, 8)
            for w in range(9):
                for i in range(5):
                    for j in range(5):
                        expected = math.comb(8, w)*local[w][i][j]
                        exact = fmpq(expected.numerator, expected.denominator)
                        value = polynomial[5*i+j][w]
                        self.assertLessEqual(value.lower().fmpq(), exact)
                        self.assertGreaterEqual(value.upper().fmpq(), exact)
        finally:
            ctx.prec = old

    def test_gv_reference_not_near_rate_ceiling(self):
        value = target.gv_comparison()
        self.assertGreater(value['gv_rate'], .84)
        self.assertGreater(value['ideal_random_code_exponent_bits_per_output'], .73)

    def test_feedthrough_branches(self):
        profile = target.feedthrough_profile(4)
        self.assertEqual(profile['birth_minima'], [8, 6, 5, 5, 3, 2, 2, 1])
        self.assertEqual(profile['return_minima'], [8, 6, 6, 4, 4, 2, 2, 1])
        self.assertFalse(profile['whole_code'])
        with self.assertRaises(ValueError):
            target.feedthrough_profile(2)

    def test_no_return_diagnostic_not_called_rigorous(self):
        profile = target.sparse_point(1, (6.58, math.log(255)), no_returns=True)
        self.assertTrue(profile['no_returns_diagnostic'])
        self.assertFalse(profile['rigorous'])


if __name__ == '__main__':
    unittest.main()
