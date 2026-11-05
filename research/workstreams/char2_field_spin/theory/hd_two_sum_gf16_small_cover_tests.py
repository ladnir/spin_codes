from fractions import Fraction as F
from math import comb, log2
import unittest

import numpy as np

from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_small_cover import (dense_moments, source_envelope,
    support_polynomial, cohort_prefix)
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_floor import column_log_coefficients


class SmallCoverTests(unittest.TestCase):
    def test_aggregate_exact(self):
        nums, den = dense_moments()
        for w in range(17):
            expected = sum((F(arity_holder_moments(a,16,16,8)[0][w],
                              arity_holder_moments(a,16,16,8)[1]) for a in range(3,17)), F(0))
            self.assertEqual(F(nums[w], den), expected)

    def test_envelope_dominates_every_full_count(self):
        raw = {(m,f): np.array([[m-f, f], [f+1, m+f]], dtype=float)
               for m in range(8) for f in range(m//2+1)}
        got = source_envelope(raw,7)
        self.assertEqual(got.shape, (16,8,2,2))
        for (m,f), matrix in raw.items():
            self.assertTrue(np.all(got[:,m] >= matrix))

    def test_pairs_polynomial_control(self):
        for n in range(5):
            caps = {kind: -0.3*(kind[0]+kind[1]) for kind in categories(n)}
            np.testing.assert_allclose(support_polynomial(caps,0,n),
                                       column_log_coefficients(np.array(list(caps.values())),n))

    def test_support_polynomial_normalization(self):
        for n1 in range(4):
            for n2 in range(4):
                caps = {kind: 0. for kind in categories(n1+n2)}
                poly = support_polynomial(caps,n1,n2)
                expected = [log2(comb(n1+2*n2,d))+d*log2(255) for d in range(n1+2*n2+1)]
                np.testing.assert_allclose(poly, expected, atol=1e-12)

    def test_once_only_prefix(self):
        z,p = F(1),F(0)
        actual = cohort_prefix(17,0,16,1,z,p)
        self.assertAlmostEqual(actual,log2(comb(8192,17)*17*120**16/255))
        with self.assertRaises(ValueError):
            cohort_prefix(17,0,15,2,z,p)


if __name__ == '__main__':
    unittest.main()
