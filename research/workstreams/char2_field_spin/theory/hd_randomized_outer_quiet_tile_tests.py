"""Checks for conditioning the quiet fractional moment on tile occupancies."""

from fractions import Fraction as F
from math import comb, log2
import unittest

from flint import arb, ctx

from hd_randomized_outer_joint_quiet import joint_moment_logs
from hd_randomized_outer_quiet_tile import (
    directed_endpoint, fractional_profile_enclosures, fractional_profile_logs,
    join, profile_moments,
)


class TileConditioningTests(unittest.TestCase):
    def test_profile_count_normalization_and_zero_atom(self):
        profiles = profile_moments(F(1, 4))
        self.assertEqual(len(profiles), 12870)
        mass = [0]*65
        for profile, active, multiplicity, quiet, free in profiles:
            self.assertEqual(active, sum(profile))
            mass[active] += multiplicity
            self.assertLessEqual(F(1, 65536**active), quiet)
            self.assertLessEqual(quiet, free)
        self.assertEqual(mass, [comb(64, j) for j in range(65)])

    def test_two_occupied_tile_direct_mask_sum(self):
        q, z = 65536, F(1, 4)
        profiles = {profile: quiet for profile, _, _, quiet, _ in profile_moments(z)}
        for a, b in ((1, 1), (2, 3), (4, 7), (8, 8)):
            direct = F(0)
            for v in range(9):
                for w in range(9):
                    d, e = max(a-v, 0), max(b-w, 0)
                    rank = join(min(4, d), min(4, e))
                    coefficient = comb(8, v)*comb(8, w)*z**(16-v-w)*(1-z)**(v+w)
                    direct += coefficient*F(q**(d+e), q**(a+b+rank))
            direct = min(direct, F(1, q**min(a+b, 4)))
            self.assertEqual(profiles[(0,)*6+(a, b)], direct)

    def test_unit_mark_exact_syndrome_and_free_moments(self):
        for _, active, _, quiet, free in profile_moments(F(1)):
            self.assertEqual(quiet, F(1, 65536**min(active, 4)))
            self.assertEqual(free, 1)

    def test_jensen_order_and_old_bound(self):
        z, alpha = F(1, 4), F(1, 32)
        quiet, free = fractional_profile_logs(z, alpha)
        mean_quiet, mean_free = fractional_profile_logs(z, F(1))
        old_quiet, old_free = joint_moment_logs(z)
        strict = False
        for j in range(65):
            self.assertLessEqual(quiet[j], float(alpha)*mean_quiet[j]+1e-12)
            self.assertLessEqual(free[j], float(alpha)*mean_free[j]+1e-12)
            self.assertLessEqual(mean_quiet[j], old_quiet[j]+1e-10)
            self.assertAlmostEqual(mean_free[j], old_free[j], places=10)
            strict |= quiet[j] < float(alpha)*mean_quiet[j]-0.01
        self.assertTrue(strict)

    def test_directed_profile_and_fixed_witness(self):
        previous = ctx.prec
        ctx.prec = 192
        try:
            for input_mark in (F(1), F(257, 256)):
                logs = fractional_profile_logs(F(1, 4), F(1, 32), input_mark=input_mark)
                bounds = fractional_profile_enclosures(F(1, 4), F(1, 32), input_mark=input_mark)
                for float_array, arb_array in zip(logs, bounds):
                    for floating, directed in zip(float_array, arb_array):
                        self.assertTrue(directed > 0)
                        self.assertAlmostEqual(floating, float(directed.log()/arb(2).log()), places=10)
            result = directed_endpoint(1536, F(1, 4), F(1, 32), exact_regional=True)
            self.assertTrue(arb(result['marker_log2_upper']) > 8584)
            self.assertTrue(arb(result['marker_log2_upper']) < 8585)
            self.assertTrue(arb(result['exact_regional_log2_upper']) > 6898)
            self.assertTrue(arb(result['exact_regional_log2_upper']) < 6899)
        finally:
            ctx.prec = previous


if __name__ == '__main__':
    unittest.main()
