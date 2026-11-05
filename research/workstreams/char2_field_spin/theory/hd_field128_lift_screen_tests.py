"""Small algebraic/replay tests, not a numerical distance gate."""

from fractions import Fraction as F
from math import log2
import unittest

import numpy as np

import hd_field128_lift_screen as lift
import hd_grouped_packet_screen as legacy_local
import hd_shifted_swap_sparse as legacy_sparse
from hd_common_refresh_packets import local_moments, matrix
from hd_grouped_quiet import quiet_cap
from hd_shifted_outer_density import log_density_float


class LiftScreenTests(unittest.TestCase):
    def test_local_exact_q16_q256_algebra(self):
        # q16 checks the formulas only. The fixed GF256 maps do not embed there.
        occ = ((0, 0, 0, 0), (0, 0, 0, 1), (0, 0, 2, 4),
               (1, 3, 7, 15), (4, 4, 4, 4), (13, 14, 15, 16), (16, 16, 16, 16))
        for q in (16, 256):
            for z in (F(7, 8), F(1, 2)):
                local = local_moments(z, q=q, h=4)
                density = F(q**16-1, (q-1)**16)
                logs = lift.local_logs(float(z), q, occ, tighten_quiet=False).reshape((-1, 4))
                for index, profile in enumerate(occ):
                    k, birth, ret, live = matrix(profile, z, local, q=q, h=4)
                    exact = (k, density*birth, ret, density*live)
                    expected = np.array([lift.log_fraction(value) for value in exact])
                    np.testing.assert_allclose(logs[index], expected, rtol=0, atol=3e-12)

    def test_local_q256_replays_legacy_all_occupancies(self):
        for z in (.5, .9, .999):
            np.testing.assert_allclose(lift.local_logs(z, 256), legacy_local.local_logs(z),
                                       rtol=0, atol=5e-12)
        self.assertEqual(legacy_local.Q, 256.)
        self.assertEqual(legacy_local.LIMIT, (1 << 21)*7//20)

    def test_quiet_q256_exact(self):
        for z in (F(7, 8), F(1, 2)):
            self.assertAlmostEqual(lift.quiet_log(float(z), 256),
                                   lift.log_fraction(quiet_cap(z, q=256)), places=12)

    def test_extension_entries_survive_below_float_underflow(self):
        q, z = 1 << 128, F(7, 8)
        occ = ((0, 0, 0, 1), (1, 1, 1, 1), (4, 4, 4, 4))
        logs = lift.local_logs(float(z), q, occ).reshape((-1, 4))
        local, density = local_moments(z, q=q, h=4), F(q**16-1, (q-1)**16)
        for index, profile in enumerate(occ):
            k, birth, ret, live = matrix(profile, z, local, q=q, h=4)
            expected = [lift.log_fraction(value) for value in (k, density*birth, ret, density*live)]
            np.testing.assert_allclose(logs[index], expected, rtol=0, atol=2e-10)
        self.assertTrue(np.all(np.isfinite(logs)))
        self.assertLess(logs[0, 2], -2048)
        self.assertGreater(np.exp2(.01*logs[0, 2]), 0)

    def test_log_polynomial_retains_tiny_paths(self):
        base = np.full((4, 2, 1), -np.inf)
        base[(0, 3), 0, 0] = 0
        base[:, 1, 0] = [0, -2048, -4096, -512]
        product = lift.polynomial_product(base, base, 2, lift.Budget(10))
        self.assertAlmostEqual(product[1, 1, 0], -2047)
        self.assertAlmostEqual(product[2, 1, 0], -4095)
        self.assertAlmostEqual(product[3, 2, 0], -1024)

    def test_tiny_mark_keeps_positive_birth_and_return(self):
        occ = ((0, 0, 0, 1), (1, 2, 3, 4))
        for q, exponent in ((256, 64), (1 << 128, 200)):
            z = F(1, 1 << exponent)
            logs = lift.local_logs(float(z), q, occ).reshape((-1, 4))
            local, density = local_moments(z, q=q, h=4), F(q**16-1, (q-1)**16)
            for index, profile in enumerate(occ):
                k, birth, ret, live = matrix(profile, z, local, q=q, h=4)
                expected = [lift.log_fraction(value) for value in (k, density*birth, ret, density*live)]
                np.testing.assert_allclose(logs[index], expected, rtol=0, atol=5e-10)
            self.assertTrue(np.all(np.isfinite(logs)))

    def test_selected_small_r_q256_replay(self):
        for r, z, alpha, paired in ((1, .999, 1., True), (2, .99, .4, True), (3, .98, .2, False)):
            old = legacy_sparse.weight_screen(r, z, alpha, paired=paired)
            new = lift.weight_screen(r, z, alpha, q=256, paired=paired,
                                     outer_measure='inherited-diagonal', max_seconds=20)
            self.assertAlmostEqual(new['log2_upper_screen'], old['log2_upper_screen'], places=7)
            self.assertEqual(new['profile_count'], old['profile_count'])
            self.assertFalse(new['directed'])
            self.assertFalse(new['whole_code_certificate'])

    def test_outer_density_and_scope_validation(self):
        for mode in lift.OUTER_MEASURES:
            self.assertAlmostEqual(lift.outer_log_density(256, mode), log_density_float(mode), places=11)
        self.assertEqual(lift.outer_log_density(1 << 128), 16384.)
        for q in (16, 32, 1024):
            with self.assertRaises(ValueError):
                lift.regional_screen(.9, .1, 1, q=q)
        with self.assertRaises(ValueError):
            lift.local_logs(.9, 16)
        with self.assertRaises(ValueError):
            lift.weight_screen(65, .9, .1)
        with self.assertRaises(ValueError):
            lift.weight_screen(1, .9, .1, distance=F(1, 2))

    def test_resource_guard(self):
        budget = lift.Budget(1e-12)
        with self.assertRaises(TimeoutError):
            budget.check()
        for maximum in (0, 65):
            with self.assertRaises(ValueError):
                lift.regional_screen(.9, .1, maximum)


if __name__ == '__main__':
    unittest.main()
