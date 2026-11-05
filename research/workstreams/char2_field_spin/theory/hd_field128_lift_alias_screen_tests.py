"""Small GF256 replays and numerical safeguards, not extension-field gates."""

import unittest
import numpy as np

import hd_field128_lift_alias_screen as alias
from hd_shifted_good_pair_screen import good_pair_coefficients
from hd_shifted_projected_screen import alias_table as legacy_table
from hd_shifted_potential_band_screen import FixedAtlasScreen
from hd_shifted_swap_screen import swap_coefficients


class AliasScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z, cls.alpha, cls.odds = .9, .125, ((-4., -3.), (-3., -2.))
        cls.coefficients = alias.AliasCoefficients(cls.z, cls.alpha, q=256)

    def test_q256_coefficients_replay_before_fft(self):
        with np.errstate(divide='ignore'):
            single = np.log2(swap_coefficients(self.z, self.alpha))
            composed = np.log2(good_pair_coefficients(self.z, self.alpha).reshape((3, 65, 4)).transpose((2, 0, 1)))
        np.testing.assert_allclose(self.coefficients.single, single, rtol=0, atol=3e-12)
        np.testing.assert_allclose(self.coefficients.composed, composed, rtol=0, atol=3e-12)

    def test_q256_full_alias_atlas_replay(self):
        table = alias.alias_table(self.z, self.alpha, self.odds[0], q=256, grid=64,
                                  _coefficients=self.coefficients)
        expected = legacy_table(self.z, self.alpha, self.odds[0], 64)
        np.testing.assert_allclose(table.logs[table.reliable], expected[table.reliable], rtol=0, atol=2e-6)
        # The below-floor cells depend on roundoff relative to a 1e-12
        # cushion. Compare them as proposal floors, not accurate coefficients.
        np.testing.assert_allclose(table.logs, expected, rtol=0, atol=.002)
        self.assertGreater(table.diagnostic['reliable_entries'], 0)
        self.assertTrue(np.all(np.isfinite(table.logs)))

    def test_q256_legacy_bank_point_replay(self):
        screen = alias.AliasBankScreen(self.z, self.alpha, self.odds, -8., q=256, grid=64,
                                      outer_measure='inherited-diagonal', floor_policy='legacy-cushion',
                                      _coefficients=self.coefficients)
        old = FixedAtlasScreen(self.z, self.alpha, self.odds, -8., grid=64)
        for weight in (1, 2, 3):
            actual = screen.evaluate(weight, method='runs')
            expected = old.evaluate(weight)
            self.assertAlmostEqual(actual['log2_upper_screen'], expected['log2_upper_screen'], places=5)
            self.assertFalse(actual['directed'])
            self.assertFalse(actual['whole_code_certificate'])

    def test_default_replaces_unreliable_entries_by_positive_bound(self):
        screen = alias.AliasBankScreen(self.z, self.alpha, (self.odds[0],), -8., q=256, grid=64,
                                      _coefficients=self.coefficients)
        weight = 3
        regions = screen.regional(weight)
        full = alias.full_counts(alias.PAIRS, weight)
        single = weight-2*full
        table = screen.tables[0]
        normalizer = (alias.LOG_FACTORIAL[alias.PAIRS]-alias.LOG_FACTORIAL[alias.PAIRS-weight+full]
                      -alias.LOG_FACTORIAL[single]-alias.LOG_FACTORIAL[full]
                      +single*self.odds[0][0]+full*self.odds[0][1])
        expected = np.where(table.reliable[:, single % 64, full % 64],
                            table.logs[:, single % 64, full % 64], table.positive_numerator[:, None])
        np.testing.assert_allclose(regions, (expected-normalizer).T.reshape((-1, 2, 2)))
        self.assertTrue(screen.regional_diagnostics['no_zero_coefficient_claims'])

    def test_perron_scaling_with_unrepresentable_unscaled_entries(self):
        matrix = np.array(((0., 2000.), (-3000., -1.)))
        rho, ratio = alias.perron_logs(matrix)
        factors = np.array(((0., ratio), (-ratio, 0.)))
        normalized = np.exp2(matrix+factors-rho)
        np.testing.assert_allclose(normalized.sum(axis=1), np.ones(2), rtol=0, atol=2e-12)
        self.assertTrue(np.isfinite(rho))
        self.assertLess(ratio, -2000)

    def test_extension_fft_and_boundary_underflow_are_explicit(self):
        coefficients = alias.AliasCoefficients(.5, 1/128, q=1 << 128)
        table = alias.alias_table(.5, 1/128, (-64., -3.), q=1 << 128, grid=64,
                                  _coefficients=coefficients)
        self.assertTrue(np.all(np.isfinite(table.logs)))
        self.assertGreater(table.diagnostic['normalized_input_underflow_terms'], 0)
        self.assertGreater(table.diagnostic['below_floor_entries'], 0)

    def test_full_bank_limit_and_memory_guard(self):
        rows = tuple((float(i), -3.) for i in range(73))
        self.assertEqual(len(alias.checked_markers(rows)), 73)
        with self.assertRaises(ValueError):
            alias.AliasBankScreen(.9, .1, rows, -8., grid=256, max_table_bytes=1)
        with self.assertRaises(ValueError):
            alias.checked_markers(tuple((float(i), 0.) for i in range(129)))


if __name__ == '__main__':
    unittest.main()
