"""Tiny width-aware coefficient/FFT tests; no selected-weight distance gate."""

from fractions import Fraction as F
import unittest

import numpy as np

import hd_field128_wide_rs_state_alias_screen as wide
import hd_field128_full_rs_state_alias_screen as narrow


class WideAliasTests(unittest.TestCase):
    def test_width64_replays_frozen_coefficients_atlas_and_regions(self):
        z, alpha, odds = .9, .1, (-4., -3.)
        for state in (16, 24, 32):
            a = wide.AliasCoefficients(z, alpha, q=256, width=64, state=state)
            b = narrow.AliasCoefficients(z, alpha, q=256, state=state)
            np.testing.assert_allclose(a.single, b.single, rtol=0, atol=1e-12)
            np.testing.assert_allclose(a.composed, b.composed, rtol=0, atol=1e-12)
            np.testing.assert_allclose(a.pair_logs(odds), b.pair_logs(odds), rtol=0, atol=1e-12)
            atlas = wide.alias_table(z, alpha, odds, q=256, width=64, state=state, grid=64, _coefficients=a)
            frozen = wide.alias.alias_table(z, alpha, odds, q=256, grid=64, _coefficients=b)
            np.testing.assert_allclose(atlas.logs[atlas.reliable], frozen.logs[atlas.reliable], rtol=0, atol=2e-6)
            np.testing.assert_allclose(atlas.logs, frozen.logs, rtol=0, atol=.002)
            np.testing.assert_allclose(atlas.positive_numerator, frozen.positive_numerator, rtol=0, atol=1e-11)
            actual = wide.AliasBankScreen(z, alpha, (odds,), -8., q=256, width=64, state=state,
                                          grid=64, _coefficients=a)
            previous = narrow.AliasBankScreen(z, alpha, (odds,), -8., q=256, state=state,
                                              grid=64, _coefficients=b)
            np.testing.assert_allclose(actual.regional(3), previous.regional(3), rtol=0, atol=2e-6)

    def test_wide_coefficients_match_direct_sparse_prefixes(self):
        for state in (48, 64):
            coefficients = wide.AliasCoefficients(.9, .1, state=state, q=256)
            single = wide.local.single_coefficients(.9, .1, 6, 256, 128, state, wide.old.Budget(60.))
            paired = wide.local.good_pair_coefficients(.9, .1, 6, 256, 128, state, wide.old.Budget(60.))
            for f in range(4):
                for s in range(7-2*f):
                    np.testing.assert_allclose(coefficients.single[:, s, f], single[:, s, f], rtol=0, atol=1e-12)
                    if s <= 2:
                        np.testing.assert_allclose(coefficients.composed[:, s, f], paired[:, s, f], rtol=0, atol=1e-12)
            self.assertEqual(coefficients.single.shape, (4, 65, 65))
            self.assertEqual(coefficients.composed.shape, (4, 3, 129))
            self.assertTrue(np.all(np.isfinite(coefficients.single[:, 1, 31])))

    def test_alpha1_wide_pair_equals_ordered_single_square(self):
        coefficients = wide.AliasCoefficients(.9, 1., state=48, q=256)
        odds = (-3., -2.)
        tilt = np.arange(65)[:, None]*odds[0]+np.arange(65)[None, :]*odds[1]
        one = np.logaddexp2.reduce((coefficients.single+tilt[None, :, :]).reshape((4, -1)), axis=1).reshape((2, 2))
        np.testing.assert_allclose(coefficients.pair_logs(odds), wide.old.multiply(one, one), rtol=0, atol=2e-12)

    def test_fft_mass_restoration_uses32_pair_powers_and64_packet_powers(self):
        z, alpha, odds = .9, .1, (-9., -9.)
        coefficients = wide.AliasCoefficients(z, alpha, q=256)
        table = wide.alias_table(z, alpha, odds, q=256, grid=64, _coefficients=coefficients)
        pair = coefficients.pair_logs(odds)
        expected = np.array(wide.old.power(pair, 32)).ravel()
        np.testing.assert_allclose(table.positive_numerator, expected, rtol=0, atol=1e-11)
        diagnostic = table.diagnostic
        self.assertEqual(diagnostic['paired_region_power'], 32)
        self.assertEqual(diagnostic['matrix_squarings'], 5)
        self.assertEqual(diagnostic['perron_pair_restoration_power'], 32)
        self.assertEqual(diagnostic['perron_packet_restoration_power'], 64)
        self.assertEqual(32*diagnostic['log2_pair_scale'], 64*diagnostic['log2_packet_scale'])
        ratio, scale = diagnostic['log2_similarity_ratio'], diagnostic['log2_pair_scale']
        similarity = np.array((0., ratio, -ratio, 0.))
        normalized = np.exp2(table.logs-32*scale+similarity[:, None, None]).sum(axis=(1, 2))
        exact_mass = np.exp2(expected-32*scale+similarity)
        np.testing.assert_allclose(normalized, exact_mass, rtol=0, atol=1e-8)
        self.assertLess(diagnostic['normalized_row_sum_error'], 1e-9)

    def test_direct_small_regional_coefficients_are_dominated(self):
        z, alpha, q, state = .9, .1, 256, 48
        bank = wide.AliasBankScreen(z, alpha, ((-9., -9.),), -8., q=q, state=state, grid=64)
        exact_prefix = wide.local.regional_screen(z, alpha, 3, q=q, width=128, state=state)
        for weight in (1, 2, 3):
            regions = bank.regional(weight)
            for full, row in enumerate(regions):
                self.assertTrue(np.all(row >= np.array(exact_prefix[weight, full])-2e-8))
        self.assertTrue(bank.regional_diagnostics['no_zero_coefficient_claims'])
        self.assertEqual(bank.metadata()['pair_slots_per_region'], 4096)
        self.assertEqual(bank.metadata()['regional_even_odd_pair_power'], 128)
        self.assertEqual(bank.metadata()['single_histogram_terms'], 2145)
        self.assertEqual(bank.metadata()['good_geometry_types'], 129)
        self.assertEqual(bank.metadata()['raw_good_pair_products'], 16641)
        self.assertEqual(bank.metadata()['maximum_packet_occupancy'], 128)
        self.assertEqual(bank.metadata()['maximum_paired_occupancy'], 256)

    def test_q128_underflow_is_reported_and_parameter_mismatch_rejected(self):
        coefficients = wide.AliasCoefficients(.5, 1/128, state=64)
        table = wide.alias_table(.5, 1/128, (-64., -3.), state=64, grid=64, _coefficients=coefficients)
        self.assertTrue(np.all(np.isfinite(table.logs)))
        self.assertGreater(table.diagnostic['normalized_input_underflow_terms'], 0)
        self.assertGreater(table.diagnostic['below_floor_entries'], 0)
        for change in (dict(state=48), dict(width=64, state=32), dict(joint=False)):
            with self.assertRaises(ValueError):
                wide.alias_table(.5, 1/128, (-4., -3.), _coefficients=coefficients, **change)
        with self.assertRaises(ValueError):
            wide.AliasBankScreen(.5, .1, ((-4., -3.),), -8., max_table_bytes=1)

    def test_source_literals_get_explicit_new_width_identity_and_pins(self):
        artifact = dict(log2_upper_rational='not inherited', fixed_witness=dict(
            pair_swaps=True, orientation_aggregation='independent-fair-swaps', grid=64,
            source_exact_witness=dict(z='7/8', alpha='1/10', potential='1/4', markers=['3/4 5/8', '1/2 1/3'])))
        proposal = wide.proposal_from_artifact(artifact, {'sha256': 'mock'}, old_source=True, r=4096)
        self.assertEqual(proposal['construction'], wide.local.construction(128, 48))
        self.assertEqual(proposal['inner_width'], 128)
        self.assertEqual(proposal['state_dimension'], 48)
        self.assertEqual(proposal['witness']['markers'], [['3/4', '5/8'], ['1/2', '1/3']])
        self.assertEqual(F(proposal['witness']['alpha']), F(str(float(F(1, 160)))))
        self.assertFalse(proposal['old_numerical_bounds_reused'])
        self.assertNotIn('log2_upper_rational', proposal)
        artifact['center'] = 4000
        with self.assertRaises(ValueError):
            wide.proposal_from_artifact(artifact, {}, old_source=True, r=4096)
        pins = wide.source_pins()
        for name in ('hd_field128_wide_rs_state_alias_screen.py', 'hd_field128_wide_rs_state_screen.py',
                     'hd_field128_lift_alias_screen.py', 'hd_field128_lift_tail_screen.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
