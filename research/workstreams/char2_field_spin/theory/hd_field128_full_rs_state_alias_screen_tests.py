"""Tiny coefficient/FFT-interface tests; no selected middle-weight gates."""

from fractions import Fraction as F
import unittest

import numpy as np

import hd_field128_full_rs_state_alias_screen as candidate


old, local, alias = candidate.old, candidate.local, candidate.alias


class StateAliasTests(unittest.TestCase):
    def test_coefficients_replay_all_public_state_sparse_prefixes(self):
        for state in local.STATES:
            coefficients = candidate.AliasCoefficients(.5, .1, state=state, q=256)
            single = local.single_coefficients(.5, .1, 6, 256, state, old.Budget(10))
            pair = local.good_pair_coefficients(.5, .1, 6, 256, state, old.Budget(10))
            for f in range(4):
                for s in range(7-2*f):
                    np.testing.assert_allclose(coefficients.single[:, s, f], single[:, s, f], rtol=0, atol=1e-12)
                    if s <= 2:
                        np.testing.assert_allclose(coefficients.composed[:, s, f], pair[:, s, f], rtol=0, atol=1e-12)
            self.assertEqual(coefficients.single.shape, (4, 33, 33))
            self.assertEqual(coefficients.composed.shape, (4, 3, 65))

    def test_alpha1_pair_is_the_ordered_full_single_square(self):
        for state in local.STATES:
            coefficients = candidate.AliasCoefficients(.9, 1., state=state, q=256)
            odds = (-3., -2.)
            tilt = np.arange(33)[:, None]*odds[0]+np.arange(33)[None, :]*odds[1]
            one = np.logaddexp2.reduce((coefficients.single+tilt[None, :, :]).reshape((4, -1)), axis=1).reshape((2, 2))
            expected = old.multiply(one, one)
            np.testing.assert_allclose(coefficients.pair_logs(odds), expected, rtol=0, atol=2e-12)

    def test_frozen_fft_engine_accepts_custom_coefficients_and_metadata(self):
        coefficients = candidate.AliasCoefficients(.9, .1, state=24, q=256)
        bank = candidate.AliasBankScreen(.9, .1, ((-4., -3.),), -8., state=24, q=256,
                                         grid=64, _coefficients=coefficients)
        self.assertEqual(bank.metadata()['construction'], local.construction(24))
        self.assertEqual(bank.witness()['state_dimension'], 24)
        self.assertEqual(bank.witness()['good_geometry_types'], 65)
        self.assertFalse(bank.witness()['whole_code_certificate'])
        self.assertFalse(bank.witness()['directed'])
        self.assertEqual(bank.tables[0].logs.shape, (4, 64, 64))
        self.assertTrue(np.all(np.isfinite(bank.tables[0].logs)))
        # Only a tiny coefficient-prefix normalization, not a distance screen.
        regions = bank.regional(3)
        self.assertEqual(regions.shape, (2, 2, 2))
        self.assertTrue(bank.regional_diagnostics['no_zero_coefficient_claims'])
        with self.assertRaises(ValueError):
            candidate.AliasBankScreen(.9, .1, ((-4., -3.),), -8., state=32, q=256,
                                       grid=64, _coefficients=coefficients)

    def test_q128_coefficients_stay_in_logs_below_float_value_range(self):
        coefficients = candidate.AliasCoefficients(2.**-200, .1, state=32)
        self.assertTrue(np.all(np.isfinite(coefficients.single[:, 1, 0])))
        self.assertLess(coefficients.single[2, 1, 0], -400)
        pair = coefficients.pair_logs((-4., -3.))
        self.assertTrue(np.all(np.isfinite(pair)))
        self.assertFalse(np.any(np.isnan(coefficients.composed)))

    def test_only_source_literals_are_reused_with_new_state_identity(self):
        source = dict(log2_upper_rational='not inherited', fixed_witness=dict(
            pair_swaps=True, orientation_aggregation='independent-fair-swaps', grid=64,
            source_exact_witness=dict(z='7/8', alpha='1/10', potential='1/4', markers=['3/4 5/8', '1/2 1/3'])))
        result = candidate.proposal_from_artifact(source, dict(path='mock', sha256='mock'),
                                                   state=32, distance=F(9, 20), old_source=True)
        self.assertEqual(result['construction'], local.construction(32))
        self.assertEqual(result['witness']['markers'], [['3/4', '5/8'], ['1/2', '1/3']])
        self.assertEqual(F(result['witness']['z']), F(str(float(F(7, 8)**16))))
        self.assertEqual(F(result['witness']['alpha']), F(str(float(F(1, 160)))))
        self.assertFalse(result['old_numerical_bounds_reused'])
        self.assertNotIn('log2_upper_rational', result)
        override = candidate.proposal_from_artifact(source, {}, state=24, distance=F(2, 5),
                                                     old_source=True, z='1/2', alpha='1/128')
        self.assertEqual(override['witness']['z'], '1/2')
        self.assertEqual(override['witness']['alpha'], '1/128')
        pins = candidate.source_pins()
        for name in ('hd_field128_full_rs_state_alias_screen.py', 'hd_field128_full_rs_state_screen.py',
                     'hd_field128_lift_alias_screen.py', 'hd_field128_lift_tail_screen.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
