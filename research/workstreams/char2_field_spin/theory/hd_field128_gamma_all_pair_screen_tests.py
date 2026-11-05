"""Tiny all-pair algebra/adapter tests; no actual K20 selected-weight gate."""

from fractions import Fraction as F
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma_all_pair_screen as screen


def product(left, right):
    return tuple(sum((left[2*i+k]*right[2*k+j] for k in range(2)), F(0))
                 for i in range(2) for j in range(2))


def toy_rows(slots):
    return tuple(tuple(F((j+2)*(i+1)+i*i+1, 31+j+3*i) for i in range(4))
                 for j in range(2*slots+1))


def direct_pairs(rows, slots, alpha, *, old_partition=False):
    """Independent scalar ordered histogram loops, with the frozen partition."""
    result = np.zeros((4, 2*slots+1, 2*slots+1))
    for s in range(slots+1):
        for f in range(slots-s+1):
            for t in range(slots+1):
                for g in range(slots-t+1):
                    count = comb(slots, s)*comb(slots-s, f)*comb(slots, t)*comb(slots-t, g)
                    left, right = rows[s+2*f], rows[t+2*g]
                    if old_partition and (s > 1 or t > 1):
                        values = product(tuple(float(x)**alpha for x in left),
                                         tuple(float(x)**alpha for x in right))
                    else:
                        values = tuple(float(x)**alpha for x in product(left, right))
                    result[:, s+t, f+g] += count*np.asarray(values, dtype=float)
    return result


class AllPairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z, cls.alpha = .5, .25
        cls.coeff = screen.AliasCoefficients(cls.z, cls.alpha, state=24)
        cls.control = screen.gamma.AliasCoefficients(cls.z, cls.alpha, state=24)
        cls.bank = screen.AliasBankScreen(cls.z, cls.alpha, ((-8., -8.),), -8., state=24,
                                         grid=64, _coefficients=cls.coeff)

    def test_exact_histogram_partition_and_full_degree(self):
        rows = screen.histogram_geometry()
        self.assertEqual(len(rows), 561)
        self.assertEqual(sum(row[2] for row in rows), 3**32)
        self.assertEqual(sorted(rows), sorted(screen.gamma.frozen.fixed16.geometry(128)))
        self.assertEqual(self.coeff.composed.shape, (4, 65, 65))
        valid = np.isfinite(self.coeff.composed[0])
        self.assertEqual(int(valid.sum()), 2145)
        for s, f in np.ndindex(65, 65):
            self.assertEqual(bool(valid[s, f]), s+f <= 64)
        self.assertTrue(np.all(np.isfinite(self.coeff.composed[:, 0, 64])))
        self.assertTrue(np.all(np.isfinite(self.coeff.composed[:, 64, 0])))
        self.assertTrue(np.all(np.isneginf(self.coeff.single)))
        # With every local entry one, alpha1 products have every entry two.
        values = screen.all_pair_coefficients(np.zeros((5, 2, 2)), 1., slots=2)
        for s, f, count in screen.histogram_geometry(4):
            np.testing.assert_allclose(values[:, s, f], log2(2*count), rtol=0, atol=2e-14)

    def test_independent_ordered_toy_products_and_alpha1(self):
        rows, slots = toy_rows(2), 2
        logs = np.log2(np.asarray(rows, dtype=float).reshape((-1, 2, 2)))
        for alpha in (1., .5, .125):
            actual = np.exp2(screen.all_pair_coefficients(logs, alpha, slots=slots))
            expected = direct_pairs(rows, slots, alpha)
            np.testing.assert_allclose(actual, expected, rtol=2e-15, atol=1e-13)
            # Both chronological orders occur, but alpha precedes their sum.
            order01, order10 = product(rows[0], rows[1]), product(rows[1], rows[0])
            self.assertNotEqual(order01, order10)
            expected10 = slots*(np.asarray(order01, dtype=float)**alpha
                                +np.asarray(order10, dtype=float)**alpha)
            np.testing.assert_allclose(actual[:, 1, 0], expected10, rtol=2e-15, atol=1e-13)

    def test_coefficientwise_old_partition_dominance_and_strict_gain(self):
        rows = toy_rows(3)
        logs = np.log2(np.asarray(rows, dtype=float).reshape((-1, 2, 2)))
        for alpha in (.25, .75, 1.):
            actual = np.exp2(screen.all_pair_coefficients(logs, alpha, slots=3))
            previous = direct_pairs(rows, 3, alpha, old_partition=True)
            self.assertTrue(np.all(actual <= previous*(1+3e-15)+1e-12))
            np.testing.assert_allclose(actual[:, :2], previous[:, :2], rtol=3e-15, atol=1e-12)
            if alpha == 1:
                np.testing.assert_allclose(actual, previous, rtol=3e-15, atol=1e-12)
            else:
                self.assertGreater(float(np.max(previous-actual)), 1.)

    def test_q128_pair_scaling_and_frozen_atlas_adapter(self):
        for odds in ((-8., -8.), (-3., -5.)):
            actual = self.coeff.pair_logs(odds)
            previous = self.control.pair_logs(odds)
            self.assertTrue(np.all(actual <= previous+2e-12))
        control_table = screen.alias.alias_table(self.z, self.alpha, (-8., -8.), grid=64,
                                                 _coefficients=self.control)
        self.assertTrue(np.all(self.bank.tables[0].positive_numerator <= control_table.positive_numerator+1e-9))
        self.assertLess(self.bank.tables[0].diagnostic['normalized_row_sum_error'], 1e-9)
        self.assertEqual(self.bank.regional(65).shape, (33, 2, 2))
        self.assertIs(screen.AliasBankScreen.regional, screen.alias.AliasBankScreen.regional)

    def test_alpha1_actual_occupancy_replay_and_untruncated_endpoint(self):
        coeff = screen.AliasCoefficients(self.z, 1., state=24)
        control = screen.gamma.AliasCoefficients(self.z, 1., state=24)
        np.testing.assert_allclose(coeff.pair_logs((-6., -6.)), control.pair_logs((-6., -6.)),
                                   rtol=0, atol=2e-12)
        raw = screen.gamma.local_logs(self.z, state=24)[64]
        expected = np.asarray(screen.old.multiply(raw, raw)).ravel()
        np.testing.assert_allclose(coeff.composed[:, 0, 64], expected, rtol=0, atol=1e-12)

    def test_metadata_source_closure_and_mock_selected_dispatch(self):
        regional = np.zeros((2, 2, 2))
        with patch.object(self.bank, 'regional', return_value=regional):
            self.bank.regional_diagnostics = dict(mock=True)
            row = self.bank.evaluate(2, method='runs', optimize_potential=False)
        self.assertEqual(row['construction'], screen.gamma.construction(24))
        self.assertTrue(row['construction_unchanged'])
        self.assertEqual(row['raw_histogram_pair_products'], 561**2)
        self.assertEqual(row['raw_occupancy_pair_products'], 65**2)
        self.assertEqual(row['maximum_paired_occupancy'], 128)
        self.assertEqual(row['composed_pairs_per_region'], 64)
        self.assertFalse(row['directed'])
        self.assertFalse(row['whole_code_certificate'])
        self.assertNotIn('good_geometry_types', row)
        self.assertEqual(self.bank.witness()['coefficient_method'], row['coefficient_method'])
        for name in ('hd_field128_gamma_all_pair_screen.py', 'hd_field128_full_rs_gamma_screen.py',
                     'hd_field128_lift_alias_screen.py', 'hd_field128_full_rs_gamma_alias_search.py',
                     'hd_field128_gamma_profile_diagnostic.py', 'hd_shifted_pairs_count.py'):
            self.assertIn(name, screen.source_pins())

    def test_mock_exact_profile_support_fresh_prefix_and_witness(self):
        regions = np.array((((-1., -7.), (-11., -4.)), ((-3., -5.), (-9., -2.))))
        expected = screen.profiles.support_sum(screen.old.PAIRS, 2, regions)
        with patch.object(self.bank, 'regional', return_value=regions):
            self.bank.regional_diagnostics = dict(mock=True)
            actual = screen.selected_screen(self.bank, 2, support_method='exact-profiles')
        prefix = self.bank.alpha*self.bank.density*2+self.bank.constant
        self.assertEqual(actual['log2_support_sum_screen'], expected['log2_support_sum_screen'])
        self.assertEqual(actual['largest_profile'], expected['largest_profile'])
        self.assertEqual(actual['log2_outer_prefix'], prefix)
        self.assertEqual(actual['log2_upper_screen'], prefix+expected['log2_support_sum_screen'])
        self.assertEqual(actual['log2_potential_ratio'], self.bank.log_ratio)
        self.assertEqual(actual['construction'], screen.gamma.construction(24))
        self.assertEqual(actual['support_method'], 'exact-profiles')
        self.assertFalse(actual['common_potential_used'])
        self.assertFalse(actual['directed'])
        self.assertFalse(actual['whole_code_certificate'])
        with patch.object(self.bank, 'evaluate', return_value=dict(mock=True)) as evaluate:
            self.assertEqual(screen.selected_screen(self.bank, 2), dict(mock=True, support_method='runs'))
            self.assertEqual(evaluate.call_args.kwargs['method'], 'runs')
        for kwargs in (dict(support_method='unknown'),
                       dict(support_method='exact-profiles', optimize_potential=True)):
            with self.assertRaises(ValueError):
                screen.selected_screen(self.bank, 2, **kwargs)

    def test_input_cache_and_budget_guards(self):
        for slots in (0, 33, 1.5):
            with self.assertRaises(ValueError):
                screen.histogram_geometry(slots)
        for logs, alpha in ((np.zeros((2, 2, 2)), .5), (np.zeros((65, 2, 2)), 0),
                            (np.full((65, 2, 2), np.nan), .5)):
            with self.assertRaises(ValueError):
                screen.all_pair_coefficients(logs, alpha)
        for kwargs in (dict(state=32), dict(q=256), dict(state=20)):
            with self.assertRaises(ValueError):
                screen.AliasCoefficients(.5, .5, **kwargs)
        for coeff, z, state in ((self.control, self.z, 24), (self.coeff, .25, 24),
                                (self.coeff, self.z, 28)):
            with self.assertRaises(ValueError):
                screen.AliasBankScreen(z, self.alpha, ((-8., -8.),), -8., state=state,
                                       _coefficients=coeff)
        with self.assertRaises(TimeoutError):
            screen.AliasCoefficients(.5, .5, _budget=screen.old.Budget(1e-12))


if __name__ == '__main__':
    unittest.main()
