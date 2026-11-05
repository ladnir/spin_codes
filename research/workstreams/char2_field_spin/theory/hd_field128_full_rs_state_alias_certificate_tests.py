"""Tiny directed full-RS alias checks; no actual GF128/K20 distance gate."""

from fractions import Fraction as F
from math import comb
import unittest
from unittest.mock import patch

from flint import arb, ctx
import numpy as np

import hd_field128_full_rs_state_alias_certificate as certificate
import hd_field128_full_rs_state_alias_screen as screen
import hd_field128_wide_rs_state_certificate as wide_control


base, tail, sparse, local = certificate.base, certificate.tail, certificate.sparse, certificate.local


def log_float(value):
    return -float('inf') if value.is_zero() else float(tail.log2arb(value))


class StateAliasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z, cls.alpha = F(7, 8), F(1, 4)
        # One small q256/grid64 table. No actual full-field witness is evaluated.
        cls.bank = certificate.AliasBank(cls.z, cls.alpha, ((F(1, 512), F(1, 512)),),
                                         q=256, state=24, grid=64, max_seconds=60)

    def test_full_geometry_and_raw_gg_prefix_against_sparse_and_float(self):
        with ctx.workprec(256):
            for state in certificate.STATES:
                single, composed = certificate.directed_parts(self.z, self.alpha, q=256, state=state,
                                                               budget=base.Budget(60))
                self.assertEqual(len(single), 561)
                self.assertEqual(len(local.fixed16.geometry(128, good=True)), 65)
                self.assertEqual(len(composed), 192)
                self.assertIn((0, 64), composed)
                self.assertIn((2, 62), composed)
                self.assertEqual(max(s + 2*f for s, f in composed), 128)
                short, pair = sparse.directed_parts(self.z, self.alpha, 3, 256, state, 80, base.Budget(10))
                for key, values in short.items():
                    self.assertEqual(str(single[key]), str(values))
                for key, values in pair.items():
                    self.assertTrue(all(a.overlaps(b) for a, b in zip(composed[key], values)))
                floats = screen.AliasCoefficients(float(self.z), float(self.alpha), q=256, state=state)
                for (s, f), values in single.items():
                    np.testing.assert_allclose([log_float(x) for x in values], floats.single[:, s, f], rtol=0, atol=1e-10)
                for (s, f), values in composed.items():
                    np.testing.assert_allclose([log_float(x) for x in values], floats.composed[:, s, f], rtol=0, atol=1e-10)

    def test_alpha1_raw_gg_equals_ordered_good_square(self):
        with ctx.workprec(256):
            single, composed = certificate.directed_parts(self.z, F(1), q=256, state=24, budget=base.Budget(60))
            expected = {}
            good = {key: row for key, row in single.items() if key[0] <= 1}
            for (s, f), left in good.items():
                for (t, g), right in good.items():
                    values = expected.setdefault((s+t, f+g), [arb(0) for unused in range(4)])
                    for i, term in enumerate(base.matrix_product(left, right)):
                        values[i] += term
            for key, values in composed.items():
                self.assertTrue(all(a.overlaps(b) for a, b in zip(values, expected[key])))
            # Multiplicities are outside alpha; degree(1,0) has both orders,
            # each with exactly32 choices for the singleton slot.
            table = sparse.packet_table(self.z, 1, q=256, state=24)
            pair = [arb(0) for unused in range(4)]
            for a, b in ((table[0], table[1]), (table[1], table[0])):
                for i, term in enumerate(base.matrix_product(a, b)):
                    pair[i] += 32*term
            self.assertTrue(all(a.overlaps(b) for a, b in zip(composed[1, 0], pair)))

    def test_atlas_dominates_sparse_carryfree_and_width64_control(self):
        with ctx.workprec(256):
            direct = sparse.regional_directed(self.z, self.alpha, 3, q=256, state=24)
            wide = wide_control.regional_directed(self.z, self.alpha, 3, q=256, width=64, state=24)
            self.assertEqual(str(direct), str(wide))
            for weight in (1, 2, 3):
                regions = self.bank.regional(weight)
                for full, row in regions.items():
                    for actual, control in zip(row, direct[weight, full]):
                        self.assertTrue(actual >= control)
                        self.assertTrue(actual.is_finite() and actual > 0)
            # The alias bank has no sparse r<=64 restriction. Inspect a small
            # collection of histograms only; do not run its selected-weight gate.
            self.assertEqual(set(self.bank.regional(65)), set(range(33)))

    def test_tiny_selected_receipt_all_profiles_and_directed_prefix(self):
        row = self.bank.weight_receipt(2, F(1, 256), method='runs', max_seconds=10)
        with ctx.workprec(256):
            direct = sparse.regional_directed(self.z, self.alpha, 2, q=256, state=24)
            regions = {full: direct[2, full] for full in tail.full_counts(base.PAIRS, 2)}
            support = tail.support_sum(base.PAIRS, 2, regions, F(1, 256), method='runs', budget=base.Budget(10))
            control = self.bank.constant + 2*self.bank.slope + tail.log2arb(support)
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= control.upper())
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= arb(row['log2_upper_enclosure']).upper())
        self.assertEqual(row['regional_histograms'], 2)
        self.assertTrue(row['covers_all_supports_at_weight'])
        self.assertTrue(row['covers_all_even_odd_splits'])
        self.assertEqual(row['state_potential_ratio'], '1/256')
        self.assertEqual(row['projective_divisor'], '255')
        self.assertEqual(row['target_distance'], '2/5')
        self.assertTrue(row['directed'])
        self.assertFalse(row['whole_code_certificate'])

    def test_new_construction_metadata_and_unchanged_frozen_dispatch(self):
        row = self.bank.metadata()
        self.assertEqual(row['construction'], local.construction(24))
        self.assertNotEqual(row['construction'], base.CONSTRUCTION)
        self.assertEqual(row['state_dimension'], 24)
        self.assertEqual(row['emission_minimum'], 41)
        self.assertEqual(row['single_histogram_terms'], 561)
        self.assertEqual(row['raw_good_pair_products'], 4225)
        self.assertEqual(row['composed_pairs_per_region'], 64)
        self.assertEqual(row['witness']['grid'], 64)
        self.assertTrue(row['no_coefficient_truncation'])
        self.assertFalse(row['old_numerical_ceilings_reused'])
        # No replacement transform or support implementation was introduced.
        self.assertIs(certificate.AliasBank.regional, certificate.interval.AliasBank.regional)
        self.assertIs(certificate.AliasBank.numerator, certificate.interval.AliasBank.numerator)
        self.assertIn('hd_shifted_alias_certificate.py', row['source_sha256'])
        self.assertIn('hd_shifted_good_pair_certificate.py', row['source_sha256'])
        self.assertIn('hd_field128_full_rs_state_screen.py', row['source_sha256'])
        self.assertIn('hd_field128_full_rs_state_alias_certificate.py', row['source_sha256'])
        self.assertNotIn('hd_field128_full_rs_state_alias_screen.py', row['source_sha256'])

    def test_exact_witness_resource_and_source_pin_guards(self):
        for kwargs in (dict(state=20), dict(q=16), dict(bits=8), dict(grid=32),
                       dict(joint=1), dict(max_table_bytes=1), dict(outer_measure='unknown')):
            with self.assertRaises(ValueError):
                certificate.AliasBank(self.z, self.alpha, ((F(1, 2), F(1, 2)),), **kwargs)
        with self.assertRaises(ValueError):
            certificate.AliasBank(.875, self.alpha, ((F(1, 2), F(1, 2)),))
        with self.assertRaises(ValueError):
            certificate.AliasBank(self.z, self.alpha, ((.5, F(1, 2)),))
        for weight in (0, 8193):
            with self.assertRaises(ValueError):
                self.bank.regional(weight)
        with self.assertRaises(ValueError):
            self.bank.weight_receipt(1, .25)
        with self.assertRaises(TimeoutError):
            certificate.directed_parts(self.z, self.alpha, budget=base.Budget(1e-12))
        with patch.object(certificate, 'source_pins', return_value={}):
            with self.assertRaises(RuntimeError):
                self.bank.check_pins()


if __name__ == '__main__':
    unittest.main()
