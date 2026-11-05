"""Tiny exact-prefix and synthetic-atlas checks; no actual K20 distance gate."""

from fractions import Fraction as F
import unittest
from unittest.mock import patch

from flint import arb, ctx
import numpy as np

import hd_field128_full_rs_gamma_alias_certificate as certificate


base, tail, sparse, local = certificate.base, certificate.tail, certificate.sparse, certificate.local


def log_float(value):
    return -float('inf') if value.is_zero() else float(tail.log2arb(value))


class GammaAliasTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.z, cls.alpha = F(7, 8), F(1, 4)
        with ctx.workprec(256):
            # Real gamma caps through occupancy3. Higher rows are synthetic
            # positive matrices, solely to exercise complete histogram/DFT code.
            prefix = sparse.packet_table(cls.z, 3, state=24)
            cls.table = {j: prefix[j] if j <= 3 else tuple(base.rational(base.round_up(value, 80))
                         for value in (F(1, j+2), F(2, j+3), F(1, j+4), F(3, j+5))) for j in range(65)}
            with patch.object(sparse, 'packet_table', return_value=cls.table):
                cls.single, cls.composed = certificate.directed_parts(cls.z, cls.alpha)
                cls.bank = certificate.AliasBank(cls.z, cls.alpha, ((F(1, 512), F(1, 512)),),
                                                 state=24, grid=64, max_seconds=30)

    def test_complete_geometry_and_exact_raw_gg_prefix(self):
        self.assertEqual(len(self.single), 561)
        self.assertEqual(len(local.frozen.fixed16.geometry(128, good=True)), 65)
        self.assertEqual(len(self.composed), 192)
        self.assertIn((0, 64), self.composed)
        self.assertIn((2, 62), self.composed)
        self.assertEqual(max(s+2*f for s, f in self.composed), 128)
        with ctx.workprec(256):
            short, pair = sparse.directed_parts(self.z, self.alpha, 3, 24, 80, base.Budget(10))
            for key, values in short.items():
                self.assertEqual(str(self.single[key]), str(values))
            for key, values in pair.items():
                self.assertTrue(all(a.overlaps(b) for a, b in zip(self.composed[key], values)))
            # The largest paired occupancy must not be dropped by sparse bounds.
            raw = base.matrix_product(self.table[64], self.table[64])
            expected = tuple(value**base.rational(self.alpha) for value in raw)
            self.assertTrue(all(a.overlaps(b) for a, b in zip(self.composed[0, 64], expected)))

    def test_real_gamma_prefix_matches_positive_log_provider_all_states(self):
        with ctx.workprec(256):
            for state in certificate.STATES:
                table = sparse.packet_table(F(1, 2), 3, state=state)
                logs = local.local_logs(.5, state=state)
                for j, row in table.items():
                    np.testing.assert_allclose([log_float(x) for x in row], logs[j].ravel(),
                                               rtol=0, atol=1e-10)
                self.assertTrue(table[1][2] > 0)

    def test_alpha1_good_products_keep_order_and_counts(self):
        with ctx.workprec(256), patch.object(sparse, 'packet_table', return_value=self.table):
            single, composed = certificate.directed_parts(self.z, F(1))
            expected = {}
            good = {key: row for key, row in single.items() if key[0] <= 1}
            for (s, f), left in good.items():
                for (t, g), right in good.items():
                    values = expected.setdefault((s+t, f+g), [arb(0) for unused in range(4)])
                    for i, term in enumerate(base.matrix_product(left, right)):
                        values[i] += term
            for key, values in composed.items():
                self.assertTrue(all(a.overlaps(b) for a, b in zip(values, expected[key])))
            ordered = [arb(0) for unused in range(4)]
            for left, right in ((self.table[0], self.table[1]), (self.table[1], self.table[0])):
                for i, term in enumerate(base.matrix_product(left, right)):
                    ordered[i] += 32*term
            self.assertTrue(all(a.overlaps(b) for a, b in zip(composed[1, 0], ordered)))

    def test_synthetic_atlas_dominates_real_low_degree_carryfree_prefix(self):
        with ctx.workprec(256):
            direct = sparse.regional_directed(self.z, self.alpha, 3, state=24)
            for weight in (1, 2, 3):
                for full, row in self.bank.regional(weight).items():
                    for actual, control in zip(row, direct[weight, full]):
                        self.assertTrue(actual >= control)
                        self.assertTrue(actual.is_finite() and actual > 0)
            # Histogram access above the sparse verifier's cutoff, not a gate.
            self.assertEqual(set(self.bank.regional(65)), set(range(33)))

    def test_mock_receipt_prefix_metadata_and_frozen_dispatch(self):
        regions = {f: (arb(1), arb(2), arb(3), arb(4)) for f in (0, 1)}
        with patch.object(self.bank, 'regional', return_value=regions):
            row = self.bank.weight_receipt(2, F(1, 256), method='runs', max_seconds=10)
        with ctx.workprec(256):
            support = tail.support_sum(base.PAIRS, 2, regions, F(1, 256), method='runs', budget=base.Budget(10))
            control = self.bank.constant+2*self.bank.slope+tail.log2arb(support)
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= control.upper())
            self.assertTrue(base.rational(F(row['log2_upper_rational'])) >= arb(row['log2_upper_enclosure']).upper())
        self.assertEqual(row['construction'], local.construction(24))
        self.assertNotEqual(row['construction'], base.CONSTRUCTION)
        self.assertEqual(row['feedback'], 'gamma A transpose')
        self.assertEqual(row['gamma_extension_degree'], 16)
        self.assertEqual(row['gamma_base_field_order'], 256)
        self.assertEqual(row['state_dimension'], 24)
        self.assertEqual(row['projective_divisor'], str((1 << 128)-1))
        self.assertEqual(row['regional_histograms'], 2)
        self.assertTrue(row['covers_all_supports_at_weight'])
        self.assertTrue(row['covers_all_even_odd_splits'])
        self.assertTrue(row['no_coefficient_truncation'])
        self.assertFalse(row['old_numerical_ceilings_reused'])
        self.assertTrue(row['directed'])
        self.assertFalse(row['whole_code_certificate'])
        self.assertIs(certificate.AliasBank.regional, certificate.interval.AliasBank.regional)
        self.assertIs(certificate.AliasBank.numerator, certificate.interval.AliasBank.numerator)
        for name in ('hd_field128_full_rs_gamma_alias_certificate.py',
                     'hd_field128_full_rs_gamma_certificate.py', 'hd_field128_full_rs_gamma_screen.py',
                     'hd_shifted_alias_certificate.py', 'hd_shifted_good_pair_certificate.py'):
            self.assertIn(name, row['source_sha256'])

    def test_exact_scope_resource_timeout_and_source_guards(self):
        for kwargs in (dict(state=32), dict(state=20), dict(q=256), dict(bits=8), dict(grid=32),
                       dict(precision=64), dict(max_table_bytes=1), dict(outer_measure='unknown')):
            with self.assertRaises(ValueError):
                certificate.AliasBank(self.z, self.alpha, ((F(1, 2), F(1, 2)),), **kwargs)
        for z, markers in ((.875, ((F(1, 2), F(1, 2)),)),
                           (self.z, ((.5, F(1, 2)),)), (self.z, ((F(0), F(1)),))):
            with self.assertRaises(ValueError):
                certificate.AliasBank(z, self.alpha, markers)
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
