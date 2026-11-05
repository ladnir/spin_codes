"""Small coefficient-prefix tests; no full-field wide source-weight gate."""

from fractions import Fraction as F
from math import comb
import unittest

from flint import arb, ctx

import hd_field128_lift_certificate as base
import hd_field128_lift_sparse_bands as bands
import hd_field128_lift_sparse_wide as wide


class WideSparseTests(unittest.TestCase):
    def assert_matrix_close(self, actual, expected):
        for a, b in zip(actual, expected):
            if a.is_zero() and b.is_zero():
                continue
            self.assertTrue(abs(a-b) <= max(a, b)*arb(2)**-180)

    def test_matches_frozen_regional_at_small_cutoff(self):
        with ctx.workprec(256):
            for alpha, paired in ((F(1, 2), True), (F(1), True), (F(1, 2), False)):
                actual = wide.regional_directed(F(7, 8), alpha, 3, q=256, paired=paired, budget=base.Budget(20))
                expected = base.regional_directed(F(7, 8), alpha, 3, q=256, paired=paired, budget=base.Budget(20))
                self.assertEqual(set(actual), set(expected))
                for index in actual:
                    self.assert_matrix_close(actual[index], expected[index])

    def test_cutoff128_retains_small_prefix_with_different_marker(self):
        # Omitted local degrees exceed3 and cannot enter any coefficient tested.
        # This exercises cutoff128 arithmetic, not a128-weight distance gate.
        with ctx.workprec(256):
            single, composed = base.directed_parts(F(7, 8), F(1, 2), 3, 256, 80, base.Budget(20), True)
            actual = wide.regional_from_parts(single, composed, 128, return_maximum=3, budget=base.Budget(20))
            expected = base.regional_directed(F(7, 8), F(1, 2), 3, q=256, budget=base.Budget(20))
            for index in actual:
                self.assert_matrix_close(actual[index], expected[index])

    def test_disjoint_weight_ranges_share_one_table(self):
        actual = wide.band_receipt(((1, 1), (3, 3)), F(7, 8), F(1, 2), maximum=5, q=256, max_seconds=20)
        expected = bands.batch_receipt(1, 3, F(7, 8), F(1, 2), maximum=5, q=256, max_seconds=20)
        self.assertEqual(actual['regional_tables_built'], 1)
        self.assertEqual(actual['pending_weight_ranges'], [])
        self.assertTrue(actual['complete_requested_ranges'])
        self.assertEqual([row['active_groups'] for row in actual['weight_receipts']], [1, 3])
        for row, old in zip(actual['weight_receipts'], expected['weight_receipts'][::2]):
            self.assertEqual(row['log2_upper_rational'], old['log2_upper_rational'])
        self.assertEqual(actual['log2_upper_rational'], bands._union_ceiling(actual['weight_receipts'], 256))
        self.assertFalse(actual['whole_code_certificate'])

    def test_cutoff256_keeps_degree128_composed_term(self):
        # Synthetic diagonal polynomial (I+y^64 I)^64. No packet bank or
        # construction-size source event is evaluated in this arithmetic test.
        with ctx.workprec(256):
            identity = tuple(arb(int(i in (0, 3))) for i in range(4))
            result = wide.regional_from_parts({}, {(0, 0): identity, (0, 64): identity}, 256,
                                               return_maximum=128, budget=base.Budget(20))
            expected = base.rational(F(64, comb(base.PAIRS, 64)))
            self.assert_matrix_close(result[128, 64], (expected, arb(0), arb(0), expected))

    def test_old_seed_plan_groups_shared_witness_and_keeps_pins(self):
        plan = wide.old_witness_plan()
        self.assertEqual(plan['input_provenance']['sha256'], wide.OLD_SHA256)
        self.assertEqual(len(plan['jobs']), 3)
        self.assertEqual([row['requested_weight_ranges'] for row in plan['jobs']],
                         [[[65, 70]], [[194, 196], [210, 225]], [[226, 256]]])
        self.assertEqual([row['coefficient_maximum'] for row in plan['jobs']], [128, 256, 256])
        for row in plan['jobs']:
            self.assertEqual(F(row['z']), F(str(float(F(row['old_z'])**16))))
            self.assertEqual(row['alpha'], '1/160')
            self.assertFalse(row['old_numerical_ceilings_reused'])
        full = wide.old_witness_plan(ranges=((65, 256),))
        self.assertEqual(len(full['jobs']), 8)
        self.assertFalse(plan['directed'])

    def test_pending_ranges_and_guards(self):
        requested = ((65, 70), (194, 196), (210, 256))
        self.assertEqual(wide.pending_ranges(requested, {65, 66, 194, 195}),
                         [[67, 70], [196, 196], [210, 256]])
        result = wide.band_receipt(((1, 1), (3, 3)), F(7, 8), F(1, 2), maximum=5,
                                   q=256, max_seconds=1e-12)
        self.assertFalse(result['complete_requested_ranges'])
        self.assertEqual(result['weight_receipts'], [])
        self.assertEqual(result['pending_weight_ranges'], [[1, 1], [3, 3]])
        self.assertIsNone(result['log2_upper_rational'])
        untouched = wide.run_plan(max_seconds=1e-12)
        self.assertEqual(untouched['completed_weight_count'], 0)
        self.assertEqual(untouched['pending_weight_ranges'], list(map(list, wide.DEFAULT_RANGES)))
        self.assertFalse(untouched['complete_requested_ranges'])
        self.assertEqual(untouched['failure']['error_type'], 'TimeoutError')
        for ranges in (((0, 3),), ((1, 257),), ((1, 3), (3, 4)), ()):
            with self.assertRaises(ValueError):
                wide.checked_ranges(ranges)
        with self.assertRaises(ValueError):
            wide.regional_directed(F(7, 8), F(1, 2), 257, q=256)
        pins = wide.source_pins()
        for name in ('hd_field128_lift_sparse_wide.py', 'hd_field128_lift_sparse_bands.py',
                     'hd_field128_lift_certificate.py'):
            self.assertIn(name, pins)


if __name__ == '__main__':
    unittest.main()
