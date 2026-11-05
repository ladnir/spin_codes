"""Validation, tiny field checks, and a mocked uniform gate only."""
from collections import Counter
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from flint import arb, ctx

import hd_lambda_extra_transvection_control as target


class VariantControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inventory = json.loads((target.DIRECTORY/target.INVENTORY_NAME).read_bytes())

    def test_complete_observed_new_inventory(self):
        moments = target.load_inventory()
        self.assertEqual(len(moments), 6885)
        self.assertEqual(set(moments), {c for b in range(6) for c in combinations(range(16), b)})
        self.assertEqual(moments[()], {2: F(1), 3: F(1), 4: F(1)})
        self.assertEqual(target.dependencies()['hd_lambda_extra_transvection.py'], target.PINS['hd_lambda_extra_transvection.py'])

    def test_incomplete_old_symmetry_wrong_variant_and_bad_moment_rejected(self):
        for key, value in [('complete', False), ('construction', target.activity.binding.SCOPE['construction']),
                           ('subgroup_size', 48), ('transvection_support', [0, 1, 4, 5]),
                           ('reviewed_source_sha256', {})]:
            bad = dict(self.inventory, **{key: value})
            with self.assertRaises(ValueError):
                target.checked_inventory(bad)
        bad = deepcopy(self.inventory)
        bad['rows'][-1] = bad['rows'][0]
        with self.assertRaises(ValueError):
            target.checked_inventory(bad)
        bad = deepcopy(self.inventory)
        bad['rows'][8]['representative_to_complement_permutation'] = list(reversed(range(16)))
        with self.assertRaises(ValueError):
            target.checked_inventory(bad)
        bad = deepcopy(self.inventory)
        bad['representatives'][-1]['moments']['4'] = '1'
        with self.assertRaises(ValueError):
            target.checked_inventory(bad)
        with patch.object(Path, 'read_bytes', return_value=b'{}'):
            with self.assertRaises(ValueError):
                target.load_inventory()

    def test_old_receipt_is_checked_as_comparison_only(self):
        baseline = target.load_baseline()
        self.assertEqual(set(baseline), set(range(198, 246)))
        old = json.loads((target.DIRECTORY/target.BASELINE_NAME).read_bytes())
        old['construction'] = target.variant.CONSTRUCTION
        with self.assertRaises(ValueError):
            target.checked_baseline(old)
        old = json.loads((target.DIRECTORY/target.BASELINE_NAME).read_bytes())
        old['selected_totals'][0]['failure_log2_upper'] = '-999999'
        with self.assertRaises(ValueError):
            target.checked_baseline(old)

    def test_all120_pair_laws_match_exact_GF4_inputs(self):
        for use_variant in (False, True):
            h = target.variant.matrix(2, 4) if use_variant else target.variant.original.matrix(2, 4)
            for source in combinations(range(16), 2):
                actual = Counter()
                for weight, mask, count in target.pair_branches(source, use_variant=use_variant, lambda_value=2, q=4):
                    actual[weight, mask] += count
                actual = +actual
                expected = Counter()
                for x, y in product(range(4), repeat=2):
                    mask = sum(bool(target.variant.original._multiply(row[source[0]], x, 4)
                                    ^ target.variant.original._multiply(row[source[1]], y, 4)) << i
                               for i, row in enumerate(h))
                    expected[int(bool(x))+int(bool(y)), mask] += 1
                self.assertEqual(actual, expected)
                self.assertEqual(sum(actual.values()), 16)

    def test_pair_profiles_raw_mass_no_collapse_and_order(self):
        profile = target.two_source_profile()
        self.assertEqual(profile['pair_count'], 120)
        self.assertEqual(profile['collapsed_pairs'], {'old': [], 'new': []})
        self.assertEqual(profile['minimum_output_weight'], {'old': 6, 'new': 6})
        self.assertEqual(profile['exceptional_ratio_weight_histograms'], {
            'old': {'6': 28, '7': 7, '8': 64, '9': 8, '10': 28, '11': 98, '13': 112, '15': 15},
            'new': {'6': 26, '7': 1, '8': 60, '9': 13, '10': 30, '11': 114, '12': 4, '13': 101, '15': 11}})
        self.assertEqual(profile['first_worst_exceptional_branch']['old']['source'], [1, 2])
        self.assertEqual(profile['first_worst_exceptional_branch']['new']['source'], [0, 7])
        self.assertTrue(all(row['exceptional_ratio_index'] == 0
                            for row in profile['first_worst_exceptional_branch'].values()))
        for row in profile['pairs']:
            for label in ('old', 'new'):
                branches = row[label+'_branches']
                self.assertEqual(tuple(count for _, _, count in branches), (1, 255, 255, 255, 255, 255, 255*252))
                self.assertEqual(sum(count for _, _, count in branches), 256**2)
                self.assertEqual(branches[1], [1, 65535, 255])
                self.assertEqual(branches[2], [1, 65535, 255])
        forward = target.pair_branches((0, 5))
        backward = target.pair_branches((5, 0))
        self.assertEqual(forward[3], backward[3])
        self.assertEqual(forward[4], backward[5])
        self.assertFalse(profile['old_lowarity_certificate_transferred'])

    def test_mock_variant_gate_one_fresh_bank_and_totals_prefix(self):
        def fake_weights(moments, p):
            # Preserving common p and six arities is all this mock tests.
            return [arb(a+p) for a in target.ARITIES], []
        with ctx.workprec(256):
            terms = [(e, o, n, arb(n)*arb(2)**-20000)
                     for (e, o), n in target.A.profile_counts(18).items()]
            uniform = target.A.log2(arb(sum(row[2] for row in terms)))-20000
            p, total = 3, 198
            theta = F(3, 10)
            expected = target.A.uniform.rational(theta)*(198*(128-8*target.A.log2(arb(255)))
                +18*target.A.log2(arb(11+p))-target.A.log2(arb(255)))
            expected += uniform-target.A.uniform.rational(F(1, 5))*734003*target.A.log2(target.A.uniform.rational(target.MARK))
            expected = str(target.A.uniform.rational_ceiling(expected))
        with patch.object(target.activity, 'scheduled_weights', side_effect=fake_weights), \
                patch.object(target.A.uniform, 'regional_directed', return_value={}) as bank, \
                patch.object(target.A, 'source_profile_terms', return_value=terms):
            receipt = target.receipt()
        bank.assert_called_once_with(F(691, 700), F(1, 5), 18, 60)
        self.assertEqual(receipt['construction'], target.variant.CONSTRUCTION)
        self.assertNotEqual(receipt['construction'], receipt['old_comparison']['construction'])
        self.assertEqual(receipt['controls'][1]['totals'][0]['failure_log2_upper'], expected)
        self.assertEqual(len(receipt['selected_totals']), 48)
        self.assertEqual(receipt['profile_count'], 99)
        self.assertEqual(receipt['actual_stabilizer_size'], 4)
        self.assertFalse(receipt['old_lowarity_certificates_transferred'])
        self.assertFalse(receipt['native_performance_claim'])
        self.assertFalse(receipt['whole_code_certificate'])
        for row in receipt['selected_totals']:
            self.assertEqual(F(row['old_minus_new_bound_bits']),
                             F(row['old_failure_log2_upper'])-F(row['failure_log2_upper']))

    def test_guards(self):
        for digest in ('', '0'*63, 'G'*64, '0'*64):
            with self.assertRaises(ValueError):
                target.load_inventory(expected_sha256=digest)
        for precision in (128, True):
            with self.assertRaises(ValueError):
                target.receipt(precision=precision)
        for budget in (0, 601, float('inf'), True):
            with self.assertRaises(ValueError):
                target.receipt(max_seconds=budget)
        for source in ((), (1,), (1, 1), (16, 0), (True, 2)):
            with self.assertRaises(ValueError):
                target.pair_branches(source)


if __name__ == '__main__':
    unittest.main()
