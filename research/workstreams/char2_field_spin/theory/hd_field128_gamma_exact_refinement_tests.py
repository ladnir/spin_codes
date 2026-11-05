"""Tiny support oracles, bounded selection, and FLOAT-only transfer checks."""

from fractions import Fraction as F
from math import log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_field128_gamma_exact_refinement as ref
from hd_shifted_pairs_count import exact_profiles
from hd_field128_gamma_near_face_certificate_tests import fixture, multiply, power


def logs(rows):
    return np.array([[log2(x) for x in row] for row in rows.values()]).reshape(-1, 2, 2)


class RefinementTests(unittest.TestCase):
    def test_all_tiny_profiles_exact_and_none_omitted(self):
        for pairs in range(1, 7):
            oracle = exact_profiles(pairs)
            for weight in range(2*pairs+1):
                counts = ref.ExactRefinedCounts(pairs, weight, exact_singletons=2*pairs)
                total = 0
                for a in counts.full:
                    for b in counts.full:
                        value = counts.count(a, b)
                        self.assertEqual(value, oracle.get((weight, a, b), 0))
                        total += value
                self.assertEqual(total, sum(v for (r, a, b), v in oracle.items() if r == weight))
                full = np.array(list(counts.full))
                actual = counts.log_counts(full[:, None], full[None, :])
                expected = np.array([[log2(counts.count(a, b)) if counts.count(a, b) else -np.inf
                                      for b in counts.full] for a in counts.full])
                np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)
                self.assertEqual(float(counts.log_counts(int(full[0]), int(full[0]))), expected[0, 0])

    def test_explicit_priority_cap_symmetry_and_auto_order(self):
        counts = ref.ExactRefinedCounts(9, 8, selected_profiles=[(1, 3), (3, 1)], max_exact_profiles=3)
        self.assertEqual(counts.selected_profiles[0], (3, 1))
        self.assertEqual(len(counts.selected_profiles), 3)
        eligible = []
        for a in counts.full:
            for b in range(a+1):
                if a*2 != 8 and (a, b) != (3, 1) and counts.base.count(a, b):
                    eligible.append((8-2*a, -counts.base.count(a, b), a, b))
        self.assertEqual(counts.selected_profiles[1:], tuple((a, b) for s, c, a, b in sorted(eligible)[:2]))
        for a in counts.full:
            for b in counts.full:
                key = max(a, b), min(a, b)
                expected = ref.exact.joint_count(9, 8, a, b) if key in counts.exact_counts else counts.base.count(a, b)
                self.assertEqual(counts.count(a, b), expected)
        with self.assertRaises(TypeError):
            counts.exact_counts[(0, 0)] = 0
        self.assertFalse(counts._keys.flags.writeable)
        self.assertTrue(counts.diagnostics()['all_profiles_retained'])

    def test_disabled_replays_frozen_and_reuse_never_recounts(self):
        rows = fixture(5, 5)
        counts = ref.ExactRefinedCounts(5, 5, max_exact_profiles=0, max_singletons=2)
        expected = ref.near.support_sum(5, 5, logs(rows), repetitions=8, max_singletons=2)
        with patch.object(ref.exact, 'joint_count', side_effect=AssertionError('no recount')):
            actual = ref.support_sum(5, 5, logs(rows), repetitions=8, counts=counts)
            repeated = ref.support_sum(5, 5, logs(rows), repetitions=8, counts=counts)
        self.assertEqual(actual['log2_support_sum_screen'], expected['log2_support_sum_screen'])
        self.assertEqual(actual, repeated)
        self.assertEqual(actual['profile_pairs_evaluated'], 9)

    def test_float_transfer_matches_independent_positive_sum(self):
        for pairs, weight in ((3, 3), (4, 4), (5, 6)):
            rows = fixture(pairs, weight)
            counts = ref.ExactRefinedCounts(pairs, weight, exact_singletons=2*pairs)
            oracle = F(0)
            for a, left in rows.items():
                for b, right in rows.items():
                    matrix = power(multiply(left, right), 4)
                    oracle += counts.count(a, b)*(matrix[0]+matrix[1])
            actual = ref.support_sum(pairs, weight, logs(rows), repetitions=8, counts=counts)
            self.assertAlmostEqual(actual['log2_support_sum_screen'], log2(oracle), places=11)

    def test_top_contribution_selector_uses_both_orders_and_bounded_heap(self):
        pairs, weight = 6, 6
        rows = fixture(pairs, weight)
        baseline = ref.near.NearFaceCounts(pairs, weight, max_singletons=2)
        expected = []
        for a, left in rows.items():
            if not 0 < weight-2*a <= 6:
                continue
            for b in range(a+1):
                count = baseline.count(a, b)
                if not count:
                    continue
                forward = power(multiply(left, rows[b]), 4)
                value = forward[0]+forward[1]
                if a != b:
                    reverse = power(multiply(rows[b], left), 4)
                    value += reverse[0]+reverse[1]
                expected.append((log2(count*value), -a, -b))
        selected = ref.select_profiles(pairs, weight, logs(rows), max_profiles=3, max_singletons=2,
                                       eligible_singletons=6, repetitions=8, chunk=2,
                                       region_source={'test_only': True})
        ranked = sorted(expected, reverse=True)[:3]
        self.assertEqual(selected['selected_profiles'], [[-a, -b] for score, a, b in ranked])
        np.testing.assert_allclose(selected['selection_scores_FLOAT'], [score for score, a, b in ranked], atol=1e-12)
        self.assertEqual(selected['baseline_profile_pairs_evaluated'], 16)
        self.assertFalse(selected['directed'])
        self.assertFalse(selected['selected_weight_certificate'])
        reused = ref.ExactRefinedCounts(pairs, weight, selected_profiles=selected['selected_profiles'], exact_singletons=0)
        self.assertEqual(list(reused.selected_profiles), [tuple(p) for p in selected['selected_profiles']])

    def test_guards_timeout_and_note_pin(self):
        for options in ({'max_exact_profiles': -1}, {'max_exact_profiles': 4097}, {'exact_singletons': True},
                        {'selected_profiles': [(3, 0)]}, {'selected_profiles': [(1, 0)], 'max_exact_profiles': 0}):
            with self.assertRaises(ValueError):
                ref.ExactRefinedCounts(4, 4, **options)
        counts = ref.ExactRefinedCounts(4, 4)
        with self.assertRaises(ValueError):
            ref.support_sum(4, 4, logs(fixture(4, 4)), counts=counts, exact_singletons=0)
        class Expired:
            def check(self):
                raise TimeoutError('test only')
        with self.assertRaises(TimeoutError):
            ref.ExactRefinedCounts(4, 4, budget=Expired())
        with self.assertRaises(TimeoutError):
            ref.select_profiles(4, 4, logs(fixture(4, 4)), budget=Expired())
        pins = ref.source_pins()
        self.assertIn('../HD_FIELD128_EXACT_PROFILE_COUNTS.md', pins)
        ref.check_source_pins(pins)
        pins['../HD_FIELD128_EXACT_PROFILE_COUNTS.md'] = 'changed'
        with self.assertRaises(RuntimeError):
            ref.check_source_pins(pins)


if __name__ == '__main__':
    unittest.main()
