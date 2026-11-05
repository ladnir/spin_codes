"""Tiny exact labelled controls for the g16 full-arity seen-count reduction."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

from hd_two_sum_gf16_full_dense import column_vectors, seen_transitions, run, Z, LOGP
from hd_two_sum_affine_dense_frozen import Z as OLD_Z, LOGP as OLD_LOGP
from hd_xor_dense_two_feature_tests import exact_tables, logs, action
from hd_xor_joint_columns import categories


def labelled_reference(tables, potential, single, double, dense_groups, q):
    """Enumerate every dense mask and every labelled sparse hole, with Fractions."""
    length = len(tables)
    result = [[F(0), F(0)] for _ in range(dense_groups+1)]
    for masks in product(range(1 << length), repeat=dense_groups):
        dense_mass = (q-1)**sum(mask.bit_count() for mask in masks)
        e = sum(mask != 0 for mask in masks)
        for holes in product(range(-1, length), repeat=double):
            mass = F(dense_mass*(q-1-length)**holes.count(-1), (q-1)**double)
            occupancy = [single+double-holes.count(i)+sum((mask >> i) & 1 for mask in masks)
                         for i in range(length)]
            vector = action(tables, occupancy, potential)
            for state in range(2):
                result[e][state] += mass*vector[state]
    return result


class Gf16FullDenseTests(unittest.TestCase):
    def test_seen_count_transitions_against_all_labelled_current_subsets(self):
        for groups in (1, 2):
            actual = {}
            for s, next_s, occupied, mass in seen_transitions(groups, 5):
                key = s, next_s, occupied
                actual[key] = actual.get(key, 0)+mass
            expected = {}
            for seen in range(groups+1):
                initial = (1 << seen)-1
                for current in range(1 << groups):
                    key = seen, (initial | current).bit_count(), current.bit_count()
                    expected[key] = expected.get(key, 0)+4**current.bit_count()
            self.assertEqual(actual, expected)

    def test_all_tiny_categories_match_exact_labelled_noncommuting_products(self):
        tables = exact_tables(length=3, maximum=4)
        potential = F(2, 3)
        for dense_groups in (1, 2):
            actual = column_vectors(logs(tables), np.log2(float(potential)), 2, dense_groups, q=7)
            for index, (single, double) in enumerate(categories(2)):
                exact = labelled_reference(tables, potential, single, double, dense_groups, 7)
                np.testing.assert_allclose(np.exp2(actual[index]), np.asarray(exact, dtype=float), rtol=4e-14)

    def test_nonempty_group_marker_and_raw_total_mass(self):
        length, q = 4, 7
        with np.errstate(divide='ignore'):
            tables = np.log2(np.tile(np.eye(2), (length, 5, 1, 1)))
        for dense_groups in (1, 2):
            actual = column_vectors(tables, -2., 2, dense_groups, q=q)
            expected = np.array([[comb(dense_groups, e)*(q**length-1)**e,
                                  comb(dense_groups, e)*(q**length-1)**e/4]
                                 for e in range(dense_groups+1)])
            for row in actual:
                np.testing.assert_allclose(np.exp2(row), expected, rtol=4e-14)
            self.assertEqual(sum(int(x) for x in expected[:, 0]), q**(length*dense_groups))

    def test_regional_order_is_not_silently_commuted(self):
        tables = logs(exact_tables(length=3, maximum=3))
        forward = column_vectors(tables, -1., 1, 2, q=7)
        backward = column_vectors(tables[::-1], -1., 1, 2, q=7)
        self.assertGreater(float(np.max(abs(forward-backward))), .01)

    def test_zero_matrices_remain_zero_without_nan(self):
        actual = column_vectors(np.full((2, 4, 2, 2), -np.inf), -1., 1, 2, q=7)
        self.assertTrue(np.all(np.isneginf(actual)))

    def test_scope_literals_budget_and_no_implicit_extra_bank(self):
        self.assertEqual((Z, LOGP), (OLD_Z, OLD_LOGP))
        with patch('hd_two_sum_gf16_full_dense.pair_regions_screen') as bank:
            result = run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertEqual(result['requested_dense_groups'], [1, 2])
        self.assertEqual(result['source_matching_counts'], [8, 0])
        self.assertFalse(result['directed'])
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['complete_selected_gates'])
        self.assertNotIn('failure_log2_screen', result)
        for invalid in (True, 0, 121, float('inf')):
            with self.assertRaises(ValueError):
                run(max_seconds=invalid)

    def test_both_gates_use_one_regional_bank_and_correct_floors(self):
        vectors = lambda regions, logp, sparse, dense, **kwargs: np.zeros((len(categories(sparse)), dense+1, 2))
        with patch('hd_two_sum_gf16_full_dense.pair_regions_screen', return_value={}) as bank, \
                patch('hd_two_sum_gf16_full_dense.ordered_regions', return_value=None), \
                patch('hd_two_sum_gf16_full_dense.column_vectors', side_effect=vectors), \
                patch('hd_two_sum_gf16_full_dense.positive_bivariate_tail', return_value=(0., {})) as tail:
            result = run()
        bank.assert_called_once_with(float(Z), 17)
        self.assertEqual([call.args[1:] for call in tail.call_args_list], [(16, 288, 9), (16, 270, 18)])
        self.assertEqual([row['sparse_normalizer_exponent'] for row in result['entries']], [256, 240])
        self.assertEqual([row['dense_per_column_normalizer_exponent'] for row in result['entries']], [8, 16])
        self.assertEqual([row['source_colorings'] for row in result['entries']], [17, 136])
        self.assertTrue(result['complete_selected_gates'])
        self.assertFalse(result['complete_source_coverage'])

    def test_invalid_tables_and_dense_counts_rejected(self):
        tables = logs(exact_tables(length=3, maximum=3))
        for invalid in (True, 0, 3):
            with self.assertRaises(ValueError):
                column_vectors(tables, -1., 1, invalid, q=7)
        for value in (float('inf'), float('nan')):
            changed = tables.copy()
            changed[0, 0, 0, 0] = value
            with self.assertRaises(ValueError):
                column_vectors(changed, -1., 1, 2, q=7)


if __name__ == '__main__':
    unittest.main()
