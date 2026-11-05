"""Exact labelled activity controls; no production regional bank or witness search."""
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_group_seen as target
import hd_two_sum_gf16_state_emptycut as cut


def multiply(a, b):
    return tuple(a[2*i]*b[j]+a[2*i+1]*b[2+j] for i in range(2) for j in range(2))


def exact_activity_words(matrices, columns, predicate):
    groups = len(matrices)-1
    total = 0
    for word in product(range(1 << groups), repeat=columns):
        if not predicate(word):
            continue
        value = (1, 0, 0, 1)
        for mask in word:
            value = multiply(value, matrices[mask.bit_count()])
        total += value[0]+value[1]
    return total


def union_masks(word):
    union = 0
    for mask in word:
        union |= mask
    return union


class GroupSeenTests(unittest.TestCase):
    def test_seen_state_dp_matches_labelled_noncommuting_words(self):
        matrices = ((1, 2, 3, 4), (2, 1, 1, 3), (3, 2, 4, 1), (2, 3, 1, 1))
        self.assertNotEqual(multiply(matrices[0], matrices[1]), multiply(matrices[1], matrices[0]))
        for groups in (1, 2, 3):
            values = np.log2(np.array(matrices[:groups+1], dtype=float).reshape(-1, 2, 2))
            for columns in (1, 2, 3):
                exact = exact_activity_words(matrices[:groups+1], columns,
                    lambda word: union_masks(word) == (1 << groups)-1)
                self.assertAlmostEqual(target.seen_all_log(values, columns=columns), log2(exact), places=12)

    def test_positive_polynomial_floor_matches_labelled_words(self):
        matrices = ((1, 2, 3, 4), (2, 1, 1, 3), (3, 2, 4, 1))
        values = np.log2(np.array(matrices, dtype=float).reshape(-1, 2, 2))
        for columns in (1, 2, 3):
            for minimum in range(columns+1):
                exact = exact_activity_words(matrices, columns,
                    lambda word: sum(mask.bit_count() for mask in word) >= minimum*2)
                self.assertAlmostEqual(target.aggregate_floor_log(values, columns=columns, minimum=minimum),
                                       log2(exact), places=12)
        # The two relaxations are genuinely different, not jointly asserted.
        self.assertNotEqual(target.seen_all_log(values, columns=3),
                            target.aggregate_floor_log(values, columns=3, minimum=1))

    def test_holder_keeps_types_fixed_before_activity_union(self):
        # One labelled local path per endpoint. Type choice remains fixed
        # through two columns in the independently enumerated lower value.
        classes = tuple((i, j, (0, 0, 0, 0), 1) for i in range(2) for j in range(2))
        raw = np.array([[[1, 2], [3, 1]], [[2, 1], [1, 4]],
                        [[1, 3], [2, 1]], [[3, 2], [1, 2]]], dtype=float)
        # raw[path, activity, fixed_type]
        branches = dict(classes=classes, path_a=np.zeros(4),
                        log_u=target.state.logsum(2*np.log2(raw), axis=2)/2)
        moments = target.column_moments(branches, 2)
        exact_seen, exact_floor = 0, 0
        for types in product(range(2), repeat=2):
            for word in product(range(4), repeat=2):
                value = (1, 0, 0, 1)
                for mask in word:
                    matrix = tuple(int(np.prod([raw[p, mask >> g & 1, types[g]] for g in range(2)]))
                                   for p in range(4))
                    value = multiply(value, matrix)
                contribution = value[0]+value[1]
                if union_masks(word) == 3:
                    exact_seen += contribution
                if sum(mask.bit_count() for mask in word) >= 2:
                    exact_floor += contribution
        self.assertGreaterEqual(target.seen_all_log(moments, columns=2)+1e-12, log2(exact_seen))
        self.assertGreaterEqual(target.aggregate_floor_log(moments, columns=2, minimum=1)+1e-12,
                                log2(exact_floor))

    def test_branch_preparation_counts_each_labelled_mask_once(self):
        curves = np.array([[0., -2., -3.], [-np.inf, -1., -2.],
                           [-np.inf, -3., -4.], [-4., -1., -2.]])
        slopes = np.array([-1., -1., -1., -1.])
        weights = {1: np.log2([2., 3., 5.]), 2: np.log2([7., 11., 13.])}
        prepared = cut.prepared_envelope(curves, slopes, weights, regions=2)
        branches = target.prepare_branches(prepared, weights, columns=2)
        for index in range(len(prepared['classes'])):
            for a in range(2):
                expected = target.state.logsum(prepared['elementary'][index]+weights[a+1])
                actual = target.state.logsum(branches['type_branch_logs'][index, :, a])
                self.assertAlmostEqual(actual, expected, places=13)
        # Branch-wise Minkowski may lose something before floor exclusion.
        separated = target.state.logsum(branches['log_u'], axis=1)
        old = target.state.logsum(2*prepared['typed'], axis=1)/2
        self.assertTrue(np.all(separated+1e-13 >= old))

    def test_path_multiplicity_stays_outside_holder(self):
        branches = dict(classes=((0, 0, (0, 0, 0, 0), 7),),
                        path_a=np.array([log2(3)]), log_u=np.log2([[2., 5.]]))
        moments = target.column_moments(branches, 3)
        for active in range(4):
            self.assertAlmostEqual(moments[active, 0, 0], log2(21*2**(3-active)*5**active))
            self.assertTrue(np.isneginf(moments[active, 0, 1]))

    def test_exact_zeros_do_not_create_nan(self):
        branches = dict(classes=((0, 0, (0, 0, 0, 0), 1),), path_a=np.array([0.]),
                        log_u=np.array([[0., -np.inf]]))
        moments = target.column_moments(branches, 2)
        self.assertFalse(np.any(np.isnan(moments)))
        self.assertTrue(np.isneginf(target.seen_all_log(moments, columns=2)))
        self.assertTrue(np.isneginf(target.aggregate_floor_log(moments, columns=2, minimum=1)))
        self.assertEqual(target.aggregate_floor_log(moments, columns=2, minimum=0), 0.)

    def test_evaluation_scope_prefix_and_minimum_not_joint_constraint(self):
        prepared = dict(regions=16)
        with patch.object(target, 'prepare_branches', return_value={}), \
                patch.object(target, 'column_moments', return_value=np.zeros((3, 2, 2))), \
                patch.object(target, 'seen_all_log', return_value=3.), \
                patch.object(target, 'aggregate_floor_log', return_value=2.):
            first = target.evaluate(prepared, {a: None for a in range(1, 17)}, 2, .99, .2, source_count=3)
            second = target.evaluate(prepared, {a: None for a in range(1, 17)}, 2, .99, .2, source_count=21)
        self.assertAlmostEqual(second['failure_log2_screen']-first['failure_log2_screen'], log2(7))
        self.assertEqual(first['failure_log2_screen'], first['aggregate_9r_failure_log2_screen'])
        self.assertFalse(first['two_constraints_enforced_jointly'])
        self.assertFalse(first['individual_nine_column_floors_enforced'])
        self.assertFalse(first['whole_code_certificate'])
        self.assertEqual(first['regional_banks_rebuilt'], 0)

    def test_guards_and_cooperative_checks(self):
        for values, columns in ((np.zeros((1, 2, 2)), 2), (np.zeros((3, 2, 2)), True),
                                (np.full((3, 2, 2), np.nan), 2)):
            with self.assertRaises(ValueError):
                target.seen_all_log(values, columns=columns)
        for minimum in (True, -1, 3):
            with self.assertRaises(ValueError):
                target.aggregate_floor_log(np.zeros((3, 2, 2)), columns=2, minimum=minimum)
        def stop():
            raise TimeoutError('stop')
        for method in (target.seen_all_log, target.aggregate_floor_log):
            with self.assertRaises(TimeoutError):
                method(np.zeros((3, 2, 2)), check_time=stop)


if __name__ == '__main__':
    unittest.main()
