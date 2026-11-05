"""Tiny fixed-type enumerations; no production regional banks or gates."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_exact_arity as target
import hd_two_sum_gf16_state_envelope as state


class ExactArityTests(unittest.TestCase):
    def test_interval_norms_equal_every_fixed_type_vector(self):
        values = np.array([[2., 3., 5.], [7., 1., 4.]])
        for columns in (1, 2, 4):
            for maximum in range(1, 5):
                intervals = {groups: (groups, min(3*groups, groups+2)) for groups in range(1, maximum+1)}
                got = target.restricted_type_norms(np.log2(values), intervals, columns=columns)
                for groups, (lower, upper) in intervals.items():
                    exact = np.zeros(2)
                    for types in product(range(3), repeat=groups):
                        if lower <= sum(a+1 for a in types) <= upper:
                            exact += np.prod(values[:, types]**columns, axis=1)
                    np.testing.assert_allclose(2**(columns*got[groups]), exact, rtol=2e-14)

    def test_full_interval_agrees_with_old_unrestricted_norm(self):
        values = np.log2(np.array([[2., 3., 5.], [7., 1., 4.]]))
        intervals = {groups: (0, 3*groups) for groups in range(1, 5)}
        got = target.restricted_type_norms(values, intervals, columns=4)
        for groups in intervals:
            expected = groups/4*state.logsum(4*values, axis=1)
            np.testing.assert_allclose(got[groups], expected, atol=1e-13)

    def test_exact_restriction_dominates_every_chernoff_witness(self):
        values = np.log2(np.array([[2., 3., 5.], [7., 1., 4.]]))
        columns = 4
        for groups in range(1, 5):
            for bound in range(groups, 3*groups+1):
                exact = target.restricted_type_norms(values, {groups: (0, bound)}, columns=columns)[groups]
                for logt in (0., -1., -4., -16.):
                    tilted = columns*values+np.arange(1, 4)*logt
                    upper = (-bound*logt+groups*state.logsum(tilted, axis=1))/columns
                    self.assertTrue(np.all(exact <= upper+1e-12))

    def test_restricted_holder_minkowski_and_noncommuting_whole_word(self):
        classes = state.path_classes(3)
        path_a = np.array([sum(n*x for n, x in zip(row[2], (0., 1., 2., 0.))) for row in classes])
        typed = np.log2(np.array([[1+row[2][0], 2+row[2][1], 1+row[2][3]] for row in classes]))
        groups, columns, bound = 2, 3, 4
        norms = target.restricted_type_norms(typed, {groups: (0, bound)}, columns=columns)[groups]
        kernel = 2**target.kernel_from_restricted_norms(classes, path_a, norms)
        exact_union = 0.
        for types in product(range(3), repeat=groups):
            if sum(a+1 for a in types) > bound:
                continue
            matrix = np.zeros((2, 2))
            for index, (initial, final, _, multiplicity) in enumerate(classes):
                matrix[initial, final] += multiplicity*2**(path_a[index]+sum(typed[index, a] for a in types))
            exact_union += np.linalg.matrix_power(matrix, columns)[0].sum()
        upper = np.linalg.matrix_power(kernel, columns)[0].sum()
        self.assertGreaterEqual(upper, exact_union*(1-1e-13))
        full_norms = target.restricted_type_norms(typed, {groups: (0, 3*groups)}, columns=columns)[groups]
        np.testing.assert_allclose(target.kernel_from_restricted_norms(classes, path_a, full_norms),
            state.kernel_from_types(classes, path_a, typed, groups, columns), atol=1e-12)

    def test_zero_types_empty_interval_and_single_exact_arity(self):
        values = np.array([[0., -np.inf, 1.]])
        got = target.restricted_type_norms(values, {2: (0, 1), 3: (3, 3), 4: (12, 12)}, columns=2)
        self.assertTrue(np.isneginf(got[2][0]))
        self.assertEqual(got[3][0], 0.)
        self.assertEqual(got[4][0], 4.)
        absent = target.restricted_type_norms(np.full((2, 3), -np.inf), {2: (0, 6)})
        self.assertTrue(np.all(np.isneginf(absent[2])))

    def test_prefix_and_no_arity_charge_or_alpha_on_type_union(self):
        classes = state.path_classes(2)
        prepared = dict(classes=classes, path_a=np.zeros(len(classes)), typed=np.zeros((len(classes), 16)))
        identity = np.array([[0., -np.inf], [-np.inf, 0.]])
        with patch.object(target, 'exact_low_kernel', return_value=identity) as kernel:
            row = target.evaluate_prepared(prepared, 18, 247, F(691, 700), F(1, 5))
        self.assertEqual(kernel.call_args.args[1:3], (18, 246))
        expected = .2*(-target.LIMIT*log2(691/700)-log2(255))+log2(comb(8192, 18))
        self.assertAlmostEqual(row['failure_log2_screen'], expected)
        self.assertFalse(row['total_arity_chernoff_marker_used'])
        self.assertFalse(row['high_side_included'])
        self.assertFalse(row['directed'])
        self.assertFalse(row['continuous_weight_band_covered'])

    def test_bounded_requests_and_deadline(self):
        for intervals in ({}, {33: (0, 33)}, {True: (0, 1)}, {2: (4, 3)}, {2: (0, 7)}, {2: (0., 4)}):
            with self.assertRaises(ValueError):
                target.restricted_type_norms(np.zeros((2, 3)), intervals)
        with self.assertRaises(ValueError):
            target.restricted_type_norms(np.zeros((2, 3)), {2: (0, 6)}, columns=17)
        def stop():
            raise TimeoutError('fixture')
        with self.assertRaises(TimeoutError):
            target.restricted_type_norms(np.zeros((2, 3)), {2: (0, 6)}, check_time=stop)


if __name__ == '__main__':
    unittest.main()
