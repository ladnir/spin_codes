"""Tiny marked-kernel and ordered-moment checks; no saved-bank control."""
from fractions import Fraction as F
from itertools import product
from math import log2, prod
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_marked_fiber_control as marked
from hd_lambda_fiber_outer import _support_enumerators
from hd_lambda_fiber_outer_tests import apply4, rank4
from hd_pair_identity_kernel import mask_table_from_entries


def logs(values):
    with np.errstate(divide='ignore'):
        return np.log2(np.asarray(values, float))


class MarkedFiberTests(unittest.TestCase):
    def test_exact_quiet_weights_and_nonquiet_surrogate_upper_tails(self):
        h = tuple(tuple(2 ^ int(i == j) for j in range(4)) for i in range(4))
        source, outside = (0, 1, 2), (3,)
        def rank_cross(rmask, cmask):
            rr = tuple(i for i in range(4) if rmask >> i & 1)
            cc = tuple(i for i in range(4) if cmask >> i & 1)
            return rank4(tuple(tuple(h[i][j] for j in cc) for i in rr), len(cc))
        C, Z = _support_enumerators(source, outside, 1, rank_cross, 4)
        atoms = (F(1, 5), F(1, 3), F(2, 3), F(3, 4))
        true = [F(0)]*4
        for values in product(range(4), repeat=3):
            y = apply4(h, source, values)
            true[sum(bool(x) for x in values)] += prod(atoms[i] for i, x in enumerate(y) if x)
            if not y[3]:
                self.assertEqual(y[:3], values)
        upper = [F(0)]*4
        for mask, n in Z:
            upper[mask.bit_count()] += n*prod(atoms[i] for i in range(4) if mask >> i & 1)
        upper[3] += sum(n*prod(atoms[i] for i in range(4) if mask >> i & 1)
                        for mask, n in C if mask)*sum(
                            n*prod(atoms[i] for i in range(4) if mask >> i & 1) for mask, n in Z)
        for floor in range(4):
            self.assertLessEqual(sum(true[floor:]), sum(upper[floor:]))
        for X in (F(1), F(2), F(5)):
            self.assertLessEqual(sum(n*X**d for d, n in enumerate(true)),
                                 sum(n*X**d for d, n in enumerate(upper)))

    def test_direct_and_factor_norms_dominate_same_noncommuting_measure(self):
        atoms = np.array((0., .25, 1.))
        matrices = np.array(([[2., 1.], [0., 3.]], [[1., 0.], [2., 1.]],
                             [[1., 3.], [1., 0.]]))
        data = dict(outside_counts={0: 1, 1: 3}, zero_fiber_counts={0: 1, 2: 3})
        for groups in (1, 2, 3):
            tables = []
            for exponent in (groups, 2*groups):
                zero = matrices.sum(axis=0)
                nonzero = (matrices*atoms[:, None, None]**exponent).sum(axis=0)
                tables.append(mask_table_from_entries(logs(zero), logs(nonzero), rows=2))
            got = marked.nonquiet_norm(data, *tables, groups)
            direct = np.zeros((2, 2))
            for a, b in product(range(3), repeat=2):
                value = (3*atoms[a])*(1+3*atoms[b])
                direct += (matrices[a] @ matrices[b])*value**groups
            for key in ('selected', 'factor_2r', 'direct_union'):
                self.assertTrue(np.all(np.exp2(groups*got[key]) >= direct*(1-2e-14)))
            np.testing.assert_array_equal(got['selected'], np.minimum(got['factor_2r'], got['direct_union']))
            self.assertEqual(got['direct_terms'], 2)

    def test_quiet_norm_degrees_and_nonunit_empty_mass(self):
        table = np.full((16, 2, 2), 7.)
        result = marked.quiet_degree_norms({0: 1, 3: 3, 7: 6}, 3, table, 2)
        np.testing.assert_allclose(result[0], 3.5)
        self.assertTrue(np.all(np.isneginf(result[1])))
        np.testing.assert_allclose(result[2], 3.5+log2(3))
        np.testing.assert_allclose(result[3], 3.5+log2(6))

    def test_chronological_columns_and_chernoff_normalization(self):
        raw = np.array(([[[1., 2.], [0., 1.]], [[2., 1.], [1., 3.]]],
                        [[[2., 0.], [3., 1.]], [[1., 4.], [2., 1.]]]))
        columns = logs(raw)
        groups, X = 2, 2.
        first = (raw[0, 0]+X*raw[0, 1])**groups
        second = (raw[1, 0]+X*raw[1, 1])**groups
        expected = float(np.log2((first @ second)[0].sum()))
        self.assertAlmostEqual(marked.word_log_moment(columns, 1., groups), expected)
        result = marked.optimize_prepared(columns, 1, 17., groups=groups,
                                          log_grid=(0., 1.), refinement_rounds=0)
        trial = result['trials'][1]
        self.assertAlmostEqual(trial['chernoff_floor_charge_log2'], -18.)
        self.assertAlmostEqual(trial['normalized_fixed_A_root_log2'], (expected-18)/2-8*log2(255))
        self.assertAlmostEqual(trial['fixed_A_repeated_r_bound_log2'], 17+expected-18-16*log2(255))

    def test_all_one_mass_and_exact_empty_nonquiet(self):
        table = np.zeros((16, 2, 2))
        data = dict(outside_counts={0: 1, 8: 3}, zero_fiber_counts={0: 1, 3: 3})
        got = marked.nonquiet_norm(data, table, table, 18)
        np.testing.assert_allclose(got['selected'], log2(12))
        quiet = marked.quiet_degree_norms(data['zero_fiber_counts'], 3, table, 18)
        quiet[3] = np.logaddexp2(quiet[3], got['selected'])
        np.testing.assert_allclose(np.logaddexp2.reduce(quiet, axis=0), 4.)
        empty = marked.nonquiet_norm(dict(outside_counts={0: 1}, zero_fiber_counts={0: 1}),
                                     table, table, 18)
        self.assertTrue(np.all(np.isneginf(empty['selected'])))

    def test_mapped_rank_guard_precedes_quiet_marking(self):
        data = dict(cross_rank=0, mapped_source_mask=1)
        with patch.object(marked, 'column_fiber_counts', return_value=data):
            with self.assertRaises(marked.UnsupportedQuietFiber):
                marked.prepare_source((0,), np.zeros((2, 2, 2)), np.zeros((2, 2, 2)))

    def test_guards_and_bounded_refinement(self):
        table = np.zeros((4, 2, 2))
        with self.assertRaises(ValueError):
            marked.nonquiet_norm(dict(outside_counts={0: 1, 1: 3},
                                       zero_fiber_counts={0: 1, 1: 3}), table, table, 2)
        with self.assertRaises(ValueError):
            marked.word_log_moment(np.zeros((1, 2, 2, 2)), -.1, 2)
        with self.assertRaises(ValueError):
            marked.run(max_seconds=121)
        columns = np.zeros((2, 2, 2, 2))
        result = marked.optimize_prepared(columns, 1, 0., groups=2)
        self.assertLessEqual(result['evaluations'], 14)
        self.assertTrue(all(0 <= row['log2_input_mark'] <= 8 for row in result['trials']))
        self.assertIsNone(marked._jsonable(float('-inf')))
        with self.assertRaises(ArithmeticError):
            marked._jsonable(float('nan'))


if __name__ == '__main__':
    unittest.main()
