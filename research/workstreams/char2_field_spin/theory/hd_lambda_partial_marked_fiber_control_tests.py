"""Tiny exact field/fiber tests only; never build a common-measure bank."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import log2, prod
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_partial_marked_fiber_control as partial
import hd_lambda_marked_fiber_control as old
from hd_lambda_fiber_outer import column_fiber_counts
from hd_lambda_outer import matrix
from hd_pair_scheduled_outer import column_permutations


def mul(a, b):
    v = (a if b & 1 else 0) ^ ((a << 1) if b & 2 else 0)
    return v ^ (7 if v & 4 else 0)


def null_basis(rows, width):
    rows, pivots = [list(row) for row in rows], []
    for col in range(width):
        pivot = next((i for i in range(len(pivots), len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        at = len(pivots)
        rows[at], rows[pivot] = rows[pivot], rows[at]
        inv = (0, 1, 3, 2)[rows[at][col]]
        rows[at] = [mul(inv, x) for x in rows[at]]
        for i in range(len(rows)):
            if i != at:
                scale = rows[i][col]
                rows[i] = [x ^ mul(scale, y) for x, y in zip(rows[i], rows[at])]
        pivots.append(col)
    result = []
    for free in range(width):
        if free not in pivots:
            v = [0]*width
            v[free] = 1
            for i, col in enumerate(pivots):
                v[col] = rows[i][free]
            result.append(v)
    return result


def linear_sum(coefficients, values):
    out = 0
    for a, b in zip(coefficients, values):
        out ^= mul(a, b)
    return out


class PartialMarkedFiberTests(unittest.TestCase):
    def test_actual_middlearity_fixed_sets(self):
        exceptional8 = tuple(i for i in range(16) if ((i & 3).bit_count() & 1) == 0)
        fixtures = ((tuple(range(8)), 4, 8), (exceptional8, 3, 0),
                    (tuple(range(7, 16)), 4, 1), (tuple(range(6, 16)), 4, 2),
                    (tuple(range(5, 16)), 4, 3), (tuple(range(4, 16)), 3, 0),
                    (tuple(i for i in range(16) if i not in (0, 1, 2, 4)), 4, 4))
        for source, rank, fixed_size in fixtures:
            fixed, got_rank = partial.fixed_coordinates(sum(1 << i for i in source))
            self.assertEqual((got_rank, fixed.bit_count()), (rank, fixed_size))

    def test_GF4_quiet_upper_tail_not_coefficientwise(self):
        h, column, schedule = matrix(2, q=4), 11, 'input-inverse-output-xor'
        before, after = column_permutations(schedule, column)
        exceptional8 = tuple(i for i in range(16) if ((i & 3).bit_count() & 1) == 0)
        for mapped in (exceptional8, tuple(range(7, 16))):
            original = tuple(before.index(i) for i in mapped)
            data = column_fiber_counts(original, column, q=4)
            geom = partial.column_geometry(data)
            outside = tuple(i for i in range(16) if i not in mapped)
            basis = null_basis(tuple(tuple(h[i][j] for j in mapped) for i in outside), len(mapped))
            self.assertEqual(len(basis), 5)  # Only 4^5 exact vectors per fixture.
            atoms = tuple(F(1+i % 3, 4) for i in range(16))
            true, upper = [F(0)]*(len(mapped)+1), [F(0)]*(len(mapped)+1)
            support_counts = Counter()
            for coefficients in product(range(4), repeat=len(basis)):
                y = [0]*16
                for j, label in enumerate(mapped):
                    y[label] = linear_sum(coefficients, [row[j] for row in basis])
                # H is an involution. This enumerates actual quiet outputs,
                # then recovers each original input rather than guessing it.
                x = tuple(linear_sum(row, y) for row in h)
                self.assertTrue(all(x[i] == 0 for i in outside))
                physical_mask = sum(bool(y[i]) << after[i] for i in mapped)
                support_counts[physical_mask] += 1
                for i in mapped:
                    if geom['mapped_fixed_mask'] >> i & 1:
                        self.assertEqual(x[i], y[i])
                degree = partial.quiet_surrogate_degree(physical_mask,
                    geom['physical_inside_mask'], geom['physical_fixed_mask'])
                weight = sum(bool(x[i]) for i in mapped)
                self.assertLessEqual(weight, degree)
                marked_weight = prod(atoms[i] for i in range(16) if physical_mask >> i & 1)
                true[weight] += marked_weight
                upper[degree] += marked_weight
            self.assertEqual(dict(support_counts), data['zero_fiber_counts'])
            self.assertEqual(true[0], upper[0])
            self.assertEqual(upper[0], 1)
            self.assertTrue(any(a > b for a, b in zip(true, upper)))
            for floor in range(len(mapped)+1):
                self.assertLessEqual(sum(true[floor:]), sum(upper[floor:]))
            for X in (F(1), F(2), F(5)):
                self.assertLessEqual(sum(v*X**d for d, v in enumerate(true)),
                                     sum(v*X**d for d, v in enumerate(upper)))

    def test_physical_and_original_fixed_labels(self):
        mapped, column = tuple(range(7, 16)), 9
        before, after = column_permutations(partial.SCHEDULE, column)
        original = tuple(before.index(i) for i in reversed(mapped))
        data = column_fiber_counts(original, column, q=4)
        geom = partial.column_geometry(data)
        self.assertEqual(geom['mapped_fixed_mask'], 1 << 7)
        self.assertEqual(geom['physical_fixed_mask'], 1 << after[7])
        self.assertEqual(geom['original_fixed_mask'], 1 << before.index(7))
        bad = dict(data, physical_inside_mask=data['mapped_source_mask'])
        with self.assertRaises(ValueError):
            partial.column_geometry(bad)

    def test_degree_bins_and_old_identity_agreement(self):
        source, schedule = tuple(range(8)), 'static'
        data = column_fiber_counts(source, 0, schedule=schedule, q=4)
        table = np.full((65536, 2, 2), 7.)
        result, geom = partial.quiet_degree_norms(data, table, 2)
        self.assertTrue(geom['quiet_pointwise_identity'])
        np.testing.assert_array_equal(result, old.quiet_degree_norms(data['zero_fiber_counts'], 8, table, 2))
        np.testing.assert_allclose(result[0], 3.5)
        exceptional = tuple(i for i in range(16) if ((i & 3).bit_count() & 1) == 0)
        data = column_fiber_counts(exceptional, 0, schedule='static', q=4)
        result, geom = partial.quiet_degree_norms(data, table, 2)
        self.assertFalse(geom['quiet_pointwise_identity'])
        self.assertTrue(np.all(np.isneginf(result[1:8])))
        np.testing.assert_allclose(result[8], 3.5+log2(4**5-1))

    def test_controller_is_scalar_only_and_preserves_prefix(self):
        saved = dict(mark='test-fixed-mark', source_count='123')
        table = np.zeros((1, 2, 2))
        columns = np.full((16, 9, 2, 2), -np.inf)
        columns[:, 0] = 0.
        with patch.object(partial, 'load_measure', return_value=(saved, {}, 23.)), \
             patch.object(partial, 'moment_table', return_value=table) as build, \
             patch.object(partial, 'prepare_source', return_value=(columns, [])):
            result = partial.run((tuple(range(8)),), log_grid=(0., 1.))
        self.assertEqual([c.args[1] for c in build.call_args_list], [18, 36])
        self.assertFalse(result['polynomial_tail_backend_used'])
        self.assertFalse(result['rowset_sum_included'])
        self.assertTrue(result['upper_tail_domination'])
        self.assertEqual(result['common_prefix_log2'], 23.)
        trial = result['controls'][0]['trials'][1]
        self.assertEqual(trial['chernoff_floor_charge_log2'], -9*8*18.)
        self.assertAlmostEqual(trial['fixed_A_repeated_r_bound_log2'],
                               23.+trial['word_log2_moment']-9*8*18-8*8*18*log2(255))

    def test_boundaries_and_representative_scope(self):
        self.assertEqual([len(s) for s in partial.SOURCES], [6, 6, 7, 7, 8, 9, 10, 11, 12])
        for arity in (0, 3, 13, 16):
            with self.assertRaises(ValueError):
                partial.fixed_coordinates((1 << arity)-1)
        for source in (tuple(range(3)), tuple(range(13)), (0, 0, 1, 2)):
            with self.assertRaises(ValueError):
                partial.prepare_source(source, None, None)
        with self.assertRaises(ValueError):
            partial.quiet_surrogate_degree(1 << 15, 255, 0)
        with self.assertRaises(ValueError):
            partial.quiet_surrogate_degree(0, 255, 256)
        with self.assertRaises(ValueError):
            partial.run((), max_seconds=1)
        with self.assertRaises(ValueError):
            partial.run((tuple(range(8)),), max_seconds=121)


if __name__ == '__main__':
    unittest.main()
