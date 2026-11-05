"""Tiny exact algebra, toy tails, and mocked controllers; no production gate."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import log2
import random
import unittest
from unittest.mock import patch

import numpy as np

import hd_lambda_double_fiber as fiber
import hd_lambda_double_middle_control as target
from hd_column_deficit_tail import column_deficit_tail


MUL = tuple(tuple(fiber.outer.core.original._multiply(i, j, 4) for j in range(4)) for i in range(4))
INV = (0, 1, 3, 2)


def rref(matrix, width):
    rows = [list(row) for row in matrix]
    pivots = []
    for column in range(width):
        pivot = next((i for i in range(len(pivots), len(rows)) if rows[i][column]), None)
        if pivot is None: continue
        index = len(pivots)
        rows[index], rows[pivot] = rows[pivot], rows[index]
        factor = INV[rows[index][column]]
        rows[index] = [MUL[factor][x] for x in rows[index]]
        for i in range(len(rows)):
            if i != index and rows[i][column]:
                factor = rows[i][column]
                rows[i] = [a ^ MUL[factor][b] for a, b in zip(rows[i], rows[index])]
        pivots.append(column)
    return rows, pivots


def kernel(matrix, width):
    rows, pivots = rref(matrix, width)
    basis = []
    for column in range(width):
        if column in pivots: continue
        vector = [0]*width
        vector[column] = 1
        for row, pivot in zip(rows, pivots):
            vector[pivot] = row[column]
        basis.append(tuple(vector))
    return basis


def apply(matrix, x):
    result = []
    for row in matrix:
        value = 0
        for a, b in zip(row, x): value ^= MUL[a][b]
        result.append(value)
    return tuple(result)


def matrix_polynomial_product(left, right):
    result = np.zeros((len(left)+len(right)-1, 2, 2), dtype=object)
    for i, a in enumerate(left):
        for j, b in enumerate(right): result[i+j] += a @ b
    return result


class DoubleMiddleTests(unittest.TestCase):
    def test_augmented_span_ranks_match_direct_GF4_including_overlap(self):
        h = fiber.outer.matrix(2, 4)
        n = tuple(tuple(x ^ int(i == j) for j, x in enumerate(row)) for i, row in enumerate(h))
        rng = random.Random(607)
        for _ in range(64):
            rm, cm = rng.randrange(65536), rng.randrange(65536)
            rows, columns = tuple(fiber._labels(rm)), tuple(fiber._labels(cm))
            actual = len(rref([[n[i][j] for j in columns] for i in rows], len(columns))[1])
            self.assertEqual(fiber.rank_update_masks(rm, cm), actual)
            self.assertEqual(fiber.rank_update_masks(cm, rm), actual)

    def test_exact_factor_masses_and_protected_coordinates(self):
        h = fiber.outer.matrix(2, 4)
        for source in (tuple(range(6)), tuple(range(8)), tuple(range(11))):
            data = fiber.fiber_counts(source, q=4)
            self.assertEqual(sum(data['outside_counts'].values()), 4**data['cross_rank'])
            self.assertEqual(sum(data['zero_fiber_counts'].values()), 4**(len(source)-data['cross_rank']))
            self.assertEqual(data['outside_counts'][0], 1)
            self.assertEqual(data['zero_fiber_counts'][0], 1)
            outside = tuple(i for i in range(16) if i not in source)
            basis = kernel([[h[i][j] for j in source] for i in outside], len(source))
            observed_protected = 0
            for position, label in enumerate(source):
                equal = all(apply(h, tuple(v[source.index(j)] if j in source else 0 for j in range(16)))[label]
                            == v[position] for v in basis)
                if equal: observed_protected |= 1 << label
            self.assertEqual(data['protected_mask'], observed_protected)
            for v in basis:
                x = tuple(v[source.index(j)] if j in source else 0 for j in range(16))
                y = apply(h, x)
                self.assertFalse(any(y[i] for i in outside))
                self.assertEqual(apply(h, y), x)

    def test_tiny_raw_coset_bound_and_input_degree_domination(self):
        source, q = tuple(range(6)), 4
        h = fiber.outer.matrix(2, q)
        data = fiber.fiber_counts(source, q=q)
        counts = Counter()
        protected_data = dict(data, physical_inside_mask=data['mapped_source_mask'],
                              physical_protected_mask=data['protected_mask'])
        restricted = tuple(tuple(row[i] for i in source) for row in h)
        for x in product(range(4), repeat=6):
            y = apply(restricted, x)
            mask = sum(bool(v) << i for i, v in enumerate(y))
            counts[mask] += 1
            if not mask & data['mapped_outside_mask']:
                self.assertLessEqual(sum(bool(v) for v in x), fiber.protected_degree(mask, protected_data))
        self.assertEqual(sum(counts.values()), 4**6)
        for t in (F(0), F(1, 3), F(1)):
            actual = sum(count*t**mask.bit_count() for mask, count in counts.items())
            bound = sum(count*t**mask.bit_count() for mask, count in data['outside_counts'].items())
            bound *= sum(count*t**mask.bit_count() for mask, count in data['zero_fiber_counts'].items())
            self.assertLessEqual(actual, bound)

    def test_schedule_transports_both_factors_and_protection(self):
        source = tuple(range(6))
        for j in (0, 7, 15):
            data = fiber.column_fiber_counts(source, j, q=4)
            before, after = target.column_permutations(fiber.SCHEDULE, j)
            raw = fiber.fiber_counts(tuple(before[i] for i in source), q=4)
            for name in ('outside_counts', 'zero_fiber_counts'):
                self.assertEqual(data[name], {fiber._permute_mask(k, after): n for k, n in raw[name].items()})
            self.assertEqual(data['physical_protected_mask'], fiber._permute_mask(raw['protected_mask'], after))
            self.assertEqual(data['original_source'], source)
            self.assertEqual(fiber.protected_degree(0, data), 0)

    def test_positive_tail_matches_exact_noncommuting_toy_and_batch(self):
        columns = np.array([[[[1, 2], [3, 1]], [[2, 1], [1, 3]]],
                            [[[2, 3], [1, 2]], [[1, 2], [4, 1]]],
                            [[[3, 1], [2, 1]], [[2, 4], [1, 1]]]], dtype=object)
        word = np.array([np.eye(2, dtype=object)])
        for column in columns:
            factor = np.zeros((3, 2, 2), dtype=object)
            for i in range(2):
                for j in range(2):
                    factor[:, i, j] = np.convolve(column[:, i, j], column[:, i, j])
            word = matrix_polynomial_product(word, factor)
        logs = np.log2(np.array(columns, dtype=float))
        for floor in range(7):
            expected = log2(int(np.sum(word[floor:, 0, :])))
            self.assertAlmostEqual(float(target.positive_deficit_tail(logs, 2, floor)), expected, places=12)
        batch = np.stack((logs, logs), axis=1)
        self.assertTrue(np.allclose(target.positive_deficit_tail(batch, 2, 4), target.positive_deficit_tail(logs, 2, 4)))

    def test_matches_frozen_lowarity_deficit_helper(self):
        logs = np.log2(np.arange(1, 16*3*4+1).reshape(16, 3, 2, 2))
        self.assertAlmostEqual(float(target.positive_deficit_tail(logs, 2, 36)),
                               float(column_deficit_tail(logs, 2, 2)), places=10)

    def test_exact_pairs_have_raw_activity_mass_and_fixed_source(self):
        table = np.zeros((65536, 2, 2))
        columns = target.exact_pair_columns((0, 7), table)
        expected = np.log2(np.array((1, 510, 65025), dtype=float))
        self.assertTrue(np.allclose(columns, expected[None, :, None, None]))
        self.assertTrue(np.allclose(columns, target.exact_pair_columns((7, 0), table)))

    def test_all_nonquiet_is_method_component_only(self):
        n = np.zeros((16, 2, 2))
        d = target.all_nonquiet_diagnostic(n, 11, 123)
        self.assertEqual(d['surrogate_degree'], 16*11*18)
        self.assertTrue(d['lower_bound_on_method_expression'])
        self.assertFalse(d['failure_probability_lower_bound'])
        self.assertAlmostEqual(d['majorant_component_log2'], 123+16-8*11*18*log2(255))

    def test_mock_controller_normalization_and_full_pair_scope(self):
        saved = dict(mark='9927/10000', source_count='123')
        with patch.object(target, 'load_measure', return_value=(saved, {}, 77.)), \
             patch.object(target, 'moment_table', return_value=None), \
             patch.object(target, 'exact_pair_columns', return_value=np.zeros((16, 3, 2, 2))), \
             patch.object(target, 'positive_deficit_tail', return_value=-1000.):
            data = target.run(mode='pairs')
        self.assertEqual(data['completed_sources'], 120)
        self.assertTrue(data['complete_all120_original_pairs'])
        self.assertEqual(len({tuple(r['original_source']) for r in data['controls']}), 120)
        expected = 77-1000-18*16*log2(255)+18*log2(120)
        self.assertAlmostEqual(data['all_pair_rowset_sum_bound_log2'], expected)
        self.assertFalse(data['old_outer_receipts_transferred'])
        self.assertEqual(data['regional_provider_runs'], 0)


if __name__ == '__main__':
    unittest.main()
