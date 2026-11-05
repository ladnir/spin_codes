from fractions import Fraction as F
from itertools import product
from math import comb, log2
from unittest import TestCase, main

import numpy as np

import hd_post_refresh_quiet_rows as gate


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def exact_column(rowsets, weights, size, constants, activity):
    """Independent labelled-mask enumeration; no path/count regrouping."""
    p, b00, b01, b10, a01, a10, a11 = constants
    answer = [[F(0), F(0)], [F(0), F(0)]]
    masks = [[sum((bits >> j & 1) << position for j, position in enumerate(row))
              for bits in range(1 << len(row))] for row in rowsets]
    for chosen in product(*masks):
        mass = F(1)
        for mask in chosen:
            mass *= weights[mask.bit_count()]*activity**mask.bit_count()
        matrix = [[F(1), F(0)], [F(0), F(1)]]
        for position in range(size):
            m = sum(mask >> position & 1 for mask in chosen)
            regional = [[(1-p if m == 0 else 0)+p*b00**m, a01*b01**m],
                        [a10*b10**m, a11]]
            matrix = multiply(matrix, regional)
        for i, j in product(range(2), repeat=2):
            answer[i][j] += mass*matrix[i][j]
    return answer


class QuietRowsTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2/3', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'),
            map(log2, self.constants)))
        self.weights = tuple(map(F, (1, 2, 3)))

    def test_all_rowsets_are_exact_gap_compositions(self):
        rows, gaps = gate.rowset_gaps(16, 8)
        self.assertEqual(len(rows), comb(16, 8))
        np.testing.assert_array_equal(gaps.sum(axis=1), 8)
        for row, gap in zip(rows, gaps):
            self.assertEqual(tuple(np.cumsum(gap[:-1])+np.arange(8)), row)
        self.assertEqual(len(set(map(tuple, gaps))), len(rows))
        self.assertFalse(gaps.flags.writeable)

    def test_noncommuting_column_matches_labelled_masks(self):
        answer = gate.quiet_column_matrices(np.log2(list(map(float, self.weights))),
            2, self.measure, log2(1.25), size=3, batch_size=2)
        for row, matrix in zip(answer['rowsets'], answer['column_log2_matrices']):
            expected = exact_column((row, row), self.weights, 3, self.constants, F(5, 4))
            np.testing.assert_allclose(2**matrix, np.array(expected, dtype=float), rtol=4e-14)
        self.assertFalse(np.allclose(answer['column_log2_matrices'][0],
                                    answer['column_log2_matrices'][-1]))
        self.assertFalse(answer['exact_empty_replacement'])
        self.assertTrue(answer['common_positive_measure'])

    def test_wholeword_group_minkowski_covers_fixed_heterogeneous_rowsets(self):
        activity = F(3, 2)
        answer = gate.quiet_row_bounds(np.log2(list(map(float, self.weights))),
            2, self.measure, log2(activity), size=3, columns=2, batch_size=2)
        actual = F(0)
        for rows in product(answer['rowsets'], repeat=2):
            column = exact_column(rows, self.weights, 3, self.constants, activity)
            word = multiply(column, column)
            actual += sum(word[0])
        self.assertGreaterEqual(answer['combined_word_log2']+1e-12, log2(actual))
        expected = 2*log2(sum((2**value)**.5 for value in answer['per_rowset_word_log2']))
        self.assertAlmostEqual(answer['combined_word_log2'], expected, places=12)
        self.assertFalse(answer['extra_rowset_factor'])

    def test_one_group_is_exact_sum_not_extra_binomial(self):
        answer = gate.quiet_row_bounds(np.log2(list(map(float, self.weights))),
            1, self.measure, 0., size=3, columns=2, batch_size=1)
        actual = F(0)
        for row in answer['rowsets']:
            column = exact_column((row,), self.weights, 3, self.constants, F(1))
            actual += sum(multiply(column, column)[0])
        self.assertAlmostEqual(answer['combined_word_log2'], log2(actual), places=12)

    def test_probability_one_and_batching(self):
        measure = dict(self.measure, log2_probability=0.)
        first = gate.quiet_row_bounds([0., -np.inf, log2(3)], 2, measure, 1.,
                                     size=4, columns=3, batch_size=1)
        second = gate.quiet_row_bounds([0., -np.inf, log2(3)], 2, measure, 1.,
                                      size=4, columns=3, batch_size=4)
        np.testing.assert_array_equal(first['column_log2_matrices'], second['column_log2_matrices'])
        np.testing.assert_array_equal(first['per_rowset_word_log2'], second['per_rowset_word_log2'])
        self.assertTrue(np.isfinite(first['combined_word_log2']))

    def test_batched_power_retains_order(self):
        matrices = np.array([[[2., 3.], [5., 7.]], [[3., 1.], [4., 2.]]])
        actual = 2**gate.log_matrix_power(np.log2(matrices), 5)
        for matrix, value in zip(matrices, actual):
            np.testing.assert_allclose(value, np.linalg.matrix_power(matrix, 5), rtol=2e-14)
        self.assertFalse(np.allclose(matrices[0]@matrices[1], matrices[1]@matrices[0]))

    def test_exact_empty_shortcut_would_break_mixed_rowset_bound(self):
        # Two groups, two positions, deterministic singleton masks. With exact
        # empty transitions, each homogeneous word is small but mixed words are not.
        epsilon = F(1, 100)
        empty, active = [[F(1), F(0)], [F(0), epsilon]], [[F(0), F(1)], [F(1), F(0)]]
        first, second = multiply(active, empty), multiply(empty, active)
        homogeneous = [sum(multiply(value, value)[0]) for value in (first, second)]
        invalid_rhs = (sum(float(value)**.5 for value in homogeneous))**2
        mixed_word = multiply(multiply(active, active), multiply(active, active))
        self.assertGreater(float(sum(mixed_word[0])), invalid_rhs)

    def test_guards_and_timeout(self):
        with self.assertRaises(ValueError):
            gate.rowset_gaps(17, 8)
        with self.assertRaises(ValueError):
            gate.moment_atoms(dict(self.measure, log2_probability=.1))
        with self.assertRaises(ValueError):
            gate.quiet_column_matrices([0., 0.], 1, self.measure, -1., size=2)
        with self.assertRaises(ValueError):
            gate.quiet_column_matrices([0., 0.], 1, self.measure, 0., size=2, batch_size=0)
        with self.assertRaises(ValueError):
            gate.run('unused', 'unused', log2_x=0., max_seconds=31)
        def stop():
            raise TimeoutError('toy budget')
        with self.assertRaises(TimeoutError):
            gate.quiet_row_bounds([0., 0.], 1, self.measure, 0., size=2, check_time=stop)


if __name__ == '__main__':
    main()
