from fractions import Fraction as F
from itertools import product
from math import log2
from unittest import TestCase, main

import numpy as np

import hd_post_refresh_full_rows as gate
from hd_post_refresh_quiet_rows_tests import multiply
from hd_two_sum_closed_kernel import branch_counts


def regional(constants, m):
    p, b0, b01, b10, a01, a10, a11 = constants
    return [[(1-p if m == 0 else 0)+p*b0**m, a01*b01**m],
            [a10*b10**m, a11]]


def mask_options(row, size, counts, activity):
    q, f, h = counts
    outside = tuple(j for j in range(size) if j not in row)
    omask = sum(1 << j for j in outside)
    answer = []
    for bits in range(1 << len(row)):
        mask = sum((bits >> j & 1) << pos for j, pos in enumerate(row))
        b = bits.bit_count()
        answer.append((mask, q[b]*activity**b))
        answer.append((mask | omask, f[b]*activity**len(row)))
        for missing in outside:
            answer.append((mask | (omask ^ (1 << missing)), h[b]*activity**len(row)))
    return answer


def exact_full_column(rows, size, counts, constants, activity):
    answer = np.zeros((2, 2), dtype=object)
    for chosen in product(*(mask_options(row, size, counts, activity) for row in rows)):
        mass, value = F(1), [[F(1), F(0)], [F(0), F(1)]]
        for _, weight in chosen:
            mass *= weight
        for position in range(size):
            m = sum(mask >> position & 1 for mask, _ in chosen)
            value = multiply(value, regional(constants, m))
        for i, j in product(range(2), repeat=2):
            answer[i, j] += mass*value[i][j]
    return answer.tolist()


def exact_marked(row, size, counts, constants, groups, hole):
    outside = tuple(j for j in range(size) if j not in row)
    answer = np.zeros((2, 2), dtype=object)
    for missing in outside if hole else (None,):
        for chosen in product(range(1 << len(row)), repeat=groups):
            mass, value = F(1), [[F(1), F(0)], [F(0), F(1)]]
            for mask in chosen:
                mass *= counts[mask.bit_count()]
            for position in range(size):
                if position in row:
                    bit = row.index(position)
                    m = sum(mask >> bit & 1 for mask in chosen)
                else:
                    m = 0 if position == missing else groups
                value = multiply(value, regional(constants, m))
            for i, j in product(range(2), repeat=2):
                answer[i, j] += mass*value[i][j]
    return answer.tolist()


class FullRowsTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2/3', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))

    def test_gap_one_hole_exact_noncommuting_order(self):
        base = np.array([[2., 3.], [5., 7.]])
        hole = np.array([[1., 4.], [2., 3.]])
        self.assertFalse(np.allclose(base@hole, hole@base))
        actual = 2**gate.gap_powers(np.log2(base), 4, hole=np.log2(hole))
        for n in range(5):
            expected = sum((np.linalg.matrix_power(base, j)@hole
                @np.linalg.matrix_power(base, n-1-j) for j in range(n)), np.zeros((2, 2)))
            np.testing.assert_allclose(actual[n, 0], np.linalg.matrix_power(base, n), rtol=1e-14)
            np.testing.assert_allclose(actual[n, 1], expected, rtol=2e-14)

    def test_full_and_one_hole_match_labelled_enumeration(self):
        counts = tuple(map(F, (2, 3, 5)))
        for hole in (False, True):
            answer = gate.marked_column_matrices(np.log2(list(map(float, counts))), 2,
                self.measure, size=4, one_hole=hole, batch_size=2)
            for row, matrix in zip(answer['rowsets'], answer['column_log2_matrices']):
                expected = exact_marked(row, 4, counts, self.constants, 2, hole)
                np.testing.assert_allclose(2**matrix, np.array(expected, dtype=float), rtol=5e-14)
            self.assertEqual(answer['costs']['hole_flag_degree'], int(hole))

    def test_group_one_full_union_is_exact(self):
        counts = tuple(tuple(map(F, row)) for row in ((1, 2), (1, 3), (2, 1)))
        logs = np.log2(np.array(counts, dtype=float))
        answer = gate.full_row_bounds(logs, 1, self.measure, log2(1.5),
                                     size=3, columns=2, batch_size=2)
        total = F(0)
        for row, matrix in zip(answer['rowsets'], answer['column_log2_matrices']):
            expected = exact_full_column((row,), 3, counts, self.constants, F(3, 2))
            np.testing.assert_allclose(2**matrix, np.array(expected, dtype=float), rtol=5e-14)
            total += sum(multiply(expected, expected)[0])
        self.assertAlmostEqual(answer['combined_word_log2'], log2(total), places=12)

    def test_young_and_wholeword_minkowski_cover_heterogeneous_fixed_rows(self):
        counts = tuple(tuple(map(F, row)) for row in ((1, 2), (1, 3), (2, 1)))
        answer = gate.full_row_bounds(np.log2(np.array(counts, dtype=float)), 2,
            self.measure, 0., size=3, columns=2, lambdas=(F(1, 3),)*3, batch_size=2)
        actual = F(0)
        for rows in product(answer['rowsets'], repeat=2):
            column = exact_full_column(rows, 3, counts, self.constants, F(1))
            actual += sum(multiply(column, column)[0])
        self.assertGreaterEqual(answer['combined_word_log2']+1e-12, log2(actual))
        self.assertFalse(answer['extra_rowset_factor'])

    def test_young_factors_and_activity_once(self):
        zeros = np.zeros((1, 2, 2))
        _, shifts = gate.combine_columns(zeros, zeros, zeros, 18, 8, 8, 2.)
        expected = [-17*log2(F(99, 100)), -17*log2(F(1, 200))+288,
                    -17*log2(F(1, 200))+51+288]
        np.testing.assert_allclose(shifts, expected, rtol=1e-15)
        # No per-b activity marker enters the marked full/hole primitive.
        first = gate.marked_column_matrices([0., 1.], 2, self.measure, size=3, batch_size=1)
        second = gate.marked_column_matrices([0., 1.], 2, self.measure, size=3, batch_size=3)
        np.testing.assert_array_equal(first['column_log2_matrices'], second['column_log2_matrices'])

    def test_g16_branch_counts_preserve_exact_uniform_total(self):
        from math import comb
        total = 0
        for b in range(9):
            quiet, hole, full = branch_counts(8, b, 16)
            total += comb(8, b)*(quiet+8*hole+full)
        self.assertEqual(total, 256**8)

    def test_zero_branch_and_guards(self):
        result = gate.marked_column_matrices([-np.inf, -np.inf], 1, self.measure,
                                            size=3, one_hole=True)
        self.assertTrue(np.isneginf(result['column_log2_matrices']).all())
        self.assertFalse(np.isnan(result['column_log2_matrices']).any())
        with self.assertRaises(ValueError):
            gate.checked_lambdas((F(1, 2),)*3)
        with self.assertRaises(ValueError):
            gate.marked_column_matrices([0., 0.], 1, self.measure, size=1)
        with self.assertRaises(ValueError):
            gate.run('unused', 'unused', log2_x=0., max_seconds=31)
        def stop():
            raise TimeoutError('toy budget')
        with self.assertRaises(TimeoutError):
            gate.full_row_bounds(np.zeros((3, 2)), 2, self.measure, 0., size=3, check_time=stop)


if __name__ == '__main__':
    main()
