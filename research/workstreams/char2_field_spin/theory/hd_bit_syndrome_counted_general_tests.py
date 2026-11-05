"""Tiny exact word tests and frozen-index overlap; no production gate."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
from unittest import TestCase, main

import numpy as np

import hd_bit_syndrome_counted_general as general
from hd_bit_syndrome_orbit_gate import counted_column_matrices as old_columns
from hd_post_refresh_quiet_rows import quiet_group_values, transition_geometry


def logs(values):
    array = np.asarray(values, dtype=object)
    return np.array([log2(value) if value else -np.inf for value in array.flat]).reshape(array.shape)


def matmul(left, right):
    return tuple(sum((left[2*i+k]*right[2*k+j] for k in range(2)), F(0))
                 for i in range(2) for j in range(2))


def elementary(markers):
    answer = [F(1)]
    for value in markers:
        answer.append(F(0))
        for j in range(len(answer)-1, 0, -1):
            answer[j] += value*answer[j-1]
    return answer


def exact_column(caps, groups, marked, constants, gaps, codes, x=F(1)):
    p, b00, b01, b10, c01, c10, c11 = constants
    masses, atoms = (1-p, p, c01, c10, c11), (F(0), b00, b01, b10, F(1))
    answer = [F(0)]*4
    for labels in product(range(5), repeat=marked):
        matrix = gaps[codes[0]]
        for position, label in enumerate(labels):
            edge = [F(0)]*4
            i, j = general.ENDPOINTS[label]
            edge[2*i+j] = masses[label]
            matrix = matmul(matmul(matrix, edge), gaps[codes[position+1]])
        e = elementary(tuple(atoms[label] for label in labels))
        scalar = sum((count*value*x**b for b, (count, value) in enumerate(zip(caps, e))), F(0))**groups
        answer = [old+entry*scalar for old, entry in zip(answer, matrix)]
    return np.asarray(answer, dtype=float).reshape(2, 2)


class CountedGeneralTests(TestCase):
    def setUp(self):
        self.constants = tuple(map(F, ('2/5', '1/3', '2', '3/5', '3/7', '1/11', '2/13')))
        self.measure = dict(zip(('log2_probability', 'log2_atom00', 'log2_atom01',
            'log2_atom10', 'log2_scale01', 'log2_scale10', 'log2_scale11'), map(log2, self.constants)))

    def test_geometry_sizes_injective_targets_and_legacy_order(self):
        for marked in range(16):
            counts, layers = general.transition_geometry_general(marked)
            self.assertEqual(len(counts), comb(marked+4, 4))
            self.assertEqual(len(layers), marked)
            self.assertEqual(tuple(sorted(counts)), counts)
            self.assertTrue(all(sum(row) == marked for row in counts))
            if marked:
                previous, _ = general.transition_geometry_general(marked-1)
                for label, target in enumerate(layers[-1][1]):
                    self.assertFalse(target.flags.writeable)
                    self.assertEqual(len(set(map(int, target))), len(previous))
                    for row, index in zip(previous, target):
                        increment = list(row)
                        increment[label] += 1
                        self.assertEqual(counts[int(index)], tuple(increment))
            if 1 <= marked <= 8:
                old_counts, old_layers = transition_geometry(marked)
                self.assertEqual(counts, old_counts)
                for (size, targets), (old_size, old_targets) in zip(layers, old_layers):
                    self.assertEqual(size, old_size)
                    for target, old in zip(targets, old_targets):
                        np.testing.assert_array_equal(target, old)
        self.assertEqual(len(general.transition_geometry_general(9)[0]), 715)
        self.assertEqual(len(general.transition_geometry_general(15)[0]), 3876)

    def test_group_values_match_exact_polynomials_and_constant(self):
        for marked in range(5):
            caps = tuple(F(2*b+1) if b != 2 else F(0) for b in range(marked+1))
            actual = general.quiet_group_values_general(logs(caps), self.measure, 1.)
            counts, _ = general.transition_geometry_general(marked)
            for row, got in zip(counts, actual):
                markers = tuple(value for value, count in zip((F(0), *self.constants[1:4], F(1)), row)
                                for _ in range(count))
                e = elementary(markers)
                expected = sum(count*value*2**b for b, (count, value) in enumerate(zip(caps, e)))
                self.assertAlmostEqual(got, log2(expected), places=12)
        self.assertEqual(general.quiet_group_values_general([3.], self.measure, 5.).tolist(), [3.])
        self.assertEqual(general.quiet_group_values_general([-np.inf], self.measure, 0.).tolist(), [-np.inf])

    def test_group_values_exact_legacy_overlap(self):
        for marked in range(1, 9):
            caps = logs(tuple(1 if b % 3 else 0 for b in range(marked+1)))
            np.testing.assert_array_equal(general.quiet_group_values_general(caps, self.measure, .75),
                                          quiet_group_values(caps, self.measure, .75))

    def test_noncommuting_exact_labelled_words(self):
        gaps = tuple(tuple(map(F, row)) for row in ((1, 0, 0, 1), (1, 2, 0, 3), (2, 0, 5, 1)))
        products = logs(gaps).reshape(3, 2, 2)
        for marked in range(5):
            caps = (tuple(F(b+1) for b in range(marked+1)),
                    tuple(F(2*b+1) if b != 1 else F(0) for b in range(marked+1)))
            terminal = np.stack([2*general.quiet_group_values_general(logs(row), self.measure, 1.) for row in caps])
            classes = np.array([0, 1, 0], dtype=np.int8)
            codes = np.array([[(position+row) % 3 for position in range(marked+1)]
                              for row in range(3)], dtype=np.int32)
            actual, meta = general.counted_column_matrices(terminal, classes, self.measure, marked,
                                                           products, codes, batch_size=2)
            for row in range(3):
                expected = exact_column(caps[int(classes[row])], 2, marked, self.constants,
                                        gaps, codes[row], F(2))
                np.testing.assert_allclose(2**actual[row], expected, rtol=8e-14, atol=1e-14)
            self.assertFalse(meta['directed'])
            self.assertFalse(meta['whole_code_certificate'])

    def test_column_overlap_and_batch_invariance(self):
        products = logs(np.array([[[1, 2], [3, 4]], [[2, 0], [1, 1]]]))
        for marked in range(1, 9):
            terminal = np.stack([2*general.quiet_group_values_general(logs(tuple(range(1, marked+2))),
                                self.measure, 0.)+shift for shift in (0., 1.)])
            classes = np.array([0, 1, 1], dtype=np.int8)
            codes = np.array([[i % 2 for i in range(marked+1)],
                              [(i+1) % 2 for i in range(marked+1)], [0]*(marked+1)], dtype=np.int32)
            actual, _ = general.counted_column_matrices(terminal, classes, self.measure, marked, products, codes, batch_size=1)
            expected, _ = old_columns(terminal, classes, self.measure, marked, products, codes, batch_size=3)
            np.testing.assert_array_equal(actual, expected)

    def test_nine_and_fifteen_constant_polynomial_one_row(self):
        p, _, _, _, c01, c10, c11 = self.constants
        empty = tuple(map(F, (1, c01, c10, c11)))
        identity = (F(1), F(0), F(0), F(1))
        for marked in (9, 15):
            counts, _ = general.transition_geometry_general(marked)
            # Terminal9 represents the entire group power of a constant3.
            terminal = np.full((1, len(counts)), log2(9))
            answer, _ = general.counted_column_matrices(terminal, np.array([0]), self.measure, marked,
                logs(identity).reshape(1, 2, 2), np.zeros((1, marked+1), dtype=np.int8), batch_size=1)
            expected = identity
            for _ in range(marked):
                expected = matmul(expected, empty)
            np.testing.assert_allclose(2**answer[0], np.asarray(expected, dtype=float).reshape(2, 2)*9, rtol=8e-14)

    def test_structural_zeros_and_probability_one(self):
        measure = dict(self.measure, log2_probability=0.)
        marked = 3
        counts, _ = general.transition_geometry_general(marked)
        result, _ = general.counted_column_matrices(np.full((1, len(counts)), -np.inf), np.array([0]),
            measure, marked, np.array([[[0., -np.inf], [-np.inf, 0.]]]),
            np.zeros((1, marked+1), dtype=np.int8))
        self.assertTrue(np.all(np.isneginf(result)))

    def test_validation_and_early_deadline(self):
        for value in (-1, 16, True, 1.5):
            with self.assertRaises(ValueError):
                general.transition_geometry_general(value)
        for coefficients in ([], [float('nan')], [float('inf')], [0.]*17):
            with self.assertRaises(ValueError):
                general.quiet_group_values_general(coefficients, self.measure, 0.)
        with self.assertRaises(ValueError):
            general.quiet_group_values_general([0.], self.measure, -1.)
        def stop():
            raise TimeoutError('tiny budget control')
        with self.assertRaises(TimeoutError):
            general.counted_column_matrices([[0.]], np.array([0]), self.measure, 0,
                np.zeros((1, 2, 2)), np.array([[0]]), check_time=stop)
        with self.assertRaises(ValueError):
            general.counted_column_matrices([[0.]], np.array([0]), self.measure, 0,
                np.zeros((1, 2, 2)), np.array([[1]]))


if __name__ == '__main__':
    main()
