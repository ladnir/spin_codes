"""Independent exact labelled-edge-word oracles; no production gate."""
from fractions import Fraction as F
from itertools import product
import unittest

import numpy as np

import hd_positive_upward as up
import hd_bit_syndrome_columns_directed as columns


def exact(value):
    return np.array([F(float(x))*F(2)**value.exponent for x in value.mantissas.flat],
                    dtype=object).reshape(value.mantissas.shape)


def mm(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(2)), F(0)) for j in range(2)] for i in range(2)]


def reference(terminal, classes, masses, marked, gaps, codes):
    """Enumerate all five-labelled words, independently of the count DP."""
    term, edge, matrices = exact(terminal), exact(masses), exact(gaps)
    counts = columns.count_geometry(marked)[0]
    indices = {row: i for i, row in enumerate(counts)}
    answer = []
    for cls, code in zip(classes, codes):
        accumulated = [[F(0), F(0)], [F(0), F(0)]]
        for word in product(range(5), repeat=marked):
            value = matrices[code[0]].tolist()
            for position, label in enumerate(word):
                matrix = [[F(0), F(0)], [F(0), F(0)]]
                i, j = columns.ENDPOINTS[label]
                matrix[i][j] = edge[label]
                value = mm(mm(value, matrix), matrices[code[position+1]])
            count = tuple(word.count(label) for label in range(5))
            scalar = term[cls, indices[count]]
            for i in range(2):
                for j in range(2):
                    accumulated[i][j] += scalar*value[i][j]
        answer.append(accumulated)
    return np.asarray(answer, dtype=object)


class DirectedColumnsTests(unittest.TestCase):
    def fixture(self, marked=3):
        count = len(columns.count_geometry(marked)[0])
        terminal = up.Scaled(np.array([[(i%7+1)/16 for i in range(count)],
                                      [(i%5+2)/32 for i in range(count)]]), 7)
        masses = up.Scaled(np.array([.5, .25, .125, .0625, .375]), -4)
        gaps = up.Scaled(np.array([[[1., 0.], [0., 1.]],
                                  [[.5, .25], [.125, .75]],
                                  [[.25, .5], [0., .125]],
                                  [[.5, 0.], [.25, .25]]]), 2)
        classes = np.array([0, 1, 0, 1], dtype=np.int8)
        codes = np.array([[(j+i)%4 for j in range(marked+1)] for i in range(4)], dtype=np.int16)
        return terminal, classes, masses, marked, gaps, codes

    def assert_close_upper(self, result, target):
        actual = exact(result)
        self.assertTrue(np.all(actual >= target))
        for got, want in zip(actual.flat, target.flat):
            if want:
                self.assertLessEqual(got/want, 1+F(1, 1 << 35))
            else:
                self.assertEqual(got, 0)

    def test_exact_labelled_noncommuting_words_and_batching(self):
        args = self.fixture()
        target = reference(*args)
        self.assertNotEqual(mm(exact(args[4])[1], exact(args[4])[2]),
                            mm(exact(args[4])[2], exact(args[4])[1]))
        for batch in (1, 3, 128):
            result, metadata = columns.counted_columns(*args, batch_size=batch)
            self.assert_close_upper(result, target)
            self.assertEqual(metadata['rowset_count'], 4)
            self.assertTrue(metadata['labelled_rows_preserved'])
            self.assertFalse(metadata['orbit_multiplicity_factor'])
            self.assertEqual(metadata['dp_passes'], 1)

    def test_order_reversal_is_not_silently_identified(self):
        args = list(self.fixture(2))
        forward, _ = columns.counted_columns(*args)
        args[-1] = args[-1][:, ::-1]
        reverse, _ = columns.counted_columns(*args)
        self.assert_close_upper(reverse, reference(*args))
        self.assertFalse(np.array_equal(exact(forward), exact(reverse)))

    def test_zero_marked_column_and_class_selection(self):
        args = self.fixture(0)
        result, metadata = columns.counted_columns(*args)
        self.assert_close_upper(result, reference(*args))
        self.assertEqual(metadata['final_atom_count_classes'], 1)
        self.assertEqual(metadata['marked_positions'], 0)

    def test_structural_zeros_and_input_preservation(self):
        args = list(self.fixture(2))
        args[2] = up.Scaled(np.array([.5, .5, 0., 0., 0.]), -2)
        args[4] = up.Scaled(np.array([[[1., 0.], [0., 1.]]]), 0)
        args[-1] = np.zeros((4, 3), dtype=np.int8)
        before = [value.mantissas.copy() for value in (args[0], args[2], args[4])]
        result, _ = columns.counted_columns(*args)
        self.assert_close_upper(result, reference(*args))
        self.assertTrue(np.all(result.mantissas[:, 0, 1] == 0))
        self.assertTrue(np.all(result.mantissas[:, 1, :] == 0))
        for wanted, value in zip(before, (args[0], args[2], args[4])):
            np.testing.assert_array_equal(wanted, value.mantissas)

    def test_large_scales_cancel_without_unscaled_overflow(self):
        args = list(self.fixture(2))
        args[0] = up.Scaled(args[0].mantissas, 5000)
        args[2] = up.Scaled(args[2].mantissas, -5000)
        args[4] = up.Scaled(args[4].mantissas, 1700)
        result, _ = columns.counted_columns(*args, batch_size=2)
        self.assert_close_upper(result, reference(*args))
        self.assertTrue(np.isfinite(result.mantissas).all())

    def test_count_geometry_matches_frozen_order_and_targets(self):
        counts, layers = columns.count_geometry(8)
        self.assertEqual(len(counts), 495)
        previous = ((0, 0, 0, 0, 0),)
        for step, (size, targets) in enumerate(layers, 1):
            current = sorted(set(tuple(x+int(i == label) for i, x in enumerate(row))
                                 for row in previous for label in range(5)))
            self.assertEqual(size, len(current))
            for label, target in enumerate(targets):
                self.assertEqual(len(set(target)), len(target))
                self.assertEqual([current[int(i)] for i in target],
                    [tuple(x+int(i == label) for i, x in enumerate(row)) for row in previous])
            previous = tuple(current)
        self.assertEqual(previous, counts)

    def test_guards_and_deadline(self):
        args = list(self.fixture(1))
        with self.assertRaises(ValueError):
            columns.counted_columns(*args, batch_size=0)
        with self.assertRaises(ValueError):
            columns.counted_columns(*args, batch_size=513)
        args[1] = np.array([0, 1, 2, 1])
        with self.assertRaises(ValueError):
            columns.counted_columns(*args)
        args = list(self.fixture(1))
        args[-1] = np.zeros((4, 1), dtype=int)
        with self.assertRaises(ValueError):
            columns.counted_columns(*args)
        def stop():
            raise TimeoutError('tiny deadline')
        with self.assertRaises(TimeoutError):
            columns.counted_columns(*self.fixture(1), check_time=stop)
        for marked in (-1, 9, True):
            with self.assertRaises(ValueError):
                columns.count_geometry(marked)


if __name__ == '__main__':
    unittest.main()
