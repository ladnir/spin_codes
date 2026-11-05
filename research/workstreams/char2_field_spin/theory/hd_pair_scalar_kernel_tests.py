"""Tiny exact tests only; no16-row production support enumeration."""
from fractions import Fraction
from itertools import product
import unittest
from unittest.mock import patch

import numpy as np

import hd_pair_scalar_kernel as kernel
from hd_blockwise_shared_transfer import multiply, rank


def matvec(matrix, values, q):
    result = []
    for row in matrix:
        value = 0
        for a, b in zip(row, values):
            value ^= multiply(a, b, q)
        result.append(value)
    return tuple(result)


def mask(values):
    return sum(1 << i for i, value in enumerate(values) if value)


def brute_counts(matrix, q):
    result = [0]*(1 << len(matrix))
    for values in product(range(q), repeat=len(matrix[0])):
        result[mask(matvec(matrix, values, q))] += 1
    return tuple(result)


def brute_pair(matrix, t, delta, q):
    """Enumerate actual inputs and fresh nonzero ratios, not support laws."""
    rows, columns = len(matrix), len(matrix[0])
    inputs = tuple(product(range(q), repeat=columns))
    outputs = tuple(matvec(matrix, values, q) for values in inputs)
    result = [Fraction(0)]*3
    for x, a in zip(inputs, outputs):
        for y, b in zip(inputs, outputs):
            degree = int(any(x))+int(any(y))
            for ratios in product(range(1, q), repeat=rows):
                weight = 0
                for i, alpha in enumerate(ratios):
                    scaled = multiply(alpha, b[i ^ delta], q)
                    weight += int((a[i] ^ scaled) != 0)
                    weight += int((a[i] ^ multiply(2, scaled, q)) != 0)
                result[degree] += t**weight
    return tuple(value/(q-1)**rows for value in result)


class PairScalarTests(unittest.TestCase):
    def test_exact_mobius_matches_small_field_inputs(self):
        examples = (
            ((1, 0), (0, 1), (1, 1), (1, 2)),
            ((1, 1), (2, 2), (0, 0), (3, 3)),
            ((1, 0, 2), (0, 1, 1), (1, 1, 3)),
        )
        for matrix in examples:
            for batch in (1, 5, 128):
                counts = kernel.support_counts_from_matrix(matrix, q=4, batch_size=batch)
                self.assertEqual(counts, brute_counts(matrix, 4))
                self.assertTrue(all(type(value) is int for value in counts))
                self.assertTrue(any(value == 0 for value in counts))

    def test_rank_batches_match_independent_elimination(self):
        for q in (4, 256):
            matrix = tuple(tuple((7*i+11*j+3*i*j) % q for j in range(4)) for i in range(6))
            ranks = kernel.zero_set_ranks(matrix, q=q, batch_size=7)
            expected = [rank(tuple(matrix[i] for i in range(6) if subset >> i & 1), q)
                        for subset in range(64)]
            np.testing.assert_array_equal(ranks, expected)

    def test_python_integer_range_and_empty_source(self):
        counts = kernel.support_counts_from_matrix(((1,)*16,), q=256)
        self.assertEqual(counts, (256**15, 255*256**15))
        self.assertGreater(counts[0], 2**63)
        empty = np.empty((4, 0), dtype=int)
        self.assertEqual(kernel.support_counts_from_matrix(empty, q=4), (1,)+(0,)*15)
        self.assertEqual(kernel.support_counts((), h=((1, 0), (0, 1)), q=4), (1, 0, 0, 0))

    def test_supplied_matrix_source_labels_retained(self):
        h = ((1, 2, 3), (2, 1, 1), (0, 3, 2), (1, 0, 3))
        selected = tuple((row[0], row[2]) for row in h)
        self.assertEqual(kernel.support_counts((0, 2), h=h, q=4), brute_counts(selected, 4))
        from hd_bit_syndrome_outer import matrix
        small = tuple(tuple(row[j] for j in (0, 3)) for row in matrix()[:4])
        counts = kernel.support_counts_from_matrix(small, q=256)
        self.assertEqual(sum(counts), 256**2)
        self.assertEqual(counts[0], 1)

    def test_scalar_pair_law_for_every_fixed_gf4_pair(self):
        t = Fraction(2, 5)
        interaction = kernel.pair_kernel(t, q=4)
        for a, b in product(range(4), repeat=2):
            exact = sum(t**(int((a ^ multiply(alpha, b, 4)) != 0)
                           +int((a ^ multiply(2, multiply(alpha, b, 4), 4)) != 0))
                        for alpha in range(1, 4))/3
            self.assertAlmostEqual(interaction[int(a != 0), int(b != 0)], float(exact))

    def test_tensor_transform_matches_labelled_bruteforce(self):
        values = np.array((0, 1, 0, 2, 3, 0, 5, 1), dtype=float)
        interaction = np.array(((0.7, 0.2), (0.4, 0.9)))
        expected = np.array([sum(values[t]*np.prod([interaction[s >> i & 1, t >> i & 1]
                                                  for i in range(3)])
                                 for t in range(8)) for s in range(8)])
        np.testing.assert_allclose(kernel.tensor_transform(values, interaction), expected, rtol=1e-14)
        for delta in range(4):
            found = kernel.xor_mask_indices(4, delta)
            expected = [sum(((s >> i) & 1) << (i ^ delta) for i in range(4)) for s in range(16)]
            np.testing.assert_array_equal(found, expected)
            np.testing.assert_array_equal(found[found], np.arange(16))

    def test_pair_degrees_match_actual_tiny_ensemble(self):
        for matrix in (((1, 0), (0, 1)), ((1, 1), (2, 2)), ((1,), (1,), (0,), (0,))):
            t = Fraction(1, 2)
            deltas = tuple(range(len(matrix)))
            counts = kernel.support_counts_from_matrix(matrix, q=4)
            data = kernel.pair_block_coefficients(counts, t, q=4, deltas=deltas)
            for delta in deltas:
                expected = brute_pair(matrix, t, delta, 4)
                np.testing.assert_allclose(data['coefficients'][delta], [float(v) for v in expected], rtol=2e-14)
            self.assertEqual(data['metadata']['activity_degree'], 'number-of-nonzero-original-input-columns')
            self.assertFalse(data['metadata']['whole_code_certificate'])

    def test_one_transform_and_separate_matching_cases(self):
        counts = [0]*256
        counts[0], counts[3], counts[240] = 1, 3, 12
        with patch.object(kernel, 'tensor_transform', wraps=kernel.tensor_transform) as transformed:
            result = kernel.pair_block_coefficients(tuple(counts), 0.4, q=4, deltas=tuple(range(8)))
            self.assertEqual(transformed.call_count, 1)
        self.assertEqual(result['metadata']['permuted_dot_products'], 8)
        self.assertGreater(np.ptp([row[2] for row in result['coefficients'].values()]), 0)
        self.assertEqual(counts[0], 1)

    def test_endpoints_and_noninjective_input_activity(self):
        counts = kernel.support_counts_from_matrix(((1, 1), (1, 1)), q=4)
        at_zero = kernel.pair_block_coefficients(counts, 0, q=4, deltas=(0, 1))
        at_one = kernel.pair_block_coefficients(counts, 1, q=4, deltas=(0, 1))
        for delta in (0, 1):
            np.testing.assert_array_equal(at_zero['coefficients'][delta], (1, 6, 9))
            np.testing.assert_array_equal(at_one['coefficients'][delta], (1, 30, 225))
        zero = kernel.pair_block_coefficients((1, 0, 0, 0), 0.2, q=4, deltas=(1,))
        np.testing.assert_array_equal(zero['coefficients'][1], (1, 0, 0))

    def test_guards_and_cooperative_budget(self):
        for value in (-1, float('nan'), 1.01):
            with self.assertRaises(ValueError):
                kernel.pair_kernel(value)
        with self.assertRaises(ValueError):
            kernel.support_counts_from_matrix(np.ones((17, 1), dtype=int))
        with self.assertRaises(ValueError):
            kernel.support_counts_from_matrix(((1,),), batch_size=8193)
        with self.assertRaises(ValueError):
            kernel.pair_block_coefficients((0, 3, 0, 0), 0.5, deltas=(1,))
        with self.assertRaises(ValueError):
            kernel.pair_block_coefficients((1, 3, 0, 0), 0.5, deltas=(1, 1))
        with self.assertRaises(ValueError):
            kernel.support_counts((1, 0), h=((1, 0), (0, 1)), q=4)
        with self.assertRaises(TimeoutError):
            kernel.support_counts_from_matrix(((1,), (2,)), q=4,
                check_time=lambda: (_ for _ in ()).throw(TimeoutError()))


if __name__ == '__main__':
    unittest.main()
