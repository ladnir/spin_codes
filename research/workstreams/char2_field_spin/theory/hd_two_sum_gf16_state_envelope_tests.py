"""Tiny exhaustive controls; never builds a production regional bank."""
from itertools import product
from math import comb
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_gf16_state_envelope as target


class StateEnvelopeTests(unittest.TestCase):
    def test_path_counts_against_every_binary_path(self):
        for length in range(1, 7):
            expected = {}
            for bits in product(range(2), repeat=length+1):
                counts = [0]*4
                for left, right in zip(bits, bits[1:]):
                    counts[2*left+right] += 1
                key = bits[0], bits[-1], tuple(counts)
                expected[key] = expected.get(key, 0)+1
            self.assertEqual({(i, j, n): m for i, j, n, m in target.path_classes(length)}, expected)
        classes = target.path_classes()
        self.assertEqual(len(classes), 274)
        self.assertEqual([sum(i == a and j == b for i, j, _, _ in classes)
                          for a, b in product(range(2), repeat=2)], [65, 72, 72, 65])
        self.assertEqual(sum(row[3] for row in classes), 2**17)

    def test_elementary_polynomials_count_each_labelled_mask_once(self):
        marks = np.array([[2., 3., 5., 7.], [11., 2., 1., 4.]])
        got = 2**target.elementary_logs(np.log2(marks))
        expected = np.zeros((2, 5))
        for bits in product(range(2), repeat=4):
            expected[:, sum(bits)] += np.prod(marks**bits, axis=1)
        np.testing.assert_allclose(got, expected, rtol=1e-14)
        np.testing.assert_allclose(target.elementary_logs(np.zeros((1, 16)))[0],
                                   np.log2([comb(16, w) for w in range(17)]))

    def test_type_weights_against_labelled_sparse_and_dense_masks(self):
        # Toy4-region mixer: pair output is full with probability2/6,
        # or misses each particular coordinate with probability1/6.
        q, group, columns, dimension, minimum = 7, 4, 5, 2, 3
        marks, floor = np.array([2., 3., 5., 7.]), 2.
        weights = {3: np.array([1., 2., 3., 4., 5.]),
                   4: np.array([2., 1., 4., 3., 6.])}
        got = 2**target.type_weights(np.log2(marks), {a: np.log2(w) for a, w in weights.items()},
            np.log2(floor), columns=columns, dimension=dimension, minimum=minimum, q=q)
        full = np.prod(marks)
        expected1 = (group/(q-1)**dimension/floor**minimum)**(1/columns)*(1+(q-1)*floor*full)
        law = (q-1-group)*full/(q-1)+sum(full/mark for mark in marks)/(q-1)
        expected2 = (comb(group, 2)/(q-1)**(2*dimension)/floor**(2*minimum))**(1/columns)
        expected2 *= 1+2*(q-1)*floor*full+(q-1)**2*floor**2*law
        expected = [expected1, expected2]
        for arity in (3, 4):
            value = 0.
            for bits in product(range(2), repeat=group):
                value += weights[arity][sum(bits)]*np.prod(marks**bits)*(floor if any(bits) else 1)
            expected.append(value/floor**(minimum/columns))
        np.testing.assert_allclose(got, expected, rtol=2e-14)

    def test_complete_curves_envelope_and_zero_occupancy_diagnostics(self):
        regions = {}
        for m in range(5):
            for f in range(m//2+1):
                regions[m, f] = np.array([[-m+f, -3.-m], [-4.-m+f, -1.-m]])
        regions[0, 0] = np.array([[0., -np.inf], [-np.inf, -1.]])
        curves = target.regional_curves(regions, 4)
        self.assertTrue(np.isneginf(curves[1, 0]))
        slopes, candidates = target.slope_candidates(curves)
        a, slack = target.affine_envelope(curves, slopes)
        self.assertTrue(np.all(slack >= -1e-14))
        self.assertTrue(np.all(np.isfinite(a)))
        self.assertTrue(all(len(row) <= 16 for row in candidates))
        # Positive birth/return envelope entries at0 are slack, not zeros.
        self.assertGreater(2**a[1], 0.)
        self.assertTrue(np.isinf(slack[1, 0]))
        with self.assertRaises(ValueError):
            target.regional_curves({k: v for k, v in regions.items() if k != (4, 2)}, 4)

    def test_fixed_type_holder_then_path_minkowski(self):
        classes = target.path_classes(3)
        counts = np.array([row[2] for row in classes])
        path_a = counts @ np.log2([2., 1., 3., 2.])
        # Arbitrary positive path-dependent type factors suffice for this
        # norm inequality. Keep each type vector fixed across ALL columns.
        typed = np.log2(np.array([[1+row[2][0], 2+row[2][3]] for row in classes]))
        groups, columns = 2, 3
        kernel = 2**target.kernel_from_types(classes, path_a, typed, groups, columns)
        matrices = []
        for types in product(range(2), repeat=groups):
            matrix = np.zeros((2, 2))
            for p, (i, j, _, multiplicity) in enumerate(classes):
                matrix[i, j] += multiplicity*2**(path_a[p]+sum(typed[p, t] for t in types))
            matrices.append(matrix)
        exact_entry_norm = sum(matrix**columns for matrix in matrices)**(1/columns)
        self.assertTrue(np.all(kernel >= exact_entry_norm*(1-1e-14)))
        exact_union = sum(np.linalg.matrix_power(matrix, columns)[0].sum() for matrix in matrices)
        upper = np.linalg.matrix_power(kernel, columns)[0].sum()
        self.assertGreaterEqual(upper, exact_union*(1-1e-14))

    def test_one_type_removes_holder_and_minkowski_loss(self):
        classes = target.path_classes(4)
        path_a = np.log2(np.arange(1, len(classes)+1, dtype=float))
        typed = np.log2(np.array([[2.+row[2][1]] for row in classes]))
        got = 2**target.kernel_from_types(classes, path_a, typed, 3, 5)
        exact = np.zeros((2, 2))
        for p, (i, j, _, multiplicity) in enumerate(classes):
            exact[i, j] += multiplicity*2**(path_a[p]+3*typed[p, 0])
        np.testing.assert_allclose(got, exact, rtol=2e-14)

    def test_matrix_power_retains_initial_state_and_noncommuting_endpoint_sum(self):
        matrix = np.array([[2., 3.], [5., 7.]])
        np.testing.assert_allclose(2**target.matrix_power_log(np.log2(matrix), 16),
                                   np.linalg.matrix_power(matrix, 16), rtol=2e-14)

    def test_bound_guards_and_no_implicit_bank(self):
        for kwargs in ({'max_seconds': 31}, {'max_seconds': True}, {'max_seconds': float('nan')},
                       {'rounds': 4}, {'curves_only': 1}):
            with patch.object(target, 'forced_nonzero_regions') as bank, self.assertRaises(ValueError):
                target.run(**kwargs)
            bank.assert_not_called()
        with patch.object(target, 'forced_nonzero_regions') as bank:
            data = target.run(max_seconds=1e-12)
        bank.assert_not_called()
        self.assertFalse(data['complete_requested_band'])
        self.assertFalse(data['directed'])


if __name__ == '__main__':
    unittest.main()
