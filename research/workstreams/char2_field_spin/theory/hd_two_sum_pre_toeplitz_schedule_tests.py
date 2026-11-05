"""Tiny exact-incidence and noncommuting tensor checks; no r18 gate."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest
from unittest.mock import patch

import numpy as np

import hd_two_sum_pre_toeplitz_schedule as target


def brute_transform(weights, slopes):
    group = len(weights).bit_length()-1
    result = np.full((2,)*(group+1), -np.inf)
    for path in product(range(2), repeat=group+1):
        terms = []
        for mask, value in enumerate(weights):
            terms.append(value+sum(slopes[path[i], path[i+1]] for i in range(group) if mask >> i & 1))
        result[path] = np.logaddexp2.reduce(terms)
    return result


class ScheduleTests(unittest.TestCase):
    def test_matching_orbit_partition(self):
        for group in (4, 8, 16):
            for weight in range(group+1):
                self.assertEqual(sum(target.orbit_size(weight, f, group)
                                     for f in range(max(0, weight-group//2), weight//2+1)), comb(group, weight))

    def test_incidence_against_direct_fixed_target_rowset_sum(self):
        group, arity, columns, alpha = 8, 2, 8, F(1, 5)
        new, old, counts = target.incidence_weights(arity, group, alpha, columns, 4)
        self.assertEqual(sum(row['input_sets'] for row in counts), comb(8, 2))
        matching = tuple((2*i, 2*i+1) for i in range(4))
        normalization = float(alpha)*4*arity*log2(255)
        for output in (0, 1, 3, 5, 15, 85, 255):
            terms = []
            for active in range(1 << group):
                if active.bit_count() != arity:
                    continue
                value = target.masks.fixed_mask_mass(active, matching, output, group=group)
                if value:
                    terms.append(float(columns*alpha)*(log2(value.numerator)-log2(value.denominator)))
            direct = (np.logaddexp2.reduce(terms)-normalization)/columns if terms else -np.inf
            self.assertAlmostEqual(new[output.bit_count(), target.full_pairs(output, matching)], direct, places=10)
        # The old fixed-target law has no matching dependence.
        for weight in range(group+1):
            values = [old[weight, f] for f in range(max(0, weight-group//2), weight//2+1)]
            self.assertTrue(np.allclose(values, values[0], atol=1e-12, rtol=0))

    def test_tensor_against_all_tiny_paths(self):
        slopes = np.array([[.2, -.7], [.6, -1.1]])
        for group in (4, 8):
            weights = np.array([-.3*mask.bit_count()+.01*(mask % 7) for mask in range(1 << group)])
            weights[0] = -np.inf
            actual = target.path_transform(weights, slopes)
            expected = brute_transform(weights, slopes)
            self.assertTrue(np.allclose(actual, expected, atol=1e-12, rtol=0))
            # Coordinate labels matter: reverse all mask bits and retain path
            # order. This fixture is intentionally not a symmetric type law.
            reflected = weights.reshape((2,)*group).transpose(tuple(reversed(range(group)))).ravel()
            self.assertGreater(np.max(np.abs(target.path_transform(reflected, slopes)-actual)), 1e-5)

    def test_kernel_and_chronological_product_against_brute(self):
        group = 4
        slopes = np.array([[.2, -.7], [.6, -1.1]])
        intercepts = np.array([.1, -.4, -.8, -.2])
        weights = np.array([-.3*mask.bit_count()+.02*mask for mask in range(1 << group)])
        empty, weights[0] = weights[0], -np.inf
        transformed = target.path_transform(weights, slopes)
        actual = target.kernel_from_transform(transformed, empty, intercepts, 2, .4, columns=4, minimum=3)
        expected = np.full((2, 2), -np.inf)
        for path in product(range(2), repeat=group+1):
            g = np.logaddexp2(empty, transformed[path]+.4)-3*.4/4
            term = 2*g+sum(intercepts[2*path[i]+path[i+1]] for i in range(group))
            expected[path[0], path[-1]] = np.logaddexp2(expected[path[0], path[-1]], term)
        self.assertTrue(np.allclose(actual, expected, atol=1e-12, rtol=0))
        other = np.array([[-.1, -.9], [-.2, -.7]])
        self.assertTrue(np.allclose(target.log_matmul(actual, other), np.log2(2**actual @ 2**other), atol=1e-12))
        self.assertFalse(np.allclose(target.log_matmul(actual, other), target.log_matmul(other, actual)))

    def test_weight_table_keeps_matching_labels(self):
        weights = np.array([[10*w+f for f in range(5)] for w in range(9)], dtype=float)
        first = tuple((2*i, 2*i+1) for i in range(4))
        second = tuple((i, i+4) for i in range(4))
        a, b = target.mask_weight_table(weights, first), target.mask_weight_table(weights, second)
        self.assertEqual(a[3], 21)
        self.assertEqual(b[3], 20)
        self.assertEqual(a[0], b[0])

    def test_path_geometry_counts_every_ordered_transition(self):
        for group in (4, 8):
            initial, final, counts = target.path_geometry(group)
            for index, path in enumerate(product(range(2), repeat=group+1)):
                self.assertEqual((initial[index], final[index]), (path[0], path[-1]))
                self.assertEqual(counts[index].tolist(),
                    [sum(2*path[i]+path[i+1] == e for i in range(group)) for e in range(4)])

    def test_zero_and_guards(self):
        self.assertTrue(np.all(np.isneginf(target.path_transform(np.full(16, -np.inf), np.zeros((2, 2))))))
        for weight, full in ((-1, 0), (3, 2), (16, 0)):
            with self.assertRaises(ValueError):
                target.orbit_size(weight, full)
        for budget in (True, 0, 31):
            with self.assertRaises(ValueError):
                target.run(budget)
        with self.assertRaises(ValueError):
            target.path_transform(np.zeros(15), np.zeros((2, 2)))
        # The budget is checked before any large incidence calculation; a
        # partial result must not claim the requested matching schedule done.
        anchor = dict(count=3, slopes=np.zeros(4), intercepts=np.zeros(4), minimum_slack=0.)
        with patch.object(target, 'load_frozen_anchor', return_value=anchor), \
                patch.object(target, 'incidence_weights', side_effect=TimeoutError('fixture')):
            result = target.run()
        self.assertEqual(result['status'], 'compute-budget-exceeded')
        self.assertFalse(result['complete_requested_controls'])
        self.assertEqual(result['completed_matching_columns'], 0)


if __name__ == '__main__':
    unittest.main()
