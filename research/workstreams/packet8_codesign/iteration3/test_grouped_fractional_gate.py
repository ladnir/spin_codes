"""Tiny exhaustive checks of coarser conditioning before fractional powers."""
from collections import defaultdict
from itertools import combinations
from math import comb
import unittest

import numpy as np

import grouped_fractional_gate as grouped


def product(matrices):
    result = np.eye(len(matrices[0]))
    for matrix in matrices:
        result = result@matrix
    return result


def scalar(matrix):
    return matrix[0].sum()


class GroupedFractionalTests(unittest.TestCase):
    def setUp(self):
        self.local = np.array([[[.8, .0], [.0, .3]],
            [[.15, .25], [.03, .08]], [[.08, .02], [.09, .04]]])

    def test_macro_matches_exhaustive_ordered_placements(self):
        for steps in (1, 2, 4):
            actual, _ = grouped.grouped_operators(self.local, steps)
            for j in range(2*steps+1):
                expected = np.zeros((2, 2))
                for positions in combinations(range(2*steps), j):
                    counts = [sum(x//2 == step for x in positions) for step in range(steps)]
                    expected += product([self.local[k] for k in counts])/comb(2*steps, j)
                np.testing.assert_allclose(actual[j], expected, rtol=2e-13, atol=1e-15)
        self.assertGreater(np.max(np.abs(self.local[0]@self.local[1]
            -self.local[1]@self.local[0])), .01)

    def test_alpha1_composition_equals_full_placement(self):
        macro, _ = grouped.grouped_operators(self.local, 2)
        grouped_logs, _ = grouped.prior.placement(macro, 8, epochs=2, windows=4)
        direct_logs, _ = grouped.prior.placement(self.local, 8, epochs=4, windows=2)
        np.testing.assert_allclose(np.exp(grouped_logs), np.exp(direct_logs),
            rtol=2e-13, atol=1e-15)

    def test_conditional_clipping_and_macro_path_bound(self):
        # Four two-slot physical steps form two two-step macros. Their total
        # marked occupancy is fixed at three; neither macro has fixed count.
        macro, _ = grouped.grouped_operators(self.local, 2)
        records = defaultdict(list)
        for positions in combinations(range(8), 3):
            fine = [sum(x//2 == step for x in positions) for step in range(4)]
            coarse = (sum(fine[:2]), sum(fine[2:]))
            records[coarse].append(scalar(product([self.local[k] for k in fine])))
        total = comb(8, 3)
        for beta in (2., 10., 40.):
            direct, conditioned = 0., 0.
            for counts, values in records.items():
                mean = scalar(macro[counts[0]]@macro[counts[1]])
                self.assertAlmostEqual(mean, sum(values)/len(values), places=14)
                direct += sum(min(1., beta*v) for v in values)/total
                conditioned += len(values)*min(1., beta*mean)/total
            self.assertLessEqual(direct, conditioned+1e-14)
            for alpha in (.3, .4, .6, 1.):
                powered = grouped.powered_operators(macro, alpha)
                path_bound = sum(len(values)*scalar(powered[c[0]]@powered[c[1]])
                    for c, values in records.items())/total
                # An independent outer-subset union must not be powered.
                union = comb(5, 2)
                self.assertLessEqual(union*conditioned, union*beta**alpha*path_bound+1e-13)

    def test_conditional_average_is_before_power(self):
        macro, _ = grouped.grouped_operators(self.local, 2)
        correct = grouped.powered_operators(macro, .4)
        wrong, _ = grouped.grouped_operators(self.local**.4, 2)
        self.assertGreater(np.max(np.abs(correct-wrong)), .01)

    def test_invalid_parameters(self):
        for group in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                grouped.grouped_operators(self.local, group)
        for alpha in (0., -1., 1.1, float('nan')):
            with self.assertRaises(ValueError):
                grouped.powered_operators(self.local, alpha)


if __name__ == '__main__':
    unittest.main()
