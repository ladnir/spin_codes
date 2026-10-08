"""Exact local and counting checks; no whole-code certificate."""

from fractions import Fraction as F
from itertools import product
import math
import unittest

import numpy as np

import cauchy_firstmoment_screen as screen


def mul(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return value


class FirstMomentTests(unittest.TestCase):
    def test_fraction_local_edges(self):
        local = screen.local_by_weight(F(1), width=4, state=2, q=4)
        self.assertEqual(local[0], [[1, 0], [0, 1]])
        for w in range(1, 5):
            self.assertEqual(local[w][0][1], 1)
            self.assertEqual(local[w][1][1], 1)
            self.assertEqual(local[w][1][0], F(1, 3**min(w, 2)))
            self.assertEqual(local[w][0][0], F(1, 9) if w > 2 else 0)

    def test_marginal_by_exhaustive_gf4_diagonals(self):
        # Three-by-two superregular matrix with rows (1,a), a=1,2,3.
        y = F(2, 5)
        local = screen.local_by_weight(y, width=3, state=2, q=4)
        for state in product(range(4), repeat=2):
            if not any(state):
                continue
            for mask in product((0, 1), repeat=3):
                total, count = F(0), 0
                for state_d in product((1, 2, 3), repeat=2):
                    refreshed = tuple(mul(a, b) for a, b in zip(state, state_d))
                    for input_values in product(*[(1, 2, 3) if bit else (0,) for bit in mask]):
                        output = [input_values[i] ^ refreshed[0] ^ mul(i+1, refreshed[1])
                                  for i in range(3)]
                        total += y**sum(value != 0 for value in output)
                        count += 1
                self.assertLessEqual(total/count, local[sum(mask)][1][1])

    def test_one_input_joint_return_bound(self):
        # This lemma needs full-support injection columns and MDS emission,
        # not CA=0.  The marginal and complete recurrence are separate claims.
        y = F(2, 5)
        for state in product(range(4), repeat=2):
            if not any(state):
                continue
            for lane in range(3):
                total, count = F(0), 0
                for d0, d1, input_d in product((1, 2, 3), repeat=3):
                    refreshed = (mul(d0, state[0]), mul(d1, state[1]))
                    following = (refreshed[0] ^ input_d,
                                 refreshed[1] ^ mul(lane+1, input_d))
                    if not any(following):
                        output = [refreshed[0] ^ mul(i+1, refreshed[1])
                                  ^ (input_d if i == lane else 0) for i in range(3)]
                        total += y**sum(value != 0 for value in output)
                    count += 1
                self.assertLessEqual(total/count, y/F(9))

    def test_aggregate_selection_count(self):
        x, y = F(2, 7), F(3, 5)
        local = screen.local_by_weight(y, width=3, state=2, q=4)
        exact = [[F(0), F(0)], [F(0), F(0)]]
        # 0=unselected; 1=selected zero; 2..4=three selected nonzero values.
        for choices in product(range(5), repeat=3):
            selected = sum(choice != 0 for choice in choices)
            weight = sum(choice >= 2 for choice in choices)
            for r in range(2):
                for k in range(2):
                    exact[r][k] += x**selected*local[weight][r][k]
        actual = screen.aggregate(float(x), float(y), width=3, state=2, q=4)
        self.assertTrue(np.allclose(actual, np.array(exact, dtype=float), rtol=2e-15))

    def test_positive_nonzero_tracker(self):
        x, y, width = .03, .7, 4
        ordinary = screen.aggregate(x, y, width=width, state=2, q=4)
        tracked = screen.aggregate(x, y, width=width, state=2, q=4, nonzero_only=True)
        for packets in (1, 2, 5, 16):
            a = np.linalg.matrix_power(ordinary, packets)[0].sum()
            b = np.linalg.matrix_power(tracked, packets)[0, 1:].sum()
            self.assertAlmostEqual(a-b, (1+x)**(width*packets), places=8)
            self.assertAlmostEqual(screen.log_terminal_power(tracked, packets, True),
                                   math.log(b), places=12)

    def test_exact_mds_support_counts(self):
        counts = {}
        for constant, slope in product(range(4), repeat=2):
            if not constant and not slope:
                continue
            word = [constant ^ mul(slope, x) for x in range(4)]
            support = tuple(i for i, value in enumerate(word) if value)
            counts[support] = counts.get(support, 0)+1
        for support, count in counts.items():
            self.assertEqual(count//3, screen.exact_mds_projective_support(len(support), 2, 4))
        self.assertEqual(screen.exact_mds_projective_support(2, 2, 4), 0)

    def test_regional_one_coefficient_matches_closed_form(self):
        for upgrade in (False, True):
            exact_coefficient = screen.regional_bound(1, .9998,
                                                      one_input_upgrade=upgrade)
            closed = screen.one_row_regional_bound(.9998, one_input_upgrade=upgrade)
            self.assertAlmostEqual(exact_coefficient, closed, places=7)


if __name__ == '__main__':
    unittest.main()
