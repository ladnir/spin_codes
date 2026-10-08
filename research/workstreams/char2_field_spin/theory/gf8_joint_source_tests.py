"""Exact toy checks for active-row/weight conditioning."""
from fractions import Fraction as F
from itertools import product
import math
import unittest

import numpy as np
from scipy.special import logsumexp

from gf8_joint_source_screen import (row_logs, source_joint_log, activation_mark_log,
                                    feasible, dominant_active, sparse_bound)
from gf8_rank_history_screen import source_row_log
from gf8_mds_geometry_screen import sparse_bound as previous_sparse


def expectation(coefficients, probabilities):
    answer = F(0)
    for mask in product((0, 1), repeat=len(probabilities)):
        mass = F(1)
        for bit, p in zip(mask, probabilities):
            mass *= p if bit else 1-p
        answer += mass*coefficients[sum(mask)]
    return answer


class JointSourceTests(unittest.TestCase):
    def test_pairwise_homogenization_exact(self):
        for beta, parity in product((F(2), F(7, 3), F(7)), (1, 2, 3)):
            empty = [F(1)]+[F(0)]*5
            occupied = [F(0) if h <= parity else beta**h/F(3)**parity for h in range(6)]
            for probabilities in ((F(0),F(1),F(1,4),F(3,4),F(1,2)),
                                  (F(1,9),F(2,9),F(3,9),F(4,9),F(5,9))):
                average = sum(probabilities)/5
                e, o = expectation(empty, probabilities), expectation(occupied, probabilities)
                eh, oh = expectation(empty, [average]*5), expectation(occupied, [average]*5)
                self.assertLessEqual(e, eh)
                self.assertLessEqual(o, oh)
                for active in range(5):
                    self.assertLessEqual(e**(4-active)*o**active, eh**(4-active)*oh**active)

    def test_activation_polynomial_and_old_source_exact_decomposition(self):
        for alpha, tau, mark in ((.13, 3., 2.), (.7, 1.1, .3)):
            for lam in (.2, 2.):
                terms = [source_joint_log(alpha, a, 8, math.log(tau), math.log(lam), 5, 3)
                         +a*math.log(mark) for a in range(9)]
                self.assertAlmostEqual(float(logsumexp(terms)),
                                       activation_mark_log(alpha, math.log(mark), 8, math.log(tau), math.log(lam), 5, 3))
                empty, occupied = row_logs(alpha, math.log(tau), math.log(lam), 5, 3)
                self.assertAlmostEqual(float(np.logaddexp(empty, occupied)),
                                       source_row_log(alpha, math.log(tau), math.log(lam), 5, 3, True))

    def test_toy_exact_column_conditioning(self):
        rows, length, parity, tau, beta = 2, 4, 1, F(3), F(3)
        groups = {}
        for bits in product((0,1), repeat=rows*length):
            weights = [sum(bits[row*length:(row+1)*length]) for row in range(rows)]
            if any(0 < w <= parity for w in weights):
                continue
            active = sum(w > 0 for w in weights)
            columns = tuple(sum(bits[row*length+j] for row in range(rows)) for j in range(length))
            weight = sum(weights)
            groups[active, columns] = groups.get((active, columns), F(0))+tau**(weight-parity*active)
        modal = F(2**rows, math.comb(rows, rows//2))
        for (active, columns), weighted_count in groups.items():
            weight = sum(columns)
            if not 0 < weight < rows*length:
                continue
            actual = weighted_count
            for h in columns:
                actual /= math.comb(rows, h)
            alpha = F(weight, rows*length)
            empty = (1-alpha)**length
            occupied = sum(math.comb(length,h)*(beta*alpha)**h*(1-alpha)**(length-h)
                           for h in range(parity+1,length+1))/tau**parity
            bound = modal**length*math.comb(rows,active)*empty**(rows-active)*occupied**active
            self.assertLessEqual(actual, bound)

    def test_sparse_weight_mark_neutral_recovers_old_bound(self):
        point = (6.1, math.log(8))
        old = previous_sparse(2, point, 8, 4, 8, 4, 8)
        for weight in (10, 12, 16):
            self.assertAlmostEqual(sparse_bound(point+(0.,), 2, weight, 8, 4, 8, 4, 8), old)

    def test_feasible_classes_and_old_middle_mode(self):
        self.assertTrue(feasible(2, 10, 8, 8, 4))
        self.assertFalse(feasible(2, 9, 8, 8, 4))
        self.assertFalse(feasible(2, 17, 8, 8, 4))
        active = dominant_active([-.5214550554611819, -.5861651606135366,
                                  1.831738468168561, .5302681560131963], 70042)
        self.assertTrue(feasible(active, 70042, 10944))


if __name__ == '__main__':
    unittest.main()
