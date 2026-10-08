"""Exact joint-rank and positive-mark checks; no global screening."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

from state8_joint_return import (joint_return_bound, joint_return_terms,
                                 minimum_return_weight, return_bounds,
                                 split_joint_return_terms)
from state_rank_cap_transfer import return_flat_factor


def mul16(a, b):
    value = 0
    while b:
        if b & 1:
            value ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19
    return value


M = tuple(tuple(mul16(a, b) for b in range(16)) for a in range(16))
INV = (0,)+tuple(next(x for x in range(1, 16) if M[a][x] == 1) for a in range(1, 16))


def matrix_rank(matrix):
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    found = 0
    for j in range(len(rows[0])):
        pivot = next((i for i in range(found, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[found], rows[pivot] = rows[pivot], rows[found]
        rows[found] = [M[INV[rows[found][j]]][x] for x in rows[found]]
        for i in range(found+1, len(rows)):
            factor = rows[i][j]
            rows[i] = [x ^ M[factor][y] for x, y in zip(rows[i], rows[found])]
        found += 1
        if found == len(rows):
            break
    return found


class JointReturnTests(unittest.TestCase):
    def test_terms_and_complete_cover_minimum(self):
        y, tau = F(2, 3), F(7)
        for width in (16, 32):
            for weight in range(1, width+1):
                pure, joint = return_bounds(weight, y, tau, width)
                self.assertEqual(pure, return_flat_factor(weight, tau, 8)
                                 *y**minimum_return_weight(weight, width)/tau**8)
                self.assertEqual(joint_return_bound(weight, y, tau, width), min(pure, joint))
                self.assertGreater(joint, 0)
                old_joint = return_flat_factor(weight, tau, 8)*sum(
                    count*y**v*(1-y)**(width-v)/tau**charge
                    for v, count, charge in joint_return_terms(weight, width))
                self.assertLessEqual(joint, old_joint)
                for v, count, charge in joint_return_terms(weight, width):
                    self.assertEqual(count, comb(width, v))
                    self.assertEqual(charge, max(8, 1+weight-v))
                for i, z, count, charge in split_joint_return_terms(weight, width):
                    self.assertEqual(count, comb(weight, i)*comb(width-weight, z))
                    self.assertLess(z, 8)
                    self.assertEqual(charge, max(8, i+z, weight-width+i+2*z+1))

    def test_zero_input_and_y_one(self):
        for width in (16, 32):
            self.assertEqual(joint_return_bound(0, F(1, 2), F(7), width), 0)
            for weight in range(1, width+1):
                pure, joint = return_bounds(weight, F(1), F(7), width)
                self.assertEqual(pure, joint)

    def test_selected_zero_identity_uses_unselected_not_exact_weight(self):
        for width in (4, 8):
            y = F(2, 5)
            for actual_weight in range(1, width+1):
                total = sum(comb(width-actual_weight, zeros)
                            *y**(width-zeros)*(1-y)**zeros
                            for zeros in range(width-actual_weight+1))
                self.assertEqual(total, y**actual_weight)

    def test_all_joint_masks_in_gf16_two_state_analogue(self):
        # Four labels in actual GF4 subset of GF16; two residue nodes form
        # an affine GF2 coset. A[u,v]=u+v+1, C=A^T, eta=2.
        labels, nodes = (0, 1, 6, 7), (2, 3)
        a = tuple(tuple(u ^ v ^ 1 for v in nodes) for u in labels)
        c = tuple(zip(*a))
        for i in range(2):
            for j in range(2):
                value = 0
                for u in range(4):
                    value ^= M[c[i][u]][a[u][j]]
                self.assertEqual(value, 0)
        checked = 0
        for hmask in range(1, 4):
            support = tuple(i for i in range(2) if hmask >> i & 1)
            h = len(support)
            for imask in range(1, 16):
                inputs = tuple(i for i in range(4) if imask >> i & 1)
                w = len(inputs)
                returns = [tuple(int(row == j) for j in support)
                           +tuple(c[row][i] for i in inputs) for row in range(2)]
                if w <= 2-h:
                    # The absent incoming-state coordinates force X=0.
                    absent = [[c[j][i] for i in inputs] for j in range(2) if j not in support]
                    self.assertEqual(matrix_rank(absent), w)
                    continue
                self.assertEqual(matrix_rank(returns), 2)
                for zmask in range(16):
                    zeros = tuple(j for j in range(4) if zmask >> j & 1)
                    equations = returns+[tuple(M[2][a[j][k]] for k in support)
                                        +tuple(int(j == i) for i in inputs) for j in zeros]
                    v = 4-len(zeros)
                    self.assertGreaterEqual(matrix_rank(equations), max(2, h+w-v))
                    inactive_zeros = tuple(j for j in zeros if j not in inputs)
                    z, i = len(inactive_zeros), len(zeros)-len(inactive_zeros)
                    if z >= h:
                        self.assertEqual(matrix_rank([[a[j][k] for k in support]
                                                      for j in inactive_zeros]), h)
                    else:
                        self.assertGreaterEqual(matrix_rank(equations),
                                                max(2, i+z, w-4+i+2*z+1))
                    checked += 1
        self.assertEqual(checked, 592)

    def test_strict_improvement_for_dense_selected_zero_return(self):
        pure, joint = return_bounds(24, F(1, 5), F(16), 32)
        self.assertLess(joint, pure/F(1000))

    def test_flat_rank_cutoff_tightens_with_inactive_zeros(self):
        # If z=7, h=8 and the only feasible old-state zero flat is empty.
        # The formula must not retain an artificial (1+1/tau)^8 multiplier.
        terms = split_joint_return_terms(8, 16)
        self.assertTrue(any(z == 7 for i, z, count, charge in terms))
        y, tau = F(1, 3), F(7)
        expected = sum(sum(F(comb(8, k), tau**k)
                           for k in range(min(7, 7-z)+1))
                       *count*y**(16-i-z)*(1-y)**(i+z)/tau**charge
                       for i, z, count, charge in terms)
        self.assertEqual(return_bounds(8, y, tau, 16)[1], expected)


if __name__ == '__main__':
    unittest.main()
