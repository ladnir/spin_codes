"""Finite algebra/composition checks; these are not full distance certificates."""

import itertools
import math
import unittest
from collections import Counter

import numpy as np

from batched_transfer import multiply
from one_route_rank import (local_operator, polynomial_product, regional_matrices,
                            regional_max_apply, regional_sum_matrix, return_patterns)


def brute_patterns(width=4, q=4):
    """Independent enumeration of the four lane roles, without XOR-DP."""
    mul = lambda a, b: multiply(a, b, 0b111)
    result = Counter()
    for roles in itertools.product(range(4), repeat=width):
        # 0: canceled active, 1: uncanceled active, 2: inactive zero,
        # 3: inactive nonzero. Every support path has exactly one encoding.
        canceled = [j for j, role in enumerate(roles) if role == 0]
        pending = [j for j, role in enumerate(roles) if role == 1]
        zeros = [j for j, role in enumerate(roles) if role == 2]
        if not canceled and not pending or len(zeros) > 1:
            continue
        f = [[1, 0], [0, 1]]
        for j in canceled:
            f[0][0] ^= 1
            f[0][1] ^= j
            f[1][0] ^= j
            f[1][1] ^= mul(j, j)
        if len(pending) >= 2:
            constraints = []
        elif pending:
            j = pending[0]
            constraints = [(mul(j, f[0][0]) ^ f[1][0], mul(j, f[0][1]) ^ f[1][1])]
        else:
            constraints = [tuple(f[0]), tuple(f[1])]
        constraints += [(1, j) for j in zeros]
        allowed = [(a, b) for a in range(q) for b in range(q) if (a or b)
                   and all(mul(c, a) ^ mul(d, b) == 0 for c, d in constraints)]
        if not allowed:
            continue
        rank = 0 if len(allowed) == q*q-1 else 1
        if rank:
            a, b = allowed[0]
            forced = canceled+[j for j, role in enumerate(roles) if role == 3]
            if any(a ^ mul(j, b) == 0 for j in forced):
                continue
            if len(pending) == 1:
                j = pending[0]
                x = mul(f[0][0], a) ^ mul(f[0][1], b)
                if x == 0 or x ^ a ^ mul(j, b) == 0:
                    continue
        key = (len(canceled)+len(pending), width-len(canceled)-len(zeros),
               len(canceled)+min(2, len(pending)), rank, width-len(pending)-len(zeros))
        result[key] += 1
    return result


class OneRouteRankTests(unittest.TestCase):
    def test_return_pattern_aggregation(self):
        expected = brute_patterns()
        actual = Counter({item[:-1]: item[-1] for item in return_patterns(4, 4)})
        self.assertEqual(actual, expected)

    def test_empty_input_and_singleton_return(self):
        z, tau, q, width = 0.83, 7, 256, 16
        fixed, ambiguous = local_operator(z, tau, width, q)
        expected = np.array([[1, 0, 0],
                             [0, (1-width/(q+1))*z**width+width/(q+1)*z**(width-1), 0],
                             [0, width/tau*z**(width-1), z**width]])
        np.testing.assert_allclose(fixed[0], expected, rtol=1e-14, atol=0)
        self.assertTrue(np.all(ambiguous[0] == 0))
        g = 1/(q-1)+(1-1/(q-1))/tau
        # Ignoring output weight, every surviving singleton return pays
        # one M-only projective-line check and one fresh diagonal pivot.
        fixed_at_one, _ = local_operator(1, tau, width, q)
        self.assertLessEqual(fixed_at_one[1, 1, 0], g/(q+1)*(1+1e-13))
        self.assertGreater(fixed_at_one[1, 1, 0], 0)

    def test_regional_normalization_against_direct_power(self):
        z, tau, rows, width, active = 0.91, 3, 8, 4, 2
        fixed, ambiguous = local_operator(z, tau, width)
        fixed[:, :, 1] += ambiguous
        fixed[:, :, 2] += ambiguous
        local = fixed[:active+1].copy()
        for w in range(len(local)):
            local[w] *= math.comb(width, w)
        direct = polynomial_product(local, local, active)
        for h in range(active+1):
            direct[h] /= math.comb(rows, h)
        normalized = regional_matrices(z, tau, active, rows, width)
        np.testing.assert_allclose(normalized, direct, rtol=1e-13, atol=1e-15)

    def test_maximum_continuation_is_below_sum(self):
        for continuation in [np.ones(3), np.array([0.2, 3.0, 0.7]), np.array([2.0, 0.0, 1.0])]:
            args = dict(z=0.91, tau=3, active_rows=2, rows=8, width=4)
            actual, scale = regional_max_apply(continuation, **args)
            matrix, sum_scale = regional_sum_matrix(**args)
            upper = matrix @ continuation * math.exp(sum_scale)
            self.assertTrue(np.all(actual*math.exp(scale) <= upper*(1+1e-13)))
        first, scale1 = regional_max_apply(np.ones(3), **args)
        second, scale2 = regional_max_apply(3*np.ones(3), **args)
        np.testing.assert_allclose(second*math.exp(scale2), 3*first*math.exp(scale1), rtol=1e-13)

    def test_reject_invalid_geometry(self):
        for width in [2, 6, 32]:
            with self.assertRaises(ValueError):
                local_operator(0.9, 2, width, q=16)
        with self.assertRaises(ValueError):
            local_operator(0.9, 0.5)


if __name__ == '__main__':
    unittest.main()
