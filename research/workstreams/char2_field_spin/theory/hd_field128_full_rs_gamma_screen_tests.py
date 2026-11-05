"""Tiny gamma algebra/local/coefficient tests; no K20 distance gate."""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
import unittest

import numpy as np

import hd_field128_full_rs_gamma_screen as gamma
from hd_blockwise_maps import matmul, mul


def subsets(values):
    values = tuple(values)
    return [item for size in range(len(values)+1) for item in combinations(values, size)]


def mul64(a, b):
    """GF8[g]/(g²+g+1): irreducible since 3 does not divide |GF8*|=7."""
    a0, a1, b0, b1 = a&7, a>>3, b&7, b>>3
    high = mul(a1, b1, 8)
    return (mul(a0, b0, 8)^high) | ((mul(a0, b1, 8)^mul(a1, b0, 8)^high)<<3)


def power64(a, exponent):
    result = 1
    while exponent:
        if exponent&1:
            result = mul64(result, a)
        a, exponent = mul64(a, a), exponent>>1
    return result


def sum64(values):
    result = 0
    for value in values:
        result ^= value
    return result


def rank64(matrix):
    a = [list(row) for row in matrix]
    rank = 0
    for col in range(len(a[0]) if a else 0):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inverse = power64(a[rank][col], 62)
        a[rank] = [mul64(x, inverse) for x in a[rank]]
        for r in range(len(a)):
            if r != rank:
                scale = a[r][col]
                a[r] = [x^mul64(scale, y) for x, y in zip(a[r], a[rank])]
        rank += 1
    return rank


class GammaTests(unittest.TestCase):
    def test_actual_new_maps_ca_zero_and_unit_gram_diagonal(self):
        for state in gamma.STATES:
            a, c, h = gamma.fixed_maps(state)
            self.assertFalse(any(x for row in matmul(c, a, 256) for x in row))
            self.assertTrue(all(sum64(row) == 1 for row in a))
            self.assertTrue(all(sum64(mul(x, x) for x in row) == 1 for row in a))
            self.assertEqual((len(h), len(h[0])), (state, state))
        with self.assertRaises(ValueError):
            gamma.construction(32)
        with self.assertRaises(ValueError):
            gamma.construction(24, 256)

    def test_toy_gamma_principal_minors_and_joint_zero_counts(self):
        a = gamma.frozen.lagrange(range(4), range(4, 6), 8)
        self.assertFalse(any(x for row in matmul(tuple(zip(*a)), a, 8) for x in row))
        self.assertNotEqual(power64(8, 8), 8)  # gamma has degree2 over GF8.
        gram = matmul(a, tuple(zip(*a)), 8)
        p = [[int(i == j)^mul64(8, gram[i][j]) for j in range(4)] for i in range(4)]
        for d in subsets(range(4)):
            self.assertEqual(rank64([[p[i][j] for j in d] for i in d]), len(d))
        for support in ((), (0,), (0, 1), (0, 2)):
            outside = tuple(i for i in range(4) if i not in support)
            joint = {(d, e): 0 for d in subsets(support) for e in subsets(outside)}
            actual_return = F(0)
            for values in product(range(64), repeat=len(support)):
                x = [0]*4
                for i, value in zip(support, values):
                    x[i] = value
                syndrome = [sum64(mul64(a[i][k], x[i]) for i in range(4)) for k in range(2)]
                if not any(syndrome):
                    continue
                y = [x[i]^mul64(8, sum64(mul64(a[i][k], syndrome[k]) for k in range(2))) for i in range(4)]
                actual_return += F(1, 2)**sum(bool(v) for v in y)
                for d, e in joint:
                    if all(y[i] == 0 for i in d+e):
                        joint[d, e] += 1
            for (d, e), count in joint.items():
                cap = gamma.joint_probability(len(support), len(d), len(e), 64, 2)
                self.assertGreaterEqual(cap, F(count, 64**len(support)))
            cap = gamma.exact_moments(len(support), F(1, 2), 64, 4, 2)[2]
            self.assertGreaterEqual(cap, actual_return/64**len(support))

    def test_joint_cap_kills_all_active_zeros_and_strengthens_old_count(self):
        q, state = 256, 24
        for j in (1, 16, 24, 32, 64):
            self.assertEqual(gamma.joint_probability(j, j, 0, q, state), 0)
            for c in range(j+1):
                for e in (0, 1, 12, 23, 24):
                    value = gamma.joint_probability(j, c, e, q, state)
                    self.assertLessEqual(value, F(q**(j-c)-1, q**j))
                    if e < state:
                        old = F((q**min(min(state, j), state-e)-1)*q**max(0, j-c-state), q**j)
                        self.assertLessEqual(value, old)

    def test_only_return_changes_and_positive_logs_match_exact(self):
        for state in (16, 24, 30):
            for z in (F(1, 2), F(1, 1<<200)):
                logs = gamma.local_logs(float(z), state=state)
                baseline = gamma.frozen.local_logs(float(z), state=state)
                np.testing.assert_array_equal(logs[:, 0, :], baseline[:, 0, :])
                np.testing.assert_array_equal(logs[:, 1, 1], baseline[:, 1, 1])
                self.assertTrue(np.all(logs[:, 1, 0] <= baseline[:, 1, 0]+1e-9))
                for j in (0, 1, 24, 64):
                    exact = gamma.exact_packet(j, z, state=state)
                    np.testing.assert_allclose(logs[j].ravel(), [gamma.old.log_fraction(v) for v in exact], atol=1e-9, rtol=0)
                self.assertTrue(np.all(np.isfinite(logs[1])))

    def test_alias_layout_counts_and_raw_composition_before_alpha(self):
        z, alpha = .5, .25
        coefficients = gamma.AliasCoefficients(z, alpha, state=24)
        logs = gamma.local_logs(z, state=24)
        self.assertEqual(coefficients.single.shape, (4, 33, 33))
        self.assertEqual(coefficients.composed.shape, (4, 3, 65))
        np.testing.assert_allclose(coefficients.single[:, 1, 0], 5+alpha*logs[1].ravel())
        terms = []
        for left, right in ((logs[0], logs[1]), (logs[1], logs[0])):
            raw = np.logaddexp2(left[:, 0, None]+right[None, 0, :], left[:, 1, None]+right[None, 1, :])
            terms.append(5+alpha*raw.ravel())
        np.testing.assert_allclose(coefficients.composed[:, 1, 0], np.logaddexp2(*terms))
        self.assertEqual(coefficients.pair_logs((0., 0.)).shape, (2, 2))

    def test_small_alias_atlas_interface_no_weight_gate(self):
        bank = gamma.AliasBankScreen(.5, .25, ((-2., -3.),), 0., state=24, grid=64, max_seconds=20)
        witness = bank.witness()
        self.assertEqual(witness['construction'], gamma.construction(24))
        self.assertEqual(witness['feedback'], 'gamma A transpose')
        self.assertFalse(witness['directed'])
        self.assertFalse(witness['whole_code_certificate'])
        self.assertEqual(bank.regional(2).shape, (2, 2, 2))
        self.assertIn('hd_field128_full_rs_gamma_screen.py', gamma.source_pins())


if __name__ == '__main__':
    unittest.main()
