"""Tiny exact/source-independent oracles. No files, bank, FFT or full-weight job."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
import unittest

import numpy as np
import hd_field128_input_weight_tilt_v39 as pilot
import hd_field128_full_rs_gamma_screen as gamma
from hd_blockwise_maps import mul


def xor_sum(values):
    result = 0
    for x in values:
        result ^= x
    return result


def emission(matrix, vector, multiply):
    return tuple(xor_sum(multiply(a, b) for a, b in zip(row, vector)) for row in matrix)


def mul64(a, b):
    # GF8[g]/(g^2+g+1), irreducible because 3 does not divide 7.
    a0, a1, b0, b1 = a & 7, a >> 3, b & 7, b >> 3
    high = mul(a1, b1, 8)
    return (mul(a0, b0, 8)^high) | ((mul(a0, b1, 8)^mul(a1, b0, 8)^high) << 3)


class Tests(unittest.TestCase):
    def test_exact_support_counts(self):
        # Independent polynomial evaluation, not the Lagrange/support formula.
        q, n = 8, 4
        observed = {}
        for a, b in product(range(q), repeat=2):
            if not (a or b):
                continue
            word = tuple(a^mul(x, b, q) for x in range(n))
            support = tuple(i for i, value in enumerate(word) if value)
            observed[support] = observed.get(support, 0)+1
        counts = pilot.support_counts(q, n, 2)
        for h in range(n+1):
            for support in combinations(range(n), h):
                self.assertEqual(observed.get(support, 0), counts[h])
        self.assertEqual(sum(comb(n, h)*counts[h] for h in range(n+1)), q**2-1)
        self.assertEqual(counts[:3], (0, 0, 0))

    def test_live_explicit_mds_enumeration_and_eta1(self):
        q, n, z = 8, 4, F(1, 5)
        a = tuple((1, i) for i in range(n))
        words = [emission(a, state, lambda x, y: mul(x, y, q))
                 for state in product(range(q), repeat=2) if any(state)]
        for eta in (F(1), F(3, 2), F(3)):
            for support in ((), (0,), (0, 1), (0, 2), (0, 1, 2, 3)):
                total = F(0)
                for word in words:
                    # Independent coordinate-wise average of every X value.
                    term = F(1)
                    for i in range(n):
                        if i in support:
                            term *= sum((eta**bool(x)*z**bool(x^word[i]) for x in range(q)), F(0))/q
                        else:
                            term *= z**bool(word[i])
                    total += term
                total /= q**2-1
                j = len(support)
                self.assertEqual(total, pilot.exact_live(j, z, eta, q=q, n=n, state=2))
                if eta == 1:
                    self.assertEqual(total, gamma.frozen.exact_moments(j, z, q, n, 2, False)[3])

    def test_return_cap_and_quiet_birth_toy_gamma(self):
        # This is the independently checked degree2 gamma/MDS4x2 toy of the
        # frozen gamma tests. Q=gamma*C*X is the unique returning input state.
        a = gamma.frozen.lagrange(range(4), range(4, 6), 8)
        c = tuple(zip(*a))
        q, z, eta, d = 64, F(1, 5), F(3), 64**2-1
        density = F(d, (q-1)**2)
        for support in ((), (0,), (0, 1), (0, 2)):
            quiet = birth = ret = F(0)
            for values in product(range(q), repeat=len(support)):
                x = [0]*4
                for i, value in zip(support, values):
                    x[i] = value
                syndrome = emission(c, x, mul64)
                input_weight = sum(bool(v) for v in x)
                mark = eta**input_weight
                if any(syndrome):
                    birth += mark*z**input_weight
                    incoming = tuple(mul64(8, v) for v in syndrome)
                    emitted = emission(a, incoming, mul64)
                    ret += mark*z**sum(bool(x[i]^emitted[i]) for i in range(4))
                else:
                    quiet += mark*z**input_weight
            mass = q**len(support)
            cap = pilot.exact_packet(len(support), z, eta, q=q, n=4, state=2)
            self.assertEqual(cap[0], quiet/mass)
            self.assertEqual(cap[1], density*birth/mass)
            self.assertGreaterEqual(cap[2], ret/(mass*d))

    def test_logs_match_exact_and_control(self):
        for eta in (F(1), F(3, 2), F(3)):
            z = F(1, 5)
            logs = pilot.tilted_logs(float(z), float(eta), q=8, n=4, state=2)
            live = pilot.log_live(float(z), float(eta), q=8, n=4, state=2)
            for j in range(5):
                expected = pilot.exact_packet(j, z, eta, q=8, n=4, state=2)
                np.testing.assert_allclose(logs[j].ravel(), [gamma.old.log_fraction(v) for v in expected],
                                           rtol=0, atol=2e-13)
                self.assertAlmostEqual(live[j], log2(pilot.exact_live(j, z, eta, q=8, n=4, state=2)), places=12)
            original = gamma.local_logs(float(z), 8, 4, 2)
            np.testing.assert_array_equal(logs[0], original[0])
            if eta == 1:
                np.testing.assert_array_equal(logs, original)

    def test_prefix_and_guards(self):
        self.assertEqual(pilot.prefix_correction(.25, 16., 5), -165.)
        self.assertEqual(pilot.prefix_correction(.25, 1., 5), 0.)
        pilot.validate_options(pilot.DEFAULT_ETAS, 110.)
        for exponents in ((1,), (0, 0), (0, 65), (0, True), tuple(range(9))):
            with self.assertRaises(ValueError):
                pilot.validate_options(exponents, 110.)
        for z, eta in ((F(1, 2), F(2)), (F(1, 2), F(1, 2)), (F(0), F(1))):
            with self.assertRaises(ValueError):
                pilot.exact_live(1, z, eta, q=8, n=4, state=2)
        with self.assertRaises(ValueError):
            pilot.exact_live(1, .25, F(2), q=8, n=4, state=2)


if __name__ == '__main__':
    unittest.main()
