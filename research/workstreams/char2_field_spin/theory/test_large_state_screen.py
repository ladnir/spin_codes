"""Independent algebra and tiny-field checks for scaled-feedback screens."""

from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np

import large_state_screen as model
import gf16_sizes
from gf16_certificate import Interval
from state_screen import birth_zero_count, zero_patterns


def gf4mul(a, b):
    product = 0
    for i in range(2):
        if (b >> i) & 1:
            product ^= a << i
    if product & 4:
        product ^= 7
    return product


def image(matrix, vector, multiply):
    return tuple(model.xor_sum(multiply(a, b) for a, b in zip(row, vector)) for row in matrix)


class LargeStateTests(unittest.TestCase):
    def test_gf16_field(self):
        for a in range(1, 16):
            self.assertEqual({model.mul(a, b) for b in range(1, 16)}, set(range(1, 16)))
        for a, b, c in itertools.product(range(16), repeat=3):
            self.assertEqual(model.mul(a, b ^ c), model.mul(a, b) ^ model.mul(a, c))

    def test_gram_and_feedback_exact_all_states(self):
        for state in range(1, 16):
            for gamma in range(1, 16):
                report = model.algebra(state, gamma)
                self.assertTrue(report['gram_anti_diagonal'])
                expected = gamma != 1 or state <= 8
                self.assertEqual(report['state_feedback_rank'] == state, expected)
                self.assertEqual(report['output_feedback_rank'] == 16, expected)

    def test_singletons_are_direct_return_outputs(self):
        for state in (4, 6, 8, 9, 10, 12, 15):
            a, _, _, b = model.matrices(state, 2)
            c = tuple(tuple(model.mul(2, a[j][i]) for j in range(16)) for i in range(state))
            weights = model.singleton_weights(state, 2)
            for j in range(16):
                x = tuple(int(i == j) for i in range(16))
                return_state = image(c, x, model.mul)
                y = tuple(u ^ v for u, v in zip(x, image(a, return_state, model.mul)))
                self.assertEqual(y, tuple(row[j] for row in b))
                self.assertEqual(sum(bool(v) for v in y), weights[j])

    def test_excluded_and_singular_parameters(self):
        for state, gamma in ((0, 2), (16, 2), (4, 0), (4, 16)):
            with self.assertRaises(ValueError):
                model.local_model(.9, state, gamma)
        with self.assertRaises(ValueError):
            model.local_model(.9, 9, 1)
        for z in (0, -1, 1.01, math.nan, math.inf):
            with self.assertRaises(ValueError):
                model.local_model(z, 10, 2)

    def test_preserved_gamma_one_envelope(self):
        for state in (1, 4, 6, 8):
            self.assertEqual(model.singleton_weights(state, 1), gf16_sizes.singleton_weights(state))
            for z in (.1, .9, .999, 1):
                np.testing.assert_array_equal(model.local_model(z, state, 1), gf16_sizes.local_model(z, state))

    def test_tiny_field_scaled_return_envelope(self):
        # Exhaust all inputs and all nonzero refreshed states for GF4, t4, s3.
        # Here CA != 0, and gamma=2 makes I+AC invertible.
        q, t, s, gamma, z = 4, 4, 3, 2, F(2, 3)
        a = tuple((1, j, gf4mul(j, j)) for j in range(q))
        c = tuple(tuple(gf4mul(gamma, a[j][i]) for j in range(q)) for i in range(s))
        emitted = [F(0) for _ in range(t+1)]
        returned = [F(0) for _ in range(t+1)]
        zero_birth = [0]*(t+1)
        input_counts = [0]*(t+1)
        singleton_weights = []
        states = [r for r in itertools.product(range(q), repeat=s) if any(r)]
        a_states = [(r, image(a, r, gf4mul)) for r in states]
        for x in itertools.product(range(q), repeat=t):
            weight = sum(bool(v) for v in x)
            input_counts[weight] += 1
            cx = image(c, x, gf4mul)
            zero_birth[weight] += not any(cx)
            for r, ar in a_states:
                y = tuple(u ^ v for u, v in zip(x, ar))
                marked = z**sum(bool(v) for v in y)
                emitted[weight] += marked
                if r == cx:
                    returned[weight] += marked
                    if weight == 1:
                        singleton_weights.append(sum(bool(v) for v in y))
        denominator = q**s-1
        patterns = zero_patterns(t, s, q)
        f = F(1, q-1)+F(q-2, q-1)*z
        for w in range(t+1):
            p0 = F(birth_zero_count(w, s, q), (q-1)**w)
            self.assertEqual(F(zero_birth[w], input_counts[w]), p0)
            g = F(0)
            for zeros, count in enumerate(patterns):
                for active_zeros in range(max(0, zeros-(t-w)), min(w, zeros)+1):
                    inactive_zeros = zeros-active_zeros
                    g += (F(count*math.comb(w, active_zeros)*math.comb(t-w, inactive_zeros), denominator)
                          *z**(t-w-inactive_zeros+active_zeros)*f**(w-active_zeros))
            self.assertEqual(emitted[w]/(input_counts[w]*denominator), g)
            j = returned[w]/(input_counts[w]*denominator)
            ret = (1-p0)/denominator
            self.assertGreaterEqual(j, ret*z**t)
            self.assertLessEqual(j, ret*z**max(0, t-s+1-w))
            inverse = sum(F(math.comb(t, h), denominator*(q-1)**max(0, h-t+w))
                          *z**(t-h)*(1-z)**h for h in range(t+1))
            self.assertLessEqual(j, inverse)
            if w == 1:
                self.assertEqual(j, sum(z**v for v in singleton_weights)/(input_counts[1]*denominator))

    def test_directed_envelope_matches_fraction_formula(self):
        # Fractions are an independent exact arithmetic oracle for each entry.
        z = F(3, 5)
        for state in (8, 9, 10, 12):
            envelope = model.interval_local_envelope(Interval.rational(z.numerator, z.denominator), state, 2)
            f, denominator = F(1, 15)+F(14, 15)*z, 16**state-1
            patterns = zero_patterns(16, state, 16)
            for w in range(17):
                p0 = F(birth_zero_count(w, state, 16), 15**w)
                g = F(0)
                for zeros, count in enumerate(patterns):
                    for az in range(max(0, zeros-(16-w)), min(w, zeros)+1):
                        iz = zeros-az
                        g += F(count*math.comb(w, az)*math.comb(16-w, iz), denominator)*z**(16-w-iz+az)*f**(w-az)
                ret = (1-p0)/denominator
                jlo, jhi = ret*z**16, ret*z**max(0, 17-state-w)
                if w == 1:
                    jlo = jhi = sum(z**v for v in model.singleton_weights(state, 2))/(16*denominator)
                elif w >= 2:
                    inverse = sum(F(math.comb(16, h), denominator*15**max(0, h-16+w))
                                  *z**(16-h)*(1-z)**h for h in range(17))
                    jhi = min(jhi, inverse)
                values = [[z**w*p0, z**w*(1-p0)], [jhi, g-jlo]]
                for i in range(2):
                    for j in range(2):
                        self.assertLessEqual(F(envelope[w][i][j].lo), values[i][j])
                        self.assertGreaterEqual(F(envelope[w][i][j].hi), values[i][j])

    def test_preserved_whole_class_screens(self):
        for active in (1, 2, 4, 8):
            self.assertAlmostEqual(model.screen_bound(18, active, .999, 6, 1),
                                   gf16_sizes.screen_bound(18, active, .999, state=6, source_rho=F(871, 1000)), places=9)
        self.assertAlmostEqual(model.middle_screen(18, 6, 256, .002, .9, gamma=1),
                               gf16_sizes.middle_screen(18, 6, 256, .002, .9), places=9)


if __name__ == '__main__':
    unittest.main()
