"""Tiny independent fixtures for the local sandwich and scalar update law."""
from itertools import combinations, product
import unittest

import numpy as np

import local_slack as slack


def mul4(a, b):
    result = 0
    for _ in range(2):
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def feedback(values):
    result = [0, 0, 0]
    for h, x in enumerate(values):
        result[0] ^= x
        result[1] ^= mul4(h, x)
        result[2] ^= mul4(mul4(h, h), x)
    return result[0] | result[1] << 2 | result[2] << 4


def mul64(a, b):
    aa, bb = [(a >> (2*i)) & 3 for i in range(3)], [(b >> (2*i)) & 3 for i in range(3)]
    c = [0]*5
    for i in range(3):
        for j in range(3):
            c[i+j] ^= mul4(aa[i], bb[j])
    for i in (4, 3):
        c[i-2] ^= c[i]
        c[i-3] ^= c[i]
    return c[0] | c[1] << 2 | c[2] << 4


class LocalSlackTests(unittest.TestCase):
    def test_exact_integer_margin_checks(self):
        value = slack.analytic_bounds()
        self.assertTrue(value['strict_1_5_bit_bound_exact_integer_check'])
        self.assertTrue(value['return_gain_below_1e_minus5_exact_check'])

    def test_potential_and_return_sandwich_in_gf4_fixture(self):
        z, offsets, states = .9, [1, 2, 3, 1], 64
        delta = (1+z)**-2
        for j in range(5):
            for support in combinations(range(4), j):
                syndrome = np.zeros(states)
                zero_input = 0.
                for labels in product(range(4), repeat=j):
                    x = [0]*4
                    for h, value in zip(support, labels):
                        x[h] = value
                    mass = z**sum((a ^ b).bit_count() for a, b in zip(x, offsets))/4**j
                    syndrome[feedback(x)] += mass
                    if not any(x):
                        zero_input += mass
                total = syndrome.sum()
                if j == 0:
                    self.assertAlmostEqual(total, zero_input)
                    continue
                self.assertLessEqual(syndrome.max(), delta**min(j, 3)*total+2e-14)
                exact_return = (total-syndrome[0])/(states-1)
                filled_return = (total-zero_input)/(states-1)
                if j <= 3:
                    self.assertAlmostEqual(syndrome[0], zero_input, places=14)
                    self.assertAlmostEqual(exact_return, filled_return, places=14)
                else:
                    self.assertGreaterEqual(exact_return+2e-14, (1-delta**3)*filled_return)
                self.assertTrue(np.all((total-syndrome[1:])/(states-1)+2e-14 >=
                    (1-delta**min(j, 3))*total/(states-1)))

    def test_shadow_preserves_births_and_zero_input(self):
        source = np.ones((9, 10, 10))
        source[0, 1:, 0] = 0.
        modified = slack.shadow_potential(source, .94)
        np.testing.assert_array_equal(modified[0], source[0])
        np.testing.assert_array_equal(modified[:, 0], source[:, 0])
        np.testing.assert_array_equal(modified[1:4, 1:, 0], source[1:4, 1:, 0])
        self.assertTrue(np.all(modified <= source))

    def test_nonzero_scalar_and_binary_adjoint_transitivity(self):
        columns = [[mul64(c, 1 << bit) for bit in range(6)] for c in range(1, 64)]
        for a in range(1, 64):
            self.assertEqual({mul64(c, a) for c in range(1, 64)}, set(range(1, 64)))
            images = {sum(((a & col).bit_count() & 1) << bit for bit, col in enumerate(matrix))
                for matrix in columns}
            self.assertEqual(images, set(range(1, 64)))


if __name__ == '__main__':
    unittest.main()
