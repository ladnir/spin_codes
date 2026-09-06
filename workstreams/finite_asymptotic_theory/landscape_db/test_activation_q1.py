"""Independent positive-arithmetic and exact GF(4) checks for the Q1 screen."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np

import activation_q1 as q1


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def exact_envelope(z):
    # A generators 1111 and 1010; spectrum {2:2,4:1}; B=A^T.
    m0 = (2 * z**2 + z**4) / 3
    m1 = (4 * z + 4 * z**3 + 4 * z**3) / 12
    return [[F(1), F(0), F(0)], [F(0), F(0), z**2], [F(0), F(0), F(3, 2) * m0]], [
        [F(0), z, F(0)], [z / 3, F(0), z], [m1 / 2, F(0), F(3, 2) * m1]]


def positive_regions(zero, one, epochs):
    rz = identity(len(zero))
    ra = [[F(0)] * len(zero) for _ in zero]
    for _ in range(epochs):
        ra = add(mul(ra, zero), mul(rz, one))
        rz = mul(rz, zero)
    return rz, [[x / epochs for x in row] for row in ra]


def exact_coefficients(zero, one, length):
    answer = []
    for weight in range(length + 1):
        total = F(0)
        for selected in itertools.combinations(range(length), weight):
            product = identity(len(zero))
            for position in range(length):
                product = mul(product, one if position in selected else zero)
            total += sum(product[0])
        answer.append(total / math.comb(length, weight))
    return answer


def gf4_mul(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def actual_epochs(z):
    maps = [0, 15, 10, 5]
    result = []
    for inputs in ([0], [1 << j for j in range(4)]):
        transfer = [[F(0)] * 4 for _ in range(4)]
        for state in range(4):
            for x in inputs:
                syndrome = ((x & 15).bit_count() % 2) | (((x & 10).bit_count() % 2) << 1)
                for alpha in (1, 2, 3):
                    destination = gf4_mul(alpha, state) ^ syndrome
                    transfer[state][destination] += z**((x ^ maps[state]).bit_count()) / (3 * len(inputs))
        result.append(transfer)
    return result


class ActivationQ1Test(unittest.TestCase):
    def test_batched_log_power_and_support_average_match_rationals(self):
        zs = [F(1, 2), F(3, 4)]
        for epochs in (1, 2, 3, 7):
            got = q1.coefficient_logs(*q1.region_logs(*q1.epoch_logs(
                4, 2, {2: 2, 4: 1}, np.array([-math.log(float(z)) for z in zs])), epochs), 4)
            for index, z in enumerate(zs):
                exact = exact_coefficients(*positive_regions(*exact_envelope(z), epochs), 4)
                np.testing.assert_allclose(got[index], [math.log(float(x)) for x in exact], atol=2e-13, rtol=2e-13)

    def test_envelope_dominates_actual_shared_setup_first_moments(self):
        for z in (F(1, 3), F(1, 2), F(3, 4)):
            for epochs in (1, 2, 3):
                actual = exact_coefficients(*positive_regions(*actual_epochs(z), epochs), 4)
                envelope = exact_coefficients(*positive_regions(*exact_envelope(z), epochs), 4)
                self.assertTrue(all(a <= b for a, b in zip(actual, envelope)))

    def test_non_native_length_rejected(self):
        with self.assertRaises(ValueError):
            q1.screen(64, 2, {32: 2, 64: 1}, 512, 256, 8, {32: 0.}, [-3.])


if __name__ == '__main__':
    unittest.main()
