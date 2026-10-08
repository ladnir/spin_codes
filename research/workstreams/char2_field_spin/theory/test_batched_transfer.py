"""Exhaustive small-field tests for the proposed batched local transfer."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import unittest

from batched_transfer import (active_output_moment, active_return_moment,
                              multiply, nonzero_sum_count, transition,
                              zero_state_probability)


class BatchedTransferTests(unittest.TestCase):
    def test_nonzero_solution_counts(self):
        q, modulus = 4, 0b111
        for n in range(q + 1):
            for columns in combinations(range(q), n):
                counts = Counter()
                for values in product(range(1, q), repeat=n):
                    h0 = h1 = 0
                    for j, value in zip(columns, values):
                        h0 ^= value
                        h1 ^= multiply(j, value, modulus)
                    counts[h0, h1] += 1
                for target in product(range(q), repeat=2):
                    self.assertEqual(nonzero_sum_count(columns, target, q, modulus), counts[target])

    def test_exact_joint_moments(self):
        q, width, modulus = 4, 4, 0b111
        states = [(a, b) for a, b in product(range(q), repeat=2) if a or b]
        for w in range(width + 1):
            for support in combinations(range(width), w):
                output_counts, return_counts = Counter(), Counter()
                for values in product(range(1, q), repeat=w):
                    x = [0] * width
                    for j, value in zip(support, values):
                        x[j] = value
                    c0 = c1 = 0
                    for j in range(width):
                        c0 ^= x[j]
                        c1 ^= multiply(j, x[j], modulus)
                    for r0, r1 in states:
                        weight = sum((x[j] ^ r0 ^ multiply(j, r1, modulus)) != 0
                                     for j in range(width))
                        output_counts[weight] += 1
                        if (c0, c1) == (r0, r1):
                            return_counts[weight] += 1
                denominator = (q*q - 1) * (q - 1)**w
                for z in [Fraction(0), Fraction(1, 3), Fraction(2, 3), Fraction(1)]:
                    direct_g = sum(count*z**weight for weight, count in output_counts.items()) / denominator
                    direct_j = sum(count*z**weight for weight, count in return_counts.items()) / denominator
                    self.assertEqual(active_output_moment(w, width, q, z), direct_g)
                    self.assertEqual(active_return_moment(support, width, q, modulus, z), direct_j)
                    matrix = transition(support, width, q, modulus, z)
                    self.assertEqual(sum(matrix[1]), direct_g)
                    if z == 1:
                        self.assertEqual(direct_j, (1-zero_state_probability(w, q))/(q*q-1))
                        self.assertEqual(tuple(map(sum, matrix)), (1, 1))

    def test_ca_zero(self):
        for width, modulus in [(4, 0b111), (8, 0b1011), (16, 0b10011)]:
            sum_j = sum_j_squared = 0
            for j in range(width):
                sum_j ^= j
                sum_j_squared ^= multiply(j, j, modulus)
            self.assertEqual((width & 1, sum_j, sum_j_squared), (0, 0, 0))


if __name__ == '__main__':
    unittest.main()
