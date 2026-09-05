#!/usr/bin/env python3

import math
import itertools
import unittest

from prime_field_ea_diagnostic import (
    gv_distance,
    regular_shell_distribution,
    shell_entropy,
    sparse_rate,
    spectral_radius,
    tail_rate,
)


class PrimeFieldEADiagnosticTests(unittest.TestCase):
    def test_binary_gv_root(self) -> None:
        self.assertAlmostEqual(gv_distance(2, 0.5), 0.1100278644383603, places=13)

    def test_sparse_rate_reduces_to_binary_formula(self) -> None:
        delta = 0.1
        expected = 1.0 - 2.0 * math.sqrt(delta * (1.0 - delta))
        self.assertAlmostEqual(sparse_rate(2, delta), expected, places=14)

    def test_dense_spectral_radius(self) -> None:
        for prime in (2, 3, 5, 17):
            s = prime - 1.0
            q = s / prime
            for z in (0.1, 0.4, 1.0):
                self.assertAlmostEqual(
                    spectral_radius(prime, q, z), (1.0 + s * z) / prime
                )

    def test_dense_tail_rate_is_qary_ball_exponent(self) -> None:
        for prime in (2, 3, 5, 7):
            delta = 0.1
            expected = math.log(prime) - shell_entropy(prime, delta)
            self.assertAlmostEqual(
                tail_rate(prime, delta, (prime - 1.0) / prime),
                expected,
                places=10,
            )

    def test_sparse_numerical_limit(self) -> None:
        for prime in (2, 3, 5):
            delta = 0.1
            q = 1e-5
            observed = tail_rate(prime, delta, q) / q
            self.assertAlmostEqual(observed, sparse_rate(prime, delta), delta=2e-4)

    def test_regular_shell_recurrence_against_enumeration(self) -> None:
        prime = 3
        length = 3
        steps = 3
        counts = [0] * (length + 1)
        choices = list(itertools.product(range(length), range(1, prime)))
        for walk in itertools.product(choices, repeat=steps):
            word = [0] * length
            for coordinate, label in walk:
                word[coordinate] = (word[coordinate] + label) % prime
            counts[sum(value != 0 for value in word)] += 1
        total = len(choices) ** steps
        exact = [count / total for count in counts]
        recurrence = regular_shell_distribution(prime, length, steps)
        for observed, expected in zip(recurrence, exact):
            self.assertAlmostEqual(observed, expected, places=15)


if __name__ == "__main__":
    unittest.main()
