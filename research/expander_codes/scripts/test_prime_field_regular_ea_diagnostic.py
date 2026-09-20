#!/usr/bin/env python3

import itertools
import math
import unittest

import numpy as np

from prime_field_regular_ea_diagnostic import (
    accumulator_input_matrices,
    dense_logterm,
    log_projective_support_count,
    log_support_grouped_bound,
    uniform_shell_transfers,
)


class PrimeFieldRegularEADiagnosticTests(unittest.TestCase):
    def test_input_matrices_are_stochastic_at_one(self) -> None:
        for prime in (2, 3, 5):
            for matrix in accumulator_input_matrices(prime, 1.0):
                np.testing.assert_allclose(matrix.sum(axis=1), np.ones(2))

    def test_uniform_shell_transfer_against_enumeration(self) -> None:
        prime = 3
        length = 3
        z = 0.37
        transfers = uniform_shell_transfers(
            prime=prime, length=length, max_weight=2, output_marker=z
        )
        for input_weight in range(3):
            total = np.zeros((2, 2), dtype=float)
            words = []
            for word in itertools.product(range(prime), repeat=length):
                if sum(value != 0 for value in word) == input_weight:
                    words.append(word)
            for start_class in range(2):
                starts = [0] if start_class == 0 else [1]
                for word in words:
                    state = starts[0]
                    output_weight = 0
                    for value in word:
                        state = (state + value) % prime
                        output_weight += state != 0
                    total[start_class, int(state != 0)] += z ** output_weight
            total /= len(words)
            np.testing.assert_allclose(transfers[input_weight], total)

    def test_dense_term_is_finite(self) -> None:
        value = dense_logterm(
            prime=3, k=50, region_length=10, region_count=10,
            cutoff=10, message_weight=25,
        )
        self.assertTrue(np.isfinite(value))

    def test_support_grouping_caps_one_support_family(self) -> None:
        prime = 101
        k = 50
        weight = 3
        log_line_bound = -math.log(2.0)
        grouped = log_support_grouped_bound(
            prime, k, weight, log_line_bound
        )
        self.assertAlmostEqual(grouped, math.log(math.comb(k, weight)))
        ungrouped = (
            log_projective_support_count(prime, k, weight) + log_line_bound
        )
        self.assertLess(grouped, ungrouped)

    def test_support_grouping_matches_first_moment_below_cap(self) -> None:
        prime = 5
        k = 40
        weight = 2
        log_line_bound = -20.0
        grouped = log_support_grouped_bound(
            prime, k, weight, log_line_bound
        )
        ungrouped = (
            log_projective_support_count(prime, k, weight) + log_line_bound
        )
        self.assertAlmostEqual(grouped, ungrouped)


if __name__ == "__main__":
    unittest.main()
