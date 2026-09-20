#!/usr/bin/env python3

from __future__ import annotations

import itertools
import math
import unittest
from collections import Counter

from regular_ec_diagnostic import (
    combine_uniform_slice_transfers,
    dense_regular_logterm,
    exact_regular_ec_logterm,
    log_poisson_conditioning_penalty,
    regular_activation_probability,
    regional_shell_distribution,
    uniform_slice_transfer_matrices,
    wrapping_input_matrices,
)


def parity_word(samples: tuple[int, ...], length: int) -> tuple[int, ...]:
    word = [0] * length
    for sample in samples:
        word[sample] ^= 1
    return tuple(word)


class RegularECDiagnosticTests(unittest.TestCase):
    def test_activation_probability(self) -> None:
        q = regular_activation_probability(7, 13)
        self.assertAlmostEqual(q, (1.0 - math.exp(-14.0 / 13.0)) / 2.0)

    def test_poisson_penalty_matches_formula(self) -> None:
        for r in range(1, 20):
            expected = math.log(math.exp(r) * math.factorial(r) / (r**r))
            self.assertAlmostEqual(log_poisson_conditioning_penalty(r), expected)

    def test_conditioned_point_probabilities_are_dominated(self) -> None:
        length = 3
        rows = 2
        exact = Counter(
            parity_word(samples, length)
            for samples in itertools.product(range(length), repeat=rows)
        )
        q = regular_activation_probability(rows, length)
        gamma = math.exp(log_poisson_conditioning_penalty(rows))
        for word, count in exact.items():
            exact_probability = count / (length**rows)
            weight = sum(word)
            bernoulli_probability = q**weight * (1.0 - q) ** (length - weight)
            self.assertLessEqual(exact_probability, gamma * bernoulli_probability + 1e-15)

    def test_dense_bound_is_finite(self) -> None:
        value = dense_regular_logterm(
            k=16,
            region_length=8,
            region_count=4,
            cutoff=3,
            message_weight=8,
        )
        self.assertTrue(math.isfinite(value))

    def test_regional_shell_distribution(self) -> None:
        distribution = regional_shell_distribution(5, 2)
        self.assertAlmostEqual(distribution[0], 1.0 / 5.0)
        self.assertAlmostEqual(distribution[2], 4.0 / 5.0)

    def test_normalized_slice_transfer_matches_direct_convolution(self) -> None:
        import numpy as np

        zero, one = wrapping_input_matrices(0.7, 3)
        slices = uniform_slice_transfer_matrices(
            length=2, max_weight=2, output_marker=0.7, memory=3
        )
        np.testing.assert_allclose(slices[0], zero @ zero)
        np.testing.assert_allclose(slices[1], (zero @ one + one @ zero) / 2.0)
        np.testing.assert_allclose(slices[2], one @ one)

    def test_exact_small_term_is_finite(self) -> None:
        term = exact_regular_ec_logterm(
            k=16,
            region_length=8,
            region_count=4,
            cutoff=3,
            memory=3,
            message_weight=2,
        )
        self.assertTrue(math.isfinite(term.log2_bound))


if __name__ == "__main__":
    unittest.main()
