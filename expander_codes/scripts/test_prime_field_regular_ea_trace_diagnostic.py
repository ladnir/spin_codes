#!/usr/bin/env python3

import itertools
import math
import unittest

import numpy as np

from prime_field_regular_ea_trace_diagnostic import (
    balanced_projective_trace_logterm,
    balanced_occupancy_size_distribution,
    dominated_trace_matrix,
    empty_set_weight,
    exact_projective_trace_logterm,
    exact_balanced_projective_trace_logterm,
    occupancy_size_distribution,
    poisson_projective_trace_logterm,
    regular_projective_trace_logterm,
    saddle_projective_trace_logterm,
    uniform_occupancy_trace_transfers,
)


class PrimeFieldRegularEATraceDiagnosticTests(unittest.TestCase):
    def test_empty_set_probability_domination(self) -> None:
        length = 4
        message_weight = 3
        beta = empty_set_weight(message_weight, length)
        endpoints = list(itertools.product(range(length), repeat=message_weight))
        for empty_size in range(length + 1):
            for empty_set in itertools.combinations(range(length), empty_size):
                probability = sum(
                    all(endpoint not in empty_set for endpoint in assignment)
                    for assignment in endpoints
                ) / len(endpoints)
                self.assertLessEqual(probability, beta ** empty_size + 1e-15)

    def test_dominated_matrix_is_positive(self) -> None:
        matrix = dominated_trace_matrix(0.2, 0.4, 0.3)
        self.assertTrue(np.all(matrix >= 0.0))

    def test_dense_term_is_finite(self) -> None:
        term = regular_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=20,
        )
        self.assertTrue(math.isfinite(term.total_log2))

    def test_poisson_term_is_finite(self) -> None:
        term = poisson_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=5,
        )
        self.assertTrue(math.isfinite(term.total_log2))

    def test_saddle_term_is_finite(self) -> None:
        term = saddle_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=5,
        )
        self.assertTrue(math.isfinite(term.total_log2))
        poisson = poisson_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=5,
        )
        self.assertLessEqual(term.total_log2, poisson.total_log2 + 1e-7)

    def test_balanced_term_is_finite(self) -> None:
        term = balanced_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=5,
        )
        self.assertTrue(math.isfinite(term.total_log2))
        full = balanced_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=30,
        )
        self.assertEqual(full.empty_set_weight, 0.0)
        self.assertTrue(math.isfinite(full.total_log2))

    def test_occupancy_size_distribution(self) -> None:
        distribution = occupancy_size_distribution(4, 2)
        np.testing.assert_allclose(distribution, (0.0, 0.25, 0.75))
        self.assertAlmostEqual(float(np.sum(distribution)), 1.0)

    def test_balanced_occupancy_size_distribution(self) -> None:
        distribution = balanced_occupancy_size_distribution(3, 2, 2)
        np.testing.assert_allclose(distribution, (0.0, 0.2, 0.8))
        self.assertAlmostEqual(float(np.sum(distribution)), 1.0)

    def test_uniform_occupancy_transfer_matches_enumeration(self) -> None:
        length = 3
        z = 0.4
        v = 0.3
        slices = uniform_occupancy_trace_transfers(
            region_length=length, max_occupied=2, z=z, v=v
        )
        empty = np.array(((1.0, 0.0), (0.0, z)))
        occupied = np.array(((v, z), (v, z)))
        brute = np.zeros((2, 2))
        for occupied_set in itertools.combinations(range(length), 2):
            matrix = np.eye(2)
            for coordinate in range(length):
                matrix = matrix @ (
                    occupied if coordinate in occupied_set else empty
                )
            brute += matrix / math.comb(length, 2)
        np.testing.assert_allclose(slices[2], brute)

    def test_exact_term_is_finite(self) -> None:
        term = exact_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=2,
        )
        self.assertTrue(math.isfinite(term.total_log2))

    def test_exact_balanced_term_is_finite(self) -> None:
        term = exact_balanced_projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            region_count=6,
            message_weight=2,
        )
        self.assertTrue(math.isfinite(term.total_log2))


if __name__ == "__main__":
    unittest.main()
