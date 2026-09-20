#!/usr/bin/env python3

import itertools
import math
import unittest

import numpy as np

from prime_field_ea_trace_diagnostic import (
    empty_probability,
    occupancy_probability,
    projective_trace_logterm,
    trace_matrix,
    trace_matrix_from_empty,
)


def brute_trace_mgf(n: int, q: float, z: float, v: float) -> float:
    total = 0.0
    for occupied in itertools.product((0, 1), repeat=n):
        event_probability = 1.0
        for event in occupied:
            event_probability *= q if event else 1.0 - q
        event_indices = [index for index, event in enumerate(occupied) if event]
        for states in itertools.product((0, 1), repeat=len(event_indices)):
            state = 0
            cursor = 0
            output_weight = 0
            zero_events = 0
            for index, event in enumerate(occupied):
                if event:
                    state = states[cursor]
                    zero_events += state == 0
                    cursor += 1
                output_weight += state
            total += event_probability * z ** output_weight * v ** zero_events
    return total


class PrimeFieldEATraceDiagnosticTests(unittest.TestCase):
    def test_occupancy_probability(self) -> None:
        self.assertAlmostEqual(occupancy_probability(0.2, 3), 1.0 - 0.8 ** 3)
        self.assertAlmostEqual(empty_probability(0.2, 3), 0.8 ** 3)

    def test_empty_and_occupied_kernels_agree(self) -> None:
        empty = 0.37
        np.testing.assert_allclose(
            trace_matrix(1.0 - empty, 0.43, 0.31),
            trace_matrix_from_empty(empty, 0.43, 0.31),
        )

    def test_trace_transfer_matches_enumeration(self) -> None:
        n = 4
        q = 0.27
        z = 0.43
        v = 0.31
        matrix_value = np.array([1.0, 0.0]) @ np.linalg.matrix_power(
            trace_matrix(q, z, v), n
        ) @ np.ones(2)
        self.assertAlmostEqual(matrix_value, brute_trace_mgf(n, q, z, v))

    def test_all_coordinates_occupied(self) -> None:
        n = 5
        z = 0.37
        v = 0.23
        matrix_value = np.array([1.0, 0.0]) @ np.linalg.matrix_power(
            trace_matrix(1.0, z, v), n
        ) @ np.ones(2)
        self.assertAlmostEqual(matrix_value, (z + v) ** n)

    def test_projective_trace_term_is_finite(self) -> None:
        term = projective_trace_logterm(
            prime=101,
            k=30,
            n=60,
            cutoff=20,
            degree=8.0,
            message_weight=3,
        )
        self.assertTrue(math.isfinite(term.total_log2))


if __name__ == "__main__":
    unittest.main()
