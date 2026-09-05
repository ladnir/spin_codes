#!/usr/bin/env python3

from __future__ import annotations

import itertools
import math
import unittest

import numpy as np

from expander_bounds import (
    ea_matrix,
    ec_enumerator_matrix,
    ec_matrix,
    log_slice_tail_bound,
    log_weight_mgf,
)


def direct_binary_mgf(code: str, n: int, m: int, q: float, z: float, wrapping: bool) -> float:
    total = 0.0
    random_taps = m - 1 if wrapping else m
    tap_choices = tuple(itertools.product((0, 1), repeat=random_taps))

    for inputs in itertools.product((0, 1), repeat=n):
        input_probability = math.prod(q if bit else 1.0 - q for bit in inputs)
        for taps_by_step in itertools.product(tap_choices, repeat=n):
            state = (0,) * m
            weight = 0
            for bit, random_part in zip(inputs, taps_by_step, strict=True):
                taps = random_part + (1,) if wrapping else random_part
                feedback = sum(a * s for a, s in zip(taps, state, strict=True)) & 1
                output = bit ^ feedback
                weight += output
                state = (output,) + state[:-1]
            total += input_probability * (z**weight) / (len(tap_choices) ** n)
    return total


def direct_binary_enumerator(n: int, m: int, x: float, z: float, wrapping: bool) -> float:
    total = 0.0
    random_taps = m - 1 if wrapping else m
    tap_choices = tuple(itertools.product((0, 1), repeat=random_taps))

    for inputs in itertools.product((0, 1), repeat=n):
        input_weight = sum(inputs)
        for taps_by_step in itertools.product(tap_choices, repeat=n):
            state = (0,) * m
            output_weight = 0
            for bit, random_part in zip(inputs, taps_by_step, strict=True):
                taps = random_part + (1,) if wrapping else random_part
                feedback = sum(a * s for a, s in zip(taps, state, strict=True)) & 1
                output = bit ^ feedback
                output_weight += output
                state = (output,) + state[:-1]
            total += x**input_weight * z**output_weight / (len(tap_choices) ** n)
    return total


class GeneratingMatrixTests(unittest.TestCase):
    def test_accumulator_matrix_matches_exhaustive_enumeration(self) -> None:
        n, q, z = 6, 0.23, 0.41
        expected = 0.0
        for inputs in itertools.product((0, 1), repeat=n):
            probability = math.prod(q if bit else 1.0 - q for bit in inputs)
            state = 0
            weight = 0
            for bit in inputs:
                state ^= bit
                weight += state
            expected += probability * z**weight

        actual = math.exp(log_weight_mgf(ea_matrix(q, z), n, 0))
        self.assertAlmostEqual(actual, expected, places=13)

    def test_nonwrapping_matrix_matches_exhaustive_enumeration(self) -> None:
        n, m, q, z = 4, 2, 0.27, 0.39
        expected = direct_binary_mgf("ec", n, m, q, z, wrapping=False)
        actual = math.exp(log_weight_mgf(ec_matrix(q, z, m, False), n, m))
        self.assertAlmostEqual(actual, expected, places=13)

    def test_wrapping_matrix_matches_exhaustive_enumeration(self) -> None:
        n, q, z = 4, 0.27, 0.39
        for m in (1, 2, 3):
            with self.subTest(memory=m):
                expected = direct_binary_mgf("ec", n, m, q, z, wrapping=True)
                actual = math.exp(log_weight_mgf(ec_matrix(q, z, m, True), n, m))
                self.assertAlmostEqual(actual, expected, places=13)

    def test_wrapping_bivariate_matrix_matches_exhaustive_enumeration(self) -> None:
        n, x, z = 4, 0.31, 0.39
        for m in (1, 2, 3):
            with self.subTest(memory=m):
                expected = direct_binary_enumerator(n, m, x, z, wrapping=True)
                matrix = ec_enumerator_matrix(x, z, m, True)
                actual = math.exp(log_weight_mgf(matrix, n, m))
                self.assertAlmostEqual(actual, expected, places=13)

    def test_wrapping_slice_bound_dominates_exhaustive_tail(self) -> None:
        n, m, u, cutoff = 5, 2, 2, 1
        inputs = tuple(bits for bits in itertools.product((0, 1), repeat=n) if sum(bits) == u)
        tap_choices = tuple(itertools.product((0, 1), repeat=m - 1))
        accepted = 0
        total = len(inputs) * len(tap_choices) ** n
        for bits in inputs:
            for taps_by_step in itertools.product(tap_choices, repeat=n):
                state = (0,) * m
                weight = 0
                for bit, random_part in zip(bits, taps_by_step, strict=True):
                    taps = random_part + (1,)
                    output = bit ^ (sum(a * s for a, s in zip(taps, state, strict=True)) & 1)
                    weight += output
                    state = (output,) + state[:-1]
                accepted += weight <= cutoff
        exact = accepted / total
        log_bound, _, _ = log_slice_tail_bound(
            input_weight=u,
            length=n,
            cutoff=cutoff,
            memory=m,
            wrapping=True,
        )
        self.assertGreaterEqual(math.exp(log_bound) + 1e-13, exact)

    def test_scaled_power_matches_direct_power(self) -> None:
        q, z, n = 0.17, 0.52, 37
        matrix = ea_matrix(q, z)
        direct = float(np.sum(np.linalg.matrix_power(matrix, n)[0, :]))
        scaled = math.exp(log_weight_mgf(matrix, n, 0))
        self.assertAlmostEqual(scaled, direct, places=13)


if __name__ == "__main__":
    unittest.main()
