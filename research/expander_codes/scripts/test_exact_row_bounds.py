#!/usr/bin/env python3

from __future__ import annotations

import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction

import numpy as np

from exact_row_bounds import (
    accumulator_slice_logtail,
    exact_row_shell_logdistribution,
    fixed_message_logtail,
)


def bits_of_weight(n: int, weight: int) -> list[int]:
    result = []
    for support in itertools.combinations(range(n), weight):
        value = 0
        for position in support:
            value |= 1 << position
        result.append(value)
    return result


def accumulator_weight(value: int, n: int) -> int:
    state = 0
    weight = 0
    for position in range(n):
        state ^= (value >> position) & 1
        weight += state
    return weight


class ExactRowBoundsTests(unittest.TestCase):
    def test_shell_chain_matches_exhaustive_row_tuples(self) -> None:
        n, d, rows = 7, 2, 3
        available_rows = bits_of_weight(n, d)
        counts: Counter[int] = Counter()
        for chosen in itertools.product(available_rows, repeat=rows):
            value = 0
            for row in chosen:
                value ^= row
            counts[value.bit_count()] += 1

        denominator = len(available_rows) ** rows
        computed = np.exp(exact_row_shell_logdistribution(n, d, rows))
        for shell in range(len(computed)):
            expected = Fraction(counts[shell], denominator)
            self.assertAlmostEqual(computed[shell], float(expected), places=13)
        self.assertAlmostEqual(float(np.sum(computed)), 1.0, places=13)

    def test_accumulator_slice_tail_matches_exhaustive_subsets(self) -> None:
        n, cutoff = 9, 3
        for u in range(n + 1):
            vectors = bits_of_weight(n, u)
            accepted = sum(accumulator_weight(value, n) <= cutoff for value in vectors)
            expected = Fraction(accepted, len(vectors))
            computed = math.exp(accumulator_slice_logtail(n, cutoff, u))
            self.assertAlmostEqual(computed, float(expected), places=13)

    def test_composed_tail_matches_exhaustive_encoder(self) -> None:
        n, d, rows, cutoff = 6, 2, 3, 2
        available_rows = bits_of_weight(n, d)
        accepted = 0
        total = len(available_rows) ** rows
        for chosen in itertools.product(available_rows, repeat=rows):
            value = 0
            for row in chosen:
                value ^= row
            accepted += accumulator_weight(value, n) <= cutoff

        shells = exact_row_shell_logdistribution(n, d, rows)
        computed = math.exp(fixed_message_logtail(n, cutoff, shells))
        self.assertAlmostEqual(computed, accepted / total, places=13)


if __name__ == "__main__":
    unittest.main()
