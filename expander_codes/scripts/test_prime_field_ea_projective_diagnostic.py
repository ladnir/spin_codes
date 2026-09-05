#!/usr/bin/env python3

import itertools
import math
import unittest

from prime_field_ea_projective_diagnostic import (
    log_gap_factor,
    projective_gap_logterm,
)


def direct_container_rank_mass(
    n: int, cutoff: int, theta: float, message_weight: int
) -> list[float]:
    result = [0.0] * (n - cutoff + 1)
    for container in itertools.combinations(range(1, n + 1), cutoff):
        container_set = set(container)
        zeros = [position for position in range(1, n + 1) if position not in container_set]
        previous = 0
        distribution = [1.0]
        for zero in zeros:
            length = zero - previous
            miss = (1.0 - theta) ** (message_weight * length)
            following = [0.0] * (len(distribution) + 1)
            for rank, probability in enumerate(distribution):
                following[rank] += probability * miss
                following[rank + 1] += probability * (1.0 - miss)
            distribution = following
            previous = zero
        for rank, probability in enumerate(distribution):
            result[rank] += probability
    return result


class PrimeFieldEAProjectiveDiagnosticTests(unittest.TestCase):
    def test_gap_factor_matches_direct_series(self) -> None:
        theta = 0.17
        weight = 3
        z = 0.23
        v = 0.41
        direct = 0.0
        for container_positions in range(80):
            miss = (1.0 - theta) ** (weight * (container_positions + 1))
            direct += z ** container_positions * (
                miss + (1.0 - miss) * v
            )
        self.assertAlmostEqual(
            math.exp(log_gap_factor(theta, weight, math.log(z), math.log(v))),
            direct,
            places=12,
        )

    def test_direct_container_mass_counts_all_containers(self) -> None:
        n = 6
        cutoff = 2
        mass = direct_container_rank_mass(n, cutoff, 0.2, 2)
        self.assertAlmostEqual(sum(mass), math.comb(n, cutoff))

    def test_projective_saddle_term_is_finite(self) -> None:
        term = projective_gap_logterm(
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
