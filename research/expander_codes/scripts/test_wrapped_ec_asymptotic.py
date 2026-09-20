import math
import unittest

import numpy as np
from scipy.optimize import minimize_scalar

from wrapped_ec_asymptotic import (
    binary_gv_distance,
    effective_eigenvalue,
    optimized_sparse_rate,
    spectral_radius,
    sparse_rate,
    wrapping_transition_matrix,
)


class WrappedECAsymptoticTests(unittest.TestCase):
    def test_rate_one_fifth_gv_distance(self) -> None:
        self.assertAlmostEqual(binary_gv_distance(0.2), 0.243003853809, places=11)

    def test_memory_one_recovers_accumulator_rate(self) -> None:
        delta = 0.05
        expected = 1.0 - 2.0 * math.sqrt(delta * (1.0 - delta))
        self.assertAlmostEqual(sparse_rate(1, delta), expected, places=14)

    def test_closed_rate_matches_effective_optimization(self) -> None:
        for memory in (1, 2, 5, 21):
            self.assertAlmostEqual(
                sparse_rate(memory, 0.05),
                optimized_sparse_rate(memory, 0.05),
                places=10,
            )

    def test_uniform_input_has_binomial_perron_root(self) -> None:
        for memory in (1, 2, 7):
            for z in (0.1, 0.5, 1.0):
                matrix = wrapping_transition_matrix(0.5, z, memory)
                expected = (1.0 + z) / 2.0
                np.testing.assert_allclose(matrix.sum(axis=1), expected)
                radius = float(np.max(np.abs(np.linalg.eigvals(matrix))))
                self.assertAlmostEqual(radius, expected, places=13)

    def test_effective_generator_matches_exact_sparse_transfer(self) -> None:
        memory = 2
        delta = 0.05
        optimum = minimize_scalar(
            lambda theta: delta * theta + effective_eigenvalue(memory, theta),
            bounds=(0.0, 16.0),
            method="bounded",
        )
        theta = float(optimum.x)
        q = 1.0e-6
        radius = spectral_radius(q, math.exp(-theta * q), memory)
        exact_scaled_rate = -delta * theta - math.log(radius) / q
        self.assertAlmostEqual(
            exact_scaled_rate,
            sparse_rate(memory, delta),
            places=5,
        )

    def test_memory_improves_sparse_rate(self) -> None:
        rates = [sparse_rate(memory, 0.05) for memory in range(1, 10)]
        self.assertTrue(all(left < right for left, right in zip(rates, rates[1:])))
        self.assertLess(rates[-1], 0.9)


if __name__ == "__main__":
    unittest.main()
