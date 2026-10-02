"""Small normalization and scope checks for the fresh tail comparison."""
from fractions import Fraction as Q
from itertools import product
from math import factorial
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb, arb_mat
import packet_tail as screen


class PacketTailTests(unittest.TestCase):
    def test_exact_rs8_shell_majorant(self):
        components, counts = screen.exact_mixture()
        self.assertEqual(len(components), 5)
        self.assertEqual(sum(counts), 2**128 - 1)
        self.assertEqual(components[0], (1, 0, 0))
        self.assertEqual([p for _, p, _ in components[1:]], [Q(75, 128), Q(45, 64), Q(105, 128), Q(15, 16)])

    def test_density_caps_match_exhaustive_integer_allocations(self):
        for n in range(1, 8):
            for a, b in product(range(n + 1), repeat=2):
                if a + b < n:
                    continue
                candidates = []
                for x in range(max(0, n - b), min(n, a) + 1):
                    value = Q(n**n, factorial(n))
                    for count in (x, n - x):
                        value *= Q(factorial(count), count**count)
                    candidates.append(value)
                self.assertEqual(screen.capped_density_loss(n, (a, b)), max(candidates))

    def test_explicit_matrix_power(self):
        matrix = np.array([[.5, .1], [.2, .4]])
        for exponent in (0, 1, 2, 7, 16):
            expected = np.log(np.linalg.matrix_power(matrix, exponent)[0].sum())
            self.assertAlmostEqual(screen.log_power(matrix, exponent), expected)

    def test_geometry_and_component_lifetime(self):
        geometry = screen.q1.Geometry(32, 2, 4)
        model = screen.TailModel(geometry, [(Q(1), Q(0), 0), (Q(3), Q(1, 2), 1)], {}, tilt=Q(1, 2))
        self.assertEqual((model.packets, model.epochs, model.threshold), (64, 2, 25))
        self.assertEqual(model.normalized_masses, (1, 3 * Q(3, 4)**2))
        self.assertEqual(model.root, (Q(1, 32), Q(1, 3)))
        with self.assertRaises(ValueError):
            model.weights((Q(0), Q(1)))

    def test_outward_bound_contains_tiny_exact_component_sum(self):
        # Constant scalar inner=1 and theta=1: the only losses are a harmless
        # Chernoff factor and shuffled-density comparison. With zero duals,
        # the bound must dominate every active-group selection counted exactly.
        geometry = screen.q1.Geometry(32, 2, 4)
        model = screen.TailModel(geometry, [(Q(1), Q(0), 0), (Q(3), Q(1, 2), 1)], {}, tilt=Q(1))
        with patch.object(screen.q1.kernel_t64, 'outward', return_value=arb_mat([[1]])):
            upper = model.outward(model.root, dict(tilt='1', parameters=['1/100', '0', '0']))
        self.assertGreater(upper, arb(4)**32)


if __name__ == '__main__':
    unittest.main()
