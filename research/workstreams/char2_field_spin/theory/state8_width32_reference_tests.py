"""Exact scalar width32/state8 checks, independent of density screening."""
from itertools import combinations
from random import Random
import unittest

import state8_width32_reference as ref
from state_capacity_reference import rank


class Width32Tests(unittest.TestCase):
    def test_all_reference_checks(self):
        self.assertEqual(ref.check(), dict(map_basis=80, singleton_returns=64,
                         pairings=2048, extended_adjoints=32, power_sums=62,
                         stream_involutions=4))

    def test_row_scaled_cauchy_identity(self):
        for phase, labels in enumerate(ref.LABELS):
            a, c = ref.matrices(phase)
            self.assertEqual(c, tuple(zip(*a)))
            for u, row in zip(labels, a):
                numerator = ref.linearized(u ^ ref.BETA)
                self.assertNotEqual(numerator, 0)
                for v, coefficient in zip(ref.NODES, row):
                    self.assertEqual(ref.mul(coefficient, u ^ v), numerator)

    def test_bounded_independent_minor_elimination(self):
        # The determinant formula proves every minor; these scalar
        # eliminations independently check representatives of every size.
        rng = Random(832)
        for phase in (0, 1):
            a, _ = ref.matrices(phase)
            for size in range(1, 9):
                for _ in range(8):
                    rows = rng.sample(range(32), size)
                    columns = rng.sample(range(8), size)
                    minor = [[a[i][j] for j in columns] for i in rows]
                    self.assertEqual(rank(minor), size)

    def test_output_weight_interface(self):
        self.assertEqual(ref.return_minimum_weight(1), 25)
        self.assertEqual(ref.return_minimum_weight(2), 23)
        self.assertEqual(ref.return_minimum_weight(8), 17)
        self.assertEqual(ref.return_minimum_weight(24), 1)
        self.assertEqual(ref.return_minimum_weight(32), 1)
        for value in (0, 33):
            with self.assertRaises(ValueError):
                ref.return_minimum_weight(value)


if __name__ == '__main__':
    unittest.main()
