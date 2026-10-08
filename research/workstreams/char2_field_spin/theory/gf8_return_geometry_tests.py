"""Independent scalar checks of the exact sparse return geometry."""
from itertools import combinations
from random import Random
import unittest

from gf8_return_geometry import return_matrix, certify_sparse_returns
from state_capacity_reference import matrices


def product(a, b):
    # Polynomial multiplication and long division, independent of the table.
    value = 0
    for bit in range(8):
        if b >> bit & 1:
            value ^= a << bit
    for bit in range(14, 7, -1):
        if value >> bit & 1:
            value ^= 0x11b << (bit-8)
    return value


def total(values):
    result = 0
    for value in values:
        result ^= value
    return result


class ReturnGeometryTests(unittest.TestCase):
    def test_exact_minima_all_phases(self):
        for state, expected in ((4, [13, 13]), (6, [15, 13]), (8, [9, 9])):
            result = certify_sparse_returns(state)
            self.assertEqual(result['minimum_return_weight'], expected)
            for phase in result['phases']:
                self.assertEqual(phase['projective_inputs_checked'], 30616)
                b = return_matrix(state, phase['phase'])
                for weight, witness in enumerate(phase['witnesses'], 1):
                    output = [total(product(row[j], x) for j, x in
                                    zip(witness['positions'], witness['values'])) for row in b]
                    self.assertEqual(sum(x != 0 for x in output),
                                     phase['minimum_return_weight'][weight-1])

    def test_direct_return_recurrence_and_ca_zero(self):
        rng = Random(138741)
        for state in (4, 6, 8):
            for phase in range(4):
                a, c = matrices(phase, True)
                for i in range(state):
                    for j in range(state):
                        self.assertEqual(total(product(c[i][h], a[h][j]) for h in range(16)), 0)
                b = return_matrix(state, phase)
                for _ in range(50):
                    x = [rng.randrange(256) for _ in range(16)]
                    refreshed = [total(product(c[h][j], x[j]) for j in range(16))
                                 for h in range(state)]
                    y = [x[j] ^ product(2, total(product(a[j][h], refreshed[h])
                                                for h in range(state))) for j in range(16)]
                    expected = [total(product(row[j], x[j]) for j in range(16)) for row in b]
                    self.assertEqual(y, expected)
                    # The map I+2AC is an involution when CA=0.
                    self.assertEqual([total(product(row[j], y[j]) for j in range(16))
                                      for row in b], x)

    def test_projective_scaling_preserves_weight(self):
        b = return_matrix(6, 1)
        for left, right in combinations(range(4), 2):
            for x, y in ((1, 1), (17, 251), (255, 3)):
                original = [product(row[left], x) ^ product(row[right], y) for row in b]
                for scalar in (2, 17, 255):
                    scaled = [product(row[left], product(scalar, x)) ^
                              product(row[right], product(scalar, y)) for row in b]
                    self.assertEqual(scaled, [product(scalar, value) for value in original])
                    self.assertEqual(sum(v != 0 for v in scaled), sum(v != 0 for v in original))


if __name__ == '__main__':
    unittest.main()
