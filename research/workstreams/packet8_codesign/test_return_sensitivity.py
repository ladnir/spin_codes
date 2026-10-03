"""Check sensitivity isolation and a small exact tilted-syndrome bound."""
from itertools import combinations, product
import unittest

import numpy as np

import maps
import return_sensitivity
import screen


def gf4(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


class SensitivityTests(unittest.TestCase):
    def test_only_return_entries_change(self):
        local = np.arange(1., 1.+3*4*4).reshape(3, 4, 4)
        changed = return_sensitivity.change_denominator(local, 24)
        expected = local.copy()
        expected[:, 1:, 0] *= 65535/(2**24-1)
        np.testing.assert_array_equal(changed, expected)
        np.testing.assert_array_equal(return_sensitivity.change_denominator(local, 16), local)
        np.testing.assert_array_equal(changed[:, 0], local[:, 0])
        np.testing.assert_array_equal(changed[:, :, 1:], local[:, :, 1:])

    def test_tilted_syndrome_atom_bounds_small_field(self):
        columns = tuple((1 << b) | (gf4(h, 1 << b) << 2) for h in range(3) for b in range(2))
        for pair in combinations(range(3), 2):
            self.assertEqual(maps.rank(columns[2*h+b] for h in pair for b in range(2)), 4)
        for z in (.73, .91):
            minimum = (1+z)**2-1
            for a in range(16):
                image = [(a & 3) ^ gf4(h, a >> 2) for h in range(3)]
                for j in (1, 2, 3):
                    for support in combinations(range(3), j):
                        masses = np.zeros(16)
                        for values in product(range(1, 4), repeat=j):
                            x = sum(value << (2*h) for h, value in zip(support, values))
                            y = x ^ sum(value << (2*h) for h, value in enumerate(image))
                            masses[maps.apply(columns, x)] += z**y.bit_count()
                        law = masses/masses.sum()
                        bound = 1/minimum**min(j, 2)
                        self.assertLessEqual(float(law.max()), bound+1e-13)
                        if j <= 2:
                            self.assertEqual(masses[0], 0)


if __name__ == '__main__':
    unittest.main()
