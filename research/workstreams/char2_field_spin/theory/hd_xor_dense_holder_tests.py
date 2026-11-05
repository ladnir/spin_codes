"""Tiny exact checks; no regional production computation or native work."""
from fractions import Fraction as F
from itertools import product
from math import prod, log2
import unittest

import numpy as np

from hd_xor_dense_holder import (uniform_mask_count, holder_moments,
    dense_feature_logs, holder_branch_caps, frozen_controls)


class DenseHolderTests(unittest.TestCase):
    def test_all_GF4_masks_against_direct_uniform_inputs(self):
        q, width = 4, 4
        for active in range(1, 1 << width):
            positions = [i for i in range(width) if active >> i & 1]
            counts = [0]*(1 << width)
            for values in product(range(q), repeat=len(positions)):
                x = [0]*width
                for i, value in zip(positions, values):
                    x[i] = value
                total = 0
                for value in x:
                    total ^= value
                target = sum((value ^ total != 0) << i for i, value in enumerate(x))
                counts[target] += 1
            self.assertEqual(counts, [uniform_mask_count(q, width, active, target) for target in range(1 << width)])

    def test_closed_moments_and_every_three_column_holder_inequality(self):
        q, width, columns, dimension = 4, 4, 3, 2
        active_sets = [a for a in range(1 << width) if a.bit_count() >= 2]
        counts = {a: [uniform_mask_count(q, width, a, t) for t in range(1 << width)] for a in active_sets}
        moments = holder_moments(q, width, columns, dimension, 2)
        for target in range(1 << width):
            expected = sum((F(counts[a][target]**columns, (q-1)**a.bit_count()) for a in active_sets), F(0))
            self.assertEqual(expected, moments[target.bit_count()])
        for word in product(range(1 << width), repeat=columns):
            left = sum((F(prod(counts[a][t] for t in word), (q-1)**a.bit_count()) for a in active_sets), F(0))
            self.assertLessEqual(left**columns, prod(moments[t.bit_count()] for t in word))

    def test_dense_feature_polynomial_counts_whole_masks(self):
        weights = [1, 2, 3, 5, 7]
        expected = {}
        for masks in product(range(16), repeat=2):
            key = (sum(mask != 0 for mask in masks), sum(mask.bit_count() for mask in masks), sum(mask == 15 for mask in masks))
            expected[key] = expected.get(key, 0)+prod(weights[mask.bit_count()] for mask in masks)
        result = dense_feature_logs(np.log2(weights), 2)
        self.assertEqual(set(result), set(expected))
        for key, value in expected.items():
            self.assertAlmostEqual(result[key], log2(value), places=12)

    def test_weighted_mask_mixture_precedes_potential_quotient(self):
        weights = [F(1, 3), F(2), F(3), F(5), F(7)]
        values = [(F(1+mask**2), F(1, 4)+((mask+7) % 16)**2) for mask in range(16)]
        potential = F(1, 4)
        empty, nonempty = holder_branch_caps(np.log2(np.array(values, dtype=float))[None, :, :], -2., np.log2(list(map(float, weights))))
        exact = [sum(weights[mask.bit_count()]*values[mask][state] for mask in range(1, 16)) for state in range(2)]
        self.assertAlmostEqual(nonempty[0], log2(max(exact[0], exact[1]/potential)), places=12)
        self.assertAlmostEqual(empty[0], log2(weights[0]*max(values[0][0], values[0][1]/potential)), places=12)

    def test_five_frozen_controls_have_no_fits(self):
        controls, reference = frozen_controls()
        self.assertEqual(len(controls), 5)
        self.assertEqual(controls[0][:3], ('inherited', F(996598217, 10**9), F(-135, 2)))
        self.assertAlmostEqual(reference, -61.000916136180436, places=10)
        self.assertTrue(all(0 < z < 1 and logp <= 0 for _, z, logp, _ in controls))


if __name__ == '__main__':
    unittest.main()
