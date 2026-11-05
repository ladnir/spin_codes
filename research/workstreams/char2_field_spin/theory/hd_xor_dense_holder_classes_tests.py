"""Exact tiny classes and mocked numerical ordering, no production replay."""
from fractions import Fraction as F
from itertools import product
from math import prod, log2
import unittest

import numpy as np

from hd_xor_dense_holder import uniform_mask_count, holder_moments
from hd_xor_dense_holder_classes import (class_moments, class_log_weights,
                                        mask_branch_caps, frozen_classes)


class HolderClassesTests(unittest.TestCase):
    def test_partition_moments_and_word_domination(self):
        q, width, columns, dimension = 4, 4, 3, 2
        classes = ((3, 7, 10, 15), (5, 6, 9, 11, 12, 13, 14))
        weights = []
        for active in classes:
            numerators, denominator = class_moments(active, q, width, columns, dimension)
            moments = [F(n, denominator) for n in numerators]
            weights.append(moments)
            for word in product(range(16), repeat=columns):
                lhs = sum((F(prod(uniform_mask_count(q, width, a, t) for t in word),
                             (q-1)**a.bit_count()) for a in active), F(0))
                self.assertLessEqual(lhs**columns, prod(moments[t] for t in word))
        total = holder_moments(q, width, columns, dimension, 2)
        for mask in range(16):
            self.assertEqual(sum(w[mask] for w in weights), total[mask.bit_count()])

    def test_non_symmetric_class_and_zero_weights(self):
        numerators, denominator = class_moments((3,), 4, 4, 3, 2)
        self.assertNotEqual(numerators[3], numerators[5])
        logs = class_log_weights(numerators, denominator, 3)
        self.assertTrue(np.isneginf(logs[4]))
        self.assertAlmostEqual(logs[0], log2(float(F(numerators[0], denominator)))/3)

    def test_mask_mixture_before_quotient(self):
        weights = [F(1, 3), F(2), F(0), F(5)]
        values = [(F(1), F(7)), (F(2), F(3)), (F(100), F(100)), (F(17), F(1))]
        logs = np.array([log2(float(w)) if w else -np.inf for w in weights])
        empty, nonempty = mask_branch_caps(np.log2(np.array(values, dtype=float))[None], -2., logs)
        self.assertAlmostEqual(empty[0], log2(float(weights[0]*max(values[0][0], 4*values[0][1]))))
        sums = [sum(weights[m]*values[m][s] for m in range(1, 4)) for s in range(2)]
        self.assertAlmostEqual(nonempty[0], log2(float(max(sums[0], 4*sums[1]))))

    def test_frozen_partition_and_guards(self):
        controls, _ = frozen_classes()
        masks = [mask for _, _, _, group, _ in controls for mask in group]
        self.assertEqual(len(controls), 5)
        self.assertEqual(sorted(masks), [m for m in range(256) if m.bit_count() >= 3])
        self.assertEqual(sorted(len(group) for _, _, _, group, _ in controls), [6, 6, 35, 68, 104])
        for active in ((), (3, 3), (0,), (16,)):
            with self.assertRaises(ValueError):
                class_moments(active, 4, 4, 3, 2)


if __name__ == '__main__':
    unittest.main()
