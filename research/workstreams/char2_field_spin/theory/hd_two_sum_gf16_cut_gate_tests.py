"""Exact finite checks of cuts, fixed-type Holder and floor normalization."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_two_sum_gf16_cut_gate import cut_logs, type_logs, source_caps, bound_at


class CutGateTests(unittest.TestCase):
    def test_exact_auxiliary_cut_identity(self):
        c = F(2, 5)
        for masks in product(range(8), repeat=3):
            union = masks[0] | masks[1] | masks[2]
            probability = sum((c**j.bit_count()*(1-c)**(3-j.bit_count())
                               for j in range(8) if all(mask & ~j == 0 for mask in masks)), F(0))
            self.assertEqual(probability, c**union.bit_count())

    def test_exact_fixed_type_holder_not_resampling(self):
        pi = (F(1, 6), F(1, 3), F(1, 2))
        types = ((F(1, 3), F(2, 5), F(3, 7)), (F(4, 9), F(5, 8), F(2, 3)))
        # Two columns, two groups: sum type assignments BEFORE any relaxation.
        lhs = sum((pi[k]*pi[j]*sum(t[k]*t[j] for t in types)**2
                   for k, j in product(range(3), repeat=2)), F(0))
        rhs = sum((pi[k]*sum(t[k]**2 for t in types) for k in range(3)), F(0))**2
        self.assertLessEqual(lhs, rhs)
        direct = sum((sum(pi[k]*types[a][k]*types[b][k] for k in range(3))**2
                      for a, b in product(range(2), repeat=2)), F(0))
        self.assertEqual(lhs, direct)

    def test_floor_marker_bound_for_heterogeneous_columns(self):
        y, floor = F(3, 2), 2
        empty, nonempty = (F(2, 3), F(1, 5), F(1, 4)), (F(1, 7), F(3, 5), F(5, 9))
        exact = sum((np.prod([nonempty[j] if bits[j] else empty[j] for j in range(3)])
                     for bits in product((0, 1), repeat=3) if sum(bits) >= floor), F(0))
        cap = np.prod([empty[j]+y*nonempty[j] for j in range(3)])/y**floor
        self.assertLessEqual(exact, cap)

    def test_once_only_sparse_counts_and_dense_floor(self):
        weights = {a: np.zeros(17) for a in range(3, 17)}
        logs = type_logs(0., 0., 0., weights)
        self.assertAlmostEqual(16*logs[0, 16], log2(16)+16*log2(256)-8*log2(255), places=11)
        self.assertAlmostEqual(16*logs[1, 16], log2(120)+32*log2(256)-16*log2(255), places=11)
        self.assertAlmostEqual(16*logs[1, 15], log2(120)+16*log2(256)-16*log2(255), places=11)
        marked = type_logs(0., 2., 3., weights)
        self.assertAlmostEqual(16*(marked[0, 0]-logs[0, 0]), -18., places=11)
        self.assertAlmostEqual(16*(marked[1, 0]-logs[1, 0]), -36., places=11)
        for row in range(2, 16):
            self.assertAlmostEqual(16*(marked[row, 0]-logs[row, 0]), -27., places=11)
        for k in range(17):
            self.assertAlmostEqual(logs[2, k], k, places=11)

    def test_cut_boundaries_and_empty_surcharge_guard(self):
        self.assertTrue(np.all(np.isneginf(cut_logs(0.)[:-1])))
        self.assertEqual(cut_logs(0.)[-1], 0.)
        self.assertAlmostEqual(float(np.exp2(cut_logs(-5.)).sum()), 1., places=13)
        with self.assertRaises(ValueError):
            cut_logs(.1)
        caps = source_caps()
        self.assertEqual(caps.shape, (33,))
        self.assertEqual(caps[0], 0.)
        caps[0] = 1e-20
        with self.assertRaises(ValueError):
            bound_at(caps, {}, 0., 0.)
        with self.assertRaises(ValueError):
            type_logs(0., -1., 0., {a: np.zeros(17) for a in range(3, 17)})


if __name__ == '__main__':
    unittest.main()
