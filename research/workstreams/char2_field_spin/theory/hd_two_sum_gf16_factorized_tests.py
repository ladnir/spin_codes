"""Exact tiny proof transitions and lightweight scope tests; no full gate."""
from fractions import Fraction as F
from itertools import product
from math import comb, log2
import unittest

import numpy as np

from hd_two_sum_gf16_factorized import (
    sparse_pair_terms, positive_binomial_tail, contraction_envelope,
    group_bounds, all_cohort_band, run)


class FactorizedTests(unittest.TestCase):
    def test_individual_sparse_row_floors_exact(self):
        counts = {}
        for word in product(range(4), repeat=4):
            if sum(bool(x & 1) for x in word) < 2 or sum(bool(x & 2) for x in word) < 2:
                continue
            key = sum(x in (1, 2) for x in word), word.count(3)
            counts[key] = counts.get(key, 0)+1
        self.assertEqual({(s, d): n for s, d, n in sparse_pair_terms(4, 2)}, counts)
        self.assertEqual(sum(counts.values()), sum(comb(4, w) for w in range(2, 5))**2)

    def test_exact_noncommuting_group_factorization(self):
        # Two groups, two labelled regions, two columns, floor1 per group.
        tables = [
            [(F(1), F(0), F(0), F(1, 3)), (F(1, 4), F(1, 3), F(1, 8), F(1, 5)),
             (F(1, 6), F(1, 8), F(1, 7), F(1, 9))],
            [(F(1), F(0), F(0), F(1, 5)), (F(1, 7), F(1, 8), F(1, 6), F(1, 4)),
             (F(1, 3), F(1, 9), F(1, 7), F(1, 8))]]
        v, b = (F(1), F(1, 2)), F(1, 2)
        intercepts = []
        for table in tables:
            intercepts.append(max((matrix[2*j]*v[0]+matrix[2*j+1]*v[1])/v[j]/b**m
                for m, matrix in enumerate(table) for j in range(2)))
        weights = (F(1, 3), F(2, 5), F(3, 7))
        exact = F(0)
        for masks in product(range(4), repeat=4):
            if not (masks[0] or masks[1]) or not (masks[2] or masks[3]):
                continue
            row, mass = (F(1), F(0)), F(1)
            for mask in masks:
                mass *= weights[mask.bit_count()]
            for column in range(2):
                for region, table in enumerate(tables):
                    m = ((masks[column] >> region) & 1)+((masks[2+column] >> region) & 1)
                    a, bb, c, d = table[m]
                    row = (row[0]*a+row[1]*c, row[0]*bb+row[1]*d)
            exact += mass*sum(row)
        empty = weights[0]
        nonempty = 2*weights[1]*b+weights[2]*b*b
        group = 2*empty*nonempty+nonempty*nonempty
        upper = (intercepts[0]*intercepts[1])**2*group**2/min(v)
        self.assertLessEqual(exact, upper)

    def test_all_cohort_polynomial_counts_once(self):
        # Distinct dense arities are whole-word terms, not a per-column sum.
        groups = (F(2, 3), F(3, 5), F(5, 7), F(7, 11))
        explicit = sum((groups[i]*groups[j]*groups[k] for i, j, k in product(range(4), repeat=3)), F(0))
        self.assertEqual(explicit, sum(groups)**3)
        self.assertNotEqual((groups[2]+groups[3])**2, groups[2]**2+groups[3]**2)

    def test_normalization_at_unit_occupancy_mark(self):
        # Artificial dense weights isolate the source-word/floor factors.
        weights = {a: np.zeros(17) for a in range(3, 17)}
        result = group_bounds(0., weights)
        row_count = sum(comb(16, w)*255**w for w in range(9, 17))/255**8
        self.assertAlmostEqual(result['singleton'], log2(16)+log2(row_count), places=11)
        self.assertAlmostEqual(result['pair'], log2(120)+2*log2(row_count), places=11)
        expected_dense = log2(sum(comb(16, j)*(2**16-1)**j for j in range(9, 17)))
        for value in result['dense_by_arity'].values():
            self.assertAlmostEqual(value, expected_dense, places=11)

    def test_complete_pair_histograms_and_scope(self):
        matrices = {(m, f): np.log2(np.array([[1., .125], [.25, .5]]))
                    for m in range(3) for f in range(m//2+1)}
        caps = contraction_envelope(matrices, 2, -1.)
        self.assertTrue(np.allclose(caps, log2(17/16)))
        del matrices[2, 1]
        with self.assertRaises(ValueError):
            contraction_envelope(matrices, 2, -1.)
        total, rows = all_cohort_band(.1, -12., 18, 20, F(99, 100), -40.)
        self.assertEqual([row['active_groups'] for row in rows], [18, 19, 20])
        self.assertGreaterEqual(total, max(row['failure_log2_screen'] for row in rows))
        with self.assertRaises(ValueError):
            run(first=17)
        self.assertAlmostEqual(positive_binomial_tail(-np.inf, 2., 3, 2), 6.)


if __name__ == '__main__':
    unittest.main()
