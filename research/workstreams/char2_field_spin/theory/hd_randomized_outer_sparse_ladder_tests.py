"""Independent checks of local replay, regional coefficients, and sparse joins."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
import unittest

from hd_randomized_outer_sparse_ladder import (
    arb_local_profiles, arb_matrix_product, arb_polynomial_product,
    directed_endpoint, float_regions, occupied_tile_counts, sparse_receipt,
)
from hd_randomized_outer_state_profile import local_profiles, one_row_receipt


def mul8(a, b):
    result = 0
    for _ in range(3):
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 8:
            a ^= 11
    return result


class SparseLadderTests(unittest.TestCase):
    def setUp(self):
        from flint import ctx
        self.previous = ctx.prec
        ctx.prec = 160

    def tearDown(self):
        from flint import ctx
        ctx.prec = self.previous

    def test_every_small_occupancy_matches_directed_local_replay(self):
        direct = arb_local_profiles(F(-3, 512), 16)
        floating = local_profiles(-3/512)
        for j in range(17):
            for row in range(3):
                for column in range(3):
                    actual = direct[j][row][column]
                    expected = floating[j][row][column]
                    if expected == float('-inf'):
                        self.assertTrue(actual.is_zero())
                    else:
                        self.assertAlmostEqual(log2(float(actual)), expected, places=10)

    def test_noncommuting_polynomial_coefficients(self):
        from flint import arb
        a = ((arb(1), arb(2), arb(0)), (arb(0), arb(1), arb(0)),
             (arb(0), arb(0), arb(1)))
        b = ((arb(0), arb(0), arb(1)), (arb(3), arb(0), arb(0)),
             (arb(0), arb(1), arb(0)))
        squared = arb_polynomial_product((a, b), (a, b), 2)
        cube = arb_polynomial_product(squared, (a, b), 2)
        terms = (arb_matrix_product(arb_matrix_product(a, b), b),
                 arb_matrix_product(arb_matrix_product(b, a), b),
                 arb_matrix_product(arb_matrix_product(b, b), a))
        for i in range(3):
            for j in range(3):
                self.assertTrue((cube[2][i][j]-sum(term[i][j] for term in terms)).is_zero())

    def test_regrouping_at_alpha_one_and_existing_one_row_receipt(self):
        import numpy as np
        from flint import arb
        one = float_regions(F(-1, 256), 1, maximum=4, steps=1)
        four = float_regions(F(-1, 256), 1, maximum=4, steps=4)
        for a, b in zip(one, four):
            np.testing.assert_allclose(a, b, rtol=2e-13, atol=1e-300)
        actual = directed_endpoint(1, F(-1, 512), steps=1)
        old = one_row_receipt()
        self.assertTrue(arb(actual['failure_log2_upper']).overlaps(arb(old['failure_log2_upper'])))

    def test_occupied_tile_normalization(self):
        for j in range(65):
            counts = occupied_tile_counts(j)
            self.assertTrue(all(count >= 0 for _, count in counts))
            self.assertEqual(sum(count for _, count in counts), comb(64, j))
        self.assertEqual(dict(occupied_tile_counts(2))[1], 8*comb(8, 2))

    def test_mds_span_membership_bound_over_gf8(self):
        columns = tuple((1, x, mul8(x, x), mul8(x, mul8(x, x))) for x in range(8))
        for j in (1, 2, 3):
            membership = Counter()
            for support in combinations(columns, j):
                span = set()
                for coefficients in product(range(8), repeat=j):
                    vector = [0]*4
                    for coefficient, column in zip(coefficients, support):
                        for coordinate in range(4):
                            vector[coordinate] ^= mul8(coefficient, column[coordinate])
                    span.add(tuple(vector))
                self.assertEqual(len(span), 8**j)
                for vector in span:
                    if any(vector):
                        membership[vector] += 1
            bound = comb(7, j-1)
            self.assertLessEqual(max(membership.values()), bound)
            self.assertTrue(all(membership[column] == bound for column in columns))

    def test_directed_contiguous_union_and_baseline(self):
        from flint import arb
        baseline = sparse_receipt('base')
        refined = sparse_receipt('membership')
        self.assertEqual([record['active_rows'] for record in refined['endpoints']], list(range(1, 6)))
        self.assertTrue(arb(baseline['union_log2_upper']) < -51)
        self.assertTrue(arb(refined['union_log2_upper']) < -56)
        self.assertTrue(arb(refined['union_log2_upper']) < arb(baseline['union_log2_upper']))
        self.assertTrue(all(record['below_2_to_minus_40'] for record in refined['endpoints']))


if __name__ == '__main__':
    unittest.main()
