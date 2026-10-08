"""Bounded actual-map rank checks and joint-moment identity tests."""

from fractions import Fraction as F
from math import comb
import unittest

from hd_randomized_outer_joint_quiet import (
    directed_endpoint, directed_interval, joint_moment_enclosures, joint_moment_logs,
    tile_residual_terms, tile_ranks,
)
from hd_randomized_outer_proof import prescribed_mds_count


class JointQuietTests(unittest.TestCase):
    def setUp(self):
        from flint import ctx
        self.previous = ctx.prec
        ctx.prec = 160

    def tearDown(self):
        from flint import ctx
        ctx.prec = self.previous

    def test_actual_beta64_tile_rank_formulas(self):
        # Two representative masks of each size, at all eight distinct tiles.
        for tile in range(8):
            for a in range(9):
                for v in range(9):
                    source = (1 << a)-1
                    zero = (1 << v)-1
                    for rotation in (0, 3):
                        mask = ((source << rotation) | (source >> (8-rotation))) & 255
                        self.assertEqual(tile_ranks(tile, mask, zero),
                                         (min(a, v), min(4, max(a-v, 0))))

    def test_local_zero_expansion_matches_mds_enumerator(self):
        z, q = F(1, 7), 256
        local = tile_residual_terms(z, q)
        for a in range(9):
            actual = sum(local[a])/comb(8, a)
            expected = sum(comb(8, w)*prescribed_mds_count(w, 8, a, q)*z**w
                           for w in range(9))
            self.assertEqual(actual, expected)

    def test_zero_small_occupancies_and_free_full_input(self):
        from flint import arb
        z, q = F(1, 7), 256
        quiet, free = joint_moment_enclosures(z, q)
        self.assertTrue((quiet[0]-1).contains(0))
        for j in range(1, 5):
            self.assertTrue((quiet[j]-arb(1)/q**j).contains(0))
        actual = (1+(q-1)*arb(1)/7)**64/q**64
        self.assertTrue((free[64]-actual).contains(0))
        for j in range(65):
            # Roundoff may overlap equality, so compare through an enclosure.
            self.assertFalse(quiet[j] < arb(1)/q**j)
            self.assertFalse(quiet[j] > free[j])

    def test_distinct_tile_residual_lines(self):
        from itertools import combinations
        from hd_blockwise_shared_maps import matrices
        from hd_blockwise_maps import mul, power
        feedthrough, _, feedback, _, _ = matrices(64)
        seen = {}
        # Bounded family: two potential inputs and one output-zero constraint.
        # The note's root-sum argument proves all one-dimensional cases.
        for tile in range(8):
            for i, j in combinations(range(8), 2):
                for v in range(8):
                    ii, jj, vv = 8*tile+i, 8*tile+j, 8*tile+v
                    vector = tuple(mul(row[ii], feedthrough[vv][jj])
                                   ^mul(row[jj], feedthrough[vv][ii]) for row in feedback)
                    inverse = power(next(x for x in vector if x), 254)
                    line = tuple(mul(x, inverse) for x in vector)
                    if line in seen:
                        self.assertEqual(seen[line], tile)
                    seen[line] = tile

    def test_refinement_and_z_one_baseline(self):
        old, _ = joint_moment_logs(F(1, 65536), distinct_lines=False)
        new, _ = joint_moment_logs(F(1, 65536), distinct_lines=True)
        self.assertTrue(all(a <= b+1e-10 for a, b in zip(new, old)))
        self.assertTrue(any(a < b-0.1 for a, b in zip(new, old)))
        quiet, free = joint_moment_logs(F(1))
        self.assertEqual(quiet, tuple(-16*min(j, 4) for j in range(65)))
        self.assertEqual(free, (0.0,)*65)

    def test_selected_directed_endpoint(self):
        from flint import arb
        result = directed_endpoint(6144, F(1, 16384), F(127, 128), F(3, 8))
        self.assertTrue(arb(result['failure_log2_upper']) < -134588)
        self.assertTrue(result['below_2_to_minus_40'])
        self.assertTrue(result['covers_only_this_occupancy'])

    def test_convex_interval_receipt(self):
        from flint import arb
        result = directed_interval()
        self.assertEqual(result['active_rows_first'], 5940)
        self.assertEqual(result['active_rows_last'], 8192)
        self.assertTrue(arb(result['union_log2_upper']) < -3832)
        self.assertTrue(result['below_2_to_minus_40'])
        self.assertEqual(result['occupancy_endpoints_evaluated'], 2)


if __name__ == '__main__':
    unittest.main()
