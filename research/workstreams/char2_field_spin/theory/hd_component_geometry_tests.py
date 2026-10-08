"""Exact small checks for joint component geometry and its convex reduction."""
from fractions import Fraction as F
from itertools import combinations
import unittest

from hd_component_geometry import (birth_exact, birth_mixture, control_factory,
                                   geometry_matrix, length_corners, rowwise_upper,
                                   run_kernel_interval, group_bins, grouped_triple_matrix)
from hd_component_rank import exact_kernels, component_moment, arb_fraction


def lengths_from_cuts(cuts, rows):
    return (cuts[0], *(b-a for a, b in zip(cuts, cuts[1:])), rows-cuts[-1])


class ComponentGeometryTests(unittest.TestCase):
    def test_birth_mixture_and_convexity(self):
        for rows in range(2, 11):
            for q, mark in ((4, F(1, 2)), (8, F(4, 5))):
                values = [birth_exact(rows, h, q, mark) for h in range(1, rows+1)]
                self.assertEqual(values, [birth_mixture(rows, h, q, mark)
                                           for h in range(1, rows+1)])
                self.assertTrue(all(a >= b for a, b in zip(values, values[1:])))
                self.assertTrue(all(a-2*b+c >= 0 for a, b, c in
                                    zip(values, values[1:], values[2:])))

    def test_uniform_kernel_domination(self):
        for rows in range(3, 9):
            q, mark = 8, F(3, 4)
            exact = exact_kernels(rows, q, mark, collapsed=True)
            control = control_factory(rows, q, mark)
            for h in range(rows+1):
                for minimum in range(1, h+1) if h else (0,):
                    bound = control(h, minimum)
                    self.assertTrue(all(exact[h][i][j] <= bound[i][j]
                                        for i in range(2) for j in range(2)))

    def test_joint_simplex_domination_all_small_cuts(self):
        for rows in range(4, 9):
            q, mark = 8, F(3, 4)
            control = control_factory(rows, q, mark)
            for active in range(1, min(4, rows)+1):
                for cuts in combinations(range(rows), active):
                    lengths = lengths_from_cuts(cuts, rows)
                    matrix = geometry_matrix(lengths, q, control)
                    zero_first = lengths[0] == 0
                    base = (0,)+(1,)*active if zero_first else (1,)*(active+1)
                    excess = rows-sum(base)
                    if not excess:
                        continue
                    indices = list(range(1, active+1) if zero_first else range(active+1))
                    vertices = length_corners(rows, active, zero_first)
                    weights = [F(lengths[i]-base[i], excess) for i in indices]
                    self.assertEqual(sum(weights), 1)
                    corner_values = [geometry_matrix(vertex, q, control) for vertex in vertices]
                    for i in range(4):
                        for j in range(4):
                            self.assertLessEqual(matrix[i][j], sum(weight*vertex[i][j]
                                for weight, vertex in zip(weights, corner_values)))

    def test_global_carry_upper_for_all_small_profiles(self):
        rows, q, mark, active = 5, 8, F(3, 4), 3
        kernels = exact_kernels(rows, q, mark, collapsed=True)
        control = control_factory(rows, q, mark)
        matrices = [geometry_matrix(vertex, q, control)
                    for first in (True, False) for vertex in length_corners(rows, active, first)]
        bound = rowwise_upper(matrices, 2)
        choices = list(combinations(range(rows), active))
        for first in choices:
            for second in choices:
                exact = component_moment((first, second), rows, q, kernels.__getitem__)
                self.assertLessEqual(exact, bound)

    def test_every_small_occupancy_interval(self):
        for rows in range(2, 9):
            q, mark = 8, F(3, 4)
            exact = exact_kernels(rows, q, mark, collapsed=True)
            for low in range(rows+1):
                for high in range(low, rows+1):
                    bound = run_kernel_interval(rows, low, high, q, mark)
                    for h in range(low, high+1):
                        for i in range(2):
                            for j in range(2):
                                self.assertLessEqual(exact[h][i][j], bound[i][j])

    def test_group_bins_cover_once(self):
        for groups in (2, 4, 8, 16, 8192):
            covered = [g for low, high in group_bins(groups) for g in range(low, high+1)]
            self.assertEqual(covered, list(range(groups)))

    def test_directed_group_bins_bound_every_small_lane_triple(self):
        from flint import ctx
        from functools import lru_cache
        previous = ctx.prec
        ctx.prec = 128
        try:
            rows, q, mark = 12, 8, F(3, 4)
            exact = exact_kernels(rows, q, mark, collapsed=True)
            @lru_cache(None)
            def kernel(low, high):
                return run_kernel_interval(rows, low, high, q, mark, directed=True)
            for low in range(3):
                for high in range(low, 3):
                    bound = grouped_triple_matrix(low, high, rows, q, mark, kernel)
                    for group in range(low, high+1):
                        for lanes in combinations(range(4), 3):
                            cuts = tuple(4*group+lane for lane in lanes)
                            matrix = geometry_matrix(lengths_from_cuts(cuts, rows), q,
                                                     lambda h, minimum=0: exact[h])
                            for i in range(4):
                                for j in range(4):
                                    self.assertGreaterEqual(bound[i][j], arb_fraction(matrix[i][j]))
        finally:
            ctx.prec = previous


if __name__ == "__main__":
    unittest.main()
