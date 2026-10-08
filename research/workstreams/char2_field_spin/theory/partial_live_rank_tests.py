from fractions import Fraction as F
from itertools import product
import unittest

from live_kernel_rank_tests import M, INV, dot, rank, subspaces, vectors
from partial_live_rank import partial_live_moment, mds_return_bound


class PartialLiveRankTests(unittest.TestCase):
    def test_exhaustive_small_field_systematic_mds(self):
        cases = 0
        for w in (1, 2):
            s = 3-w
            c = [[INV[i ^ j] for j in range(s, 3)] for i in range(s)]
            b = [tuple(c[i])+tuple(int(i == j) for j in range(s)) for i in range(s)]
            for d in (1, 2, 3):
                for basis in subspaces(3, d):
                    live = [v for v in vectors(basis) if all(v[:w])]
                    moments = {t: F(0) for t in (F(1), F(2), F(3))}
                    for diag in product(range(1, 4), repeat=3):
                        bd = [tuple(M[x][y] for x, y in zip(row, diag)) for row in b]
                        if not any(all(dot(row, v) == 0 for row in bd) for v in live):
                            continue
                        gain = rank([[dot(row, v) for v in basis] for row in bd])
                        for tau in moments:
                            moments[tau] += tau**(-gain)/27
                    for tau, actual in moments.items():
                        bound = partial_live_moment(d, w, w, F(1, 3**s), tau, 4)
                        self.assertLessEqual(actual, bound)
                        self.assertLessEqual(bound, mds_return_bound(w, s, d, tau, q=4))
                        cases += 1
        self.assertEqual(cases, 258)

    def test_kernel_cap_and_near_endpoint(self):
        self.assertEqual(partial_live_moment(3, 2, 0, F(1, 3), F(2), 4), 0)
        # The first segment is not always concave: an arbitrary interpolation
        # through adjacent points would miss the extremizer between 0 and L2.
        self.assertEqual(partial_live_moment(2, 2, 2, F(1, 5), F(3), 4), F(1, 3))
        self.assertEqual(mds_return_bound(2, 8, 10, F(6), selected_zeros=2), 0)

    def test_exhaustive_gf8_joint_output_and_return(self):
        from diagonal_refresh_rank import (MUL as mul8, INV as inv8,
                                            matrix_rank as rank8, subspaces as spaces8)
        def dot8(a, b):
            value = 0
            for x, y in zip(a, b):
                value ^= mul8[x][y]
            return value
        # Full3 Cauchy packet: two outputs/inputs and one state. Normalize
        # its bottom-right entry to one without changing the joint row space.
        matrix = [[inv8[i ^ j] for j in (4, 5, 6)] for i in range(3)]
        matrix[2] = [mul8[inv8[matrix[2][2]]][v] for v in matrix[2]]
        cases = 0
        for support in ((0,), (1,), (0, 1)):
            w = len(support)
            columns = support+(2,)
            for d in range(1, w+2):
                for basis in spaces8(w+1, d):
                    input_forms = [tuple(row[i] for row in basis) for i in range(w)]
                    for selected in ((), (0,), (1,)):
                        z = len(selected)
                        b = [tuple(matrix[i][j] for j in columns) for i in selected+(2,)]
                        histogram = [0]*(d+1)
                        for diag in product(range(1, 8), repeat=w+1):
                            bd = [tuple(mul8[x][y] for x, y in zip(row, diag)) for row in b]
                            equations = [[dot8(row, v) for v in basis] for row in bd]
                            gain = rank8(equations)
                            # Over F8, at most two proper coordinate
                            # hyperplanes cannot cover a nonzero kernel.
                            live = all(rank8(equations+[form]) > gain for form in input_forms)
                            if live:
                                histogram[gain] += 1
                        for tau in (F(1), F(2), F(7)):
                            actual = sum(F(count, 7**(w+1))*tau**(-gain)
                                         for gain, count in enumerate(histogram))
                            self.assertLessEqual(actual, mds_return_bound(w, 1, d, tau, z, q=8))
                            cases += 1
        self.assertEqual(cases, 1503)

    def test_cap_monotonicity_and_small_rank(self):
        for w in (1, 2, 4, 8):
            previous = F(0)
            for cap in range(1, 9):
                value = mds_return_bound(w, 8, cap+w, F(6))
                self.assertGreaterEqual(value, previous)
                previous = value
        self.assertLess(mds_return_bound(1, 8, 2, F(6)), F(6)**-8/F(10**10))


if __name__ == '__main__':
    unittest.main()
