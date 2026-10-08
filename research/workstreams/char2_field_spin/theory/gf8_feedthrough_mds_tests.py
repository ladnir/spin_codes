"""Exact algebra and positive-transfer checks for the full-MDS packet."""
from fractions import Fraction as Q
from itertools import combinations
import math
from random import Random
import unittest

from gf8_feedthrough_mds_screen import weighted_matrices
from gf8_feedthrough_screen import eliminate_inputs, state_profile
from gf8_mds_packet import matrices, selected_rank, fullstate_matrices, largestate_matrices
from gf8_packet_maps import MUL, determinant, rank
from gf8_rank_history_screen import cone_preserved


class MdsFeedthroughTests(unittest.TestCase):
    def test_full_selected_ranks_and_state_schur_minors(self):
        rng = Random(1296440)
        f, a, c, eta = matrices()
        for _ in range(512):
            w, z, h = rng.randrange(9), rng.randrange(9), rng.randrange(5)
            inputs = tuple(sorted(rng.sample(range(8), w)))
            outputs = tuple(sorted(rng.sample(range(8), z)))
            state = tuple(sorted(rng.sample(range(4), h)))
            self.assertEqual(selected_rank(outputs, inputs, state), min(z, w+h))
            self.assertEqual(selected_rank(outputs, inputs, state, True), min(4+z, w+h))
            b = tuple(tuple(f[j][i] for i in inputs) for j in outputs)
            emission = tuple(tuple(MUL[eta][v] for v in a[j]) for j in outputs)
            pivots, g = eliminate_inputs(b, emission)
            self.assertEqual(pivots, min(z, w))
            expected = min(4, max(0, z-w))
            self.assertEqual(state_profile(g), (expected, expected))
            for size in range(1, expected+1):
                for rows in combinations(range(len(g)), size):
                    for columns in combinations(range(4), size):
                        self.assertNotEqual(determinant(tuple(tuple(g[i][j] for j in columns) for i in rows)), 0)

    def test_transfer_cone_and_birth_caps(self):
        for y, tau in ((Q(2, 3), Q(11)), (Q(99, 100), Q(255)), (Q(1), Q(2))):
            local = weighted_matrices(y, tau)
            for matrix in local:
                self.assertTrue(cone_preserved(matrix))
            self.assertEqual(local[0][0], [1, 0, 0, 0, 0])
            self.assertEqual(local[1][0], [0, 8*y**8, 0, 0, 0])
            for w in range(1, 5):
                self.assertEqual(local[w][0][0], 0)
            for w in range(1, 9):
                self.assertTrue(all(local[w][0][c] == 0 for c in range(min(4, w)+1, 5)))

    def test_joint_and_pure_returns_have_one_diagonal_charge(self):
        y, tau = Q(3, 4), Q(13)
        old = weighted_matrices(y, tau, live_kernel=False, partial_live=False)
        new = weighted_matrices(y, tau)
        for w in range(1, 9):
            for cap in range(1, 5):
                pure = math.comb(8, w)*y**(9-w)/tau**4*sum(
                    Q(math.comb(4, k), 1)/tau**k for k in range(min(w-1, cap-1)+1))
                self.assertLessEqual(old[w][cap][0], pure)
                self.assertLessEqual(new[w][cap][0], old[w][cap][0])
        # A singleton return necessarily has all four incoming coordinates live;
        # its only selected-zero term is the empty selection.
        self.assertEqual(old[1][4][0], 8*y**8/tau**4)
        self.assertLess(new[1][1][0], old[1][1][0])

    def test_eight_state_algebra_and_cone(self):
        rng = Random(889640)
        f, a, c, eta = fullstate_matrices()
        extended = tuple(f[i]+tuple(MUL[eta][v] for v in a[i]) for i in range(8))
        extended += tuple(c[i]+tuple(int(i == j) for j in range(8)) for i in range(8))
        for _ in range(128):
            w, z, h = rng.randrange(9), rng.randrange(9), rng.randrange(9)
            inputs = tuple(sorted(rng.sample(range(8), w)))
            outputs = tuple(sorted(rng.sample(range(8), z)))
            state = tuple(sorted(rng.sample(range(8), h)))
            for returning in (False, True):
                rows = outputs+(tuple(range(8, 16)) if returning else ())
                columns = inputs+tuple(8+i for i in state)
                actual = rank(tuple(tuple(extended[i][j] for j in columns) for i in rows))
                self.assertEqual(actual, min(z+8*returning, w+h))
        for rank_caps in (False, True):
            local = weighted_matrices(Q(2, 3), Q(13), state=8, rank_caps=rank_caps)
            for w, matrix in enumerate(local):
                self.assertEqual(matrix[0][0], int(w == 0))
                self.assertTrue(cone_preserved(matrix))

    def test_sixteen_state_algebra_and_cone(self):
        rng = Random(8169640)
        f, a, c, eta = largestate_matrices()
        extended = tuple(f[i]+tuple(MUL[eta][v] for v in a[i]) for i in range(8))
        extended += tuple(c[i]+tuple(int(i == j) for j in range(16)) for i in range(16))
        for _ in range(128):
            w, z, h = rng.randrange(9), rng.randrange(9), rng.randrange(17)
            inputs = tuple(sorted(rng.sample(range(8), w)))
            outputs = tuple(sorted(rng.sample(range(8), z)))
            state = tuple(sorted(rng.sample(range(16), h)))
            for returning in (False, True):
                rows = outputs+(tuple(range(8, 24)) if returning else ())
                columns = inputs+tuple(8+i for i in state)
                actual = rank(tuple(tuple(extended[i][j] for j in columns) for i in rows))
                self.assertEqual(actual, min(z+16*returning, w+h))
        for rank_caps in (False, True):
            local = weighted_matrices(Q(2, 3), Q(13), state=16, rank_caps=rank_caps)
            for w, matrix in enumerate(local):
                self.assertEqual(matrix[0][0], int(w == 0))
                self.assertTrue(cone_preserved(matrix))


if __name__ == '__main__':
    unittest.main()
