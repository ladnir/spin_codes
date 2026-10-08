"""Exact algebra and small-field checks for the changed feedthrough proof."""
from fractions import Fraction as F
from itertools import product
import math
from random import Random
import unittest

import numpy as np

from gf8_feedthrough_screen import (BIRTH_MINIMA, RETURN_MINIMA, SUBSETS,
    eliminate_inputs, full_support_kernel, packet, rank_profiles, root_potential,
    state_profile, translate_mask, weighted_matrices, root_cap_weights, polynomial_product)
from gf8_packet_maps import (MUL, cubic_matrices, exact_return_minima,
                            matmul, rank)
from gf8_rank_history_screen import cone_preserved


GF4 = (0, 1, 188, 189)


class FeedthroughTests(unittest.TestCase):
    def test_exact_maps_symmetry_and_branch_minima(self):
        reference = None
        for phase in range(8):
            a, c = cubic_matrices(phase)
            self.assertEqual(matmul(c, a), tuple((0,)*4 for _ in range(4)))
            ac = matmul(a, c)
            if reference is None:
                reference = ac
            self.assertEqual(ac, reference)
            for shift in range(4):
                for i, j in product(range(8), range(4)):
                    self.assertEqual(a[i ^ shift][j], a[i][j ^ shift])
                    self.assertEqual(c[j][i ^ shift], c[j ^ shift][i])
        for scalar, expected in ((4, BIRTH_MINIMA), (6, RETURN_MINIMA)):
            changed = tuple(tuple(int(i == j) ^ MUL[scalar][reference[i][j]] for j in range(8))
                            for i in range(8))
            self.assertEqual(tuple(exact_return_minima(changed)), expected)
            self.assertEqual(matmul(changed, changed),
                             tuple(tuple(int(i == j) for j in range(8)) for i in range(8)))

    def test_exact_elimination_rank_identity(self):
        rng = Random(484096)
        for _ in range(512):
            a, c = cubic_matrices(rng.randrange(8))
            ac = matmul(a, c)
            inputs, zeros = SUBSETS[rng.randrange(256)], SUBSETS[rng.randrange(256)]
            left = tuple(tuple(int(z == i) ^ MUL[4][ac[z][i]] for i in inputs) for z in zeros)
            right = tuple(tuple(MUL[2][v] for v in a[z]) for z in zeros)
            r, g = eliminate_inputs(left, right)
            self.assertEqual(r, rank(left))
            self.assertEqual(rank(tuple(l+r for l, r in zip(left, right))), r+rank(g))

    def test_circuit_probability_and_live_rank_moment_over_gf4(self):
        # Varied nonsuperregular G: zero columns, parallel columns, and dense kernels.
        matrices = ((), ((1, 0, 1, 1),), ((1, 1, 1, 1),),
                    ((1, 0, 1, 1), (0, 1, 1, 188)),
                    ((1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)))
        rng = Random(44256)
        cases = 0
        diagonals = tuple(product(GF4[1:], repeat=4))
        for g in matrices:
            state_rank, charge = state_profile(g)
            for _ in range(20):
                columns = rng.randrange(1, 5)
                s = tuple(tuple(rng.choice(GF4) for _ in range(columns)) for _ in range(4))
                incoming = rank(s)
                if not incoming:
                    continue
                for tau in (F(1), F(2), F(3)):
                    total = F(0)
                    tails = [F(0)]*4
                    for diagonal in diagonals:
                        refreshed = tuple(tuple(MUL[d][v] for v in row) for d, row in zip(diagonal, s))
                        charged = rank(matmul(g, refreshed)) if g else 0
                        if charged < incoming:
                            total += tau**(-charged)
                            for threshold in range(1, incoming-charged+1):
                                tails[threshold-1] += tau**(-charged)
                    self.assertLessEqual(total/len(diagonals), root_potential(tau, state_rank, charge, q=4))
                    masses = root_cap_weights(tau, incoming, state_rank, charge, q=4)
                    for first in range(4):
                        self.assertLessEqual(tails[first]/len(diagonals), sum(masses[first:]))
                    cases += 1
        self.assertGreaterEqual(cases, 270)

    def test_support_profile_completeness_and_singleton_birth(self):
        continuation, birth = rank_profiles()
        for weight in range(9):
            self.assertEqual(sum(continuation[weight].values()), math.comb(8, weight)*256)
        self.assertEqual(dict(birth[1]), {(0, 0): 8})
        for translation in range(4):
            self.assertEqual({translate_mask(i, translation) for i in range(256)}, set(range(256)))

    def test_exact_packet_and_zero_input_return(self):
        x, y, tau = F(2, 5), F(3, 4), F(11)
        local = weighted_matrices(y, tau)
        self.assertEqual(local[0][0], [1, 0])
        self.assertEqual(local[0][1][0], 0)
        self.assertEqual(local[1][0][1], 8*y**8)
        self.assertEqual(local[1][1][0], 8*y**8/tau**4)
        for weight in range(1, 5):
            self.assertEqual(local[weight][0][0], 0)
        expected = [[sum(x**w*local[w][i][j] for w in range(9)) for j in range(2)] for i in range(2)]
        np.testing.assert_allclose(np.exp(packet(math.log(x), math.log(y), math.log(tau))),
                                   np.asarray(expected, dtype=float), rtol=2e-12)

    def test_rank_cap_successor_and_continuation_cone(self):
        rng = Random(424096)
        for _ in range(256):
            a, c = cubic_matrices(rng.randrange(8))
            ac = matmul(a, c)
            inputs, zeros = SUBSETS[rng.randrange(1, 256)], SUBSETS[rng.randrange(256)]
            b = tuple(tuple(int(z == i) ^ MUL[4][ac[z][i]] for i in inputs) for z in zeros)
            emission = tuple(tuple(MUL[2][v] for v in a[z]) for z in zeros)
            r, g = eliminate_inputs(b, emission)
            dimension = rng.randrange(1, 9)
            # These forms already include arbitrary fixed nonzero diagonal values.
            state = tuple(tuple(rng.randrange(256) for _ in range(dimension)) for _ in range(4))
            x = tuple(tuple(rng.randrange(256) for _ in range(dimension)) for _ in inputs)
            state_equations = matmul(g, state) if g else ()
            h = rank(state)-rank(state_equations)
            bx = matmul(b, x) if b else ()
            ar = matmul(emission, state) if emission else ()
            output = tuple(tuple(u ^ v for u, v in zip(left, right)) for left, right in zip(bx, ar))
            cx = matmul(tuple(tuple(row[i] for i in inputs) for row in c), x)
            successor = tuple(tuple(u ^ v for u, v in zip(left, right)) for left, right in zip(state, cx))
            actual = rank(output+successor)-rank(output)
            self.assertLessEqual(actual, min(4, h+len(inputs)-r))
        for y, tau in ((F(2, 3), F(11)), (F(99, 100), F(255)), (F(1), F(2))):
            for matrix in weighted_matrices(y, tau, rank_caps=True):
                self.assertTrue(cone_preserved(matrix))

    def test_generic_polynomial_product(self):
        rng = np.random.default_rng(5408)
        left, right = rng.uniform(.1, 1, (3, 5, 5)), rng.uniform(.1, 1, (4, 5, 5))
        expected = np.zeros((6, 5, 5))
        for i in range(3):
            for j in range(4):
                expected[i+j] += left[i]@right[j]
        np.testing.assert_allclose(np.exp(polynomial_product(np.log(left), np.log(right), 5)), expected,
                                   rtol=2e-12)


if __name__ == '__main__':
    unittest.main()
