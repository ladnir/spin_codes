"""Small exact checks for uniform-entry state refresh; no global screen."""

from fractions import Fraction as F
from itertools import combinations, product
from random import Random
import unittest

from cauchy_joint_return import MUL, dot, matrix_rank, packet_forms
from full_state_refresh import (canonical_factor, empty_nonreturn_action,
                                marginal_rank_action, return_rank_action,
                                return_terms)


def gf4_mul(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return out


G4 = [[gf4_mul(a, b) for b in range(4)] for a in range(4)]


def rank4(rows, width):
    rows = [list(row) for row in rows]
    rank = 0
    for column in range(width):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = next(x for x in range(1, 4) if G4[x][rows[rank][column]] == 1)
        rows[rank] = [G4[inv][x] for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            factor = rows[i][column]
            rows[i] = [a ^ G4[factor][b] for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def add4(a, b):
    return tuple(x ^ y for x, y in zip(a, b))


def linear4(coefficients, forms, width):
    out = (0,)*width
    for scalar, form in zip(coefficients, forms):
        out = add4(out, tuple(G4[scalar][x] for x in form))
    return out


class FullStateRefreshTests(unittest.TestCase):
    def test_canonical_factor_and_empty_return(self):
        tau, y = F(35), F(9997, 10000)
        self.assertEqual(canonical_factor(tau), sum(tau**(-j) for j in range(4)))
        self.assertLess(canonical_factor(tau), F(103, 100))
        self.assertEqual(return_rank_action(0, y, tau), canonical_factor(tau)/tau**4)
        self.assertEqual(return_rank_action(1, y, tau), canonical_factor(tau)*y**16/tau**4)
        self.assertEqual(canonical_factor(tau, 0), 0)

    def test_exact_two_input_feasible_zero_sets(self):
        # Independently enumerate all projective torus points for three
        # representative placements; arrangements are not reused here.
        for phase, inputs in ((0, (0, 1)), (1, (0, 7)), (2, (3, 12))):
            _, output = packet_forms(phase, inputs)
            feasible = {0}
            for scalar in range(1, 256):
                mask = sum((dot(row, (1, scalar)) == 0) << j for j, row in enumerate(output))
                subset = mask
                while subset:
                    feasible.add(subset)
                    subset = (subset-1) & mask
            expected = {}
            for mask in feasible:
                key = (mask.bit_count(), matrix_rank([row for j, row in enumerate(output) if mask >> j & 1], 2))
                expected[key] = expected.get(key, 0)+1
            self.assertEqual(return_terms(phase, inputs), tuple((z, rank, count) for (z, rank), count in sorted(expected.items())))

    def test_selected_zero_marginal_coefficient_ranks(self):
        from cauchy_diagonal_packet_reference import LABELS, NODES
        from rate_design_coset import power
        for inputs in ((0,), (0, 5), (1, 7, 12)):
            for phase in range(4):
                for zeros in ((0,), (0, 1, 2), (0, 1, 2, 3, 4), tuple(range(16))):
                    rows = [[MUL[2][power(LABELS[phase][j] ^ v, 3) ^ 1] for v in NODES]
                            +[int(j == i) for i in inputs] for j in zeros]
                    i = len(set(inputs).intersection(zeros))
                    z = len(zeros)-i
                    self.assertEqual(matrix_rank(rows, 4+len(inputs)), i+min(z, 4))

    def test_uniform_matrix_return_rank_action_gf4(self):
        # Exhaust all256 two-state refresh matrices and all nonzero input
        # diagonals.  Source forms may be dependent and overlap the state.
        rng, width = Random(2164401), 4
        forms = [tuple(rng.randrange(4) for _ in range(width)) for _ in range(12)]
        forms = [f for f in forms if any(f)]
        cases = [((forms[0], forms[0]), (forms[1],), ((1,),)),
                 ((forms[0], forms[1]), (forms[2],), ()),
                 ((forms[0], forms[1]), (forms[0], forms[2]), ((1, 1),)),
                 ((forms[2], forms[3]), (forms[4], forms[5]), ((1, 2),))]
        for state, inputs, selected in cases:
            r = rank4(state, width)
            b_rank = rank4(selected, len(inputs))
            histogram = {}
            for diagonal in product(range(1, 4), repeat=len(inputs)):
                values = [tuple(G4[d][x] for x in f) for d, f in zip(diagonal, inputs)]
                zero_eqs = [linear4(row, values, width) for row in selected]
                for flat in product(range(4), repeat=4):
                    refreshed = [linear4(flat[2*j:2*j+2], state, width) for j in range(2)]
                    # Any deterministic C is allowed in this joint lemma.
                    following = [add4(refreshed[j], linear4(tuple((j+i+1)%4 for i in range(len(inputs))), values, width))
                                 for j in range(2)]
                    equations = zero_eqs+following
                    rank = rank4(equations, width)
                    live_inputs = all(rank4(equations+[f], width) > rank for f in inputs)
                    live_state = any(rank4(equations+[f], width) > rank for f in state)
                    if live_inputs and live_state:
                        histogram[rank] = histogram.get(rank, 0)+1
            denominator = 256*3**len(inputs)
            for tau in (1, 2, 3):
                actual = sum(F(count, denominator*tau**rank) for rank, count in histogram.items())
                upper = sum(F(1, tau**j) for j in range(r))/tau**(2+b_rank)
                self.assertLessEqual(actual, upper)

    def test_empty_nonreturn_does_not_include_zero_image(self):
        y, tau = F(2, 3), F(7)
        self.assertLess(empty_nonreturn_action(y, tau), marginal_rank_action(0, y, tau))
        self.assertGreater(return_rank_action(0, y, tau), 0)


if __name__ == '__main__':
    unittest.main()
