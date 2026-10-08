"""Exact bounded tests for the scalar state-capacity repair, not a certificate."""
from fractions import Fraction
from itertools import combinations, product
from random import Random
import unittest

import state_capacity_reference as ref


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


G4 = tuple(tuple(gf4_mul(a, b) for b in range(4)) for a in range(4))


def rank4(matrix):
    rows = [list(row) for row in matrix]
    found = 0
    for j in range(3):
        pivot = next((i for i in range(found, len(rows)) if rows[i][j]), None)
        if pivot is None:
            continue
        rows[pivot], rows[found] = rows[found], rows[pivot]
        inverse = next(x for x in range(1, 4) if G4[x][rows[found][j]] == 1)
        rows[found] = [G4[inverse][x] for x in rows[found]]
        for i in range(found+1, len(rows)):
            factor = rows[i][j]
            rows[i] = [x ^ G4[factor][y] for x, y in zip(rows[i], rows[found])]
        found += 1
    return found


def dot4(a, b):
    return ref.xor_sum(G4[x][y] for x, y in zip(a, b))


class StateCapacityTests(unittest.TestCase):
    def test_fixed_maps_and_asymmetric_adjoint(self):
        self.assertEqual(ref.check_maps(), 160)

    def test_symmetric_maps_involution_and_adjoint(self):
        rng = Random(803822)
        for phase in range(4):
            a, c = ref.matrices(phase, True)
            for j in range(8):
                state = tuple(int(i == j) for i in range(8))
                self.assertEqual(ref.injection(ref.emission(state, phase), phase, True), (0,)*8)
            for _ in range(8):
                x = tuple(rng.randrange(256) for _ in range(16))
                state = tuple(rng.randrange(256) for _ in range(8))
                diagonal = tuple(rng.randrange(1, 256) for _ in range(8))
                z = tuple(rng.randrange(256) for _ in range(16))
                terminal = tuple(rng.randrange(256) for _ in range(8))
                y, following = ref.step(x, state, diagonal, phase, True)
                self.assertEqual(ref.injection(x, phase, True), ref.matvec(c, x))
                restored, repeated = ref.step(y, state, diagonal, phase, True)
                self.assertEqual((restored, repeated), (x, following))
                dx, ds = ref.adjoint(z, terminal, diagonal, phase, True)
                dot = lambda u, v: ref.xor_sum(ref.mul(a, b) for a, b in zip(u, v))
                self.assertEqual(dot(y, z) ^ dot(following, terminal), dot(x, dx) ^ dot(state, ds))

    def test_cauchy_identity_and_symmetric_singleton_return(self):
        for phase in range(4):
            a, c = ref.matrices(phase, True)
            for j, u in enumerate(ref.LABELS[phase]):
                state = tuple(row[j] for row in c)
                output = tuple(int(i == j) ^ ref.mul(ref.ETA, value)
                               for i, value in enumerate(ref.matvec(a, state)))
                self.assertEqual(sum(bool(x) for x in output), 9)
                for row, v in zip(c, ref.NODES):
                    self.assertEqual(ref.mul(row[j], u ^ v), ref.linearized(u ^ ref.BETA))

    def test_asymmetric_eta_gate_has_no_nonzero_solution(self):
        orbits = ref.eta_obstruction_orbits()
        self.assertEqual(tuple(map(len, orbits)), (120, 120))
        self.assertEqual(ref.subset_update_rank((0, 1, 4, 10, 12), 2), 4)
        self.assertEqual(ref.subset_update_rank((0, 3, 8, 12), 7), 3)

    def test_asymmetric_two_minors(self):
        table = ref.multiplication_table()
        for phase in range(4):
            _, c = ref.matrices(phase)
            for j, k in combinations(range(8), 2):
                for a, b in combinations(range(16), 2):
                    self.assertNotEqual(table[c[j][a]][c[k][b]] ^ table[c[j][b]][c[k][a]], 0)

    def test_pure_return_cover_and_product_counts(self):
        y, tau = Fraction(9, 10), Fraction(35)
        self.assertEqual(ref.symmetric_return_envelope(0, y, tau), 0)
        self.assertEqual(ref.symmetric_return_envelope(1, y, tau), y**9/tau**8)
        self.assertEqual(ref.symmetric_return_envelope(2, y, tau), (1+8/tau)*y**7/tau**8)
        self.assertEqual(ref.symmetric_return_envelope(3, y, tau), (1+8/tau+28/tau**2)*y**6/tau**8)
        self.assertEqual(ref.symmetric_return_envelope(2, y, tau, flat_cover=False), 9*y**7/tau**8)
        self.assertEqual(ref.symmetric_return_envelope(3, y, tau, flat_cover=False), 37*y**6/tau**8)
        self.assertEqual(ref.symmetric_return_envelope(1, y, tau, 7), 0)
        self.assertEqual(ref.symmetric_return_envelope(2, y, tau, 7), y**7/tau**8)
        counts = tuple(tuple(sum(x) for x in phase) for phase in ref.shared_product_counts())
        self.assertEqual(counts, ((30, 30), (22, 28), (29, 30), (30, 30)))

    def test_flat_cover_retains_deterministic_equation_rank(self):
        from math import comb
        f, g, h = (1, 0, 0), (0, 1, 0), (0, 0, 1)
        cases = ((f, f, f, f), (f, g, f, g),
                 (f, g, h, (1, 1, 1)), ((0, 0, 0), f, g, (1, 1, 0)))
        for states in cases:
            rank = rank4(states)
            zeros = {tuple(i for i, form in enumerate(states) if dot4(form, p) == 0)
                     for p in product(range(4), repeat=3)}
            zeros.discard(tuple(range(4)))
            # Every realized proper zero set is closed under linear span.
            for z in zeros:
                rz = rank4([states[i] for i in z])
                closure = tuple(i for i, form in enumerate(states)
                                if rank4([states[j] for j in z]+[form]) == rz)
                self.assertEqual(z, closure)
            for weight in range(1, 6):
                for tau in (1, 2, 3):
                    actual = sum(Fraction(1, tau**rank4([states[i] for i in z]))
                                 for z in zeros if len(z) < weight)
                    upper = sum(Fraction(comb(4, k), tau**k)
                                for k in range(min(weight-1, rank-1)+1))
                    self.assertLessEqual(actual, upper)

    def test_return_pivots_with_dependent_source_forms_gf4(self):
        # Two-state analogue, numeric C superregular. Enumerate actual alive
        # source witnesses, not just coefficient ranks. State and input forms
        # can coincide or have other dependencies in the common quotient.
        f, g, h = (1, 0, 0), (0, 1, 0), (0, 0, 1)
        cases = [((f, f), (f,), ((1,), (1,))),
                 ((f, g), (f, g), ((1, 1), (1, 2))),
                 ((f, g), (h, f), ((1, 1), (1, 2))),
                 ((f,), (f, g), ((1, 1), (1, 2))),
                 ((f,), (g,), ((1,), (1,)))]
        points = tuple(product(range(4), repeat=3))
        for states, inputs, c in cases:
            live_points = [p for p in points if all(dot4(form, p) for form in states+inputs)]
            histogram = {}
            for diagonal in product(range(1, 4), repeat=len(states)+len(inputs)):
                state_d = diagonal[:len(states)]
                input_d = diagonal[len(states):]
                equations = []
                for row in range(2):
                    equation = [0, 0, 0]
                    if row < len(states):
                        equation = [G4[state_d[row]][x] for x in states[row]]
                    for coefficient, d, form in zip(c[row], input_d, inputs):
                        scalar = G4[coefficient][d]
                        equation = [x ^ G4[scalar][v] for x, v in zip(equation, form)]
                    equations.append(tuple(equation))
                if any(not any(dot4(row, p) for row in equations) for p in live_points):
                    r = rank4(equations)
                    histogram[r] = histogram.get(r, 0)+1
            denom = 3**(len(states)+len(inputs))
            for tau in (1, 2, 3):
                actual = sum(Fraction(count, denom*tau**r) for r, count in histogram.items())
                expected = Fraction(int(len(inputs) > 2-len(states)), tau**2)
                self.assertLessEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
