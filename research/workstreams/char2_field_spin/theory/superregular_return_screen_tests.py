"""Independent local-cover and regional-composition checks, not a certificate.

The return path is a pure-return cover with a deterministic weight mark.
It is intentionally NOT compared with a sum of selected-zero return paths.
The latter is a different, usually more expensive cover of the same witness.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb, exp, log
from unittest.mock import patch
import unittest

import numpy as np

import superregular_return_screen as screen
from diagonal_refresh_rank import rank_potential_bound
import state_capacity_reference as maps


def matrix_product(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)] for i in range(2)]


def exact_region(active, y, tau, rows, state, flat_cover=True):
    """Sequential positive Fraction convolution, without logarithmic helpers."""
    coefficients = [[[F(1), F(0)], [F(0), F(1)]]]
    local = [screen.local_matrix(w, y, tau, state, flat_cover=flat_cover)
             for w in range(min(16, active)+1)]
    for _ in range(rows//16):
        following = [[[F(0), F(0)], [F(0), F(0)]] for _ in range(min(active+1, len(coefficients)+len(local)-1))]
        for degree, matrix in enumerate(coefficients):
            for w, packet in enumerate(local[:len(following)-degree]):
                term = matrix_product(matrix, packet)
                for i in range(2):
                    for j in range(2):
                        following[degree+w][i][j] += comb(16, w)*term[i][j]
        coefficients = following
    result = [[F(0), F(0)], [F(0), F(0)]]
    for degree, matrix in enumerate(coefficients):
        factor = F(comb(active, degree), comb(rows, degree))*tau**degree
        for i in range(2):
            for j in range(2):
                result[i][j] += factor*matrix[i][j]
    return result


def mul16(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 16:
            a ^= 19
    return result


G16 = tuple(tuple(mul16(a, b) for b in range(16)) for a in range(16))
F0, F1, F2, ZERO = (1, 0), (0, 1), (1, 1), (0, 0)


def add_forms(left, right):
    return left[0] ^ right[0], left[1] ^ right[1]


def scale_form(scalar, form):
    return G16[scalar][form[0]], G16[scalar][form[1]]


def det(left, right):
    return G16[left[0]][right[1]] ^ G16[left[1]][right[0]]


def rank2(forms):
    first = next((row for row in forms if row != ZERO), None)
    if first is None:
        return 0
    return 2 if any(det(first, row) for row in forms) else 1


def alive(forms, required, nonzero_vectors=()):
    """Exact feasibility in a two-dimensional GF16 source quotient.

    At most six nonvanishing requirements occur in these cases, so avoiding
    their proper hyperplanes is possible whenever none is forced to vanish.
    """
    rank = rank2(forms)
    if rank == 2:
        return False
    if rank == 0:
        return all(form != ZERO for form in required) and all(any(row != ZERO for row in vector)
                                                               for vector in nonzero_vectors)
    first = next(row for row in forms if row != ZERO)
    return all(det(first, form) for form in required) and all(any(det(first, row) for row in vector)
                                                             for vector in nonzero_vectors)


def small_packet_histograms(states, inputs):
    """Exhaust all fresh GF16 diagonals for a symmetric 4x2 Cauchy packet."""
    labels, nodes, eta = (0, 1, 6, 7), (2, 3), 2
    a = tuple(tuple(u ^ v ^ 1 for v in nodes) for u in labels)
    c = tuple(zip(*a))
    w = len(inputs)
    incoming = bool(rank2(states))
    nonreturn = Counter()
    returned = Counter()
    births = 0
    for diagonal in product(range(1, 16), repeat=2+w):
        r = tuple(scale_form(d, form) for d, form in zip(diagonal[:2], states))
        x = tuple(scale_form(d, form) for d, form in zip(diagonal[2:], inputs))
        following = []
        for j in range(2):
            value = r[j]
            for coefficient, form in zip(c[j], x):
                value = add_forms(value, scale_form(coefficient, form))
            following.append(value)
        output = []
        for j in range(4):
            value = add_forms(scale_form(a[j][0], r[0]), scale_form(a[j][1], r[1]))
            value = scale_form(eta, value)
            if j < w:
                value = add_forms(value, x[j])
            output.append(value)
        if not incoming:
            if alive((), inputs, (following,)):
                births += 1
            if alive(following, inputs):
                returned[rank2(following)] += 1
            continue
        # Continuing branch: selected-output-zero cover, with actual incoming
        # and next state both nonzero on a surviving source witness.
        for mask in range(16):
            equations = [output[j] for j in range(4) if mask >> j & 1]
            if alive(equations, inputs, (states, following)):
                nonreturn[rank2(equations), mask.bit_count()] += 1
        # Returning branch: exact incoming state-zero flats plus return
        # equations. No selected-output-zero branching is added here.
        for support in range(1, 4):
            equations = [states[j] for j in range(2) if not support >> j & 1]+following
            required = tuple(inputs)+tuple(states[j] for j in range(2) if support >> j & 1)
            if alive(equations, required):
                returned[rank2(equations)] += 1
                minimum = max(1, 3-w)
                # Independently check the deterministic return weight on all
                # projective source points; these do not weight the rank cover.
                for point in ((0, 1),)+tuple((1, t) for t in range(16)):
                    dot = lambda form: G16[form[0]][point[0]] ^ G16[form[1]][point[1]]
                    if all(not dot(row) for row in equations) and all(dot(row) for row in required):
                        assert sum(bool(dot(row)) for row in output) >= minimum
    return nonreturn, returned, births, 15**(2+w)


class SuperregularReturnScreenTests(unittest.TestCase):
    def test_symmetric_actual_singleton_mark(self):
        y, tau = F(9, 10), F(35)
        for phase in range(4):
            a, c = maps.matrices(phase, True)
            for j in range(16):
                state = tuple(row[j] for row in c)
                output = tuple(int(i == j) ^ maps.mul(2, value)
                               for i, value in enumerate(maps.matvec(a, state)))
                self.assertEqual(sum(bool(value) for value in output), 9)
        self.assertEqual(screen.local_matrix(1, y, tau)[1][0], y**9/tau**8)

    def test_independent_marginal_factorization(self):
        y, tau = F(7, 9), F(19, 2)
        for state in (2, 5, 8):
            roots = screen.gamma(tau, state)
            for w in range(17):
                expected = (y+(1-y)/tau)**w*sum(comb(16-w, z)*y**(16-w-z)*(1-y)**z*roots[z]
                                                     for z in range(min(state-1, 16-w)+1))
                matrix = screen.local_matrix(w, y, tau, state)
                self.assertEqual(matrix[1][1], expected)
                self.assertEqual(matrix[0][0], F(1) if not w else y**w/tau**state if w > state else F(0))
                self.assertEqual(matrix[0][1], y**w if w else 0)

    def test_raw_mask_and_weighted_flat_counts(self):
        for state in (2, 5, 8):
            for w in range(1, 17):
                actual_masks = sum(1 for mask in range(1, 1 << state) if state-mask.bit_count() < w)
                self.assertEqual(screen.return_branches(w, state), actual_masks)
                self.assertEqual(screen.return_flat_factor(w, F(1), state), actual_masks)
                for tau in (F(2), F(255)):
                    self.assertLessEqual(screen.return_flat_factor(w, tau, state), actual_masks)
                    self.assertEqual(screen.return_flat_factor(w, tau, state),
                                     sum(F(comb(state, k), tau**k) for k in range(min(w-1, state-1)+1)))

    def test_small_symmetric_maps_are_superregular_and_ca_zero(self):
        a = tuple(tuple(u ^ v ^ 1 for v in (2, 3)) for u in (0, 1, 6, 7))
        self.assertTrue(all(all(row) for row in a))
        self.assertTrue(all(det(a[i], a[j]) for i, j in combinations(range(4), 2)))
        for i, j in product(range(2), repeat=2):
            value = 0
            for row in a:
                value ^= G16[row[i]][row[j]]
            self.assertEqual(value, 0)

    def test_small_exhaustive_alive_covers_all_four_entries(self):
        cases = (((F0, F1), ()), ((F0, F1), (F0,)), ((F0, F0), (F1,)),
                 ((F0, ZERO), (F1,)), ((F0, F1), (F2, F2)),
                 ((F0, F1), (F0, F1)), ((ZERO, ZERO), (F0, F1, F2)))
        comparisons = 0
        def small_roots(tau, state, q=16):
            return tuple(max(rank_potential_bound(r, z, tau, q, state) for r in range(1, state+1))
                         for z in range(state))
        for states, inputs in cases:
            marginal, returned, births, denominator = small_packet_histograms(states, inputs)
            w, incoming = len(inputs), int(bool(rank2(states)))
            for y, tau in product((F(1, 2), F(9, 10)), (F(1), F(2), F(15))):
                with patch.object(screen, 'gamma', small_roots):
                    upper = screen.local_matrix(w, y, tau, state=2, width=4)
                if incoming:
                    actual_nonreturn = sum(F(count, denominator)*y**(4-z)*(1-y)**z/tau**r
                                           for (r, z), count in marginal.items())
                    actual_return = sum(F(count, denominator)*y**max(1, 3-w)/tau**r
                                        for r, count in returned.items())
                else:
                    actual_nonreturn = F(births, denominator)*y**w
                    actual_return = sum(F(count, denominator)*y**w/tau**r for r, count in returned.items())
                self.assertLessEqual(actual_nonreturn, upper[incoming][1], (states, inputs, y, tau, 'nonreturn'))
                self.assertLessEqual(actual_return, upper[incoming][0], (states, inputs, y, tau, 'return'))
                comparisons += 2
        self.assertEqual(comparisons, 84)

    def test_regional_log_transfer_matches_positive_fraction_convolution(self):
        y, tau = F(9, 10), F(7)
        for active, rows, flat_cover in ((0, 16, True), (1, 32, True), (3, 32, True), (2, 48, False)):
            expected = exact_region(active, y, tau, rows, 8, flat_cover)
            actual = np.exp(screen.region(active, float(y), float(tau), rows, 8, flat_cover))
            np.testing.assert_allclose(actual, np.asarray(expected, dtype=float), rtol=2e-13, atol=1e-300)

    def test_without_replacement_route_factor_exhaustively(self):
        length, active = 4, 2
        routes = tuple(permutations(range(length)))
        for mask in range(1 << length):
            occupied = [j for j in range(length) if mask >> j & 1]
            hits = sum(all(route[j] < active for j in occupied) for route in routes)
            expected = F(comb(active, len(occupied)), comb(length, len(occupied))) if len(occupied) <= active else F(0)
            self.assertEqual(F(hits, len(routes)), expected)

    def test_global_expression_matches_direct_matrix_composition(self):
        active, length, source_dimension, klog = 2, 3, 2, 4
        logtilt, logtau = 1.5, log(7)
        geometry = screen.geometry(klog, length, source_dimension)
        logy = -exp(logtilt)*active/geometry['n']
        table = np.exp(screen.region(active, exp(logy), exp(logtau), geometry['rows']))
        vector = np.asarray((1., 0.))
        for _ in range(length):
            vector = vector @ table
        expected = (log(comb(geometry['rows'], active))
                    -((length-source_dimension)*active+1)*logtau
                    -geometry['threshold']*logy+log(float(vector.sum())))/log(2)
        self.assertAlmostEqual(screen.bound(active, (logtilt, logtau), klog, length, source_dimension), expected, places=11)


if __name__ == '__main__':
    unittest.main()
