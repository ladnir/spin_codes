"""Independent math checks for an all-path FLOATING diagnostic, not a certificate.

The exact toy-field enumerations check the ingredients. Floating comparisons
check their implementation, including fractional powers and route conditioning.
No quiet-path subtraction, native encoder timing, or directed acceptance occurs.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb, exp2, fsum, isfinite, log2
from unittest.mock import patch
import unittest

import hd_randomized_outer_state_screen as screen


def field_multiply(a, b, q):
    """Independent polynomial arithmetic for the exhaustive GF4/GF8 checks."""
    modulus = {2: 3, 4: 7, 8: 11}[q]
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & q:
            a ^= modulus
    return result


def polynomial_word(coefficients, q, length):
    result = []
    for x in range(length):
        value = 0
        for coefficient in reversed(coefficients):
            value = field_multiply(value, x, q) ^ coefficient
        result.append(value)
    return tuple(result)


def exact_tilt(words, z):
    histogram = Counter(sum(x != 0 for x in word) for word in words)
    return sum((count*z**weight for weight, count in histogram.items()), F(0))/len(words)


def raw_multiply(a, b):
    return tuple(tuple(sum((a[i][k]*b[k][j] for k in range(2)), F(0))
                       for j in range(2)) for i in range(2))


def raw_power(a, n):
    result = ((F(1), F(0)), (F(0), F(1)))
    for unused in range(n):
        result = raw_multiply(result, a)
    return result


def as_logs(matrix):
    return tuple(tuple(log2(value) if value else screen.NEG for value in row)
                 for row in matrix)


def one_at_marker(local, marker, alpha=1):
    """Ordinary positive arithmetic, independent of the screen's log sums."""
    return tuple(tuple(fsum(comb(64, j)*marker**j*(1-marker)**(64-j)
                             *exp2(alpha*local[j][a][b]) for j in range(65))
                       for b in range(2)) for a in range(2))


def pair_at_marker(coefficients, marker):
    return tuple(tuple(fsum(exp2(coefficients[j][a][b])*marker**j
                             *(1-marker)**(128-j) for j in range(129))
                       for b in range(2)) for a in range(2))


class RandomizedOuterStateScreenTests(unittest.TestCase):
    def assert_log_matches(self, actual, exact, places=10):
        if not exact:
            self.assertEqual(actual, screen.NEG)
        else:
            self.assertAlmostEqual(actual, log2(exact), places=places)

    def test_mds_moments_by_exact_gf8_polynomial_enumeration(self):
        for dimension in range(4):
            words = tuple(polynomial_word(coefficients, 8, 8)
                          for coefficients in product(range(8), repeat=dimension))
            self.assertEqual(len(set(words)), 8**dimension)
            for mark in (F(1, 4), F(1, 2), F(3, 4)):
                with self.subTest(dimension=dimension, mark=mark):
                    self.assert_log_matches(screen.mds_moment(dimension, log2(mark), 8),
                                            exact_tilt(words, mark))
        # A full-dimensional tile has eight independent uniform symbols.
        for q in (8, 256, 65536):
            for mark in (F(1, 4), F(3, 4)):
                exact = ((1+(q-1)*mark)/q)**8
                self.assert_log_matches(screen.mds_moment(8, log2(mark), q), exact)

    def test_positive_fourier_and_all_gf4_affine_cosets(self):
        q, z = 4, F(1, 3)
        # The absolute trace in GF4 is x+x^2. Enumerate every character.
        for character in range(q):
            transform = F(0)
            for x in range(q):
                value = field_multiply(character, x, q)
                trace = value ^ field_multiply(value, value, q)
                self.assertIn(trace, (0, 1))
                transform += (-1)**trace*(z if x else 1)
            self.assertEqual(transform, 1+(q-1)*z if character == 0 else 1-z)
            self.assertGreaterEqual(transform, 0)
        code = tuple(polynomial_word(coefficients, q, 4)
                     for coefficients in product(range(q), repeat=2))
        zero_coset = exact_tilt(code, z)
        strict = 0
        for translate in product(range(q), repeat=4):
            coset = tuple(tuple(x ^ t for x, t in zip(word, translate)) for word in code)
            moment = exact_tilt(coset, z)
            self.assertLessEqual(moment, zero_coset)
            if translate in code:
                self.assertEqual(moment, zero_coset)
            strict += moment < zero_coset
        self.assertGreater(strict, 0)

    def test_rank_one_mixer_atoms_and_subspaces_by_enumeration(self):
        # Q=e_0. No uniform-nonzero-state approximation enters the counts.
        for q in (2, 4):
            vectors = tuple(product(range(q), repeat=4))
            products = tuple(tuple(field_multiply(a, b, q) for b in range(q))
                             for a in range(q))
            counts = Counter()
            for u in vectors[1:]:
                for v in vectors:
                    pairing = 0
                    for x, y in zip(u, v):
                        pairing ^= products[x][y]
                    if pairing == 1:
                        continue
                    successor = tuple((i == 0) ^ products[u[i]][v[0]] for i in range(4))
                    counts[successor] += 1
            denominator = (q**4-1)*(q-1)*q**3
            self.assertEqual(sum(counts.values()), denominator)
            for state in vectors:
                if not any(state):
                    expected = F(0)
                elif state == (1, 0, 0, 0):
                    expected = F(q**3, q**4-1)
                elif not any(state[1:]):
                    expected = F(1, q**4-1)
                else:
                    expected = F(q-1, q*(q**4-1))
                self.assertEqual(F(counts[state], denominator), expected)
            for rank in range(1, 5):
                containing = sum(count for state, count in counts.items()
                                 if not any(state[rank:]))
                cap = F(1, q)+(1-F(1, q))*F(q**rank-1, q**4-1)
                self.assertEqual(F(containing, denominator), cap)
                if rank < 4:
                    disjoint = sum(count for state, count in counts.items()
                                   if not state[0] and not any(state[rank+1:]))
                    self.assertEqual(F(disjoint, denominator), cap-F(1, q))

    def test_empty_and_full_packet_endpoints(self):
        z = F(1, 2)
        for q in (256, 65536):
            tau = (1+(q-1)*z)/q
            for mixer in ('rank_one', 'uniform', 'diagonal'):
                local = screen.local_operators(log2(z), mixer, q)
                self.assertEqual(local[0][0], (0.0, screen.NEG))
                self.assertEqual(local[0][1][0], screen.NEG)
                self.assert_log_matches(local[0][1][1], z**61)
                moment = tau**64
                self.assert_log_matches(local[64][0][0], min(moment, F(1, q**4)))
                self.assert_log_matches(local[64][0][1], moment-F(1, q**64))
                self.assert_log_matches(local[64][1][1], moment)
                if mixer == 'rank_one':
                    return_cap = min(F(q**3, q**4-1)*moment, F(1, q**4))
                elif mixer == 'uniform':
                    return_cap = moment/(q**4-1)
                else:
                    return_cap = moment/(q-1)**4
                self.assert_log_matches(local[64][1][0], return_cap)

    def test_one_and_two_slot_tilts_from_occupancy_counting(self):
        # One-column words have weight 8 unless zero. For two columns,
        # a nonzero degree-one polynomial has zero or one roots among 8 nodes.
        q, z = 256, F(1, 2)
        f1 = (1+(q-1)*z**8)/q
        f2 = (1+8*(q-1)*z**7+(q-1)*(q-7)*z**8)/q**2
        same_tile = F(8*comb(8, 2), comb(64, 2))
        self.assertEqual(same_tile, F(1, 9))
        zero = {1: f1, 2: same_tile*f2+(1-same_tile)*f1**2}
        live = {1: f1*z**53,
                2: same_tile*f2*z**53+(1-same_tile)*f1**2*z**45}
        local = screen.local_operators(log2(z), q=q)
        for j in (1, 2):
            self.assert_log_matches(local[j][0][0], min(zero[j], F(1, q**j)))
            self.assert_log_matches(local[j][0][1], zero[j]-F(1, q**j))
            self.assert_log_matches(local[j][1][1], live[j])
            unweighted = F(1, q**j)*(F(1, q)+(1-F(1, q))*F(q**j-1, q**4-1))
            expected = min(F(q**3, q**4-1)*live[j], unweighted)
            self.assert_log_matches(local[j][1][0], expected)

    def test_local_entries_are_nonnegative_and_have_direct_tilt_limits(self):
        q = 256
        for mixer in ('rank_one', 'uniform', 'diagonal'):
            for log_z in (-1/1024, -.25, -1, -8):
                local = screen.local_operators(log_z, mixer, q)
                for matrix in local:
                    for row in matrix:
                        for entry in row:
                            self.assertTrue(entry == screen.NEG or isfinite(entry))
                            self.assertLessEqual(entry, 1e-11)
                            self.assertGreaterEqual(exp2(entry), 0)
                # z^64 <= E z^wt(Y) <= 1 for any 64-coordinate output.
                for matrix in local:
                    self.assertGreaterEqual(matrix[1][1], 64*log_z-1e-10)
        local = screen.local_operators(-1e-10, q=q)
        for j in (1, 2, 4, 16, 64):
            rank = min(j, 4)
            self.assertAlmostEqual(exp2(local[j][0][0]), q**(-rank), delta=1e-8)
            self.assertAlmostEqual(exp2(local[j][0][1]), 1-q**(-j), delta=1e-8)
            self.assertAlmostEqual(exp2(local[j][1][1]), 1, delta=1e-8)
            expected = q**(-rank)*(1/q+(1-1/q)*(q**rank-1)/(q**4-1))
            self.assertAlmostEqual(exp2(local[j][1][0]), expected, delta=1e-10)

    def test_graph_power_matches_exact_path_enumeration(self):
        matrix = ((F(1, 3), F(1, 5)), (F(1, 7), F(2, 9)))
        for length in (0, 1, 2, 3, 7):
            actual = screen.power(as_logs(matrix), length)
            for source in range(2):
                for target in range(2):
                    if length == 0:
                        exact = F(source == target)
                    else:
                        exact = F(0)
                        for middle in product(range(2), repeat=length-1):
                            path = (source,)+middle+(target,)
                            weight = F(1)
                            for left, right in zip(path, path[1:]):
                                weight *= matrix[left][right]
                            exact += weight
                    self.assert_log_matches(actual[source][target], exact)
        triangular = ((F(1, 2), F(1, 3)), (F(0), F(1, 4)))
        actual = screen.power(as_logs(triangular), 5)
        expected = raw_power(triangular, 5)
        for a in range(2):
            for b in range(2):
                self.assert_log_matches(actual[a][b], expected[a][b])

    def test_pair_polynomial_at_alpha_one_is_one_operator_squared(self):
        q, log_z = 256, -.125
        local = screen.local_operators(log_z, q=q)
        coefficients = screen.pair_coefficients(log_z, 1, q=q)
        for marker in (.125, .4, .75):
            single = one_at_marker(local, marker)
            squared = raw_multiply(single, single)
            paired = pair_at_marker(coefficients, marker)
            for a in range(2):
                for b in range(2):
                    self.assertAlmostEqual(paired[a][b]/squared[a][b], 1, places=11)
            arguments = dict(active_rows=17, log_z=log_z, alpha=1, marker=marker,
                             q=q, rows=128, outer_n=4, outer_k=2, limit=173)
            self.assertAlmostEqual(screen.route_bound(**arguments, steps=1),
                                   screen.route_bound(**arguments, steps=2), places=10)

    def test_pair_fractional_power_is_no_larger_than_separate_steps(self):
        q, log_z, alpha = 256, -.125, .5
        local = screen.local_operators(log_z, q=q)
        coefficients = screen.pair_coefficients(log_z, alpha, q=q)
        for marker in (.125, .4):
            single = one_at_marker(local, marker, alpha)
            separate = raw_multiply(single, single)
            paired = pair_at_marker(coefficients, marker)
            for a in range(2):
                for b in range(2):
                    self.assertGreaterEqual(paired[a][b], 0)
                    self.assertLessEqual(paired[a][b], separate[a][b]*(1+1e-12))

    def test_regional_normalizer_and_projective_factor_alpha_placement(self):
        # Synthetic squared rational entries isolate the routing calculation.
        # Their alpha=1/2 powers below are EXACT rational numbers.
        roots = tuple(((F(1, 2), F(1, 4)), (F(1, 8), F(1, 2))) if j % 2
                      else ((F(1), F(1, 8)), (F(1, 4), F(1, 4))) for j in range(65))
        local = tuple(as_logs(tuple(tuple(value**2 for value in row) for row in matrix))
                      for matrix in roots)
        rows, outer_n, outer_k, active, q, limit = 128, 3, 1, 2, 4, 3
        marker, alpha, log_z = F(1, 3), F(1, 2), -1
        weights = tuple(comb(64, j)*marker**j*(1-marker)**(64-j) for j in range(65))
        transfer = tuple(tuple(sum((weights[j]*roots[j][a][b] for j in range(65)), F(0))
                               for b in range(2)) for a in range(2))
        run = raw_power(transfer, rows*outer_n//64)
        tail = run[0][0]+run[0][1]
        binomial_mass = comb(rows, active)*marker**active*(1-marker)**(rows-active)
        beta = F(q**outer_n, (q-1)**(outer_n-outer_k))
        expected = (log2(comb(rows, active))
                    +float(alpha)*(active*log2(beta)-log2(q-1)-limit*log_z)
                    +log2(tail)-outer_n*log2(binomial_mass))
        with patch.object(screen, 'local_operators', return_value=local):
            actual = screen.route_bound(active, log_z, float(alpha), float(marker),
                                        rows=rows, outer_n=outer_n, outer_k=outer_k,
                                        limit=limit, q=q)
        self.assertAlmostEqual(actual, expected, places=10)
        # Exact fixed-count coefficient versus its Bernoulli envelope.
        conditioned = [[F(0), F(0)], [F(0), F(0)]]
        for j in range(active+1):
            weight = F(comb(64, j)*comb(64, active-j), comb(rows, active))
            block = raw_multiply(roots[j], roots[active-j])
            for a in range(2):
                for b in range(2):
                    conditioned[a][b] += weight*block[a][b]
        unconditional_region = raw_multiply(transfer, transfer)
        for a in range(2):
            for b in range(2):
                self.assertLessEqual(conditioned[a][b], unconditional_region[a][b]/binomial_mass)
        conditional_run = raw_power(conditioned, outer_n)
        conditional_tail = conditional_run[0][0]+conditional_run[0][1]
        self.assertLessEqual(conditional_tail, tail/binomial_mass**outer_n)

    def test_documented_argument_domains(self):
        for log_z in (0, 1):
            with self.assertRaises(ValueError):
                screen.local_operators(log_z)
        with self.assertRaises(ValueError):
            screen.local_operators(-1, mixer='not-a-mixer')
        for change in ({'active_rows': 0}, {'active_rows': 8192}, {'alpha': 0},
                       {'alpha': 2}, {'marker': 0}, {'marker': 1}, {'steps': 3}):
            arguments = dict(active_rows=1, log_z=-1, alpha=.5, marker=.125)
            arguments.update(change)
            with self.assertRaises(ValueError):
                screen.route_bound(**arguments)


if __name__ == '__main__':
    unittest.main()
