"""Exact local regressions and continuation-cone tests; no certificate."""
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
import math
from pathlib import Path
import sys
import unittest

import numpy as np

from cauchy_firstmoment_screen import geometry, log_terminal_power
from gf8_rank_history_screen import (bound, cone_preserved, local_matrix,
                                      modal_log_penalty, packet, return_bound, source_row_log)
from higher_distance_screen import output_threshold
from state_rank_cap_transfer import local_matrix as original_rank_matrix
from state_rank_cap_transfer import rank_cap_weights, rank_tail_envelope
from superregular_return_screen import local_matrix as original_two_state
from symmetric_density_screen import bound as control_bound
from symmetric_density_screen import packet as control_packet
from symmetric_density_screen import return_bound as original_return


def rays(state):
    return [(F(1),)+(F(0),)*state]+[
        tuple(F(int(index >= start)) for index in range(state+1))
        for start in range(1, state+1)]


class RankHistoryTests(unittest.TestCase):
    def test_state8_pure_exact_backward_agreement(self):
        for weight in range(17):
            self.assertEqual(local_matrix(weight, F(3, 4), F(11, 2), 8, joint=False),
                             original_rank_matrix(weight, F(3, 4), F(11, 2)))

    def test_full_cap_recovers_two_state_local_masses(self):
        y, tau = F(9, 10), F(7)
        for state in (4, 6, 8):
            for weight in range(17):
                matrix = local_matrix(weight, y, tau, state)
                old = original_two_state(weight, y, tau, state)
                self.assertEqual(matrix[0][0], old[0][0])
                self.assertEqual(sum(matrix[0][1:]), old[0][1])
                self.assertEqual(sum(matrix[state][1:]), old[1][1])
                self.assertEqual(matrix[state][0], original_return(weight, y, tau, state)
                                 if weight else 0)

    def test_exact_packet_cone_all_generators(self):
        for state in (4, 6, 8):
            for y, tau in ((F(1), F(1)), (F(9, 10), F(7, 2)), (F(1, 3), F(255))):
                for weight in range(17):
                    matrix = np.asarray(local_matrix(weight, y, tau, state), dtype=object)
                    self.assertTrue(cone_preserved(matrix))
                    for ray in rays(state):
                        values = matrix @ np.asarray(ray, dtype=object)
                        self.assertTrue(all(value >= 0 for value in values))
                        self.assertTrue(all(values[cap] <= values[cap+1]
                                            for cap in range(1, state)))

    def test_joint_caps_tighten_uniform_return_without_losing_monotonicity(self):
        for state in (4, 6, 8):
            for weight in range(1, 17):
                values = [return_bound(weight, F(3, 4), F(7), cap, state)
                          for cap in range(1, state+1)]
                for cap, value in enumerate(values, 1):
                    self.assertLessEqual(value, return_bound(weight, F(3, 4), F(7), cap,
                                                            state, cap_joint=False))
                self.assertEqual(values[-1], original_return(weight, F(3, 4), F(7), state))
                self.assertTrue(all(a <= b for a, b in zip(values, values[1:])))

    def test_cap_tails_and_successor_bound(self):
        for state in (4, 6, 8):
            for cap in range(1, state+1):
                self.assertEqual(rank_tail_envelope(cap, 0, F(7), 256, state),
                                 (F(1),)*cap+(F(0),)*(state-cap))
                for roots in range(state+1):
                    masses = rank_cap_weights(cap, roots, F(7), 256, state)
                    self.assertTrue(all(value >= 0 for value in masses))
                    self.assertEqual(sum(masses[min(cap, state-roots):]), 0)
                for weight in (0, 1, 2):
                    matrix = local_matrix(weight, F(3, 4), F(7), state)
                    self.assertEqual(sum(matrix[cap][min(state, cap+weight)+1:]), 0)
            birth = local_matrix(1, F(3, 4), F(7), state)
            self.assertEqual(birth[0][1], F(3, 4))
            self.assertEqual(sum(birth[0][2:]), 0)

    def test_batched_aggregate_against_exact_local_sum(self):
        x, y, tau = F(4, 5), F(2, 3), F(17, 2)
        for state in (4, 6, 8):
            expected = sum((np.asarray(local_matrix(weight, y, tau, state), dtype=object)
                            *math.comb(16, weight)*x**weight for weight in range(17)),
                           np.zeros((state+1, state+1), dtype=object))
            actual = np.exp(packet(math.log(x), math.log(y), math.log(tau), state))
            # Tail differences can lose several ulps in the floating LP;
            # the exact local reference also checks this numerical error.
            np.testing.assert_allclose(actual, np.asarray(expected, dtype=float),
                                       rtol=2e-12, atol=1e-300)
            self.assertTrue(cone_preserved(expected))

    def test_provided_return_minima_preserve_cone_and_strengthen_cover(self):
        y, tau, state = F(3, 4), F(7), 6
        for weight in (1, 2, 3):
            original = local_matrix(weight, y, tau, state)
            refined = local_matrix(weight, y, tau, state, return_minima=(15, 13))
            self.assertTrue(cone_preserved(refined))
            for cap in range(state+1):
                self.assertLessEqual(refined[cap][0], original[cap][0])
                self.assertEqual(refined[cap][1:], original[cap][1:])

    def test_positive_products_and_sums_preserve_cone(self):
        for state in (4, 6, 8):
            matrices = [np.asarray(local_matrix(weight, F(2, 3), F(7), state), dtype=object)
                        for weight in (0, 1, 2, 9)]
            for left, right in zip(matrices, reversed(matrices)):
                self.assertTrue(cone_preserved(left @ right))
                self.assertTrue(cone_preserved(F(2, 5)*left+F(3, 5)*right))
            # f(0) is independent: reject a zero-column monotonicity failure.
            broken = matrices[0].copy()
            broken[1, 0] = broken[2, 0]+1
            self.assertFalse(cone_preserved(broken))

    def test_density_source_route_and_control_replacement(self):
        point = (-.23, -.72, math.log(12.3), -.09)
        g = geometry(16, 128, 86)
        weight = round(.1*g['n'])
        value, terms = bound(point, weight, 16, 128, 86, 6, components=True)
        self.assertAlmostEqual(value, sum(terms.values()), places=8)
        self.assertAlmostEqual(terms['route'], 128*math.log2(g['rows']+1), places=10)
        self.assertAlmostEqual(terms['selection'], -weight*(point[0]+point[3])/math.log(2), places=10)
        self.assertAlmostEqual(terms['output'], -output_threshold(g['n'], F(1, 10))*point[1]/math.log(2), places=10)
        self.assertAlmostEqual(terms['rank_threshold'], -point[2]/math.log(2), places=10)
        alpha = weight/g['n']
        row = (1-alpha)**128+math.exp(-42*point[2])*(1-alpha+alpha*math.exp(point[2]+point[3]))**128
        self.assertAlmostEqual(terms['source'], g['rows']*math.log2(row), places=8)
        old = control_bound(point, weight, 16, 128, 86, 6, joint=True)
        def power(logs):
            shift = float(logs.max())
            return g['n']//16*shift+log_terminal_power(np.exp(logs-shift), g['n']//16)
        difference = (power(packet(*point[:3], 6))-power(control_packet(*point[:3], 6, True)))/math.log(2)
        self.assertAlmostEqual(value-old, difference, places=8)
        self.assertLessEqual(value, old)

    def test_feasible_source_log_against_exact_positive_sum(self):
        for length, dimension in ((3, 2), (4, 2), (128, 86)):
            parity = length-dimension
            for alpha in (F(1, 10), F(1, 2), F(9, 10)):
                for tau, tilt in ((F(3), F(2, 3)), (F(7), F(3, 2))):
                    expected = (1-alpha)**length+sum(
                        math.comb(length, h)*alpha**h*(1-alpha)**(length-h)*tau**(h-parity)*tilt**h
                        for h in range(parity+1, length+1))
                    actual = source_row_log(float(alpha), math.log(tau), math.log(tilt),
                                            length, dimension, True)
                    self.assertAlmostEqual(actual, math.log(expected), places=12)
                    self.assertLessEqual(actual, source_row_log(float(alpha), math.log(tau),
                                                               math.log(tilt), length, dimension)+1e-12)

    def test_feasible_source_homogenization_and_beta_guard(self):
        for length, parity in ((3, 1), (4, 2)):
            for beta in (F(2), F(3), F(7)):
                for probabilities in product((F(0), F(1, 2), F(1)), repeat=length):
                    tau = F(3)
                    weights = [F(1)]+[F(0)]*parity+[
                        tau**(-parity)*beta**h for h in range(parity+1, length+1)]
                    heterogeneous = sum(
                        weights[sum(mask)]*math.prod(p if bit else 1-p
                                                    for p, bit in zip(probabilities, mask))
                        for mask in product((0, 1), repeat=length))
                    alpha = sum(probabilities)/length
                    homogeneous = sum(math.comb(length, h)*alpha**h*(1-alpha)**(length-h)*weights[h]
                                      for h in range(length+1))
                    self.assertLessEqual(heterogeneous, homogeneous)
        # beta=1, tau=3, p=1: densities (1,1,0) give1/3, but averaging gives23/81.
        self.assertLess(F(23, 81), F(1, 3))
        args = (2/3, math.log(3), -math.log(3), 3, 2)
        self.assertEqual(source_row_log(*args, True), source_row_log(*args, False))

    def test_modal_penalty_and_exact_source_profile_comparison(self):
        for rows in range(2, 18, 2):
            reciprocal = F(2**rows, math.comb(rows, rows//2))
            self.assertAlmostEqual(modal_log_penalty(rows), math.log(reciprocal), places=14)
            for hits in range(rows+1):
                probability = math.comb(rows, hits)*F(hits, rows)**hits*F(rows-hits, rows)**(rows-hits)
                self.assertGreaterEqual(probability, 1/reciprocal)
        rows, length, parity, tau = 2, 3, 1, F(3)
        subsets = [tuple(combinations(range(rows), h)) for h in range(rows+1)]
        reciprocal = F(2**rows, math.comb(rows, rows//2))
        for profile in product(range(rows+1), repeat=length):
            actual = sum(math.prod(1 if hits == 0 else tau**(hits-parity) if hits > parity else 0
                                   for hits in (sum(row in column for column in columns) for row in range(rows)))
                         for columns in product(*(subsets[h] for h in profile)))
            actual /= math.prod(math.comb(rows, h) for h in profile)
            alpha = F(sum(profile), rows*length)
            for tilt in (F(2, 3), F(4, 3)):
                bracket = (1-alpha)**length+sum(math.comb(length, h)*alpha**h*(1-alpha)**(length-h)
                                               *tau**(h-parity)*tilt**h for h in range(parity+1, length+1))
                upper = reciprocal**length*tilt**(-sum(profile))*bracket**rows
                self.assertLessEqual(actual, upper)

    def test_dimension80_outer_count_against_exact_lagrange(self):
        # This checks a costed scalar circuit only; no native header is generated.
        generators = Path(__file__).resolve().parents[4]/'spin/experiments/char2_field_spin'
        sys.path.insert(0, str(generators))
        try:
            from generate_cauchy256_outer import interpolation
            from generate_cauchy256_outer_fft import counts, evaluate, fft_ops, translate
            def parity_transpose(values):
                pieces = [evaluate(fft_ops(16, offset, transpose=True), values[first:first+16])
                          for first, offset in zip((0, 16, 32), (16, 32, 48))]
                combined = [pieces[0][i]^pieces[1][i]^pieces[2][i] for i in range(16)]
                tail = evaluate(fft_ops(16, 0, True, True), combined)
                return translate(tail+list(values))+tail
            expected = interpolation(80, range(80, 128))
            for column in range(48):
                self.assertEqual(parity_transpose([int(i == column) for i in range(48)]), expected[column])
            operations = (fft_ops(64, 0, True)+fft_ops(64, 64)+fft_ops(16, 0, True)
                          +sum((fft_ops(16, offset) for offset in (16, 32, 48)), ()))
            self.assertEqual(dict(counts(operations)+Counter(xor=32)), {'gfni': 434, 'xor': 978})
        finally:
            sys.path.pop(0)


if __name__ == '__main__':
    unittest.main()
