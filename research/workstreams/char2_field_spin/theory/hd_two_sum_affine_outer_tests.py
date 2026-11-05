"""Exact affine circuit, sparse laws and ordered positive polynomial checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
import random
import unittest
from unittest.mock import patch

import numpy as np

from gf8_packet_maps import MUL, matvec, matmul, xor_sum
from hd_two_sum_affine_outer import (LABELS, matrix, direct, circuit,
    uniform_mask_counts, nonzero_mask_counts, sparse_mask_law, silent_moment_count,
    column_action_exact, column_action_float, sparse_column_coefficients, screen, CONSTRUCTION)


def multiply_vector(matrix, vector):
    a, b, c, d = matrix
    x, y = vector
    return a*x+b*y, c*x+d*y


class AffineTwoSumTests(unittest.TestCase):
    def test_involution_nonzero_entries_and_all_two_column_minors(self):
        h = matrix()
        self.assertEqual(xor_sum(LABELS), 0)
        self.assertEqual(xor_sum(MUL[u][u] for u in LABELS), 0)
        self.assertEqual(matmul(h, h), tuple(tuple(int(i == j) for j in range(8)) for i in range(8)))
        self.assertTrue(all(value for row in h for value in row))
        for i, j in combinations(range(8), 2):
            for k, l in combinations(range(8), 2):
                self.assertNotEqual(MUL[h[k][i]][h[l][j]] ^ MUL[h[k][j]][h[l][i]], 0)

    def test_circuit_every_binary_basis_and_fixed_random_inputs(self):
        examples = [[(1 << bit) if i == coordinate else 0 for i in range(8)]
                    for coordinate in range(8) for bit in range(8)]
        rng = random.Random(991132)
        examples += [[rng.randrange(256) for _ in range(8)] for _ in range(100)]
        for values in examples:
            self.assertEqual(circuit(values), matvec(matrix(), values))
            self.assertEqual(direct(values), circuit(values))
            self.assertEqual(circuit(circuit(values)), tuple(values))

    def test_every_sparse_projective_law_matches_rank_kernel(self):
        for arity in (1, 2):
            expected = {mask: probability*255**arity for mask, probability in sparse_mask_law(arity).items()}
            for positions in combinations(range(8), arity):
                counts = Counter()
                for ratio in range(1, 256) if arity == 2 else (1,):
                    x = [0]*8
                    x[positions[0]] = 1
                    if arity == 2:
                        x[positions[1]] = ratio
                    target = sum(bool(value) << i for i, value in enumerate(direct(x)))
                    counts[target] += 255
                active = sum(1 << i for i in positions)
                self.assertEqual(dict(counts), expected)
                self.assertEqual(nonzero_mask_counts(active), tuple(counts[t] for t in range(256)))

    def test_silent_triples_and_uniform_boundaries(self):
        for i, j, k in combinations(range(8), 3):
            x = [0]*8
            x[i], x[j], x[k] = LABELS[j] ^ LABELS[k], LABELS[i] ^ LABELS[k], LABELS[i] ^ LABELS[j]
            self.assertEqual(direct(x), tuple(x))
            self.assertEqual(sum(bool(value) for value in x), 3)
        self.assertEqual(silent_moment_count(3), 255)
        self.assertEqual(uniform_mask_counts(0), (1,)+(0,)*255)
        self.assertEqual(uniform_mask_counts(255), tuple(255**t.bit_count() for t in range(256)))
        self.assertEqual(sum(sparse_mask_law(2).values()), 1)

    def test_ordered_hole_recurrence_against_every_labelled_assignment(self):
        single, double = 1, 2
        t = single+double
        regions = [[(F(1+i+m), F(2+i), F(1+m*m), F(3+2*i+m)) for m in range(t+1)]
                   for i in range(8)]
        potential = (F(1), F(1, 3))
        expected = [F(0), F(0)]
        for holes in product(range(-1, 8), repeat=double):
            probability = F(247**holes.count(-1), 255**double)
            occupancies = [t-holes.count(i) for i in range(8)]
            vector = potential
            for i in reversed(range(8)):
                vector = multiply_vector(regions[i][occupancies[i]], vector)
            for state in range(2):
                expected[state] += probability*vector[state]
        actual = column_action_exact(regions, potential, single, double)
        self.assertEqual(actual, tuple(expected))
        logs = np.log2(np.array(regions, dtype=float).reshape(8, t+1, 2, 2))
        screen_value = column_action_float(logs, log2(float(potential[1])), single, double)
        for state in range(2):
            self.assertAlmostEqual(screen_value[state], log2(float(expected[state])), places=12)
        reverse = column_action_exact(list(reversed(regions)), potential, single, double)
        self.assertNotEqual(reverse, actual)

    def test_constant_regional_matrices_have_unit_hole_mass(self):
        potential = (F(1), F(2, 5))
        regions = [[(F(1+i), F(2), F(3), F(4+i))]*5 for i in range(8)]
        expected = potential
        for table in reversed(regions):
            expected = multiply_vector(table[0], expected)
        for double in range(4):
            self.assertEqual(column_action_exact(regions, potential, 1, double), expected)

    def test_sparse_polynomial_matches_named_constituent_subsets(self):
        n1, n2, q = 1, 2, 8
        caps = {(s, d): F(1+s+3*d, 2) for s in range(4) for d in range(3) if s+d <= 3}
        expected = [F(0)]*(n1+2*n2+1)
        for first, pair1, pair2 in product(range(2), range(4), range(4)):
            sizes = [pair1.bit_count(), pair2.bit_count()]
            single, double = first+sizes.count(1), sizes.count(2)
            degree = first+sum(sizes)
            expected[degree] += (q-1)**degree*caps[single, double]
        self.assertEqual(sparse_column_coefficients(caps, n1, n2, q), tuple(expected))
        ones = {key: F(1) for key in caps}
        self.assertEqual(sparse_column_coefficients(ones, n1, n2, q),
                         tuple(F(comb(5, degree)*7**degree) for degree in range(6)))
        self.assertEqual(sum(comb(17, n1)*8**n1*28**(17-n1) for n1 in range(18)), 36**17)

    def test_bounded_screen_reuses_two_banks_and_has_no_geometry_strata(self):
        with patch('hd_xor_pair_route.pair_regions_screen', return_value=None) as regions, \
                patch('hd_xor_mixed_gate.ordered_region_tables', return_value=None), \
                patch('hd_two_sum_affine_outer.column_action_float', return_value=np.array([0., 0.])) as action, \
                patch('hd_xor_mixed_gate.sparse_log_count', return_value=0.):
            result = screen()
        self.assertEqual(regions.call_count, 2)
        self.assertEqual(action.call_count, 2*18*19//2)
        self.assertEqual(result['construction'], CONSTRUCTION)
        self.assertEqual(len(result['selected']), 18)
        self.assertFalse(result['pair_multiplicity_strata_used'])
        self.assertFalse(result['old_certificates_transferred'])
        self.assertFalse(result['whole_code_certificate'])
        for n1, row in enumerate(result['selected']):
            self.assertEqual(row['source_colorings'], comb(17, n1))
            self.assertEqual(row['pair_row_choices'], 28**(17-n1))
            self.assertEqual(row['sparse_support_floor'], 17*(34-n1))
            candidates = [bank['entries'][n1]['failure_log2_screen'] for bank in result['banks']]
            self.assertEqual(row['failure_log2_screen'], min(candidates))


if __name__ == '__main__':
    unittest.main()
