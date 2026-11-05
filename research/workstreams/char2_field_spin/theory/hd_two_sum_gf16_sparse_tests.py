"""Independent checks for the GF16-coset sparse proof proposal."""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb, log2
import random
import unittest
from unittest.mock import patch

import numpy as np

from gf8_packet_maps import MUL, matvec, power
from hd_two_sum_affine_outer import column_action_exact, column_action_float
from hd_two_sum_gf16_sparse import (LABELS, matrix, direct, circuit, sparse_law,
                                    exact_local_checks, ordered_regions, screen)


class GF16SparseTests(unittest.TestCase):
    def test_coset_and_every_projective_pair(self):
        self.assertEqual(tuple(x ^ 6 for x in LABELS), tuple(x for x in range(256) if power(x, 16) == x))
        self.assertEqual(exact_local_checks()['projective_pairs_checked'], 30600)
        self.assertEqual(sum(sparse_law(2).values()), 1)

    def test_all_binary_basis_and_random_circuits(self):
        rng = random.Random(1616820)
        samples = [[int(j == i)*(1 << bit) for j in range(16)] for i in range(16) for bit in range(8)]
        samples += [[rng.randrange(256) for _ in range(16)] for _ in range(100)]
        for values in samples:
            self.assertEqual(circuit(values), direct(values))
            self.assertEqual(circuit(values), matvec(matrix(), values))
            self.assertEqual(circuit(circuit(values)), tuple(values))

    def test_silent_triple_vectors_remain(self):
        for i, j, k in combinations(range(16), 3):
            values = [0]*16
            values[i], values[j], values[k] = LABELS[j] ^ LABELS[k], LABELS[i] ^ LABELS[k], LABELS[i] ^ LABELS[j]
            self.assertEqual(direct(values), tuple(values))
            self.assertEqual(sum(bool(x) for x in values), 3)

    def test_sixteen_region_hole_mixture_with_noncommuting_matrices(self):
        regions = [[(F(2+i+m), F(1+i), F(1+m), F(3+i)) for m in range(4)] for i in range(16)]
        potential = (F(1), F(1, 7))
        expected = [F(0), F(0)]
        for holes in product(range(-1, 16), repeat=2):
            probability = F(239**holes.count(-1), 255**2)
            vector = potential
            for i in reversed(range(16)):
                a, b, c, d = regions[i][3-holes.count(i)]
                vector = (a*vector[0]+b*vector[1], c*vector[0]+d*vector[1])
            for j in range(2):
                expected[j] += probability*vector[j]
        actual = column_action_exact(regions, potential, 1, 2)
        self.assertEqual(actual, tuple(expected))
        logs = np.log2(np.array(regions, dtype=float).reshape(16, 4, 2, 2))
        floated = column_action_float(logs, log2(1/7), 1, 2)
        for j in range(2):
            self.assertAlmostEqual(floated[j], log2(float(expected[j])), places=11)

    def test_region_order_and_subset_envelopes(self):
        raw = {(m, f): np.full((2, 2), 100*m+f) for m in range(18) for f in range(m//2+1)}
        tables = ordered_regions(raw, 17, 8, 0)
        self.assertEqual(tables.shape, (16, 18, 2, 2))
        self.assertTrue(np.all(tables[0::2, 8] == 804))
        self.assertTrue(np.all(tables[1::2, 8] == 800))

    def test_screen_uses_new_rs_dimensions_and_counts_once(self):
        with patch('hd_two_sum_gf16_sparse.pair_regions_screen', return_value=None) as raw, \
                patch('hd_two_sum_gf16_sparse.ordered_regions', return_value=None), \
                patch('hd_two_sum_gf16_sparse.column_action_float', return_value=np.array([0., 0.])), \
                patch('hd_two_sum_gf16_sparse.sparse_log_count', return_value=0.) as poly:
            result = screen()
        self.assertEqual(raw.call_count, 1)
        self.assertEqual(poly.call_count, 72)
        self.assertTrue(all(call.kwargs == dict(n=16, k=8) for call in poly.call_args_list))
        self.assertFalse(result['complete_source_coverage'])
        self.assertFalse(result['whole_code_certificate'])
        for source in result['sources']:
            for n1, row in enumerate(source['entries']):
                n2 = 17-n1
                self.assertEqual(row['source_colorings'], comb(17, n1))
                self.assertEqual(row['pair_row_choices'], 120**n2)
                self.assertEqual(row['singleton_row_choices'], 16**n1)
                self.assertEqual(row['sparse_support_floor'], 9*(n1+2*n2))
                self.assertEqual(row['sparse_normalizer_exponent'], 8*(n1+2*n2))


if __name__ == '__main__':
    unittest.main()
