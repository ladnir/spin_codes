"""Finite algebra only: no distance provider or native benchmark."""
from collections import Counter
from itertools import combinations, product
from math import comb
import random
import unittest

from gf8_packet_maps import determinant, matmul, matvec, rank
from hd_bit_syndrome_outer import (
    affine_matrix, affine_rank, circuit, circuit_counts, coefficient_matrix,
    direct, exact_summary, matrix, pair_mask_counts, quiet_support_counts,
    rank_histogram,
)


class BitSyndromeOuterTests(unittest.TestCase):
    def test_factorization_gram_and_involution(self):
        v, b, h = affine_matrix(), coefficient_matrix(), matrix()
        self.assertEqual(matmul(tuple(zip(*v)), v), ((0,) * 5,) * 5)
        self.assertNotEqual(determinant(b), 0)
        self.assertEqual(rank(v), 5)
        self.assertEqual(matmul(matmul(v, b), tuple(zip(*v))),
                         tuple(tuple(h[i][j] ^ int(i == j) for j in range(16)) for i in range(16)))
        self.assertEqual(matmul(h, h), tuple(tuple(int(i == j) for j in range(16)) for i in range(16)))
        self.assertTrue(all(x for row in h for x in row))

    def test_all_binary_basis_and_random_circuit(self):
        rng = random.Random(510166)
        samples = [[int(i == j) * (1 << bit) for j in range(16)]
                   for i in range(16) for bit in range(8)]
        samples += [[rng.randrange(256) for _ in range(16)] for _ in range(40)]
        for x in samples:
            self.assertEqual(circuit(x), direct(x))
            self.assertEqual(circuit(circuit(x)), tuple(x))
        self.assertEqual(circuit_counts(), {'binary_xors': 66, 'fixed_products': 8})

    def test_all_projective_pairs_and_labelled_regression(self):
        classes = Counter()
        for i, j in combinations(range(16), 2):
            law = pair_mask_counts(i, j)
            self.assertEqual(sum(count for _, count in law), 255)
            histogram = Counter()
            for mask, count in law:
                histogram[16 - mask.bit_count()] += count
            classes[tuple(sorted(histogram.items()))] += 1
        self.assertEqual(classes, {
            ((0, 239), (1, 16)): 12,
            ((0, 240), (1, 14), (2, 1)): 13,
            ((0, 241), (1, 12), (2, 2)): 2,
            ((0, 242), (1, 12), (4, 1)): 68,
            ((0, 243), (1, 10), (2, 1), (4, 1)): 12,
            ((0, 244), (1, 10), (6, 1)): 8,
            ((0, 245), (1, 6), (2, 3), (4, 1)): 1,
            ((0, 245), (1, 8), (2, 1), (6, 1)): 4,
        })
        self.assertTrue(any(mask.bit_count() == 10 for mask, _ in pair_mask_counts(3, 5)))

    def test_affine_rank_and_no_quiet_triples(self):
        v = affine_matrix()
        for arity in range(4):
            for source in combinations(range(16), arity):
                mask = sum(1 << i for i in source)
                self.assertEqual(affine_rank(mask), arity)
                self.assertEqual(rank(tuple(v[i] for i in source)), arity)
        self.assertEqual(rank_histogram(8), {4: 30, 5: 12840})
        self.assertEqual(rank_histogram(4), {3: 140, 4: 1680})
        # A parallelogram remains quiet, unlike the former quiet triples.
        x = [0] * 16
        for i in (0, 1, 2, 3):
            x[i] = 73
        self.assertEqual(direct(x), tuple(x))

    def test_small_forced_nonzero_counts_by_independent_enumeration(self):
        for q, bits in ((2, 3), (4, 2)):
            values = quiet_support_counts(q, bits)
            n = 1 << bits
            expected = [0] * (1 << n)
            for x in product(range(q), repeat=n):
                syndromes = [0] * (bits + 1)
                for i, value in enumerate(x):
                    syndromes[0] ^= value
                    for j in range(bits):
                        if i >> j & 1:
                            syndromes[j + 1] ^= value
                if not any(syndromes):
                    expected[sum(bool(value) << i for i, value in enumerate(x))] += 1
            self.assertEqual(values, tuple(expected))

    def test_arity8_no_extra_branch_confined_to_source_set(self):
        # Non-hyperplanes and their complements have rank5, so their cross
        # restriction has rank5 by the V B V^T factorization. Check the only
        # remaining30 cases directly, independently of binary elimination.
        h = matrix()
        checked = 0
        for source in combinations(range(16), 8):
            mask = sum(1 << i for i in source)
            if affine_rank(mask) != 4:
                continue
            outside = tuple(i for i in range(16) if i not in source)
            self.assertEqual(rank(tuple(tuple(h[i][j] for j in source) for i in outside)), 4)
            checked += 1
        self.assertEqual(checked, 30)

    def test_production_mobius_counts_and_scope(self):
        quiet = quiet_support_counts()
        self.assertGreaterEqual(min(quiet), 0)
        self.assertEqual(sum(quiet), 256 ** 11)
        self.assertEqual(quiet[0], 1)
        for mask, count in enumerate(quiet):
            if 0 < mask.bit_count() <= 3:
                self.assertEqual(count, 0)
            if mask.bit_count() == 4:
                self.assertEqual(count, 255 if affine_rank(mask) == 3 else 0)
        result = exact_summary()
        self.assertEqual(result['projective_pairs_checked'], 30600)
        self.assertEqual(result['maximum_pair_holes'], 6)
        self.assertEqual(result['quiet_minimum_nonempty_weight'], 4)
        self.assertEqual(sum(result['arity8_quiet_count_histogram'].values()), comb(16, 8))
        self.assertFalse(result['whole_code_certificate'])
        self.assertFalse(result['native_implementation'])

    def test_guards(self):
        for args in ((-1,), (1 << 16,), (0, 0), (0, 5)):
            with self.assertRaises(ValueError):
                affine_rank(*args)
        for q in (0, 1, 3, 255):
            with self.assertRaises(ValueError):
                quiet_support_counts(q)
        with self.assertRaises(ValueError):
            pair_mask_counts(2, 2)
        with self.assertRaises(ValueError):
            circuit([0] * 15)


if __name__ == '__main__':
    unittest.main()
