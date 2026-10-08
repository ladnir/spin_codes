"""Exact small-field checks and directed whole-class regression tests."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import unittest

from hd_randomized_outer_proof import (
    Parameters, dense_class_log2, dense_tail_receipt,
    independent_nonzero_syndrome_zero, outer_pointwise_counts,
    prescribed_mds_count, quiet_rank_marker_log2, quiet_sparse_receipt,
    quiet_weight_polynomial,
)


def mul4(a, b):
    result = 0
    for _ in range(2):
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


class RandomizedOuterProofTests(unittest.TestCase):
    def setUp(self):
        from flint import ctx
        self.previous_precision = ctx.prec
        ctx.prec = 160

    def tearDown(self):
        from flint import ctx
        ctx.prec = self.previous_precision

    def test_exact_randomized_mds_atoms(self):
        counts = Counter()
        for constant, linear in product(range(4), repeat=2):
            if constant == linear == 0:
                continue
            word = tuple(constant ^ mul4(linear, x) for x in range(3))
            for scalars in product(range(1, 4), repeat=3):
                counts[tuple(mul4(x, a) for x, a in zip(word, scalars))] += 1
        atoms = outer_pointwise_counts(3, 2, 4)
        for word in product(range(4), repeat=3):
            weight = sum(x != 0 for x in word)
            self.assertEqual(F(counts[word], 27), atoms[weight])
            self.assertLessEqual(atoms[weight], F(1, 3))
        self.assertEqual(sum(counts.values()), 27*15)

    def test_real_outer_atom_envelope_and_normalization(self):
        q, n, k = 65536, 256, 128
        atoms = outer_pointwise_counts(n, k, q)
        for atom in atoms:
            self.assertLessEqual(atom, F(1, (q-1)**128))
        total = sum(comb(n, w)*prescribed_mds_count(w, n, k, q)
                    for w in range(n+1))
        self.assertEqual(total, q**k)
        self.assertTrue(all(atoms[w] == 0 for w in range(129)))

    def test_nonzero_amplitudes_and_syndrome_law(self):
        for weight in range(5):
            hits = 0
            for values in product(range(1, 4), repeat=weight):
                first = second = 0
                for node, value in enumerate(values):
                    first ^= value
                    second ^= mul4(node, value)
                hits += first == second == 0
            self.assertEqual(F(hits, 3**weight),
                             independent_nonzero_syndrome_zero(weight, 2, 4))
        polynomial = quiet_weight_polynomial(4, 2, 4)
        self.assertEqual(polynomial, (1, 0, 0, 12, 3))
        self.assertEqual(sum(polynomial), 4**2)

    def test_uniform_potential_syndrome_rank(self):
        for size in range(5):
            for support in combinations(range(4), size):
                hits = 0
                for values in product(range(4), repeat=size):
                    first = second = 0
                    for node, value in zip(support, values):
                        first ^= value
                        second ^= mul4(node, value)
                    hits += first == second == 0
                self.assertEqual(F(hits, 4**size), F(1, 4**min(size, 2)))

    def test_regional_fractional_marker_including_final_packet(self):
        from flint import arb
        p = Parameters(field=4, rows=4, outer_n=3, outer_k=2,
                       width=2, state=1, limit=1)
        alpha, marker = F(1, 2), F(1, 3)
        for rows in range(1, 4):
            subsets = tuple(combinations(range(4), rows))
            whole = last = F(0)
            for support in subsets:
                occupied = set(i//2 for i in support)
                whole += F(1, 2**len(occupied))
                last += F(1, 2**(0 in occupied))
            whole /= len(subsets)
            last /= len(subsets)
            # q^(-alpha rank) is 2^-rank; only the final packet is free.
            moment = whole**2*last
            exact = (arb(comb(4, rows))
                     *(arb(p.beta.numerator)/p.beta.denominator)**(arb(rows)/2)
                     /arb(3).sqrt()*arb(moment.numerator)/moment.denominator)
            upper = arb(2)**quiet_rank_marker_log2(rows, marker, alpha, p)
            self.assertTrue(exact < upper)

    def test_dense_ranges_and_no_row_subset_union(self):
        from flint import arb
        all_words = dense_tail_receipt()
        quiet = dense_tail_receipt(quiet=True)
        self.assertEqual(all_words['active_rows_first'], 7548)
        self.assertEqual(quiet['active_rows_first'], 6566)
        self.assertTrue(arb(all_words['union_log2_upper']) < -326)
        self.assertTrue(arb(quiet['union_log2_upper']) < -1000)
        difference = dense_class_log2(7548)-dense_class_log2(7549)
        expected = 128*arb(65535).log()/arb(2).log()
        self.assertTrue((difference-expected).contains(0))

    def test_sparse_quiet_directed_union(self):
        from flint import arb
        receipt = quiet_sparse_receipt()
        self.assertEqual(receipt['active_rows_last'], 953)
        self.assertTrue(arb(receipt['union_log2_upper']) < -142)
        self.assertEqual(receipt['weakest_marker'], '1/8')
        self.assertEqual(receipt['weakest_alpha'], '3/32')


if __name__ == '__main__':
    unittest.main()
