"""Tiny exact checks for diagonal-XOR-diagonal state refresh; no distance gate."""
from collections import Counter, defaultdict
from itertools import product
from math import comb
import unittest

from hd_xor_column_probe import column_support_counts, column_weight_counts, xor_column


def multiply4(a, b):
    result = 0
    while b:
        if b & 1:
            result ^= a
        b >>= 1
        a <<= 1
        if a & 4:
            a ^= 7
    return result


def support(values):
    return sum(1 << i for i, value in enumerate(values) if value)


class XorStateRefreshTests(unittest.TestCase):
    def test_support_law_against_every_gf4_width4_vector(self):
        actual = defaultdict(Counter)
        for values in product(range(4), repeat=4):
            actual[support(values)][support(xor_column(values))] += 1
            self.assertEqual(xor_column(xor_column(values)), values)
        for mask in range(16):
            expected = column_support_counts(4, 4, mask)
            self.assertEqual([actual[mask][target] for target in range(16)], expected)

    def test_output_diagonal_gives_exact_conditional_product_uniformity(self):
        actual = defaultdict(Counter)
        for values in product(range(4), repeat=4):
            incoming = support(values)
            mixed = xor_column(values)
            for diagonal in product(range(1, 4), repeat=4):
                outgoing = tuple(multiply4(a, b) for a, b in zip(diagonal, mixed))
                actual[incoming][outgoing] += 1
        for incoming in range(16):
            counts = column_support_counts(4, 4, incoming)
            for outgoing in product(range(4), repeat=4):
                mask = support(outgoing)
                expected = counts[mask]*3**(4-mask.bit_count())
                self.assertEqual(actual[incoming][outgoing], expected)

    def test_refresh_lumps_exactly_by_labelled_block_counts(self):
        # Two blocks of two state coordinates, independently of an emission map.
        seen = {}
        for incoming in range(16):
            block_counts = tuple(((incoming >> (2*i)) & 3).bit_count() for i in range(2))
            actual = Counter()
            for outgoing, count in enumerate(column_support_counts(4, 4, incoming)):
                target = tuple(((outgoing >> (2*i)) & 3).bit_count() for i in range(2))
                actual[target] += count
            k = sum(block_counts)
            expected = Counter()
            if not k:
                expected[0, 0] = 1
            else:
                expected[block_counts] += (3**k+3*(-1)**k)//4
                for j0 in range(block_counts[0]+1):
                    for j1 in range(block_counts[1]+1):
                        j, remaining = j0+j1, k-j0-j1
                        residual = j % 2 if remaining == 0 else (2**remaining+(-2)**remaining)//4
                        expected[2-j0, 2-j1] += 3*comb(block_counts[0], j0)*comb(block_counts[1], j1)*residual
            self.assertEqual(+actual, +expected)
            if block_counts in seen:
                self.assertEqual(actual, seen[block_counts])
            seen[block_counts] = actual

    def test_overlapping_full_support_contributions_are_added(self):
        # Sum-zero inputs and nonzero-sum inputs with no cancelled coordinate
        # both produce full support here; neither branch may overwrite the other.
        zero_sum = (3**4+3)//4
        nonzero_sum_no_cancel = 3*(2**4+(-2)**4)//4
        self.assertEqual(column_support_counts(4, 4, 15)[15], zero_sum+nonzero_sum_no_cancel)
        self.assertEqual(zero_sum+nonzero_sum_no_cancel, 45)

    def test_width16_sparse_laws_and_all_weight_normalizations(self):
        self.assertEqual(column_weight_counts(256, 16, 0), [1]+[0]*16)
        singleton = column_weight_counts(256, 16, 1)
        self.assertEqual(singleton[15], 255)
        self.assertEqual(sum(singleton), 255)
        pair = column_weight_counts(256, 16, 2)
        self.assertEqual(pair[2], 255)
        self.assertEqual(pair[16], 255*254)
        self.assertEqual(sum(pair), 255**2)
        for active in range(1, 17):
            counts = column_weight_counts(256, 16, active)
            self.assertEqual(counts[0], 0)
            self.assertEqual(sum(counts), 255**active)
            self.assertTrue(all(count >= 0 for count in counts))


if __name__ == '__main__':
    unittest.main()
