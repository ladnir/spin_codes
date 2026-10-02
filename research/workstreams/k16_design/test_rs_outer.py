"""Exact finite checks of the MDS enumerator and independent packet transport."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
import unittest

import rs_outer
from rs_outer import (
    expected_group_support_counts,
    group_support_layers,
    mds_symbol_weight_counts,
    packet_support_numerators,
)


def evaluation_word(coefficients, q, n):
    """Evaluate every message polynomial at 0..n-1 in a prime field."""
    return tuple(sum(coefficient * pow(x, degree, q)
                     for degree, coefficient in enumerate(coefficients)) % q
                 for x in range(n))


class RsOuterTests(unittest.TestCase):
    def test_prime_field_rs_exhaustive(self):
        # Independent construction, with no MDS counting formula used.
        for q, n, k in ((2, 2, 1), (3, 3, 3), (5, 5, 2),
                        (5, 4, 3), (7, 6, 3), (7, 7, 2), (7, 4, 1)):
            counts = Counter(sum(value != 0 for value in evaluation_word(message, q, n))
                             for message in product(range(q), repeat=k))
            self.assertEqual(mds_symbol_weight_counts(q, n, k),
                             tuple(counts[w] for w in range(n + 1)))
            self.assertEqual(sum(counts.values()) - counts[0], q ** k - 1)

    def test_parallel_base_field_rows_exhaustive(self):
        # Concatenating r base-field row symbols produces the same block
        # spectrum as the MDS formula for alphabet size p**r. No extension
        # field arithmetic or weight-enumerator identity is used here.
        for p, n, k, rows in ((2, 2, 1, 4), (3, 3, 2, 2), (5, 5, 2, 2)):
            counts = Counter()
            words = tuple(evaluation_word(message, p, n)
                          for message in product(range(p), repeat=k))
            for parallel in product(words, repeat=rows):
                counts[sum(any(row[i] for row in parallel) for i in range(n))] += 1
            self.assertEqual(mds_symbol_weight_counts(p ** rows, n, k),
                             tuple(counts[w] for w in range(n + 1)))

    def test_alternative_mds_identity(self):
        # Equivalent (q-1)-factored form catches summation-endpoint mistakes.
        for q, n, k in ((4, 4, 2), (8, 8, 3), (16, 8, 4), (1 << 32, 8, 4)):
            d = n - k + 1
            alternative = [0] * (n + 1)
            alternative[0] = 1
            for w in range(d, n + 1):
                alternative[w] = comb(n, w) * (q - 1) * sum(
                    (-1) ** j * comb(w - 1, j) * q ** (w - d - j)
                    for j in range(w - d + 1))
            self.assertEqual(mds_symbol_weight_counts(q, n, k), tuple(alternative))

    def test_nonzero_block_support_exhaustive(self):
        for bits, packets in ((1, 1), (1, 4), (2, 2), (2, 4), (4, 2)):
            mask = (1 << bits) - 1
            counts = Counter(sum(bool((word >> (i * bits)) & mask) for i in range(packets))
                             for word in range(1, 1 << (bits * packets)))
            self.assertEqual(packet_support_numerators(bits, packets),
                             tuple(counts[w] for w in range(packets + 1)))

    def test_fixed_gl_setups_and_messages_exhaustive(self):
        # The systematic [3,2] code (u,v,u+v) over GF(4) is MDS. Addition
        # is bitwise XOR; multiplication is unnecessary for this generator.
        # Enumerate all independent GL(2,2) setup matrices, then all messages
        # with each setup held fixed. Each bit is a one-bit packet here.
        matrices = tuple((a, b) for a in range(1, 4) for b in range(1, 4) if a != b)
        self.assertEqual(len(matrices), 6)
        counts = [0] * 7
        for setup in product(matrices, repeat=3):
            for u, v in product(range(4), repeat=2):
                if u == v == 0:
                    continue
                support = sum((mask & word).bit_count() & 1
                              for matrix, word in zip(setup, (u, v, u ^ v))
                              for mask in matrix)
                counts[support] += 1
        self.assertEqual(expected_group_support_counts(3, 2, 1, 2),
                         tuple(Fraction(count, len(matrices) ** 3) for count in counts))

    def test_full_dimension_transport(self):
        # A full-space code remains the full binary space under every GL
        # setup. Its packet enumerator therefore needs no averaging.
        for n, bits, packets in ((1, 1, 1), (2, 1, 2), (3, 2, 2)):
            count = n * packets
            exact = (Fraction(0),) + tuple(Fraction(comb(count, w) * ((1 << bits) - 1) ** w)
                                          for w in range(1, count + 1))
            self.assertEqual(expected_group_support_counts(n, n, bits, packets), exact)

    def test_target_size_exact_counts_and_layers(self):
        q = 1 << 32
        symbol_counts = mds_symbol_weight_counts()
        self.assertEqual(symbol_counts[:5], (1, 0, 0, 0, 0))
        self.assertEqual(symbol_counts[5], comb(8, 5) * (q - 1))
        self.assertEqual(sum(symbol_counts[1:]), (1 << 128) - 1)
        self.assertEqual(packet_support_numerators(),
                         (0,) + tuple(comb(8, w) * 15 ** w for w in range(1, 9)))
        layers = group_support_layers()
        self.assertEqual(tuple(layer.symbol_weight for layer in layers), (5, 6, 7, 8))
        for layer in layers:
            h = layer.symbol_weight
            self.assertEqual(len(layer.numerators), 65)
            self.assertEqual(layer.denominator, (q - 1) ** h)
            self.assertEqual(Fraction(sum(layer.numerators), layer.denominator), symbol_counts[h])
            self.assertTrue(all(value == 0 for value in layer.numerators[:h]))
            self.assertTrue(all(value == 0 for value in layer.numerators[8 * h + 1:]))
        counts = expected_group_support_counts()
        self.assertEqual(len(counts), 65)
        self.assertEqual(sum(counts), (1 << 128) - 1)
        self.assertEqual(counts[:5], (Fraction(0),) * 5)
        self.assertTrue(all(value > 0 for value in counts[5:]))
        self.assertEqual(counts[5], Fraction(comb(8, 5) * (8 * 15) ** 5, (q - 1) ** 4))

    def test_k16_geometry(self):
        self.assertEqual(rs_outer.GROUP_DIMENSION,
                         rs_outer.RS_DIMENSION * rs_outer.PACKET_BITS * rs_outer.PACKETS_PER_SYMBOL)
        self.assertEqual(rs_outer.GROUP_OUTPUT_BITS, 2 * rs_outer.GROUP_DIMENSION)
        self.assertEqual(rs_outer.REGION_COUNT,
                         rs_outer.RS_LENGTH * rs_outer.PACKETS_PER_SYMBOL)
        self.assertEqual(rs_outer.GROUP_COUNT * rs_outer.GROUP_DIMENSION, 1 << 16)
        self.assertEqual(rs_outer.GROUP_COUNT * rs_outer.GROUP_OUTPUT_BITS, 2 * (1 << 16))
        self.assertEqual(rs_outer.GROUP_COUNT * rs_outer.PACKET_BITS,
                         rs_outer.PHYSICAL_STEPS_PER_REGION * rs_outer.PHYSICAL_T)
        self.assertEqual(rs_outer.PHYSICAL_STEPS_PER_REGION * rs_outer.PHYSICAL_T,
                         rs_outer.MACROS_PER_REGION * rs_outer.MACRO_T)
        self.assertEqual(rs_outer.REGION_COUNT * rs_outer.MACROS_PER_REGION * rs_outer.MACRO_T,
                         2 * rs_outer.MESSAGE_BITS)

    def test_invalid_parameters(self):
        for invalid in (0, -1, 1.5, True):
            for index in range(3):
                args = [16, 8, 4]
                args[index] = invalid
                with self.assertRaises(ValueError):
                    mds_symbol_weight_counts(*args)
        for args in ((1, 1, 1), (6, 4, 2), (12, 4, 2), (4, 5, 2), (7, 4, 5)):
            with self.assertRaises(ValueError):
                mds_symbol_weight_counts(*args)
        for invalid in (0, -1, 1.5, True):
            for kwargs in ({"packet_bits": invalid}, {"packets_per_symbol": invalid},
                           {"n": invalid}, {"k": invalid}):
                with self.assertRaises(ValueError):
                    expected_group_support_counts(**kwargs)
            for args in ((invalid, 2), (2, invalid)):
                with self.assertRaises(ValueError):
                    packet_support_numerators(*args)


if __name__ == "__main__":
    unittest.main()
