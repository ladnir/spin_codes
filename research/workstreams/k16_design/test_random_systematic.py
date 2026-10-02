"""Exact finite checks for the random-systematic expected-count identities."""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
import unittest

from random_systematic import (
    active_row_support_numerators,
    expected_group_support_counts,
    expected_row_weight_counts,
    group_support_layers,
    systematic_union_counts,
)


def union_weight(values):
    union = 0
    for value in values:
        union |= value
    return union.bit_count()


class RandomSystematicTests(unittest.TestCase):
    def test_systematic_messages_exhaustive(self):
        for n in range(1, 4):
            for h in range(4):
                counts = Counter(union_weight(messages)
                                 for messages in product(range(1, 1 << n), repeat=h))
                self.assertEqual(systematic_union_counts(n, h),
                                 tuple(counts[a] for a in range(n + 1)))

    def test_message_and_parity_pairs_exhaustive(self):
        # No probability formula is used by the enumeration: list every
        # nonzero message row and every possible independent parity row.
        for n in range(1, 4):
            for h in range(4):
                counts = Counter()
                parities = tuple(product(range(1 << n), repeat=h))
                for messages in product(range(1, 1 << n), repeat=h):
                    systematic = union_weight(messages)
                    for parity in parities:
                        counts[systematic + union_weight(parity)] += 1
                self.assertEqual(active_row_support_numerators(n, h),
                                 tuple(counts[v] for v in range(2 * n + 1)))

    def test_all_row_subsets_exhaustive(self):
        for n, rows in ((1, 1), (1, 4), (2, 2), (2, 3)):
            # An inactive message row has parity zero. Every active message
            # row has 2^n parity possibilities of equal probability.
            choices = [(0, 0)] + list(product(range(1, 1 << n), range(1 << n)))
            expected = [Fraction(0)] * (2 * n + 1)
            by_layer = [Counter() for _ in range(rows + 1)]
            for group in product(choices, repeat=rows):
                messages, parities = zip(*group)
                h = sum(value != 0 for value in messages)
                if not h:
                    continue
                support = union_weight(messages) + union_weight(parities)
                expected[support] += Fraction(1, 1 << (h * n))
                by_layer[h][support] += 1
            self.assertEqual(expected_group_support_counts(n, rows), tuple(expected))
            for layer in group_support_layers(n, rows):
                self.assertEqual(layer.denominator, 1 << (layer.active_rows * n))
                self.assertEqual(layer.numerators,
                                 tuple(by_layer[layer.active_rows][v] for v in range(2 * n + 1)))

    def test_sampled_matrices_exhaustive(self):
        # Average every binary R independently for each row. Matrices are
        # fixed while the complete message space is enumerated, as in setup.
        n, rows = 2, 2
        matrices = tuple(product(range(1 << n), repeat=n))
        counts = [0] * (2 * n + 1)
        for setup in product(matrices, repeat=rows):
            for messages in product(range(1 << n), repeat=rows):
                if not any(messages):
                    continue
                parity = tuple(sum(((mask & message).bit_count() & 1) << bit
                                   for bit, mask in enumerate(matrix))
                               for matrix, message in zip(setup, messages))
                counts[union_weight(messages) + union_weight(parity)] += 1
        self.assertEqual(expected_group_support_counts(n, rows),
                         tuple(Fraction(count, len(matrices) ** rows) for count in counts))

    def test_one_row_formula(self):
        for n in (1, 2, 3, 8, 128):
            direct = expected_row_weight_counts(n)
            self.assertEqual(direct, expected_group_support_counts(n, 1))
            self.assertEqual(sum(direct), (1 << n) - 1)
            self.assertEqual(direct[0], 0)
            self.assertEqual(direct[1], Fraction(n, 1 << n))

    def test_target_size_exact_counts(self):
        n, rows = 128, 4
        layers = group_support_layers(n, rows)
        counts = expected_group_support_counts(n, rows)
        self.assertEqual(len(counts), 257)
        self.assertEqual(sum(counts), (1 << 512) - 1)
        self.assertEqual(counts[0], 0)
        self.assertTrue(all(value > 0 for value in counts[1:]))
        # Support one requires a single common systematic position and zero
        # parity. Every active row then has the same one-bit message.
        self.assertEqual(counts[1],
                         sum((Fraction(comb(rows, h) * n, 1 << (h * n))
                              for h in range(1, rows + 1)), Fraction(0)))
        for layer in layers:
            h = layer.active_rows
            self.assertEqual(Fraction(sum(layer.numerators), layer.denominator),
                             comb(rows, h) * ((1 << n) - 1) ** h)
        cumulative = Fraction(0)
        for value in counts:
            cumulative += value
            self.assertGreaterEqual(cumulative, value)
        self.assertEqual(cumulative, (1 << 512) - 1)

    def test_invalid_parameters(self):
        for invalid in (0, -1, 1.5, True):
            for function in (expected_group_support_counts, expected_row_weight_counts,
                             group_support_layers, systematic_union_counts,
                             active_row_support_numerators):
                with self.assertRaises(ValueError):
                    function(invalid, 1) if function in (systematic_union_counts,
                                                       active_row_support_numerators) else function(invalid)
        for invalid in (0, -1, 1.5, True):
            with self.assertRaises(ValueError):
                expected_group_support_counts(2, invalid)
        for invalid in (-1, 1.5, True):
            with self.assertRaises(ValueError):
                systematic_union_counts(2, invalid)


if __name__ == "__main__":
    unittest.main()
