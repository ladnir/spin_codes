import itertools
import unittest
from fractions import Fraction as Q

import canonical_counts as counts
import local_models as model


def binary_rank(rows):
    basis = {}
    for row in rows:
        while row:
            high = row.bit_length()-1
            if high not in basis:
                basis[high] = row
                break
            row ^= basis[high]
    return len(basis)


class CanonicalCounts(unittest.TestCase):
    def setUp(self):
        self.length, self.width, self.rows = 8, 2, 2
        generators = [0b11000011, 0b00111100, 0b10101010]
        self.code = [0]
        for row in generators:
            self.code += [word ^ row for word in self.code]
        self.dimensions = []
        for size in range(9):
            largest = 0
            for positions in itertools.combinations(range(8), size):
                mask = sum(1 << i for i in positions)
                mass = sum((word & ~mask) == 0 for word in self.code)
                self.assertEqual(mass & (mass-1), 0)
                largest = max(largest, mass.bit_length()-1)
            self.dimensions.append(largest)

    def test_rank_counts(self):
        for dimension in range(4):
            actual = [0]*3
            for pair in itertools.product(range(1 << dimension), repeat=2):
                actual[binary_rank(pair)] += 1
            self.assertEqual(actual, [counts.rank_total(dimension, 2, r) for r in range(3)])

    def test_all_tuple_block_counts_are_dominated(self):
        actual = [[0]*5 for _ in range(2)]
        original = [[0]*9 for _ in range(2)]
        for words in itertools.product(self.code, repeat=2):
            r = binary_rank(words)
            if not r:
                continue
            union = words[0] | words[1]
            h = sum(bool(union & (3 << (2*i))) for i in range(4))
            for bound in range(h, 5):
                actual[r-1][bound] += 1
            for bound in range(union.bit_count(), 9):
                original[r-1][bound] += 1
        for original_arg in (None, original):
            caps = counts.canonical_rank_caps(self.dimensions, 2, 2, original_arg)
            for row, upper in zip(actual, caps):
                self.assertTrue(all(a <= b for a, b in zip(row, upper)))
                self.assertEqual(row[-1], upper[-1])

    def test_transport_against_exact_profile_spectrum(self):
        shells = [0]*5
        for words in itertools.product(self.code, repeat=2):
            union = words[0] | words[1]
            if union:
                shells[sum(bool(union & (3 << (2*i))) for i in range(4))] += 1
        actual_cdf = list(itertools.accumulate(shells))
        caps = counts.total_caps(counts.canonical_rank_caps(self.dimensions, 2, 2))
        local = model.full_block(2, 2)
        exact = counts.transport_cdf(actual_cdf, local)
        upper = counts.transport_cdf(caps, local)
        self.assertTrue(all(a <= b for a, b in zip(exact, upper)))
        self.assertEqual(exact[-1], len(self.code)**2-1)
        # Independent direct enumeration of each nonzero4bit block's output.
        direct = [Q(0)]*9
        for h, mass in enumerate(shells):
            if not mass:
                continue
            for blocks in itertools.product(range(1, 16), repeat=h):
                w = sum(bool(block & 3) + bool(block & 12) for block in blocks)
                direct[w] += Q(mass, 15**h)
        self.assertEqual(exact, tuple(itertools.accumulate(direct)))

    def test_invalid_interfaces(self):
        with self.assertRaises(ValueError):
            counts.canonical_rank_caps([0, 0, 1], width=3)
        with self.assertRaises(ValueError):
            counts.canonical_rank_caps([0, 1, 0], width=1)
        with self.assertRaises(ValueError):
            counts.transport_cdf([0, 3, 2], model.full_block(2))
        with self.assertRaises(ValueError):
            counts.transport_cdf([0, 1], (Q(1), Q(0), Q(0)))

    def test_width_one_is_support_preserving(self):
        caps = counts.total_caps(counts.canonical_rank_caps(self.dimensions, 1, 2))
        self.assertEqual(model.full_block(1, 2), (Q(0), Q(1)))
        self.assertEqual(counts.transport_cdf(caps, model.full_block(1, 2)), caps)

    def test_width_compare_uses_same_table(self):
        import width_compare
        original = self.dimensions[:]
        result = width_compare.compare(self.dimensions, (1, 2, 4, 8))
        self.assertEqual(self.dimensions, original)
        for width, row in result.items():
            self.assertEqual(len(row['output_cdf']), 9)
            self.assertEqual(row['output_cdf'][-1], (1 << (4*self.dimensions[-1]))-1)
            self.assertEqual(row['output_cdf'], counts.transport_cdf(
                counts.total_caps(counts.canonical_rank_caps(self.dimensions, width)), model.full_block(width)))


if __name__ == '__main__':
    unittest.main()
