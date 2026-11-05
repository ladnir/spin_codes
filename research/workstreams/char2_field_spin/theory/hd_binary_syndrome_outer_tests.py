from collections import Counter
from itertools import combinations, product
import random
import unittest

from hd_binary_syndrome_outer import (
    matrix, column_mask, direct, circuit, primitive_counts,
    two_row_joint_counts, two_row_support_counts)
from hd_parity_outer import support_counts
from hd_bit_syndrome_outer import affine_rank


class BinarySyndromeTests(unittest.TestCase):
    def test_gram_zero_involution_and_circuit(self):
        v = tuple((1,)+tuple(i >> k & 1 for k in range(4)) for i in range(16))
        self.assertTrue(all(sum(row[a]*row[b] for row in v) % 2 == 0 for a in range(5) for b in range(5)))
        for i in range(16):
            x = tuple(int(j == i) for j in range(16))
            self.assertEqual(direct(direct(x)), x)
            self.assertEqual(circuit(x), direct(x))
        rng = random.Random(1729)
        for _ in range(24):
            x = tuple(rng.randrange(256) for _ in range(16))
            self.assertEqual(circuit(x), direct(x))
        self.assertEqual(primitive_counts()['binary_xors'], 57)

    def test_rank_five_quiet_minimum_four(self):
        self.assertEqual(affine_rank(65535), 5)
        for size in range(1, 4):
            for rows in combinations(range(16), size):
                self.assertEqual(affine_rank(sum(1 << i for i in rows)), size)
        plane = tuple(int(i < 4) for i in range(16))
        self.assertEqual(direct(plane), plane)
        self.assertEqual(affine_rank(15), 3)

    def test_every_two_source_GF4_value_for_both_mixers(self):
        for mixer in ('binary_syndrome', 'parity'):
            for first, second in combinations(range(16), 2):
                counts = Counter()
                for a, b in product(range(4), repeat=2):
                    x = [0]*16
                    x[first], x[second] = a, b
                    output = sum(bool(y) << i for i, y in enumerate(direct(x, mixer=mixer)))
                    counts[int(bool(a))+int(bool(b)), output] += 1
                self.assertEqual(two_row_joint_counts(first, second, mixer=mixer, q=4), counts)
                self.assertEqual(sum(two_row_joint_counts(first, second, mixer=mixer).values()), 256**2)

    def test_parity_interface_and_input_degree(self):
        for first, second in ((0, 1), (0, 8), (5, 14)):
            mask = (1 << first) | (1 << second)
            self.assertEqual(two_row_support_counts(first, second, mixer='parity'), support_counts(mask))
            by_weight = Counter()
            for (weight, _), count in two_row_joint_counts(first, second).items():
                by_weight[weight] += count
            self.assertEqual(by_weight, {0: 1, 1: 2*255, 2: 255**2})

    def test_singleton_tradeoff_and_distinct_construction(self):
        self.assertEqual(Counter(column_mask(i).bit_count() for i in range(16)), {7: 7, 9: 8, 15: 1})
        self.assertNotEqual(matrix(), matrix('parity'))
        with self.assertRaises(ValueError):
            two_row_joint_counts(3, 3)


if __name__ == '__main__':
    unittest.main()
