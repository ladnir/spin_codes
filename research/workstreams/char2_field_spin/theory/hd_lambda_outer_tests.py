from collections import Counter
from itertools import combinations, product
import random
import unittest

from hd_lambda_outer import (
    FULL, matrix, direct, circuit, primitive_counts,
    exceptional_ratios, column_branches, joint_counts)
from hd_bit_syndrome_outer import affine_rank


class LambdaOuterTests(unittest.TestCase):
    def test_involution_nonzero_entries_and_circuit(self):
        for q, lam in ((4, 2), (4, 3), (256, 2)):
            self.assertEqual(set(x for row in matrix(lam, q) for x in row), {lam, lam ^ 1})
            for i in range(16):
                x = tuple(int(j == i) for j in range(16))
                self.assertEqual(direct(direct(x, lam, q), lam, q), x)
                self.assertEqual(circuit(x, lam, q), direct(x, lam, q))
            rng = random.Random(20261009)
            for _ in range(12):
                x = tuple(rng.randrange(q) for _ in range(16))
                self.assertEqual(circuit(x, lam, q), direct(x, lam, q))
        self.assertEqual(primitive_counts()['binary_xors'], 57)
        self.assertEqual(primitive_counts()['fixed_products'], 1)

    def test_quiet_kernel_and_minimum_four(self):
        self.assertEqual(affine_rank(FULL), 5)
        for a in range(1, 4):
            for sources in combinations(range(16), a):
                self.assertEqual(affine_rank(sum(1 << i for i in sources)), a)
        plane = tuple(int(i < 4) for i in range(16))
        self.assertEqual(direct(plane), plane)
        for sources in combinations(range(4), 3):
            for values in product(range(1, 4), repeat=3):
                x = [0]*16
                for i, value in zip(sources, values):
                    x[i] = value
                self.assertNotEqual(direct(x, q=4), tuple(x))

    def test_all_singletons_and_ordered_controller_schema(self):
        for i in range(16):
            self.assertEqual(column_branches((i,)), ((0, 0, 1), (1, FULL, 255)))
        for a, b in combinations(range(16), 2):
            forward = column_branches((a, b))
            backward = column_branches((b, a))
            self.assertEqual(len(forward), 7)
            self.assertEqual(tuple(count for _, _, count in forward), (1, 255, 255, 255, 255, 255, 255*252))
            self.assertEqual(forward[3], backward[3])
            self.assertEqual(forward[4], backward[5])
            self.assertEqual(forward[5], backward[4])
            self.assertEqual(joint_counts((a, b)), joint_counts((b, a)))

    def test_all_GF4_values_for_all120_pairs(self):
        for a, b in combinations(range(16), 2):
            expected = Counter()
            for x, y in product(range(4), repeat=2):
                values = [0]*16
                values[a], values[b] = x, y
                mask = sum(bool(v) << i for i, v in enumerate(direct(values, q=4)))
                expected[int(bool(x))+int(bool(y)), mask] += 1
            self.assertEqual(joint_counts((a, b), q=4), expected)
            self.assertEqual(column_branches((a, b), q=4)[-1], (2, FULL, 0))

    def test_field_ratios_and_validation(self):
        for q in (4, 256):
            ratios = exceptional_ratios(q=q)
            self.assertEqual(len(set(ratios)), 3)
            self.assertNotIn(0, ratios)
        for sources in ((), (1, 1), (0, 1, 2), (16,)):
            with self.assertRaises(ValueError):
                column_branches(sources)
        with self.assertRaises(ValueError):
            matrix(lambda_value=1)
        with self.assertRaises(ValueError):
            matrix(q=16)


if __name__ == '__main__':
    unittest.main()
