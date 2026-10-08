"""Bounded exact tests for the new diagonal-MDS-diagonal state mixer."""

from fractions import Fraction as F
from itertools import combinations, product
import unittest

import hd_blockwise_diagonal_state as proof
from hd_blockwise_shared_transfer import all_subspaces, dot, rank


class DiagonalStateTests(unittest.TestCase):
    def test_actual_fixed_map_is_superregular(self):
        h = proof.cauchy_state_map()
        checked = 0
        for size in range(1, 5):
            for rows in combinations(range(4), size):
                for columns in combinations(range(4), size):
                    minor = tuple(tuple(h[i][j] for j in columns) for i in rows)
                    self.assertEqual(rank(minor, 256), size)
                    checked += 1
        self.assertEqual(checked, 69)

    def test_all_small_fixed_point_probabilities(self):
        mixers = tuple(proof.small_mixers())
        vectors = tuple(v for v in product(range(4), repeat=2) if any(v))
        self.assertEqual(len(mixers), 81)
        for h in vectors:
            counts = {c: 0 for c in vectors}
            for mixer in mixers:
                self.assertEqual(rank(mixer, 4), 2)
                c = tuple(dot(row, h, 4) for row in mixer)
                counts[c] += 1
            for c, count in counts.items():
                expected = proof.point_probability(sum(x != 0 for x in h),
                                                   sum(x != 0 for x in c), 2, 4)
                self.assertEqual(F(count, 81), expected)
                self.assertLessEqual(expected, F(1, 9))

    def test_all_small_joint_subspaces(self):
        mixers = tuple(proof.small_mixers())
        checked = 0
        for joint in all_subspaces(4, 4):
            dimension = len(joint[0])
            state, feedback = joint[:2], joint[2:]
            hits = 0
            for mixer in mixers:
                successor = tuple(tuple(dot(mixer[i], state_column, 4) ^ feedback[i][j]
                                        for j, state_column in enumerate(zip(*state)))
                                  for i in range(2))
                hits += int(rank(successor, 4) < dimension)
            self.assertLessEqual(F(hits, 81), proof.joint_return_bound(dimension, 2, 4))
            checked += 1
        self.assertEqual(checked, 529)

    def test_directed_k20_receipt(self):
        receipt = proof.k20_receipt()
        self.assertEqual(receipt['exact_rational_upper'], str(F(65537, 65535**4)))
        self.assertIn('47.9998', receipt['failure_bits_lower'])
        self.assertTrue(receipt['quiet_subcode_unchanged'])
        self.assertFalse(receipt['whole_code_certificate'])


if __name__ == '__main__':
    unittest.main()
