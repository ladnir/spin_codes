"""Bounded exact algebra tests, no full-distance computation or files."""
from fractions import Fraction as F
from itertools import combinations
import unittest
import numpy as np
import hd_field128_dense_feedback_return_v42 as candidate
import hd_field128_feedback_rank_probe_v41 as gf


class Tests(unittest.TestCase):
    def test_probabilities(self):
        for j in range(65):
            for c in range(j+1):
                for e in range(65-j):
                    u = min(32, j-c)
                    rank = candidate.rank_floor(j, c, e)
                    self.assertTrue(max(0, min(e, 32)+u-32) <= rank <= min(e, u))
                    self.assertGreaterEqual(candidate.probability(j, c, e, 1 << 128), 0)
        self.assertEqual(candidate.rank_floor(1, 0, 1), 1)
        self.assertEqual(candidate.probability(1, 0, 1, 1 << 128), 0)
        # The remaining one-rank loss must not be silently replaced by full rank.
        self.assertEqual(candidate.rank_floor(32, 31, 31), 0)

    def test_toy_actual_ranks(self):
        # GF16 is the unique order16 subfield of AES GF256. The outer field
        # has degree2, so the same bordered-minor proof uses degree=2.
        table = np.array([[gf.mul(a, b) for b in range(256)] for a in range(256)], dtype=np.uint8)
        a = gf.power(3, 17)
        base = {0} | {gf.power(a, k) for k in range(15)}
        self.assertEqual(len(base), 16)
        points = [0, 1, a, a ^ 1]
        b, d = gf.mul(a, a), gf.power(a, 3)
        self.assertTrue(len(set(points+[b, b ^ 1, d, d ^ 1])) == 8)
        def interpolation(nodes):
            return np.array([[gf.mul(x ^ nodes[1], gf.power(nodes[0] ^ nodes[1], 254)),
                              gf.mul(x ^ nodes[0], gf.power(nodes[0] ^ nodes[1], 254))]
                             for x in points], dtype=np.uint8)
        emission, feedback = interpolation([b, b ^ 1]), interpolation([d, d ^ 1])
        p = gf.product(emission, feedback.T, table)
        self.assertFalse(np.any(gf.product(feedback.T, emission, table)))
        tested = 0
        subsets = [tuple(i for i in range(4) if mask >> i & 1) for mask in range(16)]
        for gamma in [x for x in range(256) if x not in base][:8]:
            l = np.eye(4, dtype=np.uint8) ^ table[gamma, p]
            if any(gf.rank(l[np.ix_(s, s)], table) != len(s) for s in subsets):
                continue
            for support in subsets:
                outside = [x for x in range(4) if x not in support]
                for c in range(len(support)+1):
                    for active in combinations(support, c):
                        for e in range(len(outside)+1):
                            for inactive in combinations(outside, e):
                                j, u = len(support), min(2, len(support)-c)
                                observed = gf.rank(l[np.ix_(active+inactive, support)], table)
                                exact = F(256**(j-observed)-256**(j-c-u), 256**j)
                                cap = candidate.probability(j, c, e, 256, state=2, degree=2)
                                self.assertLessEqual(exact, cap)
                                tested += 1
        self.assertGreater(tested, 1000)


if __name__ == '__main__':
    unittest.main()
