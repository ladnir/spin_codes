"""Exact small-field checks for marked local input/output moments."""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
import random
import unittest

from hd_blockwise_shared_transfer import multiply
from hd_randomized_outer_joint_input import split_enumerator, marked_pattern, profile_moments
from hd_randomized_outer_quiet_tile import profile_moments as old_profiles


def tile(values):
    # Degree-below-two interpolation on GF4 nodes0,1, evaluated at2,3.
    a, b = values
    slope = a ^ b
    return (a ^ multiply(2, slope, 4), a ^ multiply(3, slope, 4))


class JointInputTests(unittest.TestCase):
    def test_exact_split_mds_enumerator(self):
        for a in range(3):
            for support in combinations(range(2), a):
                actual = Counter()
                for values in product(range(4), repeat=a):
                    x = [0, 0]
                    for i, value in zip(support, values):
                        x[i] = value
                    actual[sum(value != 0 for value in x), sum(value != 0 for value in tile(x))] += 1
                counts = split_enumerator(a, 2, 4)
                self.assertEqual(sum(map(sum, counts)), 4**a)
                for i, row in enumerate(counts):
                    for w, count in enumerate(row):
                        self.assertEqual(count, actual[i, w])

    def test_quiet_bound_and_exact_free_moment(self):
        z, v = F(1, 3), F(3, 2)
        for profile in product(range(3), repeat=2):
            for first in combinations(range(2), profile[0]):
                for second in combinations(range(2), profile[1]):
                    support = first+tuple(2+i for i in second)
                    quiet, free = F(0), F(0)
                    for values in product(range(4), repeat=len(support)):
                        x = [0]*4
                        for i, value in zip(support, values):
                            x[i] = value
                        y = tile(x[:2])+tile(x[2:])
                        weight = v**sum(value != 0 for value in x)*z**sum(value != 0 for value in y)
                        free += weight
                        if x[0] ^ x[1] ^ x[2] ^ x[3] == 0:
                            quiet += weight
                    bound, terminal = marked_pattern(profile, z, v, width=2, state=1,
                                                       field=4, distinct=False)
                    self.assertEqual(terminal, free/4**len(support))
                    self.assertGreaterEqual(bound, quiet/4**len(support))

    def test_residual_hyperplanes_from_distinct_tiles(self):
        from hd_blockwise_shared_maps import matrices, rank
        feedthrough, _, feedback, _, _ = matrices(64)
        rng = random.Random(78623)
        for _ in range(48):
            left, right = rng.sample(range(8), 2)
            a, b = rng.randrange(3, 9), rng.randrange(3, 9)
            columns = sorted(rng.sample(range(8*left, 8*left+8), a)
                             +rng.sample(range(8*right, 8*right+8), b))
            rows = rng.sample(range(8*left, 8*left+8), a-3)+rng.sample(range(8*right, 8*right+8), b-3)
            output = tuple(tuple(feedthrough[i][j] for j in columns) for i in rows)
            combined = output+tuple(tuple(row[j] for j in columns) for row in feedback)
            self.assertEqual(rank(combined)-rank(output), 4)

    def test_production_profile_normalization_and_caps(self):
        from math import comb
        z, v = F(1, 8), F(17, 16)
        original = {p:k for p, j, mult, k, t in old_profiles(z)}
        mass = [0]*65
        for profile, j, mult, quiet, free in profile_moments(z, v):
            mass[j] += mult
            zero = F(1, 65536**j)
            self.assertGreaterEqual(quiet, zero)
            self.assertLessEqual(quiet, free)
            self.assertLessEqual(quiet, zero+v**j*(original[profile]-zero))
        self.assertEqual(mass, [comb(64, j) for j in range(65)])


if __name__ == '__main__':
    unittest.main()
