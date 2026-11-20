"""Exact conditioned-character square identity for a nonuniform slice subset."""
import itertools
import math
from fractions import Fraction as F
import unittest
from certify_tail_parity_mixing import kraw_values


class ConditionedPairTests(unittest.TestCase):
    def test_conditional_square_and_cauchy(self):
        n = 6
        family = [(1 << j) | (1 << ((j+1) % n)) for j in range(n)]
        conditioned = [[word >> 1 for word in family if word & 1 == bit] for bit in [0, 1]]
        permutations = list(itertools.permutations(range(n-1)))
        means = [[], []]
        for bit in [0, 1]:
            words = conditioned[bit]
            for perm in permutations:
                row = []
                for char in range(1 << (n-1)):
                    mapped = sum(((char >> j) & 1) << perm[j] for j in range(n-1))
                    row.append(sum(F((-1)**((word & mapped).bit_count()), len(words)) for word in words))
                means[bit].append(row)
            for char in range(1, (1 << (n-1))-1):
                h = char.bit_count()
                square = sum(row[char]**2 for row in means[bit])/len(permutations)
                distance = sum(F(kraw_values(n-1, (u ^ v).bit_count())[h], math.comb(n-1, h))
                               for u, v in itertools.product(words, repeat=2))/len(words)**2
                self.assertEqual(square, distance)
        for a, b in [(1, 2), (3, 5), (7, 14), (15, 23)]:
            mixed = sum(x[a]*y[b] for x, y in zip(*means))/len(permutations)
            first = sum(x[a]**2 for x in means[0])/len(permutations)
            second = sum(y[b]**2 for y in means[1])/len(permutations)
            self.assertLessEqual(mixed**2, first*second)


if __name__ == '__main__':
    unittest.main()
