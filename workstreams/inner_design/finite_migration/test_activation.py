from collections import Counter
from fractions import Fraction as F
import itertools
import math
import unittest
from flint import arb, ctx
import activation_density as density


class ActivationTests(unittest.TestCase):
    def test_pointwise_domination_exhaustive_toy(self):
        ctx.prec = 128
        feedback = (1, 2, 3, 4, 5, 7)
        expansion = (1, 2, 3, 4, 6, 7)
        weights = {q: sum((q & a).bit_count() % 2 for a in expansion) for q in range(1, 8)}
        spectrum = dict(Counter(weights.values()))
        levels = sorted(spectrum)
        for j in range(7):
            counts = Counter()
            for subset in itertools.combinations(range(6), j):
                state = 0
                for i in subset:
                    state ^= feedback[i]
                counts[state] += 1
            total = math.comb(6, j)
            cap = max([counts[q] for q in range(1, 8)])
            row = density.zero_row(spectrum, counts[0], total, cap, arb(1) / 2**j)
            self.assertEqual(row[1], 0)
            for q in range(8):
                exact = arb(counts[q]) / (total * 2**j)
                upper = row[0] if q == 0 else row[levels.index(weights[q]) + 2] / spectrum[weights[q]]
                # The reconstructed uniform-shell measure dominates every
                # individual syndrome, not merely the total mass.
                self.assertLessEqual(exact.upper(), upper.upper())

    def test_instance_identity_unchanged(self):
        self.assertEqual(density.Engine(16).identity(), density.ladder.Engine(16).identity())


if __name__ == '__main__':
    unittest.main()
