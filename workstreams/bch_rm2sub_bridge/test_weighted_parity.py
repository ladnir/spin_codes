"""Exhaustive marked-row Fourier checks, retaining shared permutations."""
from fractions import Fraction as F
import itertools
import math
import unittest


class WeightedParityTests(unittest.TestCase):
    def test_two_row_shared_permutation_experiment(self):
        n, q = 4, 2
        # A transitive shell, but not the full weight-two slice.
        word = (1, 1, 0, 0)
        family = [word[j:]+word[:j] for j in range(n)]
        masks = [(1 << j) | (1 << ((j+1) % n)) for j in range(n)]
        chars = range(1 << n)
        trivial = {0, (1 << n)-1}
        char_sign = lambda v, x: (-1)**((v & x).bit_count())
        # Aggregate exactly the law of one independent row pair and shared pi.
        law = {}
        for perm in itertools.permutations(range(n)):
            for u, v in itertools.product(masks, repeat=2):
                a = sum(((u >> perm[j]) & 1) << j for j in range(n))
                b = sum(((v >> perm[j]) & 1) << j for j in range(n))
                law[a, b] = law.get((a, b), F(0))+F(1, math.factorial(n)*len(family)**2)
        single = max(abs(sum(prob*char_sign(c, a) for (a, b), prob in law.items()))
                     for c in chars if c not in trivial)
        paired = max(abs(sum(prob*char_sign(c, a)*char_sign(d, b)
                             for (a, b), prob in law.items()))
                     for c in chars if c not in trivial for d in chars if d not in trivial)
        self.assertLessEqual(single*single, paired)
        # The monomial (J_0/q) is expanded into two marked-row terms.
        for c, d in itertools.product(chars, repeat=2):
            weighted = F(0)
            for ((a, b), prob), ((e, f), prob2) in itertools.product(law.items(), repeat=2):
                weight = F((a & 1)+(e & 1), q)
                weighted += prob*prob2*weight*char_sign(c, a ^ e)*char_sign(d, b ^ f)
            if c in trivial and d in trivial:
                self.assertEqual(weighted, F(1, 2))
            elif c in trivial or d in trivial:
                self.assertLessEqual(abs(weighted), single**(q-1))
            else:
                self.assertLessEqual(abs(weighted), paired**(q-1))

    def test_remaining_support_inequality(self):
        a, b = F(3, 8), F(23, 64)
        for q1 in range(1, 10):
            for q2 in range(1, 10):
                for overlap in range(min(q1, q2)+1):
                    self.assertLessEqual(a**(q1+q2-2*overlap)*b**overlap,
                                         b**min(q1, q2))


if __name__ == '__main__':
    unittest.main()
