"""Exhaustive shared-setup checks of the inverse-probability witness.

Row randomness changes which candidates are retained; region randomness
is shared by all candidates. Conditional probabilities are nonconstant.
"""
from fractions import Fraction as F
import itertools
import unittest


def permute(word, pi):
    return sum(((word >> pi[j]) & 1) << j for j in range(4))


def passes(word):
    return (word & 3).bit_count() % 2 == 0 and (word >> 2).bit_count() % 2 == 0


def row_swap(pair, row):
    if row is None:
        return pair
    a, b = pair
    if ((a ^ b) >> row) & 1:
        a ^= 1 << row
        b ^= 1 << row
    return a, b


class WeightedWitnessTests(unittest.TestCase):
    def test_full_shared_setup(self):
        candidates = [(3, 3), (5, 15), (15, 5)]
        permutations = list(itertools.permutations(range(4)))
        total, square, nonzero = F(0), F(0), F(0)
        formula_mean, formula_second = F(0), F(0)
        for row in [None, 1]:
            words = [row_swap(pair, row) for pair in candidates]
            keep = [all(word.bit_count() in [2, 4] for word in pair) for pair in words]
            probabilities = []
            for pair, retained in zip(words, keep):
                if retained:
                    probability = F(1)
                    for word in pair:
                        probability *= F(sum(passes(permute(word, pi)) for pi in permutations), len(permutations))
                    probabilities.append(probability)
                else:
                    probabilities.append(F(0))
            formula_mean += sum(keep)/F(2)
            for i, j in itertools.product(range(len(words)), repeat=2):
                if not (keep[i] and keep[j]):
                    continue
                joint = F(1)
                for a, b in zip(words[i], words[j]):
                    joint *= F(sum(passes(permute(a, pi)) and passes(permute(b, pi))
                                   for pi in permutations), len(permutations))
                formula_second += joint/(2*probabilities[i]*probabilities[j])
            for first, second in itertools.product(permutations, repeat=2):
                value = sum((F(passes(permute(pair[0], first)) and passes(permute(pair[1], second)))/probability
                             if probability else F(0)) for pair, probability in zip(words, probabilities))
                total += value
                square += value*value
                nonzero += value > 0
        setups = 2*len(permutations)**2
        mean, second, probability = total/setups, square/setups, nonzero/setups
        self.assertEqual(mean, F(2))
        self.assertEqual(mean, formula_mean)
        self.assertEqual(second, formula_second)
        self.assertLessEqual(mean*mean/second, probability)
        self.assertEqual(probability, F(7, 18))
        self.assertGreater(second, mean*mean)


if __name__ == '__main__':
    unittest.main()
