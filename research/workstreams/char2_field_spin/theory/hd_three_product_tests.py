"""Exact bounded checks, including independent finite-field evaluation."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest

from hd_additive_two_state import support_law as original_law
from hd_rational_post_tests import mul
from hd_three_product import collapsed_majorant, exact_probability, support_law


def encode(source, setup):
    first = second = 0
    out = []
    for x, (d1, d2, f) in zip(source, setup):
        first ^= mul(d1, x)
        out.append(first ^ mul(f, second))
        second ^= mul(d2, x)
    return tuple(out)


class ThreeProductTests(unittest.TestCase):
    def test_all_gf4_local_values_and_local_invertibility(self):
        tuples = list(product(range(1, 4), repeat=3))
        for setup in tuples:
            d1, d2, f = setup
            images = set()
            for s1, s2, x in product(range(4), repeat=3):
                a, b = s1 ^ mul(d1, x), s2 ^ mul(d2, x)
                images.add((a ^ mul(f, s2), a, b))
            self.assertEqual(len(images), 64)
        for s1, s2, x in product(range(4), repeat=3):
            count = Counter()
            for d1, d2, f in tuples:
                a, b = s1 ^ mul(d1, x), s2 ^ mul(d2, x)
                y = a ^ mul(f, s2)
                count[2*int(bool(a))+int(bool(b)), int(bool(y))] += 1
            law = {key: F(value, len(tuples)) for key, value in count.items()}
            self.assertEqual(law, support_law(2*int(bool(s1))+int(bool(s2)), int(bool(x)), 4))

    def test_all_length_two_setups_and_messages(self):
        tuples = list(product(range(1, 4), repeat=3))
        counts = [[0]*3 for _ in range(3)]
        for setup in product(tuples, repeat=2):
            images = set()
            for source in product(range(4), repeat=2):
                output = encode(source, setup)
                images.add(output)
                counts[sum(bool(x) for x in source)][sum(bool(y) for y in output)] += 1
            self.assertEqual(len(images), 16)
        for weight in range(3):
            denominator = len(tuples)**2 * comb(2, weight) * 3**weight
            for cutoff in range(3):
                self.assertEqual(F(sum(counts[weight][:cutoff+1]), denominator),
                                 exact_probability(2, weight, cutoff, 4))

    def test_probability_mass_and_changed_law(self):
        for q, old, active in product((4, 8, 16, 256, 65536), range(4), range(2)):
            law = support_law(old, active, q)
            self.assertEqual(sum(law.values()), 1)
            self.assertTrue(all(value > 0 for value in law.values()))
        self.assertNotEqual(support_law(2, 1, 4), original_law(2, 1, 4))

    def test_common_positive_envelope_both_ensembles(self):
        for q, mark in product((4, 8, 256, 65536), (F(1, 7), F(1, 2), F(1))):
            bound = collapsed_majorant(mark, q)
            for law_function, old, active, target in product((support_law, original_law), range(4), range(2), range(2)):
                actual = sum((mass*mark**out for (new, out), mass in law_function(old, active, q).items()
                              if int(bool(new)) == target), F(0))
                self.assertLessEqual(actual, bound[active][int(bool(old))][target])

    def test_inverse_tail_has_rank_two(self):
        # Deterministic bounded fixture search; values and supports are exact.
        tuples = list(product(range(1, 4), repeat=3))
        found = False
        for index in range(50):
            setup = tuple(tuples[(index*(2*i+1)+7*i) % len(tuples)] for i in range(4))
            columns = []
            for j in range(2):
                first = second = 0
                source = []
                for i, (d1, d2, f) in enumerate(setup):
                    x = mul(mul(d1, d1), int(i == j)^first^mul(f, second))
                    source.append(x)
                    first ^= mul(d1, x)
                    second ^= mul(d2, x)
                self.assertEqual(encode(source, setup), tuple(int(i == j) for i in range(4)))
                columns.append(source)
            if mul(columns[0][2], columns[1][3]) ^ mul(columns[0][3], columns[1][2]):
                found = True
                break
        self.assertTrue(found)


if __name__ == '__main__':
    unittest.main()
