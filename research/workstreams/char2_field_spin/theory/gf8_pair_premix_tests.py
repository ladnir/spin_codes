"""Exact regressions for the deterministic pair-mix zero-flat cover."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
import math
from random import Random
import unittest

import numpy as np

from gf8_pair_premix_screen import mixed_cover, pair_patterns, weighted_matrices, packet
from gf8_high_rate_outer_screen import RETURN_MINIMA
from gf8_packet_width_screen import local_matrix
from rate_design_coset import mul, xor_sum
from state_capacity_reference import rank


GF4 = (0, 1, 188, 189)


def classify(forms):
    """Enumerate exactly all feasible mixed zero flats in the GF4 model."""
    dimension = len(forms[0])
    mixed = []
    singles = doubles = 0
    eligible = []
    for start in range(0, len(forms), 2):
        first, second = forms[start:start+2]
        occupied = int(any(first))+int(any(second))
        singles += occupied == 1
        doubles += occupied == 2
        pair = [tuple(a^b for a, b in zip(first, second)),
                tuple(a^mul(188, b) for a, b in zip(first, second))]
        if occupied:
            eligible.extend((start, start+1))
        mixed.extend(pair)
    flats = set()
    for witness in product(GF4, repeat=dimension):
        evaluate = lambda form: xor_sum(mul(a, b) for a, b in zip(form, witness))
        if any(any(form) and not evaluate(form) for form in forms):
            continue
        zeros = frozenset(index for index in eligible if not evaluate(mixed[index]))
        flats.add(zeros)
    return singles, doubles, [(len(zeros), rank([mixed[index] for index in sorted(zeros)]))
                              for zeros in flats]


class PairPremixTests(unittest.TestCase):
    def test_support_multiplicities(self):
        by_weight = Counter()
        for singles, doubles, weight, multiplicity in pair_patterns():
            self.assertEqual(weight, singles+2*doubles)
            by_weight[weight] += multiplicity
        self.assertEqual(dict(by_weight), {w: math.comb(8, w) for w in range(9)})

    def test_exact_singleton_and_empty(self):
        y, tau = F(2, 3), F(11)
        for singles in range(5):
            self.assertEqual(mixed_cover(singles, 0, y, tau),
                             local_matrix(2*singles, y, tau, 8, 4, RETURN_MINIMA))

    def test_exact_small_field_zero_flat_count_and_cover(self):
        rng = Random(21496)
        cases = [tuple((v,) for v in row) for row in product(GF4, repeat=4)]
        two_forms = tuple(product(GF4, repeat=2))
        cases += [(f, g, (0, 0), (0, 0)) for f, g in product(two_forms, repeat=2)]
        cases += [tuple(rng.choice(two_forms) for _ in range(4)) for _ in range(256)]
        y, tau = F(2, 3), F(7)
        locals_ = [local_matrix(w, y, tau, 8, 4, RETURN_MINIMA) for w in range(9)]
        for forms in cases:
            singles, doubles, flats = classify(forms)
            counts = Counter(k for _, k in flats)
            for k, count in counts.items():
                self.assertLessEqual(count, math.comb(doubles, k)*2**k)
            for z, k in flats:
                self.assertLessEqual(k, z)
                self.assertLessEqual(z, doubles)
            upper = mixed_cover(singles, doubles, y, tau, locals_)
            for row, column in product(range(2), repeat=2):
                actual = sum(locals_[2*singles+2*doubles-z][row][column]/tau**k for z, k in flats)
                self.assertLessEqual(actual, upper[row][column])
        self.assertEqual(len(cases), 768)

    def test_deterministic_cancellation_can_have_rank_zero(self):
        # Both raw coordinates are required live, but f+f is identically zero.
        singles, doubles, flats = classify(((1,), (1,), (0,), (0,)))
        self.assertEqual((singles, doubles), (0, 1))
        self.assertEqual(flats, [(1, 0)])

    def test_packet_matches_fraction_sum(self):
        x, y, tau = F(2, 5), F(3, 4), F(13)
        local = weighted_matrices(y, tau)
        expected = [[sum(x**w*local[w][r][c] for w in range(9)) for c in range(2)] for r in range(2)]
        actual = np.exp(packet(math.log(x), math.log(y), math.log(tau)))
        np.testing.assert_allclose(actual, np.asarray(expected, dtype=float), rtol=2e-12)


if __name__ == '__main__':
    unittest.main()
