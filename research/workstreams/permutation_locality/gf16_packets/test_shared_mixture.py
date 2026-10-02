import unittest
from fractions import Fraction as Q
from itertools import product
from math import comb
from collections import Counter

import shared_mixture as shared
from test_shared_support import exact_law
from joint_support import span


class SharedMixtureTests(unittest.TestCase):
    def test_exact_shell_domination(self):
        caps = [0, 2, 7, 5, 1]
        mixture = shared.envelope(caps, [Q(1, 4), Q(1, 2), Q(3, 4), Q(1)])
        shared.verify(caps, mixture)
        with self.assertRaises(ArithmeticError):
            shared.verify(caps, [(c/100, p) for c, p in mixture])

    def test_toy_code_labeled_measure(self):
        # Four rows in the binary code {00,11}. The 15 nonzero message
        # tuples have full union support, regardless of their binary rank.
        words = span([0b11])
        actual = Counter()
        shells = [0, 0, 0]
        for rows in product(words, repeat=4):
            packets = tuple(sum(((row >> j) & 1) << i for i, row in enumerate(rows)) for j in range(2))
            if not any(packets):
                continue
            shells[sum(bool(v) for v in packets)] += 1
            for image, probability in exact_law(packets).items():
                actual[image] += probability
        cdf = [sum(shells[:u+1]) for u in range(3)]
        mixture = shared.envelope(cdf, [Q(1, 2), Q(1)])
        for image, mass in actual.items():
            u = sum(bool(v) for v in image)
            bound = sum((c*p**u*(1-p)**(2-u) for c, p in mixture), Q(0)) / 15**u
            self.assertLessEqual(mass, bound)

    def test_cdf_differences_are_not_shell_bounds(self):
        # Valid CDF caps [0, 4, 4] allow all four words in shell 2.
        # The differences [0, 4, 0] would incorrectly remove that shell.
        caps = [0, 4, 4]
        mixture = shared.envelope(caps, [Q(1, 2)])
        self.assertGreaterEqual(sum(c*p*p for c, p in mixture), 4)

    def test_active_comparison_zero_remains_active(self):
        rows = shared.as_components([(Q(7), Q(1, 4))])
        self.assertEqual(rows[0], ('zero', Q(1), (Q(1), Q(0), Q(0), Q(0), Q(0)), 0))
        self.assertEqual(rows[1][3], 1)
        self.assertEqual(rows[1][2][0], Q(3, 4))
        self.assertEqual(sum(rows[1][2]), 1)

    def test_invalid_majorants(self):
        for caps, centers in (([1, 1], [Q(1, 2)]), ([0, -1], [Q(1, 2)]),
                              ([0, 1], [0]), ([0, 1, 0], [1]),
                              ([0, 1], [Q(1, 2), Q(1, 2)])):
            with self.assertRaises(ValueError):
                shared.envelope(caps, centers)

    def test_empty_active_mass_budget(self):
        caps = [0, 2, 7, 5, 1]
        centers = [Q(1, 4), Q(1, 2), Q(3, 4), Q(255, 256), Q(1)]
        for bits in (0, 2, 8):
            mixture = shared.envelope(caps, centers, zero_bits=bits)
            shared.verify(caps, mixture)
            self.assertTrue(all(c*(1-p)**4 <= Q(2)**-bits for c, p in mixture))

    def test_impossible_empty_budget_rejected(self):
        with self.assertRaises(ValueError):
            shared.envelope([0, 2, 7], [Q(1, 2)], zero_bits=128)

    def test_tilted_assignment_preserves_domination(self):
        for tilt in (Q(1, 8), Q(1, 4), Q(1)):
            caps = [0, 2, 7, 5, 1]
            mixture = shared.envelope(caps, [Q(1, 4), Q(1, 2), Q(3, 4), Q(1)], cost_tilt=tilt)
            shared.verify(caps, mixture)

    def test_component_mass_budget(self):
        caps = [0, 2, 7, 5, 1]
        centers = [Q(1, 4), Q(1, 2), Q(3, 4), Q(255, 256), Q(1)]
        for tilt in (Q(1, 8), Q(1, 4), Q(1)):
            mixture = shared.envelope(caps, centers, zero_bits=2, cost_tilt=tilt, mass_bits=8)
            shared.verify(caps, mixture)
            self.assertTrue(all(c <= 256 and c*(1-p)**4 <= Q(1, 4) for c, p in mixture))
        with self.assertRaises(ValueError):
            shared.envelope(caps, centers, mass_bits=1)
        for invalid in (-1, True, Q(1, 2), 2049):
            with self.assertRaises(ValueError):
                shared.envelope(caps, centers, mass_bits=invalid)

    def test_comparison_empty_shell_preserved(self):
        from monotone_comparison import thinned_shells
        caps = thinned_shells([0, 5, 9], q=16)
        mixture = shared.envelope(caps, [Q(1, 4), Q(1, 2), Q(3, 4), Q(1)], allow_empty=True)
        shared.verify(caps, mixture, allow_empty=True)
        self.assertGreaterEqual(sum(c*(1-p)**2 for c, p in mixture), caps[0])
        self.assertTrue(all(row[3] == 1 for row in shared.as_components(mixture)[1:]))
        with self.assertRaises(ValueError):
            shared.verify(caps, mixture)
        tiny = [Q(1, 1024), Q(1)]
        fine = shared.envelope(tiny, [Q(1024, 1025)], allow_empty=True)
        shared.verify(tiny, fine, allow_empty=True)
        self.assertEqual(sum(c*(1-p) for c, p in fine), tiny[0])


if __name__ == '__main__':
    unittest.main()
