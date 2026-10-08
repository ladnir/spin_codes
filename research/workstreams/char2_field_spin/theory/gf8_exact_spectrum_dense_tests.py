"""Independent exact spectrum, profile-averaging, and bound checks."""
from fractions import Fraction as F
from itertools import product
from math import comb, prod
import random
import unittest

from gf8_exact_spectrum_dense import (bracket, convex_tilt_floor, convexity_slacks,
    excess_bracket, positive_excess, support_counts, check_ball_exponent)


def direct_mds_supports(q, length, dimension):
    counts = [0]*(1 << length)
    for message in product(range(q), repeat=dimension):
        word = [sum(a*pow(x, j, q) for j, a in enumerate(message)) % q for x in range(length)]
        support = sum(1 << j for j, value in enumerate(word) if value)
        counts[support] += 1
    return counts


def row_average(probabilities, tilt, support_counts_by_mask, q):
    result = F(0)
    for mask, count in enumerate(support_counts_by_mask):
        weight = mask.bit_count()
        probability = prod(a if mask >> j & 1 else 1-a for j, a in enumerate(probabilities))
        result += probability*tilt**weight*F(count, (q-1)**weight)
    return result


def exact_support_weight_law(probabilities):
    """Independent direct Poisson-binomial convolution, no homogenization."""
    result = [F(1)]
    for probability in probabilities:
        following = [F(0)]*(len(result)+1)
        for w, mass in enumerate(result):
            following[w] += mass*(1-probability)
            following[w+1] += mass*probability
        result = following
    return result


class ExactSpectrumTests(unittest.TestCase):
    def test_exhaustive_prime_field_spectra(self):
        for q, length, dimension in ((5, 4, 2), (5, 4, 3), (7, 5, 3)):
            exact = support_counts(length, dimension, q)
            for mask, count in enumerate(direct_mds_supports(q, length, dimension)):
                self.assertEqual(count, exact[mask.bit_count()])

    def test_full_code_cardinality(self):
        for dimension in (1, 2, 67, 74, 80, 84, 94, 127):
            counts = support_counts(128, dimension, 256)
            self.assertEqual(sum(comb(128, w)*c for w, c in enumerate(counts)), 256**dimension)

    def test_exhaustive_unequal_profile_excess_bound(self):
        counts = direct_mds_supports(5, 4, 2)
        grid = (F(0), F(1, 4), F(1, 2), F(3, 4), F(1))
        checked = 0
        for profile in product(grid, repeat=4):
            alpha = sum(profile)/4
            for tilt in (F(1), F(3, 2), F(4)):
                self.assertLessEqual(row_average(profile, tilt, counts, 5),
                                     excess_bracket(alpha, alpha, tilt, 4, 2, 5))
                checked += 1
        self.assertEqual(checked, 1875)

    def test_exhaustive_convex_profile_bound(self):
        counts = direct_mds_supports(5, 4, 2)
        tilt = convex_tilt_floor(4, 2, 5)
        self.assertTrue(all(x >= 0 for x in convexity_slacks(tilt, 4, 2, 5)))
        for profile in product((F(0), F(1, 3), F(2, 3), F(1)), repeat=4):
            alpha = sum(profile)/4
            self.assertLessEqual(row_average(profile, tilt, counts, 5),
                                 bracket(alpha, alpha, tilt, 4, 2, 5))
        with self.assertRaises(ValueError):
            bracket(F(1, 2), F(1, 2), F(1), 4, 2, 5)

    def test_endpoint_uniformity(self):
        for tilt in (F(1), F(5, 4), F(4)):
            upper = excess_bracket(F(3, 4), F(1), tilt, 128, 74, 256)
            for step in range(11):
                alpha = F(3, 4)+F(step, 40)
                self.assertLessEqual(excess_bracket(alpha, alpha, tilt, 128, 74, 256), upper)

    def test_full_support_removes_diagonal_overcount(self):
        parity = 54
        f = F(support_counts(128, 74, 256)[128], 255**128)
        refined = excess_bracket(F(1), F(1), F(1), 128, 74, 256)
        self.assertEqual(refined, max(F(1, 256**parity), f))
        old = F(256, 255)**128/F(256**parity)
        self.assertLess(refined, old)
        self.assertTrue(all(v >= 0 for v in positive_excess(128, 74, 256)))

    def test_unequal_gf256_profiles_by_exact_convolution(self):
        rng = random.Random(20261005)
        grid = (F(0), F(1), F(1, 100), F(99, 100), F(1, 2), F(31, 32))
        profiles = [(F(0),)*128, (F(1),)*128,
                    (F(0),)*16+(F(1),)*112, (F(1, 2),)*64+(F(1),)*64]
        profiles += [tuple(rng.choice(grid) for _ in range(128)) for _ in range(8)]
        counts = support_counts(128, 74, 256)
        checked = 0
        for profile in profiles:
            law = exact_support_weight_law(profile)
            self.assertEqual(sum(law), 1)
            alpha = sum(profile)/128
            for tilt in (F(1), F(5, 4), F(4)):
                direct = sum(mass*count*(tilt/F(255))**w
                             for w, (mass, count) in enumerate(zip(law, counts)))
                self.assertLessEqual(direct, excess_bracket(alpha, alpha, tilt, 128, 74, 256))
                checked += 1
        self.assertEqual(checked, 36)

    def test_exact_input_and_coefficient_checks(self):
        for dimension in (67, 74, 80, 84, 94):
            self.assertTrue(all(c >= 0 for c in support_counts(128, dimension, 256)))
            floor = convex_tilt_floor(128, dimension, 256)
            self.assertTrue(all(c >= 0 for c in convexity_slacks(floor, 128, dimension, 256)))
        with self.assertRaises(ValueError):
            excess_bracket(F(1, 2), F(1), F(-1))
        with self.assertRaises(ValueError):
            support_counts(128, 129, 256)
        self.assertEqual(check_ball_exponent(F(12685, 10000)), F(12685, 10000))
        with self.assertRaises(ValueError):
            check_ball_exponent(F(1268, 1000))


if __name__ == '__main__':
    unittest.main()
