"""Small exact checks for the unsigned swap-inheritance lemma.

These test coefficient algebra, not whole-code coverage. Fractional powers
have already been applied in the synthetic one-step coefficient fixture.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, factorial
import unittest

from hd_common_refresh_packets import local_moments, matrix


ZERO = (F(0),) * 4
IDENTITY = (F(1), F(0), F(0), F(1))


def add(*rows):
    return tuple(sum(entries) for entries in zip(*rows))


def scale(row, value):
    return tuple(value * entry for entry in row)


def multiply(left, right):
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def poly_add(*polynomials):
    result = {}
    for polynomial in polynomials:
        for degree, coefficient in polynomial.items():
            result[degree] = add(result.get(degree, ZERO), coefficient)
    return result


def poly_multiply(left, right):
    result = {}
    for (s, f), a in left.items():
        for (t, g), b in right.items():
            degree = s+t, f+g
            result[degree] = add(result.get(degree, ZERO), multiply(a, b))
    return result


def poly_power(polynomial, exponent):
    result = {(0, 0): IDENTITY}
    for _ in range(exponent):
        result = poly_multiply(result, polynomial)
    return result


class SwapInheritanceTests(unittest.TestCase):
    def assert_dominated(self, left, right):
        for degree, coefficient in left.items():
            self.assertTrue(all(a <= b for a, b in zip(coefficient, right[degree])))

    def test_fair_swap_coefficients_and_unsigned_multiplicity(self):
        # States 0,1,2,3 denote empty, first singleton, second singleton, full.
        oriented = {}
        unsigned_arrangements = {}
        for slots in product(range(4), repeat=3):
            u, v, full = slots.count(1), slots.count(2), slots.count(3)
            # Positive asymmetric matrices retain orientation and position.
            coefficient = tuple(F(1+sum((i+1)*(state+1)**(j+1)
                                       for i, state in enumerate(slots)))
                                for j in range(4))
            key = u, v, full
            oriented[key] = add(oriented.get(key, ZERO), coefficient)
            unsigned = tuple(1 if state == 2 else state for state in slots)
            unsigned_arrangements.setdefault(unsigned, []).append(coefficient)
        for singleton in range(4):
            for full in range(4-singleton):
                coefficients = [oriented[(u, singleton-u, full)]
                                for u in range(singleton+1)]
                averaged = scale(add(*coefficients), F(1, 2**singleton))
                maximum = tuple(max(row[j]/comb(singleton, u)
                                    for u, row in enumerate(coefficients))
                                for j in range(4))
                independent = ZERO
                count = 0
                for slots, rows in unsigned_arrangements.items():
                    if slots.count(1) == singleton and slots.count(3) == full:
                        self.assertEqual(len(rows), 2**singleton)
                        independent = add(independent, scale(add(*rows), F(1, len(rows))))
                        count += 1
                self.assertEqual(averaged, independent)
                self.assertTrue(all(a <= b for a, b in zip(averaged, maximum)))
                self.assertEqual(count, factorial(3)//(
                    factorial(3-singleton-full)*factorial(singleton)*factorial(full)))

    def test_constant_caps_preserve_unsigned_normalizer(self):
        cap = (F(2), F(3), F(5), F(7))
        for slots in range(1, 7):
            for singleton in range(slots+1):
                for full in range(slots-singleton+1):
                    count = comb(slots, full)*comb(slots-full, singleton)
                    oriented = [scale(cap, count*comb(singleton, u))
                                for u in range(singleton+1)]
                    self.assertEqual(scale(add(*oriented), F(1, 2**singleton)),
                                     scale(cap, count))

    def test_actual_good_caps_and_ordered_good_good_are_orientation_invariant(self):
        z = F(1, 3)
        local = local_moments(z, q=256, h=4)
        fixtures = []
        for first_full, second_full in ((0, 0), (3, 7), (15, 16), (16, 9)):
            for half in range(2):
                if (first_full, second_full)[half] == 16:
                    continue
                caps = []
                for orientation in range(2):
                    occupancy = [first_full, first_full, second_full, second_full]
                    occupancy[2*half+orientation] += 1
                    caps.append(matrix(tuple(sorted(occupancy)), z, local, q=256, h=4))
                self.assertEqual(caps[0], caps[1])
                fixtures.append(caps)
        # No commutation is used: all four ordered orientation products agree.
        first, second = fixtures[1], fixtures[-1]
        expected = multiply(first[0], second[0])
        self.assertNotEqual(expected, multiply(second[0], first[0]))
        for left, right in product(first, second):
            self.assertEqual(multiply(left, right), expected)

    def test_order_survives_good_dirty_composition_powers_and_positive_aliases(self):
        good = {(0, 0): (F(1), F(2), F(3), F(1)),
                (1, 0): (F(2), F(1), F(4), F(3))}
        dirty_swap = {(2, 0): (F(2), F(3), F(1), F(4)),
                      (2, 1): (F(3), F(1), F(2), F(2))}
        dirty_max = {degree: add(row, (F(1), F(2), F(0), F(3)))
                     for degree, row in dirty_swap.items()}
        self.assertNotEqual(poly_multiply(good, dirty_swap),
                            poly_multiply(dirty_swap, good))
        # The independently established good-good polynomial is common to both.
        composed = {degree: scale(row, F(1, 2))
                    for degree, row in poly_multiply(good, good).items()}
        def paired(dirty):
            return poly_add(composed, poly_multiply(good, dirty),
                            poly_multiply(dirty, good), poly_multiply(dirty, dirty))
        lower, upper = paired(dirty_swap), paired(dirty_max)
        self.assert_dominated(lower, upper)
        lower, upper = poly_power(lower, 3), poly_power(upper, 3)
        self.assert_dominated(lower, upper)
        def positive_alias(polynomial):
            folded = {}
            for (singleton, full), row in polynomial.items():
                residue = singleton % 3, full % 2
                weighted = scale(row, F(2, 5)**singleton*F(3, 7)**full)
                folded[residue] = add(folded.get(residue, ZERO), weighted)
            return folded
        self.assert_dominated(positive_alias(lower), positive_alias(upper))

    def test_one_fixed_orientation_need_not_dominate_swaps(self):
        # Therefore the old oriented receipt cannot be inherited from this lemma.
        first, second = (F(1),)*4, (F(9),)*4
        swapped = scale(add(first, second), F(1, 2))
        self.assertTrue(all(a > b for a, b in zip(swapped, first)))


if __name__ == '__main__':
    unittest.main()
