"""Exact checks of the positive-coefficient regional Chernoff identity."""

from fractions import Fraction as F
import math
import unittest

from coupled_chernoff import log_terminal_power, interval_bound
from subfield_bound_audit import (
    ZERO, IDENTITY, matadd, matscale, matmul, matpower, terminal_from_zero,
)
from subfield_certificate import Interval


def product(left, right):
    out = [ZERO]*(len(left)+len(right)-1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j] = matadd(out[i+j], matmul(a, b))
    return out


def power(polynomial, exponent):
    out = [IDENTITY]
    for _ in range(exponent):
        out = product(out, polynomial)
    return out


class CoupledChernoffTest(unittest.TestCase):
    def setUp(self):
        # Deliberately noncommuting positive matrices: no eigenbasis assumption.
        self.local = [((F(1, 2), F(1, 4)), (F(1, 8), F(1, 3))),
                      ((F(1, 3), F(1, 5)), (F(1, 7), F(1, 2))),
                      ((F(1, 6), F(1, 3)), (F(1, 4), F(1, 5)))]

    def test_exact_coefficient_identity(self):
        width, rows, q = 2, 6, 4
        base = [matscale(matrix, math.comb(width, w))
                for w, matrix in enumerate(self.local)]
        regional = power(base, rows//width)
        hpoly = [ZERO]*(width+1)
        for w, matrix in enumerate(self.local):
            for extra in range(width-w+1):
                factor = math.comb(width, w)*q**w*math.comb(width-w, extra)
                hpoly[w+extra] = matadd(hpoly[w+extra], matscale(matrix, factor))
        hpower = power(hpoly, rows//width)
        for a in range(rows+1):
            direct = ZERO
            for h in range(a+1):
                factor = F(math.comb(a, h)*q**h, math.comb(rows, h))
                direct = matadd(direct, matscale(regional[h], factor))
            self.assertEqual(direct, matscale(hpower[a], F(1, math.comb(rows, a))))

    def test_positive_coefficient_majorization(self):
        coefficients = power(self.local, 3)
        x = F(2, 5)
        value = ZERO
        for a, matrix in enumerate(coefficients):
            value = matadd(value, matscale(matrix, x**a))
        for a, matrix in enumerate(coefficients):
            for i in range(2):
                for j in range(2):
                    self.assertLessEqual(matrix[i][j], value[i][j]/x**a)

    def test_interval_bound_matches_exact_formula(self):
        x, z = F(1, 5), F(2, 3)
        width, rows, length, dimension, q, a = 2, 6, 4, 2, 4, 3
        h = ZERO
        for w, matrix in enumerate(self.local):
            factor = math.comb(width, w)*(q*x)**w*(1+x)**(width-w)
            h = matadd(h, matscale(matrix, factor))
        threshold = 5  # ceil((1/4)*6*4)-1
        exact = (terminal_from_zero(matpower(h, rows*length//width))
                 / (math.comb(rows, a)**(length-1)*q**((length-dimension)*a)
                    *(q-1)*z**threshold*x**(length*a)))

        def provider(_z):
            return [[[Interval.rational(entry.numerator, entry.denominator)
                      for entry in row] for row in matrix] for matrix in self.local]

        enclosed = interval_bound(a, x, z, provider, rows=rows, q=q, width=width,
                                  length=length, dimension=dimension, distance=F(1, 4))
        self.assertLessEqual(F(enclosed.lo), exact)
        self.assertGreaterEqual(F(enclosed.hi), exact)

    def test_endpoint_bound(self):
        rows, power_value, geometric = 23, 15, F(3, 7)
        values = [geometric**a/F(math.comb(rows, a)**power_value) for a in range(rows+1)]
        for first in range(rows+1):
            for last in range(first, rows+1):
                self.assertLessEqual(max(values[first:last+1]), max(values[first], values[last]))

    def test_scaled_power(self):
        import numpy as np
        matrix = self.local[0]
        floating = np.array(matrix, dtype=float)
        for exponent in [0, 1, 2, 7, 16]:
            expected = math.log(float(terminal_from_zero(matpower(matrix, exponent))))
            self.assertAlmostEqual(log_terminal_power(floating, exponent), expected, places=12)


if __name__ == '__main__':
    unittest.main()
