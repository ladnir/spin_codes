"""Independent exact tests of the Arb sparse-class replay."""

from fractions import Fraction as F
import math
import unittest

from flint import arb,ctx,fmpq

from state8_sparse_certificate import (ball,integer_bits,matrix_product,
    polynomial_power,regional,source_maximum,terminal_power,rational_marks)
from superregular_return_screen import gamma,local_matrix


class SparseReplayTests(unittest.TestCase):
    def assert_contains_rational(self,value,expected):
        # contains(fmpq) first rounds the rational to a ball. That enclosing
        # ball can be wider than value, even when the exact rational lies
        # inside value. Compare exported directed dyadic endpoints instead.
        exact = fmpq(expected.numerator,expected.denominator)
        self.assertLessEqual(value.lower().fmpq(),exact)
        self.assertGreaterEqual(value.upper().fmpq(),exact)

    def setUp(self):
        self.precision = ctx.prec
        ctx.prec = 160

    def tearDown(self):
        ctx.prec = self.precision

    def test_exact_polynomial_coefficients(self):
        # Use a synthetic positive matrix polynomial, with binomial factors
        # explicitly undone so the production input interface is exercised.
        coeff = [[[F(1),F(0)],[F(0),F(1)]],
                 [[F(1,2),F(1,3)],[F(2,3),F(1,5)]],
                 [[F(1,7),F(1,11)],[F(1,13),F(1,17)]]]
        local = [[[value/math.comb(16,w) for value in row] for row in matrix]
                 for w,matrix in enumerate(coeff)]
        local += [[[F(0),F(0)],[F(0),F(0)]] for _ in range(14)]
        actual = [[[F(1),F(0)],[F(0),F(1)]]]
        for _ in range(3):
            following = [[[F(0),F(0)],[F(0),F(0)]] for _ in range(min(5,len(actual)+2))]
            for a,left in enumerate(actual):
                for b,right in enumerate(coeff):
                    if a+b<len(following):
                        for i in range(2):
                            for j in range(2):
                                following[a+b][i][j] += sum(left[i][k]*right[k][j] for k in range(2))
            actual = following
        result = polynomial_power(local,3,4)
        for w,matrix in enumerate(actual):
            for slot,(i,j) in enumerate(((0,0),(0,1),(1,0),(1,1))):
                expected = matrix[i][j]
                self.assert_contains_rational(result[slot][w],expected)

    def test_regional_weight_and_terminal(self):
        local = [local_matrix(w,F(2,3),F(7)) for w in range(17)]
        polys = polynomial_power(local,1,3)
        got = regional(polys,3,16,F(7))
        exact = tuple(sum(F(math.comb(3,w))*7**w*local[w][i][j] for w in range(4))
                      for i,j in ((0,0),(0,1),(1,0),(1,1)))
        for value,expected in zip(got,exact):
            self.assert_contains_rational(value,expected)
        expected = exact
        for _ in range(4):
            expected = matrix_product(expected,exact)
        result = terminal_power(got,5)
        value = expected[0]+expected[1]
        self.assert_contains_rational(result,value)

    def test_source_maximum_and_rational_marks(self):
        for tau in (F(1),F(3,2),F(255)):
            expected = max(F(math.comb(12,a))/tau**(3*a+1) for a in range(2,7))
            result = source_maximum(2,6,12,3,tau)
            self.assertTrue(result.upper()>=ball(expected).lower())
        y,tau = rational_marks(1,(7.268964250212447,math.log(255)),1540096)
        self.assertIsInstance(y,F)
        self.assertEqual(tau,255)

    def test_no_float_cache_contamination(self):
        gamma.cache_clear()
        gamma(255.,8)
        self.assertTrue(all(isinstance(value,F) for value in gamma(F(255),8)))

    def test_integer_margin_uses_definite_inequality(self):
        self.assertEqual(integer_bits(arb(2)**(-119)),119)
        self.assertEqual(integer_bits(ball(F(3,4))*arb(2)**(-119)),119)


if __name__=='__main__':
    unittest.main()
