from fractions import Fraction as F
import math
import unittest

from flint import ctx
from gf8_periodic_climb_sparse import periodic_polynomial_power,certify
from gf8_punctured_screen import local_matrices


def exact_product(left,right,degree):
    out=[[[F(0),F(0)],[F(0),F(0)]] for _ in range(degree+1)]
    for a in range(min(degree,len(left)-1)+1):
        for b in range(min(degree-a,len(right)-1)+1):
            for i in range(2):
                for j in range(2):
                    out[a+b][i][j]+=sum(left[a][i][k]*right[b][k][j] for k in range(2))
    return out


class PeriodicSparseTests(unittest.TestCase):
    def test_ordered_polynomial_against_exact_small_product(self):
        y,tau=F(2,3),F(7)
        degree=3
        result=[[[F(1),F(0)],[F(0),F(1)]]]
        for outputs in (6,7,6,7):
            local=local_matrices(y,tau,8,outputs,8)
            weighted=[[[x*math.comb(8,w) for x in row] for row in m]
                      for w,m in enumerate(local)]
            result=exact_product(result,weighted,degree)
        previous=ctx.prec
        try:
            ctx.prec=512
            computed=periodic_polynomial_power(y,tau,(6,7),4,degree)
            for index,(i,j) in enumerate(((0,0),(0,1),(1,0),(1,1))):
                for w in range(degree+1):
                    value=computed[index][w]
                    self.assertLessEqual(F(str(value.lower().fmpq())),result[w][i][j])
                    self.assertGreaterEqual(F(str(value.upper().fmpq())),result[w][i][j])
        finally:
            ctx.prec=previous

    def test_reject_incomplete_period(self):
        with self.assertRaises(ValueError):
            periodic_polynomial_power(F(1,2),F(8),(6,7),3,5)

    def test_low_range84(self):
        receipt=certify(last=16,precision=1024)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertFalse(receipt['whole_code'])
        self.assertEqual(receipt['output_pattern'],[6,6,6,7])
        self.assertEqual(receipt['nominal_rate'],'21/25')
        self.assertEqual(receipt['input_n'],1605632)
        self.assertEqual(receipt['output_n'],1254400)
        self.assertGreaterEqual(receipt['margin_bits_lower'],230)

    def test_low_range74_period_two(self):
        receipt=certify(last=16,dimension=74,pattern=(5,6),precision=1024)
        self.assertTrue(receipt['complete'] and receipt['rigorous'])
        self.assertEqual(receipt['output_pattern'],[5,6])
        self.assertEqual(receipt['nominal_rate'],'37/44')
        self.assertEqual(receipt['input_n'],1818624)
        self.assertEqual(receipt['output_n'],1250304)
        self.assertEqual(receipt['output_threshold'],125030)
        self.assertGreaterEqual(receipt['margin_bits_lower'],312)


if __name__=='__main__':
    unittest.main()
