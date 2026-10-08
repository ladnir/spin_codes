from fractions import Fraction as F
import random
import unittest

from flint import arb,arb_poly,ctx

from hd_sumrank_paired_certificate import round_up,grouped_packets,pair_coefficients,poly_product,scalar_power


class PairedCertificateTests(unittest.TestCase):
    def test_directed_dyadic_rounding(self):
        rng=random.Random(53192)
        for _ in range(400):
            x=F(rng.randrange(1,1<<200),rng.randrange(1,1<<250))
            y=round_up(x,12)
            self.assertGreaterEqual(y,x)
            self.assertLessEqual(y,x*(1+F(1,2048)))
            self.assertEqual(y.denominator&(y.denominator-1),0)
        self.assertEqual(round_up(F(0)),0)

    def test_polynomial_matrix_product(self):
        with ctx.workprec(128):
            a=tuple(arb_poly(list(row)) for row in ((1,2),(3,4),(0,1),(2,1)))
            b=tuple(arb_poly(list(row)) for row in ((0,2),(1,0),(3,1),(4,5)))
            actual=poly_product(a,b,2)
            for r in range(2):
                for s in range(2):
                    for j in range(3):
                        expected=sum(a[2*r+h][i]*b[2*h+s][j-i] for h in range(2) for i in range(j+1))
                        self.assertTrue(actual[2*r+s][j].contains(expected))

    def test_pairs_at_integer_exponent(self):
        with ctx.workprec(128):
            groups=grouped_packets(F(3,4),state=4,cutoff=2,maximum=3)
            local=[arb_poly([]) for _ in range(4)]
            for j,count,entries in groups:
                for h,x in enumerate(entries):
                    local[h]+=arb_poly([0]*j+[count*arb(x.numerator)/x.denominator])
            expected=poly_product(tuple(local),tuple(local),6)
            actual,_=pair_coefficients(F(3,4),F(1),state=4,cutoff=2,maximum=3)
            for h in range(4):
                for j in range(7):
                    self.assertTrue(actual[h][j].overlaps(expected[h][j]))
            self.assertEqual(actual[0][0],1)
            self.assertEqual(actual[1][0],0)
            self.assertEqual(actual[2][0],0)
            self.assertTrue(actual[3][0] > 0)
            self.assertTrue(actual[3][0] < 1)

    def test_scalar_power_matches_repetition(self):
        matrix=(arb(1),arb(2),arb(3),arb(4))
        for n in range(9):
            expected=(arb(1),arb(0),arb(0),arb(1))
            for _ in range(n):
                expected=tuple(expected[2*r]*matrix[s]+expected[2*r+1]*matrix[2+s] for r in range(2) for s in range(2))
            self.assertEqual(scalar_power(matrix,n),expected)


if __name__=='__main__':
    unittest.main()
