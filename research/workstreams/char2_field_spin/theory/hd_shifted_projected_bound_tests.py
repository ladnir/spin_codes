"""Projected histogram normalization and chronology, independently checked."""
from fractions import Fraction as F
from itertools import product
import unittest

from flint import arb,ctx

from hd_shifted_alias_certificate import rational
from hd_shifted_projected_bound import ProjectedBound
from hd_shifted_unoriented_certificate import matrix_product
from hd_sumrank_paired_certificate import scalar_power


class Region:
    def __init__(self,parity):
        self.parity = parity
    def numerator(self,ranges):
        (s,sh),(f,fh) = ranges
        assert s==sh and f==fh
        return tuple(arb(1+s*(i+1)+f*(4-i)+self.parity*(i*i+1)) for i in range(4))


class Witness:
    precision,grid,uniform_aliases = 256,128,False
    construction = 'synthetic-algebra-test-only'
    def __init__(self):
        with ctx.workprec(self.precision):
            self.slope,self.constant = arb(3),arb(-7)
            self.regions = Region(0),Region(1)
            self.log_markers = tuple(tuple(rational(F(m)).log()/arb(2).log() for m in row)
                                     for row in ((2,3),(5,7)))
    def log_factorial(self,n):
        return (arb(n)+1).lgamma()/arb(2).log()
    def witness(self):
        return dict(test_only=True)


class ProjectedTests(unittest.TestCase):
    def test_normalization_and_distributed_tilt(self):
        w = Witness()
        p = ProjectedBound(w,F(2,3),F(3,5))
        r,a,b = 6,1,2
        with ctx.workprec(256):
            normalized = []
            for parity,f in enumerate((a,b)):
                hist = 4096-r+f,r-2*f,f
                log_m = w.log_factorial(4096)-sum((w.log_factorial(n) for n in hist),arb(0))
                log_m += sum((n*x for n,x in zip(hist[1:],w.log_markers[parity])),arb(0))
                normalized.append(tuple(v/(arb(2)**log_m)
                                        for v in w.regions[parity].numerator(((hist[1],hist[1]),(f,f)))))
            direct = scalar_power(matrix_product(*normalized),128)
            actual = scalar_power(matrix_product(p.regional(0,r,a),p.regional(1,r,b)),128)
            tilt = arb(2)**(w.slope*r)/rational(p.x**r*p.y**(r-a-b))
            for observed,expected in zip(actual,direct):
                expected *= tilt
                self.assertGreaterEqual(observed.upper(),expected.lower())
                self.assertLessEqual(observed.upper(),expected.upper()*(1+arb(2)**-220))

    def test_envelopes_dominate_every_projected_histogram(self):
        p = ProjectedBound(Witness(),F(2,3),F(3,5))
        ranges = ((4,8),(0,3),(1,3))
        with ctx.workprec(256):
            envelopes,counts = p.envelopes(ranges)
            for parity in (0,1):
                count = 0
                for r,f in product(range(4,9),range(ranges[parity+1][0],4)):
                    if f<=r//2:
                        count += 1
                        for v,cap in zip(p.regional(parity,r,f),envelopes[parity]):
                            self.assertLessEqual(v,cap)
                self.assertEqual(counts[parity],count)
            whole = scalar_power(matrix_product(*envelopes),128)
            for r,a,b in product(range(4,9),range(4),range(1,4)):
                if max(a,b)<=r//2:
                    point = scalar_power(matrix_product(p.regional(0,r,a),p.regional(1,r,b)),128)
                    for value,cap in zip(point,whole):
                        self.assertLessEqual(value.upper(),cap.upper()*(1+arb(2)**-220))

    def test_receipt_never_claims_complete_and_retains_identity(self):
        p = ProjectedBound(Witness(),F(2,3),F(3,5))
        receipt = p.receipt(((4,5),(0,1),(0,1)))
        self.assertTrue(receipt['directed'])
        self.assertFalse(receipt['whole_code_certificate'])
        self.assertEqual(receipt['construction'],'synthetic-algebra-test-only')
        self.assertEqual(receipt['support_count_method'],'projected-normalization-with-cyclic-run-tilt')

    def test_invalid_markers_histograms_and_empty_geometry(self):
        with self.assertRaises(TypeError):
            ProjectedBound(Witness(),.5,F(1,2))
        p = ProjectedBound(Witness(),F(1,2),F(1,2))
        with self.assertRaises(ValueError):
            p.regional(0,1,1)
        with self.assertRaises(ValueError):
            p.receipt(((4,4),(2,2),(2,2)))


if __name__=='__main__':
    unittest.main()
