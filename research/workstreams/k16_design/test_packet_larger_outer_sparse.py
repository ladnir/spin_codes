import unittest
from fractions import Fraction as Q
from math import comb
from flint import arb, arb_mat, ctx
from packet_larger_outer_sparse import outer, occupancy_bound


class LargerOuterSparseTests(unittest.TestCase):
    def test_geometry(self):
        self.assertEqual((outer(2).message_bits,outer(2).regions),(256,128))
        self.assertEqual((outer(4).message_bits,outer(4).regions),(512,256))
        for bad in (1,3,True):
            with self.assertRaises(ValueError):
                outer(bad)

    def test_scalar_exact(self):
        previous=ctx.prec
        try:
            ctx.prec=192
            regional=[arb_mat([[1]]),arb_mat([[arb(1)/2]]),arb_mat([[arb(1)/4]])]
            value=occupancy_bound(regional,q=2,groups=32,regions=3,beta=Q(5),tilt=Q(1,10),cutoff=7)
            expected=comb(32,2)*25*(arb(7)/10).exp()*(arb(17)/32)**6
            self.assertGreaterEqual(value,expected.lower())
            self.assertLess(value,expected.upper()*(1+arb(2)**-150))
        finally:
            ctx.prec=previous

    def test_continuous_state(self):
        M=arb_mat([[arb(4)/5,arb(1)/10],[arb(1)/50,arb(3)/5]])
        value=occupancy_bound([M,M],q=1,groups=32,regions=3,beta=Q(1),tilt=Q(1,10),cutoff=0)
        direct=32*sum((M**3)[0,j] for j in range(2))
        reset=32*sum(M[0,j] for j in range(2))**3
        self.assertGreaterEqual(value,direct.lower())
        self.assertLess(value,reset.lower())


if __name__=='__main__':
    unittest.main()
