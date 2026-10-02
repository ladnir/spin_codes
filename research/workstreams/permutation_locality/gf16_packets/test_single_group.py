import unittest
from itertools import product
from math import comb
from fractions import Fraction as Q
from flint import fmpq, fmpq_mat, arb_mat, arb, ctx
import single_group as sg


class ExactSupportTests(unittest.TestCase):
    def test_noncommuting_products(self):
        matrices=[fmpq_mat([[1,2],[0,1]]),fmpq_mat([[1,0],[3,2]])]
        for n in range(7):
            values=sg.support_moments(*matrices,n,matrix=fmpq_mat,scalar=fmpq,rounding=lambda x:x)
            totals=[fmpq(0)]*(n+1)
            for bits in product((0,1),repeat=n):
                row=fmpq_mat([[1,0]])
                for bit in bits:row=row*matrices[bit]
                totals[sum(bits)]+=sum((row[0,k] for k in range(2)),fmpq(0))
            self.assertEqual(values,[value/comb(n,j) for j,value in enumerate(totals)])

    def test_outward_contains_exact_products(self):
        exact=[fmpq_mat([[fmpq(1,3),fmpq(2,7)],[0,1]]),fmpq_mat([[1,0],[fmpq(3,5),fmpq(2,3)]])]
        ctx.prec=192
        floating=[arb_mat(m) for m in exact]
        wanted=sg.support_moments(*exact,8,matrix=fmpq_mat,scalar=fmpq,rounding=lambda x:x)
        actual=sg.support_moments(*floating,8)
        for upper,lower in zip(actual,wanted):self.assertGreaterEqual(upper,arb(lower))

    def test_cdf_fold_uses_tail_majorant(self):
        for shells in product(range(3),repeat=3):
            counts=[0]+[sum(shells[:j])+j for j in range(1,4)]
            for weights in ((0,1,5,2),(0,5,3,1),(0,2,3,6)):
                upper=sg.fold_cdf(counts,list(map(Q,weights)),scalar=Q,rounding=lambda x:x)
                self.assertGreaterEqual(upper,sum(a*b for a,b in zip(shells,weights[1:])))
        self.assertEqual(sg.fold_cdf([0,2,4],[Q(0),Q(1),Q(3)],scalar=Q,rounding=lambda x:x),12)

    def test_invalid_inputs(self):
        with self.assertRaises(ValueError):sg.support_moments(arb_mat(2,3),arb_mat(2,2),1)
        with self.assertRaises(ValueError):sg.support_moments(arb_mat(2,2),arb_mat(2,2),-1)
        with self.assertRaises(ValueError):sg.fold_cdf([1,2],[arb(1),arb(1)])
        with self.assertRaises(ValueError):sg.fold_cdf([0,3,2],[arb(1)]*3)


if __name__=='__main__':unittest.main()
