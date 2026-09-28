from fractions import Fraction as Q
from math import comb
import unittest

import model
from flint import fmpq,fmpq_mat,arb,ctx
from occupancy_model import placement
from dense_cover import interval_bound


class DenseTests(unittest.TestCase):
    def test_exact_noncommuting_thinning_identity(self):
        slots,epochs=2,3;G=slots*epochs;a,p=fmpq(2,7),fmpq(3,5)
        ops=[fmpq_mat([[1,1],[0,1]]),fmpq_mat([[1,0],[1,1]]),fmpq_mat([[2,1],[1,3]])]
        regions=placement(ops,epochs=epochs,windows=slots,matrix=fmpq_mat,rounding=lambda x:x,maximum_groups=G)
        total=fmpq_mat(2,2);components=[]
        for q in range(G+1):
            regional=sum((regions[r]*comb(q,r)*p**r*(1-p)**(q-r) for r in range(q+1)),fmpq_mat(2,2))
            beta=comb(G,q)*a**q*(1-a)**(G-q)
            components.append(regional*beta);total+=components[-1]
        local=sum((ops[j]*comb(slots,j)*(a*p)**j*(1-a*p)**(slots-j) for j in range(slots+1)),fmpq_mat(2,2))
        self.assertEqual(total,local**epochs)
        for part in components:
            self.assertTrue(all(part[i,j]<=total[i,j] for i in range(2) for j in range(2)))

    def test_frozen_witness_log_convexity(self):
        for G in (6,12):
            for regions in (1,3,7):
                a,outer=Q(2,7),Q(11,3)
                values=[Q(comb(G,q))*outer**q/(comb(G,q)*a**q*(1-a)**(G-q))**regions for q in range(G+1)]
                for q in range(1,G):self.assertLessEqual(values[q]**2,values[q-1]*values[q+1])
                for lo in range(G+1):
                    for hi in range(lo,G+1):
                        self.assertLessEqual(max(values[lo:hi+1]),max(values[lo],values[hi]))

    def test_interval_includes_its_length(self):
        ctx.prec=192
        values={q:arb(2)**(-60+(q-4)**2) for q in range(2,7)}
        bound=interval_bound(values.__getitem__,2,6)
        self.assertTrue(bound>=sum(values.values(),arb(0)))
        self.assertEqual(bound,5*max(values[2],values[6]))


if __name__=='__main__':unittest.main()
