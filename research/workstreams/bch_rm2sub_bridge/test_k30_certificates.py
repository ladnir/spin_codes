"""Small exact checks for the large-length certificate kernels."""
from fractions import Fraction as F
import itertools
import math
import unittest

from flint import arb,ctx
import certify_refresh_q1 as matrix
from certify_k30_q1 import polynomial_region
import k30_sparse as sparse
import k30_kernel as general
from k30_dense_screen import upper_hull


class K30Test(unittest.TestCase):
    def test_shuffled_poisson_binomial_density_bound(self):
        for length in range(1,6):
            for ps in itertools.product((F(0),F(1,4),F(3,4),F(1)),repeat=length):
                distribution=[F(1)]
                for p in ps:
                    updated=[F(0)]*(len(distribution)+1)
                    for j,value in enumerate(distribution):
                        updated[j]+=value*(1-p);updated[j+1]+=value*p
                    distribution=updated
                r=sum(ps,F(0))/length
                for j,value in enumerate(distribution):
                    reference=math.comb(length,j)*r**j*(1-r)**(length-j)
                    self.assertLessEqual(value,(length+1)*reference)

    def test_upper_cost_hull(self):
        points=[(0.,0.),(.2,.1),(.4,.8),(.7,.2),(1.,0.)]
        hull=upper_hull(points)
        self.assertEqual(hull,[(0.,0.),(.4,.8),(1.,0.)])
        for x,y in points:
            left,right=next((a,b) for a,b in zip(hull,hull[1:]) if a[0]<=x<=b[0])
            bound=left[1]+(right[1]-left[1])*(x-left[0])/(right[0]-left[0])
            self.assertLessEqual(y,bound+1e-15)

    def test_kernel_split_for_arbitrary_uniform_and_missing_one_states(self):
        generators=(255,170,204,240)
        words=[0]
        for g in generators:words += [w^g for w in words]
        z=F(2,3);m=15
        entries=general.epochs(8,4,{4:14,8:1},{0:1,4:14,8:1},z,8,F)
        for j,entry in enumerate(entries):
            supports=[sum(1<<i for i in v) for v in itertools.combinations(range(8),j)]
            classes=[[u for u in supports if any((u&g).bit_count()%2 for g in generators)==b]
                     for b in (False,True)]
            beta=F(len(classes[0]),len(supports))
            self.assertEqual(entry[0],beta*z**j)
            self.assertEqual(entry[1],(1-beta)*z**j)
            for kind,inputs in enumerate(classes):
                moments=[sum((z**((u^a).bit_count()) for u in inputs),F(0))/len(supports) for a in words[1:]]
                for state,actual in ((1,max(moments)),(2,sum(moments,F(0))/m),
                                     (3,max((sum(moments,F(0))-v)/(m-1) for v in moments))):
                    if kind==0:
                        self.assertGreaterEqual(entry[4*state+2],actual)
                    else:
                        self.assertGreaterEqual(entry[4*state],actual/m)
                        self.assertGreaterEqual(entry[4*state+3],actual*(m-1)/m)

    def test_general_transfer_agrees_below_kernel_distance(self):
        ctx.prec=256
        for maximum in (1,2,3):
            a=sparse.epochs(8,4,{4:14,8:1},{0:1,4:14,8:1},F(2,3),maximum,F)
            b=general.epochs(8,4,{4:14,8:1},{0:1,4:14,8:1},F(2,3),maximum,F)
            # The general entry can tighten L by the arbitrary-state bound.
            self.assertTrue(all(y<=x for ma,mb in zip(a,b) for x,y in zip(ma,mb)))

    def test_q1_polynomial_power_encloses_exact_region(self):
        ctx.prec=256
        for count in (1,2,3,17):
            exact=matrix.region(*matrix.epoch(4,2,{2:2,4:1},F(2,3),F),count,F,linear=True)
            got=polynomial_region(*matrix.epoch(4,2,{2:2,4:1},arb(2)/3,arb),count)
            for wanted,actual in zip(exact,got):
                for x,y in zip(wanted,actual):
                    self.assertLessEqual(sparse.rational(y.lower()),x)
                    self.assertGreaterEqual(sparse.rational(y.upper()),x)

    def test_multibit_moments_and_kernel_exclusion(self):
        generators=(255,170,204,240)
        words=[0]
        for g in generators:words += [w^g for w in words]
        spectrum={w:sum(x.bit_count()==w for x in words[1:]) for w in (4,8)}
        kernel={0:1,4:14,8:1}
        z=F(2,3);m=15
        for j,entry in enumerate(sparse.epochs(8,4,spectrum,kernel,z,3,F)):
            if j==0:continue
            supports=[sum(1<<i for i in support) for support in itertools.combinations(range(8),j)]
            # A^T u is nonzero for every tested input; this premise is essential.
            self.assertTrue(all(any((u&g).bit_count()%2 for g in generators) for u in supports))
            moments=[sum((z**((u^x).bit_count()) for u in supports),F(0))/len(supports) for x in words[1:]]
            self.assertEqual(entry[8]+entry[11],sum(moments,F(0))/m)
            self.assertGreaterEqual(entry[4]+entry[7],max(moments))
            self.assertEqual(entry[1],z**j)
        with self.assertRaises(AssertionError):
            sparse.epochs(8,4,spectrum,kernel,z,4,F)

    def test_region_polynomial_against_exact_epoch_enumeration(self):
        ctx.prec=256
        spectrum={4:14,8:1};kernel={0:1,4:14,8:1};maximum=3
        epoch=sparse.epochs(8,4,spectrum,kernel,F(2,3),maximum,F)
        want=[(F(int(i==j)) for i in range(4) for j in range(4))]
        want=[tuple(want[0])]
        for step in range(2):
            updated=[(F(0),)*16 for _ in range(maximum+1)]
            for i,a in enumerate(want):
                for j,b in enumerate(epoch[:maximum+1-i]):
                    updated[i+j]=matrix.add(updated[i+j],tuple(x*math.comb(8,j) for x in matrix.product(a,b)))
            want=updated
        got=sparse.regions(8,4,spectrum,kernel,arb(2)/3,maximum,length=16)
        for degree,row in enumerate(want):
            for value,upper in zip(row,got[degree]):
                exact=value/math.comb(16,degree)
                self.assertGreaterEqual(sparse.rational(upper),exact)


if __name__=='__main__':unittest.main()
