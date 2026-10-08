"""Exact small-field checks of the global-emission return cap."""
from fractions import Fraction as F
from itertools import product
from math import comb
import unittest
from hd_blockwise_shared_transfer import dot,mm,rank
from hd_common_emission_alternate_screen import check_global_return_invertible


class CommonGlobalEmissionTests(unittest.TestCase):
    def test_concrete_return_is_invertible(self):
        check_global_return_invertible()

    def test_subspace_bound_on_all_small_input_supports(self):
        q=4;s=2;n=4;a=tuple((1,x) for x in range(q));c=tuple(zip(*a))
        self.assertFalse(any(x for row in mm(c,a,q) for x in row))
        for mask in range(16):
            u=[i for i in range(n) if mask>>i&1];j=len(u)
            rho=rank(tuple(tuple(row[i] for i in u) for row in c),q)
            for z in (F(1,2),F(3,4)):
                actual=F(0);quiet=F(0)
                for values in product(range(q),repeat=j):
                    x=[0]*n
                    for i,value in zip(u,values):x[i]=value
                    r=tuple(dot(row,x,q) for row in c)
                    if any(r):actual+=z**sum(x[i]^dot(row,r,q)!=0 for i,row in enumerate(a))
                    else:quiet+=z**sum(value!=0 for value in x)
                actual/=q**j;quiet/=q**j
                inactive=sum(F(comb(n-j,e))*z**(n-j-e)*(1-z)**e*
                             (F(1,q**max(0,rho+e-s))-F(1,q**rho))
                             for e in range(min(s,n-j+1)))
                cap=min(inactive*(z+(1-z)/q)**(j-rho),
                        (1-F(1,q**rho))*z**max(0,n-s+1-j),
                        (z+(1-z)/q)**j-quiet)
                self.assertLessEqual(actual,cap)


if __name__=='__main__':unittest.main()
