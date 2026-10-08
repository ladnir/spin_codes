"""Exact checks of the fresh-scalar determinant bound, separate from timing."""
from fractions import Fraction as F
import unittest
import numpy as np
from hd_blockwise_shared_transfer import multiply,mm,rank
from hd_scattered_feedback import maps
from hd_scattered_return_rank import rank_data
from hd_scalar_scatter import return_caps


class ScalarScatterTests(unittest.TestCase):
    def test_all_small_field_input_and_output_zero_sets(self):
        q=4;s=2;a=((1,0),(1,0),(0,1),(0,1))
        c=mm(((3,2),(2,3)),tuple(zip(*a)),q);ac=mm(a,c,q)
        for umask in range(16):
            u=[i for i in range(4) if umask>>i&1]
            rho=rank(tuple(tuple(row[j] for j in u) for row in c),q)
            for zmask in range(16):
                d=[i for i in u if zmask>>i&1]
                e=[i for i in range(4) if zmask>>i&1 and i not in u]
                rest=[i for i in u if i not in d]
                remainder_rank=rank(tuple(tuple(row[j] for j in rest) for row in c),q)
                full_rank=rank(tuple(tuple(ac[i][j] for j in u) for i in e),q)
                residual_rank=rank(tuple(tuple(ac[i][j] for j in rest) for i in e),q)
                generic=F(1,q**(len(d)+residual_rank))-F(1,q**(len(d)+remainder_rank))
                inactive=F(1,q**full_rank)-F(1,q**rho)
                cap=generic+F(min(len(d),s),q-1)*inactive
                actual=F(0);bad=0
                for gamma in range(1,q):
                    p=tuple(tuple(int(i==j)^multiply(gamma,ac[i][j],q) for j in u)
                            for i in d+e)
                    r=rank(p,q)
                    bad+=r<len(d)+residual_rank
                    joint=rank(p+tuple(tuple(row[j] for j in u) for row in c),q)
                    actual+=F(1,q**r)-F(1,q**joint)
                self.assertLessEqual(bad,min(len(d),s))
                self.assertLessEqual(actual/(q-1),cap)

    def test_single_input_average_over_all_nonzero_scalars(self):
        a,c=maps();ac=mm(a,c,65536);occ,_,_,_=rank_data();q=65536
        for z in (.07,.5,.7):
            for order in (256,65536):
                exact=0.
                for j in range(64):
                    weight=sum(bool(ac[i][j]) or i==j for i in range(64))
                    exact+=((order-2)*z**weight+z**(weight-1))/(order-1) if ac[j][j] else z**weight
                exact*=F(q-1,64*q)
                bound=return_caps(z,order)[np.all(occ==(0,0,0,1),axis=1)][0]
                self.assertGreaterEqual(bound,float(exact)*(1-1e-12))


if __name__=='__main__':unittest.main()
