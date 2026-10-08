from fractions import Fraction as F
from itertools import combinations, product
import math
import unittest
import numpy as np
import gf8_fullstate_screen as old
import gf8_fullstate_partial_live as new
from live_kernel_rank_tests import M, INV, dot, rank


class PartialLiveFullStateTests(unittest.TestCase):
    def test_local_entrywise_domination(self):
        for state in (8,16):
            for weight in range(9):
                for y,tau in ((F(9,10),F(6)),(F(3,4),F(255))):
                    before=old.local_matrix(weight,y,tau,state)
                    after=new.local_matrix(weight,y,tau,state)
                    for i in range(2):
                        for j in range(2):
                            self.assertLessEqual(after[i][j],before[i][j])
                            if (i,j)!=(1,0):
                                self.assertEqual(after[i][j],before[i][j])

    def test_aggregate_matches_local(self):
        x,y,tau=.7,.9,6.
        for state in (8,16):
            expected=sum(math.comb(8,w)*x**w*np.asarray(new.local_matrix(w,y,tau,state)) for w in range(9))
            actual=new.packet(math.log(x),math.log(y),math.log(tau),state)
            np.testing.assert_allclose(actual,expected,rtol=3e-13,atol=1e-20)

    def test_exact_gf4_joint_fixed_vector_hypotheses(self):
        # Two inputs/two state symbols permit an exhaustive fixed-vector
        # check, including vectors whose incoming state has zero coordinates.
        width=state=2
        a=[[INV[i^v] for v in (2,3)] for i in range(2)]
        c=[list(row) for row in zip(*a)]
        eta=next(e for e in (1,2,3) if rank([
            [int(i==j)^M[e][dot(a[i],[c[k][j] for k in range(2)])]
             for j in range(2)] for i in range(2)])==2)
        cases=0
        for w in (1,2):
            for support in combinations(range(width),w):
                for mask in range(1<<width):
                    selected=[z for z in range(width) if mask>>z&1]
                    active=len(set(selected)&set(support))
                    inactive=len(selected)-active
                    charge=max(state,active+inactive,w-width+active+2*inactive+1)
                    b=[[c[r][j] for j in support]+[int(r==j) for j in range(state)]
                       for r in range(state)]
                    b += [[int(r==j) for j in support]+[M[eta][v] for v in a[r]]
                          for r in selected]
                    self.assertLessEqual(w+state-rank(b),min(w,width-len(selected)))
                    for inputs in product((1,2,3),repeat=w):
                        for previous in product(range(4),repeat=state):
                            v=inputs+previous
                            hits=0
                            for diagonal in product((1,2,3),repeat=w+state):
                                dv=tuple(M[x][d] for x,d in zip(v,diagonal))
                                hits += all(dot(row,dv)==0 for row in b)
                            self.assertLessEqual(F(hits,3**(w+state)),F(1,3**charge))
                            cases += 1
        self.assertEqual(cases,960)


if __name__=='__main__':
    unittest.main()
