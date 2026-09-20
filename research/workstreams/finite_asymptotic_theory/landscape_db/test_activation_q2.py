"""Exact support enumeration and Q1-axis checks for the native Q2 recurrence."""
import itertools
import math
import unittest
from fractions import Fraction as F

import numpy as np
import activation_q1 as q1
import activation_q2 as q2
from test_activation_q1 import mul, identity


def exact_epochs(z):
    t, s, counts = 8, 4, {4:14, 8:1}
    result = []
    for j in range(3):
        moment = sum((F(c*math.comb(w,a)*math.comb(t-w,j-a),15*math.comb(t,j))*z**(w+j-2*a)
                      for w,c in counts.items() for a in range(j+1) if a<=w and j-a<=t-w),F(0))
        if j == 0:
            result.append([[F(1),F(0),F(0)],[F(0),F(0),z**4],[F(0),F(0),F(15,14)*moment]])
        else:
            result.append([[F(0),z**j,F(0)],[z**(4-j)/15,F(0),z**(4-j)],
                           [moment/14,F(0),F(15,14)*moment]])
    return result


def exact_regions(z, epochs=2):
    kernels = exact_epochs(z)
    result=[]
    for j in range(3):
        total=[[F(0)]*3 for _ in range(3)]
        for positions in itertools.combinations(range(8*epochs),j):
            product=identity(3)
            for e in range(epochs):
                count=sum(e*8<=p<(e+1)*8 for p in positions)
                product=mul(product,kernels[count])
            total=[[a+b for a,b in zip(x,y)] for x,y in zip(total,product)]
        result.append([[x/math.comb(8*epochs,j) for x in row] for row in total])
    return result


class ActivationQ2Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.kernel=q2.PairKernel()

    def test_regions_against_all_placements(self):
        for epochs in (1,2,3):
            z=F(3,4)
            exact=exact_regions(z,epochs)
            got=q2.region_logs(8,4,{4:14,8:1},[-math.log(float(z))],8*epochs)[0]
            np.testing.assert_allclose(np.exp(got),np.array(exact,dtype=float),atol=2e-14,rtol=2e-13)

    def test_pair_dp_against_exact_support_pairs(self):
        exact=exact_regions(F(3,4))
        logs=np.full((3,3,3),-np.inf)
        for j,i,k in itertools.product(range(3),repeat=3):
            if exact[j][i][k]: logs[j,i,k]=math.log(float(exact[j][i][k]))
        got=self.kernel.coefficients(logs,4)
        for a,b in itertools.product(range(5),repeat=2):
            total=F(0)
            for first in itertools.combinations(range(4),a):
                for second in itertools.combinations(range(4),b):
                    product=identity(3)
                    for p in range(4):
                        product=mul(product,exact[int(p in first)+int(p in second)])
                    total+=sum(product[0])
            expected=total/(math.comb(4,a)*math.comb(4,b))
            self.assertAlmostEqual(got[a,b],math.log(float(expected)),places=12)

    def test_q1_axes_and_exchange_symmetry(self):
        lam=np.array([.0001,.03,.7])
        for length in (16,56):
            regions=q2.region_logs(8,4,{4:14,8:1},lam,length)
            rz,ra=q1.region_logs(*q1.epoch_logs(8,4,{4:14,8:1},lam),length//8)
            axis=q1.coefficient_logs(rz,ra,32)
            for i,region in enumerate(regions):
                got=self.kernel.coefficients(region,32)
                np.testing.assert_allclose(got[0],axis[i],atol=2e-12,rtol=2e-12)
                np.testing.assert_allclose(got,got.T,atol=1e-12,rtol=1e-12)

    def test_actual_weight_two_activation_and_following_epochs(self):
        # RM(1,3) state code and GF(16) modulo x^4+x+1.
        generators=[255,170,204,240]
        words=[0]*16
        for state in range(16):
            for i,g in enumerate(generators):
                if (state>>i)&1: words[state]^=g
        def field_mul(a,b):
            out=0
            while b:
                if b&1: out^=a
                b>>=1; a<<=1
                if a&16: a^=19
            return out
        z=F(3,4)
        actual=[]
        for weight in range(3):
            matrix=[[F(0)]*16 for _ in range(16)]
            for selected in itertools.combinations(range(8),weight):
                x=sum(1<<j for j in selected)
                syndrome=sum(((x&g).bit_count()%2)<<i for i,g in enumerate(generators))
                if weight: self.assertNotEqual(syndrome,0)
                for state in range(16):
                    factor=z**((x^words[state]).bit_count())/(15*math.comb(8,weight))
                    for alpha in range(1,16):
                        matrix[state][field_mul(alpha,state)^syndrome]+=factor
            actual.append(matrix)
        envelope=exact_epochs(z)
        starts=[([F(q==j) for q in range(16)],int(j!=0)) for j in range(16)]
        starts.append(([F(0)]+[F(1,15)]*15,2))
        for distribution,kind in starts:
            for first,second in itertools.product(range(3),repeat=2):
                exact=mul(mul([distribution],actual[first]),actual[second])
                bound=mul(mul([[F(i==kind) for i in range(3)]],envelope[first]),envelope[second])
                self.assertLessEqual(sum(exact[0]),sum(bound[0]))


if __name__=='__main__':
    unittest.main()
