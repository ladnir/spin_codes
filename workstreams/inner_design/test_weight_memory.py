from fractions import Fraction as F
import itertools
import math
import unittest

import numpy as np
import transvection as tv
import weight_memory as wm


def exact_operator(columns,s,j,rounds,z):
    """Enumerated rational state operator, independent of the class envelope."""
    size=1<<s;m=size-1;epsilon=F(1,1<<rounds)
    matrix=[[F(0) for _ in range(size)] for _ in range(size)]
    inputs=[0] if j==0 else [1<<p for p in range(len(columns))]
    for q in range(size):
        emitted=sum(tv.dot(q,col)<<p for p,col in enumerate(columns))
        for x in inputs:
            syndrome=0
            for p,col in enumerate(columns):
                if x>>p&1:syndrome^=col
            weight=z**((emitted^x).bit_count())/len(inputs)
            if q==0:matrix[q][syndrome]+=weight;continue
            matrix[q][q^syndrome]+=epsilon*weight
            for fresh in range(1,size):matrix[q][fresh^syndrome]+=(1-epsilon)*weight/m
    return matrix


class WeightMemoryTests(unittest.TestCase):
    def test_sequence_domination(self):
        for columns,s in (([1,2,3],2),([1,2,3,4,5],3)):
            record=wm.build_neighbors(columns,s)
            weights=wm.all_weights(columns,s)
            levels=record['levels'];spectrum=record['spectrum']
            for rounds,mode in itertools.product((1,2),('forget','neighbors')):
                for z in (F(1,2),F(9,10),F(1)):
                    logs=wm.transfers(record,np.array([-math.log(float(z))]),rounds,mode)
                    envelopes=[np.exp(m[0]) for m in logs]
                    exact=[exact_operator(columns,s,j,rounds,z) for j in (0,1)]
                    initial=[[F(q==0) for q in range(1<<s)]]
                    coordinates=[[1.,0.]+[0.]*len(levels)]
                    for q in range(1,1<<s):
                        initial.append([F(i==q) for i in range(1<<s)])
                        coordinates.append([0.,1.]+[0.]*len(levels))
                    for i,w in enumerate(levels):
                        initial.append([F(int(weights[q])==w,spectrum[w]) for q in range(1<<s)])
                        coordinates.append([0.,0.]+[float(k==i) for k in range(len(levels))])
                    for law,coords in zip(initial,coordinates):
                        for pattern in itertools.product((0,1),repeat=3):
                            actual=law.copy();bound=np.array(coords)
                            for j in pattern:
                                actual=[sum(actual[q]*exact[j][q][y] for q in range(1<<s)) for y in range(1<<s)]
                                bound=bound@envelopes[j]
                                self.assertLessEqual(float(actual[0]),bound[0]+1e-11)
                                caps={w:bound[i+2]/spectrum[w] for i,w in enumerate(levels)}
                                residual=sum(max(0.,float(actual[q])-caps[int(weights[q])]) for q in range(1,1<<s))
                                self.assertLessEqual(residual,bound[1]+1e-11)

    def test_generic_three_state_helpers(self):
        spectrum={2:3};columns=[1,2,3];lam=np.array([.01,.1,.5])
        zero,one=tv.epoch_logs(spectrum,columns,lam,2)
        left=tv.q1.region_logs(zero,one,8)
        rz,ra=wm.regions(zero,one,8);right=(rz,ra-math.log(8))
        for a,b in zip(left,right):np.testing.assert_allclose(a,b,rtol=1e-13,atol=1e-13)
        np.testing.assert_allclose(tv.q1.coefficient_logs(*left,7),wm.coefficients(*right,7),rtol=1e-13,atol=1e-13)


if __name__=='__main__':unittest.main()
