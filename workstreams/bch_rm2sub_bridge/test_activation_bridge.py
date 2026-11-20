"""Exact finite-field toy checks and independent region/support checks."""
import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction as F

import numpy as np
import bridge as base
import activation_bridge as model


def word(q,rows):
    value=0
    for i,row in enumerate(rows):
        if (q>>i)&1:
            value^=row
    return value


def field_mul(a,b):
    result=0
    while b:
        if b&1:
            result^=a
        b>>=1; a<<=1
        if a&16:
            a^=0x13
    return result


def exact_kernels(rows,z):
    t,s=8,4
    columns=[sum(((r>>j)&1)<<i for i,r in enumerate(rows)) for j in range(t)]
    kernels=[]
    for active in (False,True):
        mat=[[F(0) for _ in range(1<<s)] for _ in range(1<<s)]
        inputs=[(1<<j,columns[j]) for j in range(t)] if active else [(0,0)]
        for q in range(1<<s):
            for x,b in inputs:
                weighted=z**(x^word(q,rows)).bit_count()/((1<<s)-1)/len(inputs)
                for alpha in range(1,1<<s):
                    mat[q][field_mul(alpha,q)^b]+=weighted
        kernels.append(mat)
    return kernels


class BridgeTests(unittest.TestCase):
    def test_toy_exact_paths_and_all_entering_classes(self):
        # RM(1,3): self-orthogonal A, B=A^T, distinct nonzero B columns.
        rows=[0xff,0xaa,0xcc,0xf0]
        spectrum=Counter(word(q,rows).bit_count() for q in range(1,16))
        checked=0
        for z in (F(1,3),F(3,4),F(99,100)):
            kernels=exact_kernels(rows,z)
            envelope=model.matrices(8,4,spectrum,z,F)
            # All point masses in D, zero, uniform-live, all punctured laws.
            starts=[([F(int(q==0)) for q in range(16)],0)]
            starts += [([F(int(q==j)) for q in range(16)],1) for j in range(1,16)]
            starts += [([F(0)]+[F(1,15)]*15,2)]
            starts += [([F(0)]+[F(int(q!=omit),14) for q in range(1,16)],2) for omit in range(1,16)]
            for exact,kind in starts:
                approx=tuple(F(int(j==kind)) for j in range(3))
                def visit(v,u,depth):
                    nonlocal checked
                    self.assertLessEqual(sum(v),sum(u)); checked+=1
                    if depth==0:
                        return
                    for a in (0,1):
                        nxt=[sum((v[q]*kernels[a][q][j] for q in range(16)),F(0)) for j in range(16)]
                        visit(nxt,model.positive_vector_mul(u,envelope[a]),depth-1)
                visit(exact,approx,4)
        self.assertEqual(checked,2976)

    def test_region_position_and_support_averaging(self):
        t,s,spectrum=8,4,{4:14,8:1}; z=F(3,4); epochs=3; length=5
        zero,active=model.matrices(t,s,spectrum,z,F)
        rz,ra=model.positive_regions(t,s,spectrum,z,F,epochs)
        identity=tuple(F(int(i==j)) for i in range(3) for j in range(3))
        summed=[F(0)]*9
        for position in range(epochs):
            product=identity
            for j in range(epochs):
                product=model.positive_mul(product,active if j==position else zero)
            summed=[x+y/epochs for x,y in zip(summed,product)]
        self.assertEqual(tuple(summed),ra)
        coefficients=model.positive_coefficients(rz,ra,F,length)
        for weight in range(length+1):
            total=F(0)
            for support in itertools.combinations(range(length),weight):
                current=(F(1),F(0),F(0))
                for j in range(length):
                    current=model.positive_vector_mul(current,ra if j in support else rz)
                total+=sum(current)
            self.assertEqual(total/math.comb(length,weight),coefficients[weight])
        logs=model.log_coefficients(np.array([[math.log(float(x)) if x else -np.inf for x in rz[i:i+3]] for i in (0,3,6)]),
                                    np.array([[math.log(float(x)) if x else -np.inf for x in ra[i:i+3]] for i in (0,3,6)]),length)
        for actual,expected in zip(logs,coefficients):
            self.assertAlmostEqual(actual,math.log(float(expected)),places=12)

    def test_actual_maps_and_log_region_implementation(self):
        for name in base.CONFIGS:
            t,s,spectrum=base.load_map(name)
            lam=0.001
            positive=model.positive_regions(t,s,spectrum,math.exp(-lam),float,base.ROWS//t)
            logs=model.log_regions(t,s,spectrum,lam)
            for a,b in zip(positive,logs):
                np.testing.assert_allclose(np.array(a).reshape(3,3),np.exp(b),rtol=2e-11,atol=1e-15)


def activation_diagnostics():
    rows=[]
    for name in base.CONFIGS:
        t,s,spectrum=base.load_map(name)
        selected=base.read(base.HERE/'inputs'/f'{name}_selection.json')['selected']
        generators=[int(a,16) for a in selected['A_generator_words_hex']]
        hist=Counter(word(int(b,16),generators).bit_count() for b in selected['B_columns_hex'])
        ratios={}
        for lam in (0.0001,0.001,0.01,0.1):
            exact=sum(c*math.exp(-lam*w) for w,c in hist.items())/t
            old=sum(c*math.exp(-lam*w) for w,c in spectrum.items())/((1<<s)-2)
            ratios[str(lam)]=exact/old
        rows.append(dict(configuration=name,activation_output_spectrum=dict(sorted(hist.items())),
                         two_epoch_actual_over_old_bound=ratios))
    print(rows)


if __name__=='__main__':
    activation_diagnostics()
    unittest.main()
