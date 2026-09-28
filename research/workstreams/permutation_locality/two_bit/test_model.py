"""Small exhaustive references for the isolated two-bit proof model."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
import unittest

import numpy as np
from flint import arb,arb_mat,ctx

import model
import probe
from second import moments


def words(windows):
    result={s:[] for s in model.shapes(windows)}
    for masks in product(range(4),repeat=windows):
        shape=(sum(m.bit_count()==1 for m in masks),sum(m==3 for m in masks))
        result[shape].append(sum(m<<(2*w) for w,m in enumerate(masks)))
    return result


def linear(columns,x):
    value=0
    for i,c in enumerate(columns):
        if x>>i&1:value^=c
    return value


class ModelTests(unittest.TestCase):
    def setUp(self):ctx.prec=192

    def test_characters_and_inverse_exhaustive(self):
        columns=[1,2,3,5,4,7];reference=words(3)
        shapes,table,denominators,_,inverse=model.characters(columns,3,3)
        for shape,row in zip(shapes,table):
            selected=reference[shape]
            self.assertEqual(len(selected),denominators[shape])
            for a in range(8):
                expected=sum((-1)**((a&linear(columns,x)).bit_count()%2) for x in selected)
                self.assertEqual(expected,int(row[inverse[a]]))
            counts=model.inverse_counts(row[inverse],len(selected))
            expected=Counter(linear(columns,x) for x in selected)
            self.assertEqual(list(map(int,counts)),[expected[a] for a in range(8)])

    def test_output_polynomial_exhaustive(self):
        reference=words(3)
        for tilt in ('0','.17'):
            for y in range(64):
                result=model.output_coefficients(model.pattern_histogram(y,3),3,tilt)
                for shape,selected in reference.items():
                    ctx.prec=384
                    lam=arb(Q(tilt).numerator)/Q(tilt).denominator
                    exact=sum(((-lam*(x^y).bit_count()).exp() for x in selected),arb(0))/len(selected)
                    # The tested quantity is a directed upper endpoint, not a ball.
                    self.assertTrue(result[shape]>=exact.lower())
                    self.assertTrue(abs(result[shape]-exact)<arb(2)**-160)
                    ctx.prec=192

    def test_second_moments_exhaustive(self):
        feedback=[1,2,3,5,4,7];expansion=[0b010101,0b101010,0b110001]
        selected,table,_,_,inverse=model.characters(feedback,3,3)
        actual=moments(expansion,feedback,3,width=2,character_data=(selected,table,inverse))
        for shape,inputs in words(3).items():
            a,b=shape
            if not a+b:continue
            key=(1,)*a+(2,)*b;W=a+2*b
            values=[]
            for target in range(8):
                values.append(Q(sum((2*(x&linear(expansion,target^linear(feedback,x))).bit_count()-W)**2
                                    for x in inputs),len(inputs)))
            D,exact,upper=actual[key]
            self.assertEqual(D,len(inputs));self.assertEqual(exact,values[0])
            self.assertGreaterEqual(upper,max(values))

    def test_mass_column_bounds_for_every_small_target(self):
        expansion=[0b010101,0b101010,0b110001]
        for feedback in ([1,2,3,5,4,7],[0]*6):
            for shape,inputs in words(3).items():
                if not sum(shape):continue
                law=Counter(linear(feedback,x) for x in inputs)
                atom=Q(max(law.values()),len(inputs))
                nonzero_atom=Q(max((n for s,n in law.items() if s),default=0),len(inputs))
                z=Q(3,5)
                tilted={(s,x):z**(x^linear(expansion,s)).bit_count() for s in range(1,8) for x in inputs}
                raw=max(tilted.values())
                point=max(sum(tilted[s,x] for x in inputs)/len(inputs) for s in range(1,8))
                for s in range(1,8):
                    for target in range(8):
                        actual=sum((tilted[s,x] for x in inputs if s^linear(feedback,x)==target),Q(0))/len(inputs)
                        bound=min(point,raw*(nonzero_atom if target==0 else atom))
                        self.assertLessEqual(actual,bound)

    def test_normalized_region_matches_outward_polynomial(self):
        ops=[arb_mat([[1,.125],[0,.75]]),arb_mat([[.5,.25],[.125,.75]]),
             arb_mat([[.25,.5],[.25,.5]])]
        exact=model.region(ops,2);floating=model.float_region(ops,2)
        for r in range(3):
            expected=np.array([[float(exact[r][i,j]) for j in range(2)] for i in range(2)])
            np.testing.assert_allclose(floating[r],expected,rtol=2e-13,atol=1e-25)

    def test_invalid_inputs(self):
        for updates in (0,-1,33,1.5,True):
            with self.assertRaises(ValueError):model.epoch_operators({},'.1',updates=updates)
        for columns,bits,degree in (([1],2,1),([1,4],2,1),([1,2],2,2)):
            with self.assertRaises(ValueError):model.characters(columns,bits,degree)
        for hist,degree,tilt in (((1,1),1,'.1'),((1,-1,1),1,'.1'),((1,1,1),4,'.1'),((1,1,1),1,'-.1')):
            with self.assertRaises(ValueError):model.output_coefficients(hist,degree,tilt)

    def test_dense_float_placement_closed_form(self):
        # T_j=t^j makes every placement of r packets have the same value.
        t=arb(999)/1000
        operators=[arb_mat([[t**j]]) for j in range(65)]
        values=model.float_region(operators,512)
        np.testing.assert_allclose(values[:,0,0],.999**np.arange(513),rtol=2e-12,atol=0)

    def test_tiny_matrix_proposal_normalization(self):
        # Squaring before normalization underflows for this perfectly valid matrix.
        value=model.log_power_moment(np.eye(9)*2.**-900,256,model.terminal())
        self.assertAlmostEqual(value/np.log(2),-900*256,places=8)
        self.assertEqual(model.log_power_moment(np.zeros((9,9))),float('inf'))

    def test_one_group_coefficients(self):
        # Diagonal matrices commute, so each conditional support coefficient
        # has an independent closed form, including both endpoint supports.
        zero=arb_mat([[arb(3)/4 if i==j else 0 for j in range(9)] for i in range(9)])
        active=arb_mat([[arb(1)/4 if i==j else 0 for j in range(9)] for i in range(9)])
        actual=probe.one_group_moments([zero,active])
        ctx.prec=768
        for u,value in enumerate(actual):
            expected=(arb(3)/4)**(256-u)*(arb(1)/4)**u
            self.assertTrue(value>=expected)
            self.assertTrue(abs(value/expected-1)<arb(2)**-170)

    def test_cdf_and_shell_fold_dominate(self):
        # Alternating support costs prevent accidental monotonicity assumptions.
        shells=[Q((u%5)+1) for u in range(257)]
        cdf=[];running=Q(0)
        for u,a in enumerate(shells):
            running+=a;cdf.append(running+u+1)
        caps=[a+1 for a in shells]
        costs=[arb((u%7)+1)/11 for u in range(257)]
        actual=sum((probe.aq(a)*p for a,p in zip(shells,costs)),arb(0))
        self.assertTrue(probe.fold(cdf,caps,costs)>=actual)


if __name__=='__main__':unittest.main()
