"""Small independent checks; no regional provider or selected screen."""
from itertools import product
import unittest
import numpy as np

from hd_scheduled_column_control import column_word_polynomial, column_branches
from hd_pair_state_tail_general import word_polynomial


class ColumnControlTests(unittest.TestCase):
    def test_single_column_is_exact_interface_overlap(self):
        rng=np.random.default_rng(615)
        logs=rng.uniform(-10,0,(3,2,2))
        counts=np.array([0.,2.,3.]);weights=np.array([0,1,2])
        np.testing.assert_array_equal(column_word_polynomial([logs],counts,weights,3),
                                      word_polynomial([logs],counts,weights,3))

    def test_cartesian_pair_is_no_tighter_coefficientwise(self):
        rng=np.random.default_rng(419)
        for groups in (1,2,4):
            left=rng.uniform(-5,1,(2,2,2));right=rng.uniform(-5,1,(2,2,2))
            counts=np.array([0.,3.]);weights=np.array([0,1])
            paired=[];pc=[];pw=[]
            for i,j in product(range(2),repeat=2):
                paired.append(np.log2(np.exp2(left[i]) @ np.exp2(right[j])))
                pc.append(counts[i]+counts[j]);pw.append(weights[i]+weights[j])
            separate=column_word_polynomial([left,right],counts,weights,groups)
            grouped=word_polynomial([np.array(paired)],pc,np.array(pw),groups)
            self.assertTrue(np.all(separate<=grouped+1e-12))
            if groups==1:np.testing.assert_allclose(separate,grouped,atol=1e-12)

    def test_sixteen_columns_retain_chronology(self):
        a=np.array([[1.,2.],[3.,1.]]);b=np.array([[2.,1.],[1.,4.]])
        columns=[np.log2(x)[None,...] for x in [a,b]*8]
        observed=column_word_polynomial(columns,[0.],np.array([1]),1)
        expected=np.eye(2)
        for x in [a,b]*8:expected=expected@x
        np.testing.assert_allclose(np.exp2(observed[16]),expected,rtol=1e-13)
        self.assertTrue(np.all(np.isneginf(observed[:16])))

    def test_gf4_lambda_preserves_small_support_counts(self):
        for source in ((0,),(2,),(0,7),(1,14)):
            for j in range(16):
                for schedule in ('static','input-xor','input-inverse'):
                    self.assertEqual(column_branches(source,schedule,j,'lambda',lambda_value=2),
                                     column_branches(source,schedule,j,'lambda',lambda_value=0xbc))


if __name__=='__main__':unittest.main()
