import unittest
from math import comb,log2
import numpy as np
from hd_column_deficit_tail import column_deficit_tail
from hd_pair_state_tail_general import entrywise_power,matrix_polynomial_product,tail_from_polynomial


class DeficitTests(unittest.TestCase):
    def test_agrees_with_untruncated_activity_tail(self):
        rng=np.random.default_rng(917)
        columns=[np.log2(rng.uniform(.01,2,(3,2,2))) for _ in range(16)]
        word=np.full((1,2,2),-np.inf)
        word[0,0,0]=word[0,1,1]=0.
        for column in columns:
            word=matrix_polynomial_product(word,entrywise_power(column,3))
        expected=tail_from_polynomial(word,9*2*3)
        self.assertAlmostEqual(float(column_deficit_tail(columns,2,3)),float(expected),places=10)

    def test_arity5_exceeds_old_total_degree_but_truncates_safely(self):
        a,r=5,18
        column=np.full((a+1,2,2),-np.inf)
        for w in range(a+1):
            column[w,0,0]=column[w,1,1]=log2(comb(a,w))+w*log2(255)
        actual=float(column_deficit_tail([column]*16,a,r))
        total=16*a*r
        expected=log2(sum(comb(total,d)*255**(total-d) for d in range(7*a*r+1)))
        self.assertAlmostEqual(actual,expected,places=9)

    def test_bad_shape_and_capacity_rejected(self):
        columns=[np.zeros((3,2,2))]*16
        for bad in (0,9):
            with self.assertRaises(ValueError): column_deficit_tail(columns,bad,18)
        with self.assertRaises(ValueError): column_deficit_tail(columns,3,18)


if __name__=='__main__': unittest.main()
