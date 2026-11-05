from fractions import Fraction as F
from itertools import combinations,product
import unittest

import numpy as np

from hd_blockwise_shared_transfer import multiply,dot,rank,mm
from hd_xor_state_local import (batch_ranks,current_emission,local_counts,
                                marked_kernel,checked_emission,field_tables)


def brute_counts(a,U,S,q):
    n,h=len(a),len(a[0])
    state=[i for i in range(h) if S>>i&1]
    counts=[[0]*(n+1) for _ in range(1<<h)]
    for state_values in product(range(1,q),repeat=len(state)):
        old=[0]*h
        for i,x in zip(state,state_values):old[i]=x
        for values in product(range(q),repeat=len(U)):
            x=[0]*n
            for i,v in zip(U,values):x[i]=v
            y=[x[i]^dot(a[i],old,q) for i in range(n)]
            w=[old[j]^dot(tuple(a[i][j] for i in range(n)),x,q) for j in range(h)]
            target=sum((v!=0)<<j for j,v in enumerate(w))
            counts[target][sum(v!=0 for v in y)]+=1
    return tuple(tuple(row) for row in counts)


class XorStateLocalTests(unittest.TestCase):
    def test_every_tiny_support_matches_direct_enumeration(self):
        a=tuple((1,x) for x in range(4))
        for size in range(3):
            for U in combinations(range(4),size):
                for S in range(4):
                    data=local_counts(U,S,emission=a,field_order=4)
                    self.assertEqual(data['counts'],brute_counts(a,U,S,4))
                    self.assertEqual(sum(sum(row) for row in data['counts']),4**size*3**S.bit_count())
                    z=F(2,3)
                    positive=tuple(sum((F(c,data['denominator'])*z**(4-k)*(1-z)**k
                                         for k,c in enumerate(row)),F(0))
                                   for row in data['positive_zero_coefficients'])
                    self.assertEqual(marked_kernel(data,z),positive)

    def test_native_four_point_butterflies_reproduce_current_basis(self):
        a=current_emission()
        inverse6=next(i for i in range(1,256) if multiply(6,i,256)==1)
        def pair(x,y,c):
            x^=multiply(c,y,256)
            return x,x^y
        for coefficient in range(4):
            outputs=[]
            for offset in (0,4,8,12):
                q=[int(i==coefficient) for i in range(4)]
                c2=multiply(multiply(offset,offset^1,256),inverse6,256)
                q[0],q[2]=pair(q[0],q[2],c2)
                q[1],q[3]=pair(q[1],q[3],c2)
                q[0],q[1]=pair(q[0],q[1],offset)
                q[2],q[3]=pair(q[2],q[3],offset+2)
                outputs.extend(q)
            self.assertEqual(outputs,[row[coefficient] for row in a])
        self.assertEqual(mm(tuple(zip(*a)),a,256),tuple((0,)*4 for _ in range(4)))

    def test_selected_coefficient_columns_are_not_assumed_mds(self):
        a=current_emission()
        self.assertEqual(rank(tuple(tuple(a[i][j] for j in (0,2)) for i in (0,1)),256),1)
        self.assertEqual(checked_emission(a,256),(16,4))
        first=local_counts((),1)
        third=local_counts((),4)
        self.assertEqual(first['counts'][1][16],255)
        self.assertEqual(third['counts'][4][14],255)
        self.assertNotEqual(sum(marked_kernel(first,F(1,2))),sum(marked_kernel(third,F(1,2))))

    def test_batch_ranks_against_independent_elimination(self):
        for q in (4,256):
            matrices=np.array([[[((i+1)*(j+2)+3*k+ij)%q for ij in range(6)]
                                 for j in range(7)] for k,i in enumerate(range(37))],dtype=np.uint8)
            expected=[rank(tuple(tuple(int(v) for v in row) for row in matrix),q) for matrix in matrices]
            np.testing.assert_array_equal(batch_ranks(matrices,q),expected)
        np.testing.assert_array_equal(batch_ranks(np.zeros((3,7,0),dtype=np.uint8),256),0)

    def test_actual_one_input_exact_counts_and_denominator(self):
        a=current_emission()
        data=local_counts((0,),2)
        self.assertEqual(data['counts'],brute_counts(a,(0,),2,256))
        self.assertEqual(data['denominator'],255*256)
        self.assertEqual(sum(marked_kernel(data,F(1))),1)
        self.assertTrue(data['exact_local'])
        self.assertFalse(data['whole_code_certificate'])

    def test_full_state_two_input_return_fiber_by_direct_enumeration(self):
        a=np.array(current_emission(),dtype=np.uint8)
        table,_=field_tables(256)
        values=np.array(tuple(product(range(256),repeat=2)),dtype=np.uint8)
        # W=0 uniquely fixes Q=C_U X. Enumerate just this fiber, not255^4*256^2 points.
        state=table[values[:,0,None],a[0][None,:]]^table[values[:,1,None],a[5][None,:]]
        keep=np.all(state!=0,axis=1)
        state=state[keep];values=values[keep]
        y=np.zeros((len(state),16),dtype=np.uint8)
        for j in range(4):y^=table[state[:,j,None],a[None,:,j]]
        y[:,0]^=values[:,0];y[:,5]^=values[:,1]
        histogram=np.bincount((y!=0).sum(axis=1),minlength=17)
        data=local_counts((0,5),15)
        np.testing.assert_array_equal(data['counts'][0],histogram)
        self.assertEqual(int(histogram.sum()),255**2)
        self.assertEqual(data['denominator'],255**4*256**2)

    def test_budget_and_full_emission_guards(self):
        for U,S in (((0,1,2),1),((0,0),1),((16,),1),((0,),16)):
            with self.assertRaises(ValueError):local_counts(U,S)
        with self.assertRaises(ValueError):local_counts((),1,emission=((1,0),(1,0)),field_order=4)
        with self.assertRaises(ValueError):batch_ranks(np.zeros((8193,1,1),dtype=np.uint8),256)


if __name__=='__main__':unittest.main()
