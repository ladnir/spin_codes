"""Exact small-domain checks of the vectorized joint census."""
import unittest
from itertools import combinations,product
import numpy as np

from joint_four import census,atoms,pair_cache,batch,shape_ids,direct
from group_moment import maps


class JointFour(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs=maps()

    def test_shapes_and_mask_order(self):
        for j in range(1,5):
            shapes,ids=shape_ids(j)
            for index,masks in enumerate(product(range(1,16),repeat=j)):
                self.assertEqual(shapes[int(ids[index])],tuple(sorted(m.bit_count() for m in masks)))

    def test_each_small_domain_matches_scalar_enumeration(self):
        images,columns,spectrum=self.inputs
        _,membership=atoms(images,columns)
        positions=(0,9,17,31)
        for j in range(1,5):
            shapes,_=shape_ids(j); levels=[0]+sorted(spectrum)
            expected=np.zeros((len(shapes),len(levels),129),dtype=np.int64)
            expected_fresh=np.zeros((4,len(shapes),129),dtype=np.int64)
            for windows in combinations(positions,j):
                counts,fresh=direct(windows,images,columns,levels,shapes,membership)
                expected+=counts; expected_fresh+=fresh
            actual,_=census(j,positions,self.inputs,check_samples=False)
            for i,shape in enumerate(shapes):
                rows,fresh,denominator=actual[shape]
                self.assertTrue(np.array_equal([rows[v] for v in levels],expected[i]))
                self.assertTrue(np.array_equal(fresh,expected_fresh[:,i]))

    def test_invalid_domains(self):
        for j,positions in ((0,(0,)),(5,tuple(range(6))),(2,(0,0)),(2,(2,1)),(1,(32,)),(True,(0,))):
            with self.assertRaises(ValueError): census(j,positions,self.inputs)


if __name__=='__main__':
    unittest.main()
