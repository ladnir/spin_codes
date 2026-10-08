from fractions import Fraction as F
from itertools import combinations
import unittest
from high_distance_fields import local_matrices as old_local, bound as old_bound, rectangular_cauchy_maps
from high_distance_narrow import local_matrices, bound, profile
from gf8_packet_maps import rank, matmul
from gf8_mds_packet import identity, add


class NarrowTests(unittest.TestCase):
    def test_width_four_agreement(self):
        for bits in (8,16,32):
            for state in (1,2,4):
                for outputs in (1,2,4):
                    self.assertEqual(local_matrices(F(3,4),F(17),4,outputs,state,bits),old_local(F(3,4),F(17),4,outputs,state,bits))
        for sparse,value,point in ((True,16,(5.,2.)),(False,100000,(-2.,-.3,2.,0.))):
            self.assertAlmostEqual(bound(point,value,sparse,width=4,outputs=2,state=4,dimension=16),
                old_bound(point,value,sparse,bits=32,distance=F(2,5),width=4,outputs=2,state=4,dimension=16),places=7)

    def test_narrow_full_cauchy_ranks(self):
        for width,state in ((1,1),(1,2),(1,4),(2,1),(2,2),(2,4)):
            f,a,c,raw,e,ei=rectangular_cauchy_maps(width,state)
            transformed=tuple(f[i]+a[i] for i in range(width))+tuple(c[i]+identity(state)[i] for i in range(state))
            # Normalizing all state rows preserves the rank of selected output
            # rows together with either all state rows or no state rows.
            for nr in range(width+1):
                for rows in combinations(range(width),nr):
                    for include_state in (False,True):
                        selected=rows+tuple(range(width,width+state)) if include_state else rows
                        for nc in range(width+state+1):
                            for cols in combinations(range(width+state),nc):
                                self.assertEqual(rank(tuple(tuple(transformed[i][j] for j in cols) for i in selected)),min(len(selected),len(cols)))
            self.assertEqual(rank(add(f,matmul(a,c))),width)

    def test_width_one_explicit_transfer(self):
        y,t=F(3,5),F(7)
        self.assertEqual(local_matrices(y,t,1,1,1,8),[[[1,0],[0,y]],[[0,y],[y/t,y+(1-y)/t]]])

    def test_geometry_and_endpoint(self):
        g=profile(1,1,64,32)
        self.assertEqual(g['input_n'],2**21)
        self.assertEqual(g['output_n'],2**21)
        g=profile(2,1,64,16)
        self.assertEqual(g['input_n'],2**22)
        self.assertEqual(g['output_n'],2**21)
        with self.assertRaises(ValueError):bound((-2.,-.3,2.,0.),2**21,False)


if __name__=='__main__':unittest.main()
