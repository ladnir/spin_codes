from fractions import Fraction as F
from itertools import combinations,product
import unittest

import numpy as np

from hd_rankone_mask_screen import (packet_matrix,mask_matrices,pair_coefficients,
                                    log_product,log_power,exact_sparse_regions,project_flat_to_masks)
from hd_rankone_common_screen import packet_matrix as three_state_matrix
import hd_rankone_common_screen_tests as common_tests


class RankOneMaskTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        common_tests.RankOneCommonTests.setUpClass()
        cls.states=common_tests.RankOneCommonTests.states
        cls.index=common_tests.RankOneCommonTests.index
        cls.kernel=common_tests.RankOneCommonTests.kernel

    def test_exact_toy_measures_fit_every_mask_and_the_flat_row(self):
        z=F(2,3)
        masks=np.array([sum((value!=0)<<i for i,value in enumerate(state)) for state in self.states])
        for support in ((),(0,),(0,1),(0,4),(0,1,4)):
            occ=(sum(i<4 for i in support),sum(i>=4 for i in support))
            cap=packet_matrix(occ,float(z),q=4,h=1)
            measures=[]
            for old in self.states:
                measure=[F(0) for _ in self.states]
                for values in product(range(4),repeat=len(support)):
                    x=[0]*8
                    for i,value in zip(support,values):x[i]=value
                    syndrome=[0,0]
                    for i,value in enumerate(x):syndrome[i//4]^=value
                    before=tuple(a^b for a,b in zip(old,syndrome))
                    weight=sum((value^old[i//4])!=0 for i,value in enumerate(x))
                    scale=z**weight/F(4**len(support))
                    for j,p in enumerate(self.kernel[self.index[before]]):measure[j]+=scale*p
                measures.append(np.array([float(v) for v in measure]))
            cases=list(zip(measures,masks))+[(sum(measures[1:])/15,4)]
            for measure,row in cases:
                self.assertLessEqual(measure[0],cap[row,0]+2e-14)
                for target in range(1,4):
                    excess=np.maximum(0.,measure[masks==target]-cap[row,4]/15).sum()
                    self.assertLessEqual(excess,cap[row,target]+2e-14)
            projected=project_flat_to_masks(cap,q=4,h=1)
            for measure,row in zip(measures,masks):
                for target in range(4):
                    self.assertLessEqual(measure[masks==target].sum(),projected[row,target]+2e-14)

    def test_flat_zero_return_agrees_with_existing_return_bound(self):
        for occ in ((0,0,0,0),(0,0,0,1),(0,2,5,16),(16,16,16,16)):
            matrix=packet_matrix(occ,.73)
            old=three_state_matrix(occ,.73)
            self.assertAlmostEqual(matrix[16,0],old[2,0],delta=2e-14*max(old[2,0],1e-150))
            self.assertEqual(matrix.shape,(17,17))

    def test_label_permutation_conjugates_state_masks(self):
        occ=(0,1,3,7)
        permutation=(2,0,3,1)
        old=packet_matrix(occ,.7)
        new=packet_matrix(tuple(occ[i] for i in permutation),.7)
        remap=[sum(((mask>>i)&1)<<permutation[i] for i in range(4)) for mask in range(16)]+[16]
        np.testing.assert_allclose(new,old[np.ix_(remap,remap)],rtol=1e-14,atol=0)
        self.assertNotEqual(packet_matrix((1,0,0,0),.7)[1,0],
                            packet_matrix((0,1,0,0),.7)[1,0])

    def test_pair_average_matches_every_labeled_support_of_weight_at_most_two(self):
        z=.87
        coefficients=pair_coefficients(z,1.,2)
        expected={key:np.zeros((17,17)) for key in coefficients}
        partner={j:j+16 for j in range(16)}|{j+16:j for j in range(16)}
        partner|={32+j:48+j for j in range(16)}|{48+j:32+j for j in range(16)}
        profiles={}
        for degree in range(3):
            for support in combinations(range(64),degree):
                full=sum(partner[j] in support for j in support)//2
                single=degree-2*full
                occ=tuple(sum(j//16==block for j in support) for block in range(4))
                if occ not in profiles:profiles[occ]=packet_matrix(occ,z)
                expected[single,full]+=profiles[occ]/2**single
        for key,logs in coefficients.items():
            np.testing.assert_allclose(np.exp2(logs),expected[key],rtol=1e-12,atol=0)

    def test_noncommuting_product_and_complete_sparse_extraction(self):
        a=np.arange(1,26,dtype=float).reshape(5,5)/40
        b=np.flip(a,axis=0)+.1
        np.testing.assert_allclose(np.exp2(log_product(np.log2(a),np.log2(b))),a@b)
        np.testing.assert_allclose(np.exp2(log_power(np.log2(a),5)),np.linalg.matrix_power(a,5))
        self.assertEqual(set(exact_sparse_regions(.99,1.,2)),{(0,0),(1,0),(2,0),(2,1)})

    def test_work_and_input_guards(self):
        for occ in ([],[[0]*5],[[1.5]*4],[[17]*4],[[0]*4]*4097):
            with self.assertRaises(ValueError):mask_matrices(occ,.7)
        with self.assertRaises(ValueError):pair_coefficients(.7,.1,65)
        with self.assertRaises(ValueError):exact_sparse_regions(.7,.1,3)

    def test_projection_never_sends_flat_mass_to_zero(self):
        matrix=np.zeros((17,17))
        matrix[:16,16]=1.
        projected=project_flat_to_masks(matrix)
        np.testing.assert_array_equal(projected[:,0],0.)
        np.testing.assert_array_equal(projected[:,16],0.)
        np.testing.assert_allclose(projected[:16,1:16].sum(axis=1),1.)
        # Projection precedes the fractional power, preserving addition of
        # the two refresh branches instead of summing their alpha powers.
        matrix[0,15]=.25
        self.assertLess(project_flat_to_masks(matrix)[0,15]**.1,
                        project_flat_to_masks(matrix**.1)[0,15])


if __name__=='__main__':unittest.main()
