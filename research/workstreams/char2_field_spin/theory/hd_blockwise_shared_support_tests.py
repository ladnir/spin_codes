import unittest
from itertools import combinations

from hd_blockwise_shared_maps import rank,identity
from hd_blockwise_shared_support import (tile_generator,paired_masks,rank_profile,
    state_residual,stacked_output_rank,output_constraints,kernel,plane_comparison)


class SharedSupportTests(unittest.TestCase):
    def test_tile_is_mds_16_12(self):
        for beta in (8,64):
            generator=tile_generator(0,beta)
            for selected in combinations(range(16),12):
                self.assertEqual(rank(tuple(generator[i] for i in selected)),12)

    def test_local_and_stacked_rank(self):
        for beta in (8,64):
            for w in range(9):
                for z in range(9):
                    for tile in (0,1,7):
                        local=state_residual(tile,range(w),range(z),beta)
                        self.assertEqual(local['input_pivots'],min(w,z))
                        self.assertEqual(local['rank'],min(4,max(0,z-w)))
            s=(0,2,5,8,11,19,20,33,48,49,50)
            z=(0,1,2,4,5,8,9,10,11,12,16,17,34,35,48,49,50,51)
            self.assertEqual(stacked_output_rank(s,z,beta)['output_rank'],rank(output_constraints(s,z,beta)))

    def test_paired_defect(self):
        for weight in range(5,9):
            support,zeros=paired_masks(tuple(range(weight//2))+tuple(range(8,8+weight-weight//2)))
            profile=rank_profile(support,zeros)
            self.assertEqual(profile['output_nullity'],weight-4)
            self.assertEqual(profile['state_image_rank'],weight-4)
            self.assertEqual(profile['feedback_image_rank'],0)
            self.assertTrue(profile['full_input_support_possible'])
            self.assertFalse(profile['live_input_and_return_possible'])
        support,zeros=paired_masks((0,1,2,8,9,10))
        result=stacked_output_rank(support,zeros)
        self.assertEqual(result['local_state_ranks'][:2],[2,2])
        self.assertEqual(result['state_rank'],2)

    def test_prescribed_source_dependencies(self):
        support,zeros=paired_masks((0,1,2,8,9,10))
        constraints=output_constraints(support,zeros)
        embedding=kernel(constraints,10)
        profile=rank_profile(support,zeros,source_embedding=embedding)
        self.assertEqual(profile['source_dimension'],2)
        self.assertEqual(profile['output_rank'],0)
        self.assertEqual(profile['state_image_rank'],2)
        self.assertEqual(profile['feedback_image_rank'],0)

    def test_disjoint_shift_and_remaining_defect(self):
        support,zeros=paired_masks((0,1,2,8,9,10))
        fixed=rank_profile(support,zeros,64)
        self.assertEqual(fixed['output_rank'],10)
        self.assertFalse(fixed['full_input_support_possible'])
        support=(0,1,2,9,12,15)
        zeros=(0,1,2,3,4,8,9,13,14,15)
        profile=rank_profile(support,zeros,64)
        self.assertEqual(profile['output_rank'],9)
        self.assertEqual(profile['output_nullity'],1)
        self.assertEqual(profile['state_image_rank'],1)
        self.assertEqual(profile['feedback_image_rank'],1)
        self.assertTrue(profile['full_input_support_possible'])
        self.assertEqual(stacked_output_rank(support,zeros,64)['state_rank'],3)

    def test_bounded_plane_catalog(self):
        result=plane_comparison()
        self.assertEqual(result['inventories']['8']['unique_planes'],[1408,1408])
        self.assertEqual(result['inventories']['8']['identical_planes'],1408)
        self.assertEqual(result['inventories']['64']['unique_planes'],[1408,1408])
        self.assertEqual(result['inventories']['64']['identical_planes'],0)
        self.assertEqual(result['inventories']['64']['plane_pairs'],1982464)
        self.assertEqual(result['inventories']['64']['rank_three_pairs'],4704)
        self.assertFalse(result['uniform_probability_claim'])


if __name__=='__main__':
    unittest.main()
