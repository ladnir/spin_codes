"""FLOAT projected proposals retain scope and agree with fixed replay controls."""
from fractions import Fraction as F
from math import log2
import unittest

import numpy as np

from hd_shifted_projected_screen import ProjectedScreen,alias_table
from hd_shifted_swap_screen import CONSTRUCTION


class ProjectedScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        row = tuple(log2(F(q)) for q in ('21153727/526426284','62479461/686625565'))
        cls.screen = ProjectedScreen(F('6231346319206105/9007199254740992'),
                                    F('3529422830148241/36028797018963968'),
                                    (row,row),('375/473','49/11825'))

    def test_float_cushion_agrees_with_frozen_directed_control(self):
        value = self.screen.evaluate(((500,524),(188,212),(188,212)))
        # Numerical agreement is a regression diagnostic, NOT certification.
        self.assertAlmostEqual(value['log2_upper_screen'],-385.442619507605,delta=.001)
        self.assertEqual(value['projected_histogram_counts'],[625,625])
        self.assertFalse(value['directed'])
        self.assertFalse(value['whole_code_certificate'])
        self.assertEqual(value['construction'],CONSTRUCTION)

    def test_box_envelope_dominates_selected_singletons(self):
        outer = self.screen.evaluate(((508,516),(196,204),(196,204)))['log2_upper_screen']
        for r,a,b in ((508,196,204),(512,200,200),(516,204,196)):
            point = self.screen.evaluate(((r,r),(a,a),(b,b)))['log2_upper_screen']
            self.assertGreaterEqual(outer,point)

    def test_scalar_and_vector_regional_entries_agree(self):
        vector = self.screen.regional(0,np.array([512,514]),np.array([200,201]))
        np.testing.assert_allclose(vector[:,0],self.screen.regional(0,512,200))
        np.testing.assert_allclose(vector[:,1],self.screen.regional(0,514,201))

    def test_cached_atlas_is_immutable(self):
        value = alias_table(self.screen.z,self.screen.alpha,self.screen.odds[0],128)
        self.assertEqual(value.shape,(4,128,128))
        self.assertFalse(value.flags.writeable)
        self.assertIs(value,self.screen.tables[0])

    def test_identity_and_exact_count_markers_survive_serialization(self):
        witness = self.screen.serializable()
        self.assertEqual(witness['count_markers'],['375/473','49/11825'])
        self.assertTrue(witness['projected_normalization'])
        self.assertTrue(witness['pair_swaps'])
        self.assertEqual(witness['orientation_aggregation'],'independent-fair-swaps')
        self.assertEqual(witness['construction'],CONSTRUCTION)

    def test_invalid_geometry_and_inputs_rejected(self):
        for parity,r,f in ((2,1,0),(0,1,1),(0,512.5,200),(0,-1,0)):
            with self.assertRaises(ValueError):
                self.screen.regional(parity,r,f)
        with self.assertRaises(ValueError):
            self.screen.evaluate(((4,4),(2,2),(2,2)))
        with self.assertRaises(ValueError):
            alias_table(.5,.1,(-2.,-2.),96)


if __name__=='__main__':
    unittest.main()
