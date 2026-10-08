from fractions import Fraction as F
import unittest

import numpy as np

from high_distance_expansion_screen import bound,period_matrix,profile
from high_distance_fields import bound as fixed_bound


class ExpansionScreenTests(unittest.TestCase):
    def test_fixed_profile_agreement(self):
        for active in (1,8,16):
            self.assertAlmostEqual(bound((5.,2.),active,length=64,dimension=16,
                pattern=(2,),distance=F(9,25)),fixed_bound((5.,2.),active,bits=16),places=7)
        point=(-2.,-.3,2.,0.)
        self.assertAlmostEqual(bound(point,100000,False,length=64,dimension=16,
            pattern=(2,),distance=F(9,25)),fixed_bound(point,100000,False,bits=16),places=7)

    def test_threefold_geometry(self):
        g=profile()
        self.assertEqual(g['rate'],'1/2')
        self.assertEqual(g['internal_expansion'],'3')
        self.assertEqual(g['rows']%192,0)
        self.assertEqual(g['input_n'],3*g['padded_k'])
        self.assertEqual(g['output_n'],2*g['padded_k'])
        self.assertGreaterEqual(g['padded_k'],g['requested_k'])
        self.assertEqual(g['output_threshold'],(3*g['output_n']+9)//10-1)

    def test_period_repetition_agreement(self):
        for sparse,value,point in ((True,8,(5.,2.)),(False,100000,(-2.,-.3,2.,0.))):
            self.assertAlmostEqual(bound(point,value,sparse,length=64,dimension=16,pattern=(2,)),
                bound(point,value,sparse,length=64,dimension=16,pattern=(2,2)),places=7)

    def test_chronological_not_averaged(self):
        point=(-2.,-.3,2.)
        p,scale=period_matrix(point,(3,3,2),4,4,16)
        a,sa=period_matrix(point,(3,),4,4,16)
        b,sb=period_matrix(point,(2,),4,4,16)
        np.testing.assert_allclose(p*np.exp(scale),(a@a@b)*np.exp(2*sa+sb),rtol=1e-13)
        reverse,sr=period_matrix(point,(2,3,3),4,4,16)
        self.assertGreater(np.max(abs(p*np.exp(scale)-reverse*np.exp(sr))),1e-8)

    def test_invalid_geometry(self):
        with self.assertRaises(ValueError):profile(pattern=())
        with self.assertRaises(ValueError):profile(pattern=(0,))
        with self.assertRaises(ValueError):profile(length=300)
        with self.assertRaises(ValueError):bound((-2.,-.3,2.,0.),0,False)

    def test_promising_fourfold_fixed_witnesses(self):
        options=dict(bits=16,length=256,dimension=64,pattern=(2,),distance=F(2,5))
        self.assertAlmostEqual(-bound((6.306439155378674,2.6479784410681093),64,**options),
            828.4568807307,places=4)
        self.assertGreater(-bound((-2.624817162721374,-.15292178839610845,
            2.5482573881579826,2.663796147433451),125829,False,**options),1600)
        self.assertGreater(-bound((5.2830720593937714,-11.487465733567783,
            11.090339630053647,-7.349901065339367),4152361,False,**options),1400000)
        g=profile(length=256,dimension=64,pattern=(2,),distance=F(2,5))
        first=(3*g['input_n']+99)//100
        self.assertEqual((first-1)//193,651)
        self.assertGreaterEqual(652*193,first)

    def test_larger_outer_needs_a_larger_fixed_field(self):
        with self.assertRaises(ValueError):profile(length=512,dimension=256,pattern=(4,))
        g=profile(length=512,dimension=256,pattern=(4,),outer_bits=16)
        self.assertEqual(g['outer_field_bits'],16)
        self.assertEqual(g['outer_max_nodes'],65536)
        self.assertEqual(g['fixed_inner_field_bits'],8)
        self.assertEqual(g['input_n'],2**21)
        self.assertEqual(g['rate'],'1/2')
        with self.assertRaises(ValueError):bound((5.,2.),1,bits=8,length=512,
            dimension=256,pattern=(4,),outer_bits=16)


if __name__=='__main__':unittest.main()
