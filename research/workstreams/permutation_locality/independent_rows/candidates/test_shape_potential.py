"""Independent small-geometry checks of the continuation relaxation."""
from itertools import combinations, product
from math import comb, exp
import unittest
from unittest.mock import patch

import numpy as np
from flint import arb_mat

import shape_potential as candidate


class ShapePotentialTests(unittest.TestCase):
    def setUp(self):
        self.families=[np.array([[[.8,.1],[.2,.6]]]),
                       np.array([[[.6,.4],[.1,.7]],[[.2,.8],[.5,.3]]]),
                       np.array([[[.1,.2],[.4,.1]],[[.25,.05],[.1,.4]]])]

    def test_one_shape_matches_matrix_placement(self):
        families=[x[:1] for x in self.families]
        v=np.array([.3,.7])
        for epochs in (1,2,3,5):
            weights=candidate.placement_weights(2*epochs,epochs,2)
            actual=candidate.region_action(families,v,weights)
            expected=candidate.screen.float_placement(np.array([x[0] for x in families]),
                                                     2*epochs,epochs,2)@v
            np.testing.assert_allclose(actual,expected,rtol=2e-14,atol=2e-14)

    def test_bounds_every_small_shape_sequence(self):
        v=np.array([.3,.7])
        for epochs in (1,2,3):
            actual=candidate.region_action(self.families,v,
                       candidate.placement_weights(2*epochs,epochs,2))
            for r in range(2*epochs+1):
                total=np.zeros(2)
                for slots in combinations(range(2*epochs),r):
                    occupancies=[sum(i//2==epoch for i in slots) for epoch in range(epochs)]
                    outcomes=[]
                    for choices in product(*(self.families[j] for j in occupancies)):
                        matrix=np.eye(2)
                        for a in choices:matrix=matrix@a
                        outcomes.append(matrix@v)
                    total+=np.max(outcomes,axis=0)
                total/=comb(2*epochs,r)
                self.assertTrue(np.all(actual[r]>=total-2e-14))

    def test_no_larger_than_entrywise_envelope(self):
        v=np.array([.3,.7])
        envelope=np.array([a.max(axis=0) for a in self.families])
        for epochs in (1,2,3,5):
            actual=candidate.region_action(self.families,v,
                       candidate.placement_weights(2*epochs,epochs,2))
            loose=candidate.screen.float_placement(envelope,2*epochs,epochs,2)@v
            self.assertTrue(np.all(actual<=loose+2e-14))
        self.assertLess(actual[3,0],loose[3,0])

    def test_placement_weights_sum_to_one(self):
        for limit,records in candidate.placement_weights(13,epochs=7,windows=3):
            mass=np.zeros(limit+1)
            for _,rs,_,p in records:mass[rs]+=p
            np.testing.assert_allclose(mass,1,rtol=2e-14,atol=2e-14)

    def test_invalid_placement_geometry(self):
        for args in ((-1,2,2),(5,2,2),(1,0,2),(1,2,0),(True,2,2)):
            with self.assertRaises(ValueError):candidate.placement_weights(*args)

    def test_complete_shape_catalog_required(self):
        base=[arb_mat([[int(i==j) for j in range(11)] for i in range(11)]) for _ in range(33)]
        shaped={2:{s:base[2]*1 for s in candidate.expected_shapes(2)}}
        families=candidate.as_families(base,shaped)
        self.assertEqual(families[2].shape,(10,11,11))
        shaped[2].pop((4,4))
        with self.assertRaises(ValueError):candidate.as_families(base,shaped)
        with self.assertRaises(ValueError):candidate.as_families(base,{33:{}})

    def test_shape_metadata_must_match(self):
        base=[arb_mat(11,11) for _ in range(33)]
        for params in ((candidate.Q('.06'),candidate.Q('.9'),2),
                       (candidate.Q('.056'),candidate.Q('.8'),2),
                       (candidate.Q('.056'),candidate.Q('.9'),3)):
            with self.assertRaisesRegex(ValueError,'must match'):
                candidate.shape_matrices(base,{}, {'parameters':params},'.056','.9',2)

    def test_single_packet_activation_retains_its_weight(self):
        # Isolate the exact Z->F contribution from the expensive, independently
        # tested census refinements. Other entries are irrelevant to this test.
        base=[arb_mat([[1]*11 for _ in range(11)]) for _ in range(33)]
        shapes=candidate.expected_shapes(1)
        spectrum={v:1 for v in (48,56,64,72,80)}
        details={'histograms':None,'spectrum':spectrum,
                 'moments':{tuple(s.count(b) for b in range(1,5)):None for s in shapes},
                 'cancellations':({s:None for s in shapes},spectrum),
                 'density':{s:{'density':0,'uniform':{v:0 for v in spectrum}} for s in shapes}}
        def unchanged(matrices,*args,**kwargs):return matrices
        with patch.object(candidate.screen.window_histogram,'refine',side_effect=unchanged), \
             patch.object(candidate.screen.full_feedback_refinement,'refine',side_effect=unchanged), \
             patch.object(candidate.screen.cancellation_joint,'refine',side_effect=unchanged):
            for rounds in (1,2,3):
                details['parameters']=(candidate.Q('.056'),candidate.Q('.9'),rounds)
                result=candidate.shape_matrices(base,{s:None for s in shapes},details,
                                                '.056','.9',rounds,maximum=1)
                for (b,),matrix in result[1].items():
                    actual=float(matrix[candidate.Z,candidate.F])
                    expected=exp(-.056*b)*(.9 if b==4 else 1)
                    self.assertAlmostEqual(actual,expected,places=14)

    def test_common_potential_bounds_finite_linear_moment(self):
        matrix=.01*np.ones((11,11))+.5*np.eye(11)
        families=[matrix[None,:,:],matrix[None,:,:]]
        weights=candidate.placement_weights(2,epochs=2,windows=1)
        masses=np.array([.2,.5,.3])
        result=candidate.common_potential(families,weights,masses,np.ones(11),regions=7,iterations=3)
        exact=(np.linalg.matrix_power(matrix,14)@candidate.screen.baseline.TAIL_TERMINAL)[0]
        self.assertGreaterEqual(exp(result['log_moment']),exact-1e-14)

    def test_zero_auxiliary_potential_is_rejected(self):
        v=np.ones(11);v[3]=0
        with self.assertRaises(ValueError):
            candidate.common_potential([],[],np.array([1.]),v)

    def test_invalid_potential_iteration_counts(self):
        for values in ({'regions':0},{'iterations':0},{'regions':True}):
            with self.assertRaises(ValueError):
                candidate.common_potential([],[],np.array([1.]),np.ones(11),**values)


if __name__=='__main__':unittest.main()
