"""Combinatorial checks for the h4 shifted-pair FLOAT screen."""
from math import log,log2
import unittest

import numpy as np
from scipy.special import logsumexp

from hd_bucket_feedback import profiles
from hd_grouped_packet_screen import OCC, coefficients, local_logs
from hd_shifted_pairs_screen import half_weights,mixed_step,mixed_coefficients,positive_alias_power,proposal


class ShiftedPairScreenTests(unittest.TestCase):
    def test_half_polynomial_normalizes(self):
        for singleton,full in ((-2.,-4.),(-3.,-.5),(1.,2.)):
            actual = logsumexp(half_weights(singleton,full)*log(2))/log(2)
            expected = 16*log2(1+2*2**singleton+2**full)
            self.assertAlmostEqual(actual,expected,places=11)

    def test_full_only_limit_matches_unshifted_profile_operator(self):
        z,alpha,odds = .87,.6,-2.
        degree,coeff = coefficients(z,alpha,(2,),h=4,route_width=2)
        expected = logsumexp((coeff+degree[:,0,None,None]*odds)*log(2),axis=0)/log(2)
        actual = mixed_step(z,alpha,-100.,odds)
        np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=0)

    def test_independent_lane_markers_match_all_scalar_profiles(self):
        z,alpha,odds = .91,.7,-2.
        items = profiles(16,64,4)
        self.assertEqual(tuple(tuple(row) for row in OCC),tuple(row for row,_ in items))
        weights = np.array([log2(count)+sum(row)*odds for row,count in items])
        expected = logsumexp((weights[:,None,None]+alpha*local_logs(z))*log(2),axis=0)/log(2)
        actual = mixed_step(z,alpha,odds,2*odds)
        np.testing.assert_allclose(actual,expected,atol=2e-12,rtol=0)

    def test_shifted_singleton_orientations_are_balanced(self):
        for word in range(1,255):
            contracted = [(word>>i)&1 for i in range(8)]
            active = sum(contracted)
            runs = sum(contracted[i] and not contracted[(i-1)%8] for i in range(8))
            counts = {(a,b):0 for a in (0,1) for b in (0,1)}
            for i in range(8):
                counts[contracted[i],contracted[(i+1)%8]] += 1
            self.assertEqual(counts[1,0],runs)
            self.assertEqual(counts[0,1],runs)
            self.assertEqual(counts[1,1],active-runs)

    def test_boundary_initialization_with_zero_empty_pairs(self):
        for active,runs in ((2048,2048),(3072,1024)):
            result = proposal(active,runs,iterations=1)
            self.assertTrue(np.isfinite(result['log2_upper_screen']))

    def test_three_variable_polynomial_matches_positive_evaluation(self):
        z,alpha,singleton,full = .91,.7,-2.,-3.
        polynomial = mixed_coefficients(z,alpha,singleton,full)
        expected = np.exp2(mixed_step(z,alpha,singleton,full))
        np.testing.assert_allclose(polynomial.sum(axis=(1,2,3)).reshape((2,2)),expected,rtol=2e-13,atol=0)

    def test_alias_extraction_matches_unaliased_single_step(self):
        polynomial = mixed_coefficients(.48,.095,-4.3,-2.4)
        target = (1,1,0)
        logs,diagnostics = positive_alias_power(polynomial,target,64,steps=1,cushion=0.)
        expected = polynomial[(slice(None),)+target].reshape((2,2))
        np.testing.assert_allclose(np.exp2(logs),expected,rtol=2e-10,atol=0)
        self.assertLess(diagnostics['maximum_target_imaginary'],1e-12)

    def test_alias_rejects_ill_conditioned_tiny_target_coefficient(self):
        # This fixture lost 3.47e-5 relatively in its tiny zero-to-zero
        # coefficient despite an absolute difference of only 1.11e-19.
        polynomial = mixed_coefficients(.91,.7,-2.,-3.)
        with self.assertRaisesRegex(ArithmeticError,'cancellation floor'):
            positive_alias_power(polynomial,(3,4,5),64,steps=1,cushion=0.)


if __name__=='__main__':
    unittest.main()
