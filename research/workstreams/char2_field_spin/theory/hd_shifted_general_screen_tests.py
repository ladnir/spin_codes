import unittest
import numpy as np
from hd_shifted_general_screen import general_step,general_coefficients,bound,regional,homogeneous_coefficients
from hd_shifted_pairs_screen import mixed_step,mixed_coefficients,regional_matrices
from hd_shifted_general_count import matching_histograms


class GeneralShiftedScreenTests(unittest.TestCase):
    def test_balanced_marker_reduces_to_existing(self):
        for z,alpha,single,full in ((.91,.7,-2.,-3.),(.6,.1,-4.,-2.)):
            np.testing.assert_allclose(general_step(z,alpha,(single,single,full)),
                                       mixed_step(z,alpha,single,full),atol=2e-12,rtol=0)

    def test_oriented_polynomial_matches_evaluation(self):
        for odds in ((-2.,-3.,-4.),(-1.,-4.,-2.)):
            polynomial = general_coefficients(.91,.7,odds)
            expected = np.exp2(general_step(.91,.7,odds))
            np.testing.assert_allclose(polynomial.sum(axis=(1,2,3)).reshape((2,2)),
                                       expected,rtol=2e-13,atol=0)

    def test_balanced_coefficients_match_existing(self):
        np.testing.assert_allclose(general_coefficients(.87,.6,(-2.,-2.,-3.)),
                                   mixed_coefficients(.87,.6,-2.,-3.),atol=0,rtol=0)

    def test_orientation_swap_preserves_local_operator(self):
        np.testing.assert_allclose(general_step(.81,.35,(-2.,-4.,-3.)),
                                   general_step(.81,.35,(-4.,-2.,-3.)),atol=2e-12,rtol=0)

    def test_full_pair_global_region_matches_existing(self):
        hist = matching_histograms(4096,512,512,512,384)[0]
        actual,_ = regional(.7,.08,hist,(-40.,-40.,-2.8))
        expected,_ = regional_matrices(512,128,.7,.08,(-2.8,-4.,-3.))
        # Absent categories have a finite tiny marker in the global bound,
        # so the new nonnegative coefficient upper bound is slightly looser.
        np.testing.assert_allclose(actual,expected,atol=2e-7,rtol=0)

    def test_paired_full_coefficients_match_existing(self):
        from hd_grouped_packet_screen import coefficients
        degrees,logs = coefficients(.7,.08,(2,),h=4,route_width=2,steps=2)
        expected = np.exp2(logs)
        np.testing.assert_allclose(homogeneous_coefficients(.7,.08,2),expected,rtol=1e-12,atol=0)

    def test_homogeneous_composition_at_alpha_one(self):
        for category in (0,1,2):
            single = homogeneous_coefficients(.87,1.,category,1)
            expected = np.zeros((65,2,2))
            for i,first in enumerate(single):
                for j,second in enumerate(single):
                    expected[i+j] += first@second
            actual = homogeneous_coefficients(.87,1.,category,2)
            np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=0)


if __name__=='__main__':
    unittest.main()
