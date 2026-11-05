from fractions import Fraction as F
import unittest
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp

from hd_rankone_mask_screen import pair_coefficients
from hd_rankone_mask_proposal import (marked_matrix,refresh_power,perron_diagnostic,
                                     bounded_proposal)
import hd_rankone_common_screen_tests as common_tests


class RankOneProposalTests(unittest.TestCase):
    def test_direct_contraction_matches_small_polynomial(self):
        z,alpha,odds=.83,.37,(-2.3,-1.7)
        coefficients=pair_coefficients(z,alpha,2)
        logs=np.array(tuple(coefficients.values()))
        degrees=np.array(tuple(coefficients))
        expected=logsumexp((logs+(degrees@odds)[:,None,None])*np.log(2),axis=0)/np.log(2)
        np.testing.assert_allclose(marked_matrix(z,alpha,odds,maximum=2),expected,
                                   rtol=1e-13,atol=2e-13)

    def test_second_refresh_preserves_zero_and_scales_whole_branches(self):
        cap=np.arange(1,290,dtype=float).reshape(17,17)
        second=refresh_power(cap,2)
        np.testing.assert_array_equal(second[:,0],cap[:,0])
        np.testing.assert_allclose(second[:,1:16],cap[:,1:16]/256)
        np.testing.assert_allclose(second[:,16],cap[:,16]*257/256)
        np.testing.assert_array_equal(refresh_power(cap,1),cap)

    def test_exact_two_refresh_kernel_retains_the_scalar_line(self):
        common_tests.RankOneCommonTests.setUpClass()
        states=common_tests.RankOneCommonTests.states
        kernel=common_tests.RankOneCommonTests.kernel
        a,eta=F(3,15),F(1,4)
        for i,x in enumerate(states):
            line={tuple(common_tests.mul(c,v) for v in x) for c in range(1,4)}
            for j,y in enumerate(states):
                actual=sum(kernel[i][k]*kernel[k][j] for k in range(16))
                expected=(F(int(i==j)) if i==0 else F(0) if j==0 else
                    a*a*int(i==j)+(eta*eta-a*a)*F(int(y in line),3)+(1-eta*eta)/15)
                self.assertEqual(actual,expected)

    def test_perron_weights_for_a_noncommuting_positive_fixture(self):
        matrix=np.arange(1,290,dtype=float).reshape(17,17)/300
        data=perron_diagnostic(np.log2(matrix))
        self.assertAlmostEqual(sum(data['state_weights']),1.)
        self.assertAlmostEqual(sum(data['live_mask_weight_by_size'])+data['flat_weight'],1.)
        self.assertAlmostEqual(data['flat_to_nonflat'],data['nonflat_to_flat'])

    def test_search_budget_is_hard_and_flags_stay_float(self):
        def fixture(z,alpha,markers):
            return dict(log2_bound_FLOAT=z+alpha+sum(v*v for v in markers),
                        z=z,alpha=alpha,marker_log2=list(markers),directed=False,
                        whole_code_certificate=False)
        with patch('hd_rankone_mask_proposal.fixed_gate',side_effect=fixture) as call:
            result=bounded_proposal(max_evaluations=12)
            self.assertLessEqual(call.call_count,12)
            self.assertEqual(call.call_count,result['evaluations'])
            self.assertFalse(result['directed'])
            self.assertFalse(result['whole_code_certificate'])
        with self.assertRaises(ValueError):bounded_proposal(max_evaluations=49)
        with self.assertRaises(ValueError):marked_matrix(.5,.1,(0.,0.),maximum=65)
        with self.assertRaises(ValueError):refresh_power(np.eye(17),3)


if __name__=='__main__':unittest.main()
