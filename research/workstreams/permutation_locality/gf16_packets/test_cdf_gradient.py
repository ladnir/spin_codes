from math import comb
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch
import unittest
import numpy as np
from scipy.special import expit
from flint import arb,arb_mat,ctx
import cdf_gradient as gradient
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from occupancy_cdf_cover import fold_log,cover
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment


class CdfGradientTests(unittest.TestCase):
    def test_optional_gradient_proposals_still_run_outward_cover(self):
        ctx.prec=192
        exact=[arb_mat([[arb(1)/4]]) for _ in range(3)]
        operators={('.001','1'):(exact,[np.array([[.25]]) for _ in exact])}
        args=SimpleNamespace(groups=2,joint_witness=True,joint_top=1,probe_supports=[],probe_vector=[],
            retain_parents=True,target_bits=48,max_splits=0,check_interval=1,screen_only=False,precision=192)
        counts={'1':[int(u>=38) for u in range(257)]}
        for enabled in (False,True):
            args.analytic_gradient=enabled
            with patch('gf16_packets.cdf_gradient.JointObjective',wraps=gradient.JointObjective) as objective:
                result=cover(args,operators,counts,np.ones(1),cutoff=0)
                self.assertEqual(bool(objective.call_count),enabled)
            self.assertIsNotNone(result)
            self.assertTrue(0<result<arb(2)**-48)

    def test_grouped_count_derivatives(self):
        logits=np.array([-.8,.3,1.1]);counts=[1,4,7];step=1e-6
        masses,slopes=gradient.count_distribution(expit(logits),counts)
        expected=np.array([1.])
        for p,n in zip(expit(logits),counts):
            expected=np.convolve(expected,[comb(n,j)*p**j*(1-p)**(n-j) for j in range(n+1)])
        np.testing.assert_allclose(masses,expected,rtol=2e-14,atol=1e-15)
        np.testing.assert_allclose(slopes.sum(axis=1),0,atol=1e-14)
        for j in range(len(logits)):
            direction=np.eye(len(logits))[j]*step
            plus=gradient.count_distribution(expit(logits+direction),counts)[0]
            minus=gradient.count_distribution(expit(logits-direction),counts)[0]
            np.testing.assert_allclose(slopes[j],(plus-minus)/(2*step),rtol=1e-6,atol=5e-11)

    def test_scaled_power_derivatives(self):
        matrix=np.array([[.1,.3],[.2,.25]]);slopes=np.array([[[.01,-.02],[.03,.01]],[[.03,.01],[-.01,.02]]])
        terminal=np.array([1.,.5]);step=1e-6
        for length in (1,3,16,256):
            value,actual=gradient.log_power_gradient(matrix,slopes,length,terminal)
            self.assertAlmostEqual(value,log_power_moment(matrix,length,terminal),places=11)
            for i,slope in enumerate(slopes):
                expected=(log_power_moment(matrix+step*slope,length,terminal)
                          -log_power_moment(matrix-step*slope,length,terminal))/(2*step)
                self.assertAlmostEqual(actual[i],expected,delta=1e-7)
        value,actual=gradient.log_power_gradient(np.array([[1e-200]]),np.array([[[2e-200]]]),256,np.ones(1))
        self.assertAlmostEqual(value,256*np.log(1e-200),places=8)
        self.assertAlmostEqual(actual[0],512,places=10)

    def test_joint_value_and_gradient_match_existing_objective(self):
        box=((38,70),)*3+((90,145),)*2+((155,230),)*4
        counts=[u*u*u+1 for u in range(257)]
        region=np.array([[[.7-.03*j,.06+.003*j],[.1+.01*j,.4-.01*j]] for j in range(len(box)+1)])
        terminal=np.ones(2);offset=18.;objective=gradient.JointObjective(region,box,counts,terminal,offset)
        def reference(logits):
            ps=expit(logits);probabilities=[ps[objective.intervals.index(interval)] for interval in box]
            return (log_power_moment(matrix_for_probabilities(region,probabilities),256,terminal)
                    +sum(fold_log(counts,*interval,p) for interval,p in zip(box,probabilities))+offset)
        for logits in (np.array([-.91,.13,1.37]),np.array([-.3,-.6,.8])):
            value,slopes=objective(logits)
            self.assertAlmostEqual(value,reference(logits),places=10)
            for j in range(len(logits)):
                direction=np.eye(len(logits))[j]*1e-6
                difference=(reference(logits+direction)-reference(logits-direction))/(2e-6)
                self.assertAlmostEqual(slopes[j],difference,delta=2e-6)


if __name__=='__main__':unittest.main()
