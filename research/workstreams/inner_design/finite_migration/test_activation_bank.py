import itertools
import math
from types import SimpleNamespace
import unittest

import numpy as np
import parameter_activation_bank as bank
import parameter_activation_partition as partition


class ActivationBankTests(unittest.TestCase):
    def test_fixed_probability_proposal(self):
        checker = SimpleNamespace(length=12,block=128,log_gammas=np.array([0.,3.,0.]),
                                  cutoff=20,ps=np.array([0.,.5,1.]))
        for box in (dict(lower=[3,4,0],upper=[8,9,5]),
                    dict(lower=[10,2,0],upper=[10,2,0])):
            corners = bank.base.typed.vertices(box['lower'],box['upper'],12)
            for theta in (.1,.35,.5,.8,.99):
                tilt,term = -.5,-70.
                value,witness = bank.proposal_bound(checker,box,theta,tilt,term)
                proposal = np.array(witness['proposal'])
                self.assertTrue(np.all(proposal > 0))
                self.assertAlmostEqual(float(sum(proposal)),1.,places=14)
                self.assertAlmostEqual(float(proposal@checker.ps),theta,places=14)
                penalty = bank.base.typed.lattice_log_count(box['lower'],box['upper'])
                recomputed = float(max(bank.base.typed.point_logs(corners,12,128,checker.log_gammas,
                    proposal,term,20,math.exp(tilt))))+penalty
                self.assertAlmostEqual(value,recomputed,places=11)
                lo,hi = max(0.,2*theta-1),theta
                for u in np.linspace(lo,hi,102)[1:-1]:
                    trial = np.array([1-2*theta+u,2*(theta-u),u])
                    bound = float(max(bank.base.typed.point_logs(corners,12,128,checker.log_gammas,
                        trial,term,20,math.exp(tilt))))+penalty
                    self.assertLessEqual(value,bound+1e-8)

    def test_partition_exact_small_domains(self):
        for total in (3,7,12):
            for box in (dict(lower=[0,0,0],upper=[total]*3),
                        dict(lower=[total-2,1,0],upper=[total,total,total]),
                        dict(lower=[total-1,1,0],upper=[total-1,total,total])):
                points = {p for p in itertools.product(range(total+1),repeat=3)
                          if sum(p)==total and all(a<=v<=b for a,v,b in zip(box['lower'],p,box['upper']))}
                children = partition.split_box(box,total)
                if len(points) == 1:
                    self.assertEqual(children,[])
                    continue
                self.assertEqual(len(children),2)
                groups = [{p for p in points if all(a<=v<=b for a,v,b in zip(c['lower'],p,c['upper']))}
                          for c in children]
                self.assertFalse(groups[0]&groups[1])
                self.assertEqual(groups[0]|groups[1],points)


if __name__ == '__main__':
    unittest.main()
