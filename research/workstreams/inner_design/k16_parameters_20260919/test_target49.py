"""Regression and exact arithmetic checks for the 49-bit exploration adapters."""
from fractions import Fraction as F
import unittest
import numpy as np
from flint import ctx
import model
import subspace
import target49_q1
import target49_rounds_model
import target49_rounds_q1
import target49_mixing


class Target49(unittest.TestCase):
    def setUp(self): ctx.prec=256

    def test_sharp_baseline(self):
        new=target49_q1.Engine(128,19)
        old=model.study.core.mixing.Engine(16,'1')
        self.assertEqual(new.epoch(F(-5)),old.epoch(F(-5)))

    def test_small_step_float_transfer(self):
        new=target49_q1.Engine(32,15)
        rz,ra=model.study.core.epoch_logs(new.record,np.array([-5.]),sharp=True)
        for a,b in zip(new.epoch(F(-5)),(rz[0],ra[0])):
            np.testing.assert_allclose(np.array([float(x) for x in a]).reshape(new.n,new.n),np.exp(b),rtol=2e-13,atol=1e-300)

    def test_bernoulli_bridge_ratio(self):
        for m in (3,7,4095,8191):
            f=F(4*m+3,4*m+2)
            for cap in (F(1),F(3,2),F(100),F(m+1)):
                old=cap/F(2*(m+1))+F(1,2*m)
                new=cap/F(4*(m+1))+F(3,4*m)
                self.assertLessEqual(new,f*old)
                if cap==1: self.assertEqual(new,f*old)

    def test_two_round_identity_and_zero_source(self):
        record=subspace.search(12)
        old=subspace.Engine(record)
        new=target49_rounds_model.Engine(record)
        q1=target49_rounds_q1.Engine(record)
        self.assertEqual(new.identity(),q1.identity())
        self.assertEqual(new.identity()['inner']['transvection_rounds'],2)
        for a,b in zip(old.epoch(F(-4)),new.epoch(F(-4))):
            self.assertEqual(a[:old.n],b[:new.n])
            self.assertTrue(all(v>=0 for v in b))

    def test_diagnostic_does_not_mutate_historical_functions(self):
        old=model.study.core.epoch_logs
        target49_mixing.evaluate(model.study.core.inner(16,10),2)
        self.assertIs(old,model.study.core.epoch_logs)


if __name__=='__main__': unittest.main()
