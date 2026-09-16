"""Regression checks for the discovery evaluator, not proof certification."""
import math
import unittest
import numpy as np
import screen
import certify_imt_dense as certify_dense
import sparse
from fractions import Fraction as F


class ScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = screen.Model()

    def test_occupation_matches_original(self):
        model = self.model
        e = model.engine
        for beta,lam in ((.0001,.0002),(.01,.03),(.2,.8),(.8,2.)):
            epochs = screen.candidate.model.independent.g.epochs(e.spectrum,e.kernel,e.columns,e.caps,e.low,lam)
            probabilities = model.logcombs+model.j*math.log(beta)+(128-model.j)*math.log1p(-beta)
            reference = np.exp(np.logaddexp.reduce(epochs+probabilities[:,None,None],axis=0))
            np.testing.assert_allclose(model.occupation(beta,lam),reference,rtol=1e-12,atol=1e-100)

    def test_fourier_matches_original(self):
        model = self.model
        e = model.engine
        for beta,lam in ((.0001,.0002),(.01,.03),(.2,.8),(.8,2.)):
            reference = np.exp(screen.screen_dense.bernoulli(e.spectrum,e.b_spectrum,e.kernel,beta,lam))
            np.testing.assert_allclose(model.fourier(beta,lam),reference,rtol=1e-12,atol=1e-100)

    def test_exact_mean_expansion(self):
        e = self.model.engine
        self.assertTrue(all(e.a_columns))
        self.assertEqual(sum(w*n for w,n in e.spectrum.items()),128*(1<<18))

    def test_uniform_input_identity(self):
        for lam in (.001,.1,1.,2.0907410969):
            self.assertAlmostEqual(self.model.rate(.5,lam,'statefree'),math.log((1+math.exp(-lam))/2),places=14)
        row = screen.optimize(self.model,1.,.5,math.log(2)/2,.11)
        entropy = -.11*math.log(.11)-.89*math.log(.89)
        self.assertLess(row['exponent'],0)
        self.assertAlmostEqual(row['exponent'],entropy-math.log(2)/2,places=10)

    def test_sparse_first_order(self):
        v,h,mean = sparse.witness(self.model)
        alpha = 1e-8
        matrix = self.model.occupation(.8*alpha,-math.log1p(-1.6*alpha))
        vec = np.r_[1.,float(v)*(1+alpha*np.array(list(map(float,h))))]
        residual = (matrix @ vec-(1-96*alpha)*vec)/(alpha*vec)
        expected = [-6.3,float(-h[0]/2-F(8,5)*48+F(4,5)*(F(1,2)+F(64,self.model.m))/v+96)]
        expected += [float(mean+96)]*len(self.model.levels)
        np.testing.assert_allclose(residual,expected,rtol=0,atol=.01)

    def test_geometry_rejects_gap_and_overlap(self):
        from types import SimpleNamespace
        segments = [SimpleNamespace(index=0,lower=F(13,125),upper=F(112,125))]
        full = dict(segment=0,alpha=['1/10000','1'],row_density=['13/125','112/125'])
        certify_dense.check_geometry([full],segments)
        with self.assertRaises(AssertionError):
            certify_dense.check_geometry([dict(full,alpha=['1/1000','1'])],segments)
        with self.assertRaises(AssertionError):
            certify_dense.check_geometry([full,full],segments)


if __name__=='__main__':
    unittest.main()
