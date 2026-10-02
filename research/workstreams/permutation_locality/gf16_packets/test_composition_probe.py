from fractions import Fraction as Q
from math import comb
import unittest
import numpy as np
from scipy.special import logsumexp
from composition_probe import round_counts,poisson_binomial_logs,zero_path_logs


class CompositionTests(unittest.TestCase):
    def test_rounding_preserves_total_and_one_unit_error(self):
        for p in ([0.,.1,.9],[.2,.3,.5],[1.,0.],[1/3]*3):
            for total in (1,7,2048):
                counts=round_counts(p,total)
                self.assertEqual(sum(counts),total)
                self.assertTrue(all(abs(n-total*x)<=1 for n,x in zip(counts,p)))
        for p,total in (([-.1,1.1],3),([.1,.8],3),([1.],0)):
            with self.assertRaises(ValueError):round_counts(p,total)

    def test_direct_count_logs_match_rational_polynomial(self):
        for features,counts in (([Q(0),Q(1)],[3,4]),([Q(1,5),Q(3,4)],[3,4]),
                ([Q(0),Q(1),Q(1,2),Q(3,5)],[2,1,4,2]),([Q(0)],[4]),([],[])):
            exact=[Q(1)]
            for p,n in zip(features,counts):
                for _ in range(n):
                    row=[Q(0)]*(len(exact)+1)
                    for j,x in enumerate(exact):row[j]+=x*(1-p);row[j+1]+=x*p
                    exact=row
            logs=poisson_binomial_logs(features,counts)
            np.testing.assert_allclose(np.exp(logs),list(map(float,exact)),rtol=1e-13,atol=1e-15)
        logs=poisson_binomial_logs([Q(1,1000),Q(999,1000)],[1024,1024])
        self.assertTrue(np.isfinite(logs).all())
        self.assertLess(logs[0],-700)
        self.assertAlmostEqual(logsumexp(logs),0.,places=10)
        for features,counts in (([Q(1,2)],[-1]),([Q(2)],[1]),([Q(0)],[True])):
            with self.assertRaises(ValueError):poisson_binomial_logs(features,counts)

    def test_zero_path_matches_exact_conditional_polynomial(self):
        for local in ((Q(1),Q(1,20),Q(1,30)),(Q(1),Q(0),Q(1,300))):
            for epochs in (1,2,4):
                polynomial=[Q(1)]
                for _ in range(epochs):
                    result=[Q(0)]*(len(polynomial)+2)
                    for k,c in enumerate(polynomial):
                        for j,q in enumerate(local):result[k+j]+=c*comb(2,j)*q
                    polynomial=result
                exact=[c/comb(2*epochs,k) for k,c in enumerate(polynomial)]
                np.testing.assert_allclose(np.exp(zero_path_logs(list(map(float,local)),epochs)),
                    list(map(float,exact)),rtol=1e-12,atol=1e-15)
        logs=zero_path_logs([1e-100,1e-120],64)
        self.assertTrue(np.isfinite(logs).all())
        self.assertAlmostEqual(logs[0],64*np.log(1e-100),places=8)
        for values in ([-1,1],[1,np.nan],[1,np.inf],[[1,2]]):
            with self.assertRaises(ValueError):zero_path_logs(values,2)


if __name__=='__main__':unittest.main()
