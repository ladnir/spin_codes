from itertools import combinations
from math import comb
import unittest
import numpy as np
import regional_count_probe as diagnostic
import scalar_cover as sc


class RegionalCountDiagnosticTests(unittest.TestCase):
    def test_ordered_conditional_products_and_iid_mixture(self):
        local=np.array([[[.8,.1],[.2,.7]],[[.2,.5],[.1,.4]],[[.5,.1],[.3,.2]]])
        for epochs in (1,2,3):
            region=diagnostic.placement(local,epochs);slots=2*epochs
            for q in range(slots+1):
                total=np.zeros((2,2))
                for selected in combinations(range(slots),q):
                    matrix=np.eye(2)
                    for step in range(epochs):matrix=matrix@local[sum(v//2==step for v in selected)]
                    total+=matrix
                np.testing.assert_allclose(region[q],total/comb(slots,q),rtol=2e-14,atol=1e-15)
            for p in (.1,.5,.9):
                weights=np.array([comb(slots,q)*p**q*(1-p)**(slots-q) for q in range(slots+1)])
                mixed=np.tensordot(weights,region,axes=1)
                step=sum(comb(2,q)*p**q*(1-p)**(2-q)*local[q] for q in range(3))
                np.testing.assert_allclose(mixed,np.linalg.matrix_power(step,epochs),rtol=2e-14,atol=1e-15)

    def test_count_caps_cover_small_actual_families(self):
        features=[0.,.2,.5,.8,1.];active=[0,1,1,1,1]
        for probabilities in ([.2]*4,[0,.2,.8,1],[.2,.5,.5,.8],[0,0,1,1]):
            mean=np.mean(probabilities);variance=np.mean([p*(1-p) for p in probabilities])
            ratios,success=diagnostic.count_ratio(features,active,mean,(variance,variance),
                sum(p>0 for p in probabilities),4,10.)
            self.assertGreater(success,0)
            exact=np.array([1.])
            for p in probabilities:exact=np.convolve(exact,[1-p,p])
            binomial=np.array([comb(4,q)*mean**q*(1-mean)**(4-q) for q in range(5)])
            self.assertTrue(np.all(np.exp(ratios)*binomial+1e-12>=exact))

    def test_tiny_regional_moments_are_scaled_before_squaring(self):
        region=np.array([np.eye(2)*1e-200,np.eye(2)*2e-200])
        matrix,shift=diagnostic.weighted_region(region,np.log([.25,.75]))
        value=sc.log_power(matrix,256)+256*shift
        self.assertTrue(np.isfinite(value))
        self.assertAlmostEqual(value,256*np.log(1.75e-200),places=8)

    def test_countwise_scales_survive_underflow_during_placement(self):
        local=np.full((3,1,1),1e-200);epochs=3;slots=2*epochs
        old=diagnostic.placement(local,epochs)
        self.assertTrue(np.all(old==0))
        region,scales=diagnostic.scaled_placement(local,epochs)
        weights=np.array([comb(slots,q)*.25**q*.75**(slots-q) for q in range(slots+1)])
        matrix,shift=diagnostic.weighted_region(region,np.log(weights)+scales)
        self.assertAlmostEqual(sc.log_power(matrix,5)+5*shift,5*epochs*np.log(1e-200),places=9)
        # A high-weight but structurally zero count cannot set the scale.
        local[1:]=0
        region,scales=diagnostic.scaled_placement(local,epochs)
        self.assertTrue(np.all(np.isneginf(scales[1:])))
        matrix,shift=diagnostic.weighted_region(region,np.log(weights)+scales)
        self.assertAlmostEqual(sc.log_power(matrix,5)+5*shift,
                               5*(epochs*np.log(1e-200)+slots*np.log(.75)),places=9)


if __name__=='__main__':unittest.main()
