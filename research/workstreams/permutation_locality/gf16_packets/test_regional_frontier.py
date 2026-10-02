from math import comb
import unittest
import numpy as np
import regional_frontier_probe as probe


class RegionalFrontierTests(unittest.TestCase):
    def test_binomial_weighted_placement_matches_direct_product(self):
        local=np.array([[[1.,.2],[0.,.9]],[[.5,0.],[.3,.7]]])
        for epochs in (1,2,3,5):
            for p in (.2,.5,.8):
                weights=np.array([comb(epochs,j)*p**j*(1-p)**(epochs-j) for j in range(epochs+1)])
                for regions in (1,3):
                    actual=probe.evaluate(local,[(np.log(weights),.75)],.5,epochs,regions)
                    mean=(1-p)*local[0]+p*local[1]
                    exact=(1.25+np.log(np.linalg.matrix_power(mean,epochs*regions)[0].sum()))/np.log(2)
                    self.assertAlmostEqual(actual,exact,places=12)

    def test_counterfactuals_include_class_rows_and_keep_baseline(self):
        original=np.arange(2*8*8,dtype=float).reshape(2,8,8)+1
        old=original.copy();variants=probe.variants(original,{'bits':3})
        np.testing.assert_array_equal(original,old)
        np.testing.assert_array_equal(variants['baseline'],old)
        self.assertTrue(np.all(variants['no_returns'][:,1:,0]==0))
        self.assertTrue(np.all(variants['no_births_from_zero'][:,0,1:]==0))
        self.assertTrue(np.all(variants['no_returns_or_lazy_mass_to_arbitrary'][:,1:,:2]==0))
        self.assertEqual(np.count_nonzero(variants['entirely_empty_input_path']),1)
        self.assertEqual(variants['entirely_empty_input_path'][0,0,0],1)

    def test_count_scaling_retains_tiny_products_and_structural_zeros(self):
        local=np.array([[[1.,.2],[0.,.9]],[[.5,0.],[.3,.7]]])
        weights=np.log(np.array([1,5,10,10,5,1])/32)
        basic=probe.evaluate(local,[(weights,0)],0,5,3)
        tiny=probe.evaluate(local*1e-80,[(weights,0)],0,5,3)
        self.assertAlmostEqual(tiny,basic+15*np.log2(1e-80),places=10)
        zeros=np.zeros_like(local);zeros[0,0,0]=1
        empty=probe.evaluate(zeros,[(weights,0)],0,5,3)
        self.assertAlmostEqual(empty,-15,places=12)


if __name__=='__main__':unittest.main()
