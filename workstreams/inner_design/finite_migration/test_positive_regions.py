"""Independent log-domain checks for scaled matrix-polynomial products."""
import unittest
import numpy as np
import positive_regions
import parameter_full


class RegionTests(unittest.TestCase):
    def test_positive_toy_regions(self):
        rng = np.random.default_rng(389)
        for states in (3,7):
            epoch = np.log(rng.random((9,states,states))/states)
            actual = positive_regions.regions(epoch,8,64,8)
            expected = parameter_full.g.regions(epoch,8,64,8)
            np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-10)

    def test_real_imt_regions(self):
        record,bs,kernel,caps,low = parameter_full.prepare(64,10)
        for lam in (0.0003,0.003,0.03):
            epoch = parameter_full.g.epochs(record['spectrum'],kernel,record['feedback_columns'],
                                            caps,low,lam,maximum=8)
            actual = positive_regions.regions(epoch,64,1024,8)
            expected = parameter_full.g.regions(epoch,64,1024,8)
            np.testing.assert_allclose(actual,expected,rtol=1e-11,atol=1e-9)


if __name__ == '__main__':
    unittest.main()
