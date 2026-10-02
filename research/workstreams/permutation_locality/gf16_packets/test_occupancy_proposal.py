import unittest
import numpy as np
import sparse_cover  # Load the existing occupancy proof modules.
import occupancy_screen as screen


class OccupancyProposalTests(unittest.TestCase):
    def test_grouped_positive_products_match_sequential_convolution(self):
        rng=np.random.default_rng(280926)
        for n in (0,1,7,8,32,96,160,512):
            region=rng.random((n+1,3,3))
            for probabilities in ([0.]*n,[1.]*n,[.5]*n,[.999999]*n,
                    [(.2,.7,.95)[i%3] for i in range(n)],list(rng.random(n))):
                masses=np.array([1.])
                for p in probabilities:masses=np.convolve(masses,[1-p,p])
                expected=sum((mass*matrix for mass,matrix in zip(masses,region)),np.zeros((3,3)))
                for source in (region,list(region)):
                    np.testing.assert_allclose(screen.matrix_for_probabilities(source,probabilities),expected,rtol=3e-13,atol=2e-15)


if __name__=='__main__':unittest.main()
