import unittest
from itertools import product
import numpy as np
import row_mixing_screen as candidate


class RowMixingTests(unittest.TestCase):
    def test_shape_average_matches_labeled_packets(self):
        theta=[.1,.2,.3,.4]
        families=[np.ones((1,2,2))]
        for j in range(1,5):
            families.append(np.array([np.full((2,2),1+sum(s)**2)
                                      for s in sorted(candidate.shape.expected_shapes(j))]))
        averaged=candidate.average(families,theta)
        for j in range(1,5):
            exact=sum((1+sum(s)**2)*np.prod([theta[b-1] for b in s])
                      for s in product(range(1,5),repeat=j))
            np.testing.assert_allclose(averaged[j],exact,rtol=2e-14)

    def test_remove_penalty_without_reweighting_packet_law(self):
        rho=.9
        families=[np.ones((1,1,1))]
        for j in range(1,4):
            families.append(np.array([[[rho**s.count(4)]]
                                      for s in sorted(candidate.shape.expected_shapes(j))]))
        restored=candidate.unpenalized(families,rho)
        for family in restored:np.testing.assert_allclose(family,1,rtol=1e-15)

    def test_reject_invalid_probability(self):
        with self.assertRaises(ValueError):candidate.average([], [.1]*4)


if __name__=='__main__':unittest.main()
