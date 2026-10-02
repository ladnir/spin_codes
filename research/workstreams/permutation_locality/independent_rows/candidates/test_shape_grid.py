import unittest
import numpy as np
from flint import arb_mat, ctx

import shape_grid as candidate
from mixing_attack import retarget_matrix


class ShapeGridTests(unittest.TestCase):
    def test_float_update_conversion_matches_outward_formula(self):
        ctx.prec=192
        spectrum={48:3,56:7,64:11,72:13,80:17}
        old=np.arange(1,122,dtype=float).reshape(11,11)/128
        for i in (2,9,10):old[i,3]=0
        for occupancy in (0,1):
            families=[old[None,:,:]]*(occupancy+1)
            for rounds in (2,3,8,None):
                actual=candidate.retarget(families,spectrum,'.072',rounds)[occupancy][0]
                exact=retarget_matrix(arb_mat(old.tolist()),occupancy,spectrum,'.072',2,rounds)
                reference=np.array([[float(exact[i,j]) for j in range(11)] for i in range(11)])
                np.testing.assert_allclose(actual,reference,rtol=2e-15,atol=1e-17)

    def test_shape_penalty_and_unresolved_shape_fallback(self):
        families=[np.ones((1,1,1)),np.ones((4,1,1)),np.ones((1,1,1))]
        lower=candidate.penalize(families,'.45')
        np.testing.assert_allclose(lower[1].ravel(),[1,1,1,.5])
        self.assertEqual(lower[2][0,0,0],1)
        higher=candidate.penalize(families,'1')
        self.assertAlmostEqual(higher[2][0,0,0],(1/.9)**2)
        with self.assertRaises(ValueError):candidate.penalize([families[0],np.ones((3,1,1))],'.9')


if __name__=='__main__':unittest.main()
