"""Small, noncommuting derivative checks for a diagnostic, not a proof."""
import unittest

import numpy as np

from local_sensitivity import moment_gradient, sensitivity, placement_tape
from mass_density_screen import float_placement


class Adjoint(unittest.TestCase):
    def test_matrix_gradient_matches_finite_differences(self):
        matrix=np.array([[.9,.2],[.1,.7]])
        terminal=np.array([1.,.3])
        value,gradient=moment_gradient(matrix,7,terminal)
        self.assertAlmostEqual(value,np.log((np.linalg.matrix_power(matrix,7)@terminal)[0]),places=12)
        for s,t in np.ndindex(matrix.shape):
            plus=matrix.copy(); minus=matrix.copy()
            plus[s,t]+=1e-6; minus[s,t]-=1e-6
            derivative=(moment_gradient(plus,7,terminal)[0]-moment_gradient(minus,7,terminal)[0])/2e-6
            self.assertAlmostEqual(gradient[s,t],derivative,places=7)

    def test_local_noncommuting_gradient(self):
        ops=np.array([[[.8,.1],[.03,.6]],[[.4,.3],[.2,.5]],[[.2,.5],[.7,.1]]])
        terminal=np.array([1.,.7])
        for q in (1,3,6):
            self.assertTrue(np.allclose(placement_tape(ops,q,3,2)[-1],float_placement(ops,q,3,2),rtol=1e-14,atol=0))
            value,elastic=sensitivity(ops,q,.63,terminal,epochs=3,windows=2,regions=5)
            for j,s,t in np.ndindex(ops.shape):
                plus=ops.copy(); minus=ops.copy()
                plus[j,s,t]*=np.exp(1e-5); minus[j,s,t]*=np.exp(-1e-5)
                upper=sensitivity(plus,q,.63,terminal,epochs=3,windows=2,regions=5)[0]
                lower=sensitivity(minus,q,.63,terminal,epochs=3,windows=2,regions=5)[0]
                self.assertAlmostEqual(elastic[j,s,t],(upper-lower)/2e-5,places=7)
            self.assertAlmostEqual(elastic.sum(),15.,places=10)

    def test_rejects_invalid_data(self):
        with self.assertRaises(ValueError): moment_gradient(np.eye(2),0,np.ones(2))
        with self.assertRaises(ValueError): placement_tape(np.ones((1,2,2)),2,2,2)
        with self.assertRaises(ValueError): sensitivity(np.ones((3,2,2)),2,1.,np.ones(2),epochs=2,windows=2)


if __name__=='__main__':
    unittest.main()
