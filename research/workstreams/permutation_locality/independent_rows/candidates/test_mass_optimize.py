"""Check the proposal optimizer's derivatives and coupled family assembly."""
import unittest
import numpy as np

from mass_optimize import family, objective


class FixedMixtures(unittest.TestCase):
    def setUp(self):
        self.base=np.arange(3*4*4,dtype=float).reshape(3,4,4)/100+.2
        self.zero=self.base.copy(); self.zero[1:,:,0]*=.7
        self.density=self.base.copy(); self.density[1:,:,3]*=1.3

    def test_column_endpoints_and_scope(self):
        result=family(self.base,self.zero,self.density,np.array([[1.,0.],[0.,1.]]))
        expected=self.base.copy(); expected[1,:,0]=self.zero[1,:,0]; expected[2,:,3]=self.density[2,:,3]
        self.assertTrue(np.array_equal(result,expected))
        self.assertTrue(np.array_equal(family(self.base,self.zero,self.density,np.zeros((2,2))),self.base))
        with self.assertRaises(ValueError): family(self.base,self.zero,self.density,np.ones((2,2))*1.1)

    def test_fraction_derivative(self):
        fractions=np.full((2,2),.4)
        args=(self.base,self.zero,self.density)
        terminal=np.ones(4)
        _,gradient=objective(*args,fractions,3,.6,terminal,epochs=2,windows=2,regions=3)
        for i,j in np.ndindex(fractions.shape):
            high=fractions.copy(); low=fractions.copy(); high[i,j]+=1e-5; low[i,j]-=1e-5
            actual=(objective(*args,high,3,.6,terminal,epochs=2,windows=2,regions=3)[0]
                    -objective(*args,low,3,.6,terminal,epochs=2,windows=2,regions=3)[0])/2e-5
            self.assertAlmostEqual(actual,gradient[i,j],places=7)


if __name__=='__main__':
    unittest.main()
