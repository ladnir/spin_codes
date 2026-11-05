import unittest
from math import log2
import numpy as np
from hd_lambda_fiber_control import active_norm,support_norm
from hd_lambda_fiber_outer import column_fiber_counts


class FiberControlTests(unittest.TestCase):
    def test_all_one_atoms_preserve_exact_mass(self):
        table=np.zeros((65536,2,2))
        for source in ((0,1,2,3,4),(0,1,2,4,8),tuple(range(8))):
            data=column_fiber_counts(source,7)
            actual=active_norm(data,table,table,18)
            np.testing.assert_allclose(actual,log2(256**len(source)-1),rtol=0,atol=3e-14)

    def test_empty_active_factor_is_exact_zero(self):
        table=np.zeros((1,2,2))
        self.assertTrue(np.all(np.isneginf(support_norm({0:1},table,18,omit_zero=True))))
        self.assertTrue(np.all(active_norm(dict(outside_counts={0:1},
            zero_fiber_counts={0:1}),table,table,18)==-np.inf))

    def test_nonunit_empty_mass_not_dropped(self):
        # Constant monomial norm on a finite measure of mass16 is 16**(1/p).
        table=np.full((1,2,2),4.)
        np.testing.assert_allclose(support_norm({0:1},table,8),0.5)
        with self.assertRaises(ValueError):
            active_norm(dict(outside_counts={0:2},zero_fiber_counts={0:1}),table,table,4)


if __name__=='__main__': unittest.main()
