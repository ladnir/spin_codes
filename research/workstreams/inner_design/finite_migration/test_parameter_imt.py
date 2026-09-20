"""Validate IMT grid geometry, independent maps, and scaled coefficients."""
import math
import unittest

import numpy as np
import parameter_q1 as q1
import zero_constant_dense
from fractions import Fraction as F
from flint import arb,ctx


class ParameterTests(unittest.TestCase):
    def test_complete_geometry(self):
        rows = q1.geometries()
        self.assertEqual(len(rows),130)
        self.assertEqual(len(set(rows)),130)
        self.assertTrue(all(((1 << m)//(b//2)) % t == 0 for b,t,s,m in rows))

    def test_nested_independent_maps(self):
        for t in (64,128,256):
            small,large = q1.inner(t,t.bit_length()),q1.inner(t,t.bit_length()+1)
            mask = (1 << small['s'])-1
            for field in ('expansion_columns','feedback_columns'):
                self.assertEqual(small[field],[c & mask for c in large[field]])
            self.assertNotEqual(small['expansion_columns'],small['feedback_columns'])
            self.assertEqual(sum(small['spectrum'].values()),(1 << small['s'])-1)

    def test_scaled_coefficients_against_log_domain(self):
        rng = np.random.default_rng(937)
        for states in (3,7,11):
            rz,ra = (np.log(rng.random((5,states,states))/states) for _ in range(2))
            actual = q1.coefficients(rz,ra,16)
            expected = q1.wm.coefficients(rz,ra,16)
            np.testing.assert_allclose(actual,expected,rtol=1e-12,atol=1e-11)

    def test_actual_imt_region_coefficients(self):
        for t,s in ((64,7),(128,10),(256,12)):
            record = q1.inner(t,s)
            zero,one = q1.wm.transfers(record,np.exp(np.array([-12.,-8.,-4.])),1)
            rz,ra = q1.wm.regions(zero,one,64)
            ra -= math.log(64)
            np.testing.assert_allclose(q1.coefficients(rz,ra,32),q1.wm.coefficients(rz,ra,32),
                                       rtol=1e-11,atol=1e-9)

    def test_epoch_agrees_with_independent_map_general_transfer(self):
        search = q1.ladder.model.independent.search
        record = q1.inner(64,10)
        a,b = record['expansion_columns'],record['feedback_columns']
        low = search.low_cancellation(a,b,10,maximum=1)
        # Only degrees zero and one are evaluated. Distinct, nonzero B
        # columns make their zero-syndrome counts exactly 1 and 0.
        kernel = [1,0]
        caps = [dict(cap='0'),dict(cap='1')]
        for lam in (0.0001,0.01,0.3):
            generic = search.g.epochs(record['spectrum'],kernel,b,caps,low,lam,maximum=1)
            zero,one = q1.wm.transfers(record,np.array([lam]),1)
            np.testing.assert_allclose(zero[0],generic[0],rtol=1e-12,atol=1e-12)
            np.testing.assert_allclose(one[0],generic[1],rtol=1e-12,atol=1e-12)

    def test_zero_constant_point(self):
        import constant_point
        import retune_fixed_input
        ctx.prec = 256
        c = zero_constant_dense.Checker(16)
        q,coordinate = 206,F(49,256)
        witness = retune_fixed_input.proposal(c,q,coordinate,-12,-0.3)
        old = c.fixed_input_bound(q,q,coordinate,coordinate,witness)
        nu = c.density(coordinate,coordinate)[0]
        ratio = constant_point.constant_density.factor(c.rows,q,q,nu,nu,c.ps,F(0))
        original = c.comparison(q,q,coordinate,coordinate,ctx.prec)
        ordinary = old+256*(ratio.log()-original.log())
        exceptional = constant_point.evaluate(c,q,coordinate,witness,F(1,q),F(-3335,16))
        # The zero-row-count implementation also removes the all-one term
        # from its ordinary scalar, so it can only improve this comparator.
        comparator = (ordinary.exp()+exceptional.exp()).log()
        self.assertTrue(c.split_bound(q,q,coordinate,coordinate,witness) < comparator+arb('1e-50'))


if __name__ == '__main__':
    unittest.main()
