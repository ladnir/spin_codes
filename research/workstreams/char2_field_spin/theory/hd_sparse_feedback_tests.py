"""Checks for ordered positive-matrix composition used by the schedule screen."""
import unittest
import numpy as np
from hd_sparse_feedback_screen import compose


class SparseFeedbackTests(unittest.TestCase):
    def test_directed_period_evaluation_matches_pair_products(self):
        from fractions import Fraction as F
        from flint import arb,ctx
        from hd_sparse_feedback_certificate import period_coefficients
        from hd_sumrank_paired_certificate import pair_coefficients
        from flint import arb_poly
        with ctx.workprec(192):
            actual,_=period_coefficients(F(3,4),F(1),cutoff=4,maximum=4)
            left,_=pair_coefficients(F(3,4),F(1),cutoff=4,maximum=4,feedthrough='bucket4')
            right,_=pair_coefficients(F(3,4),F(1),cutoff=4,maximum=4,feedthrough='bucket4',right_feedthrough='bucket8')
            # Evaluate independently as a product of two scalar matrices.
            marker=arb(1)/3
            a=[arb_poly(row)(marker) for row in left]
            b=[arb_poly(row)(marker) for row in right]
            for r in range(2):
                for s in range(2):
                    expected=a[2*r]*b[s]+a[2*r+1]*b[2+s]
                    self.assertTrue(actual[2*r+s](marker).overlaps(expected))

    def test_direct_ordered_polynomial_product(self):
        rng=np.random.default_rng(4418)
        a=rng.integers(1,8,size=(4,2,2));b=rng.integers(1,8,size=(3,2,2))
        direct=np.zeros((6,2,2))
        for j in range(4):
            for k in range(3):direct[j+k]+=a[j]@b[k]
        np.testing.assert_allclose(np.exp2(compose(np.log2(a),np.log2(b))),direct,rtol=2e-14)
        self.assertFalse(np.allclose(compose(np.log2(a),np.log2(b)),compose(np.log2(b),np.log2(a))))

    def test_zero_entries_and_scalar_rescaling(self):
        a=np.array([[[0.,2.],[3.,0.]],[[1.,0.],[0.,4.]]])
        b=np.array([[[5.,0.],[0.,6.]],[[0.,7.],[8.,0.]]])
        with np.errstate(divide='ignore'):logs_a=np.log2(a);logs_b=np.log2(b)
        expected=compose(logs_a,logs_b)
        actual=compose(logs_a-700,logs_b+500)
        np.testing.assert_allclose(actual,expected-200,atol=2e-13)

    def test_fractional_product_majorant(self):
        rng=np.random.default_rng(4148)
        for alpha in (.03,.25,.5,1.):
            for _ in range(100):
                a=rng.uniform(0,1,(2,2));b=rng.uniform(0,1,(2,2))
                self.assertTrue(np.all((a@b)**alpha<=(a**alpha)@(b**alpha)+1e-14))


if __name__=='__main__':unittest.main()
