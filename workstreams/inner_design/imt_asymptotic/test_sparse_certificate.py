"""Independent regression checks for the exact polynomial sparse envelope."""
from fractions import Fraction as F
import math
import unittest

from flint import fmpq,fmpq_poly
import numpy as np
import screen
import certify_imt_sparse as exact


class SparseCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model=screen.Model()
        cls.vector,cls.rows,cls.residuals,cls.dominance=exact.build(cls.model)

    def test_all_rows_exact(self):
        for p in self.residuals:
            self.assertEqual(p[0],0)
            self.assertLess(F(exact.sign_certificate(p,strict=True)['upper']),0)

    def test_upper_envelope_matches_numeric_transfer(self):
        for x in (fmpq(1,100),fmpq(1,3),fmpq(1)):
            alpha=float(x)/10000
            vec=np.array([float(p(x)) for p in self.vector])
            actual=self.model.occupation(.8*alpha,-math.log1p(-1.6*alpha)) @ vec
            upper=np.array([float(p(x)) for p in self.rows])
            self.assertTrue(np.all(upper+1e-14>=actual))
            self.assertTrue(np.all(upper<=(1-96*alpha)*vec))

    def test_wrong_gamma_rejected(self):
        # The same positive vector does not support contraction 1-110 alpha.
        alpha=fmpq_poly([0,fmpq(1,10000)])
        with self.assertRaises(AssertionError):
            exact.sign_certificate(self.rows[0]-(1-110*alpha)*self.vector[0],strict=True)

    def test_sign_checks_are_not_point_tests(self):
        exact.sign_certificate(fmpq_poly([0,-1,1]),strict=False)
        with self.assertRaises(AssertionError):
            exact.sign_certificate(fmpq_poly([0,-1,1]),strict=True)
        with self.assertRaises(AssertionError):
            exact.sign_certificate(fmpq_poly([0,-1,2]),strict=False)


if __name__=='__main__':unittest.main()
