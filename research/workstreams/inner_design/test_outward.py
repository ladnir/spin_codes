"""Regression checks between independently implemented envelope arithmetic.

These detect drift; the certificate, not binary64 agreement, verifies bounds.
"""
from fractions import Fraction as F
import json
import math
from pathlib import Path
import unittest

from flint import ctx
import numpy as np
import certify_no_constant as cert
import general_occupancies as g

HERE=Path(__file__).resolve().parent


class OutwardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.previous_precision=ctx.prec
        ctx.prec=256
        cls.engine=cert.Engine(json.loads((HERE/'NO_CONSTANT_MAP.json').read_text()))

    @classmethod
    def tearDownClass(cls):ctx.prec=cls.previous_precision

    def test_fixed_weight_matrices(self):
        e=self.engine
        for log_lam in (F(-2),F(-8)):
            actual=np.array([[float(x) for x in row] for row in e.epoch(log_lam)]).reshape(129,e.n,e.n)
            reference=np.exp(g.epochs(e.spectrum,e.kernel,e.columns,e.caps,e.low,math.exp(float(log_lam)),1,128))
            np.testing.assert_allclose(actual,reference,rtol=1e-10,atol=1e-14)

    def test_bernoulli_matrices(self):
        e=self.engine
        for theta in (F(3,100),F(1,10),F(1,2),F(24,25)):
            lam=cert.old.number(F(1,5))
            actual=np.array([float(x) for x in e.bernoulli(cert.old.number(theta),lam)]).reshape(e.n,e.n)
            reference=np.exp(g.bernoulli_epoch(e.spectrum,e.kernel,float(theta),.2))
            np.testing.assert_allclose(actual,reference,rtol=1e-10,atol=1e-14)

    def test_integer_coverage(self):
        for suffix,qmin in (('33_200',65),('19_100',129)):
            cover=json.loads((HERE/f'NO_CONSTANT_DENSE_COVER_{suffix}.json').read_text())
            g.check_coverage(cover['selected_boxes'],32768,qmin)
            self.assertEqual(cover['unresolved_leaves'],0)


if __name__=='__main__':unittest.main()
