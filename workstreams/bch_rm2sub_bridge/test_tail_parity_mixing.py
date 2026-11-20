import unittest
from fractions import Fraction as F
from certify_tail_parity_mixing import kraw_values
from christoffel_caps import kraw


class MixingTests(unittest.TestCase):
    def test_integer_recurrence(self):
        for n in range(2,18):
            for w in range(n+1):
                self.assertEqual(kraw_values(n,w),[kraw(n,j,w) for j in range(n+1)])

    def test_uniform_overlap_bound(self):
        rho=F(2,5);tau=F(3,8)
        for q in range(1,30):
            for r in range(q+1):
                self.assertLessEqual(rho**(2*(q-r))*tau**r,tau**q)


if __name__=='__main__':unittest.main()
