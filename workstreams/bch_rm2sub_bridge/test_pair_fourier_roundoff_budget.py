import unittest
from fractions import Fraction as F


class RoundoffTests(unittest.TestCase):
    def test_exact_error_budgets(self):
        u=F(1,1<<53)
        arithmetic=(1+8*u)**130*(1+4*u)**156-1
        uniform=(1+16*u)**128*(1+8*u)**130*(1+4*u)**156-1+F(1,1<<1000)
        self.assertLess(arithmetic,F(1,1<<40))
        self.assertLess(uniform,F(1,1<<40))
        self.assertLess(F(1)/(1-2047*u),1+F(1,1<<38))


if __name__=='__main__':unittest.main()
