"""Small exact-reference checks for the floating smaller-outer screen."""
from fractions import Fraction as F
from math import comb, log
import unittest

import numpy as np

from joint_outer_screen import bound, generic_six_density_bound, log_coefficients, packet, shortening_coefficients, shortening_cost
from rate_design_coset import mul, power, xor_sum
from return_cover_certificate import Witness
from sparse_active_interval_screen import regional_from_coefficients
from state8_density_certificate import transfer
from state8_joint_return import joint_return_bound


def multiply(left, right):
    return [[sum(left[i][k]*right[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


class JointOuterTests(unittest.TestCase):
    def test_exact_local_matches_directed_density_operator(self):
        x, y, tau = F(2, 7), F(3, 5), F(9)
        expected = [[F(0), F(0)], [F(0), F(0)]]
        for w in range(17):
            value = packet(w, y, tau)
            self.assertEqual(value[1][0], joint_return_bound(w, y, tau))
            for i in range(2):
                for j in range(2):
                    expected[i][j] += comb(16, w)*x**w*value[i][j]
        enclosed = transfer(Witness(x, y, tau, F(1)))
        for i in range(2):
            for j in range(2):
                self.assertLessEqual(F(enclosed[i][j].lo), expected[i][j])
                self.assertGreaterEqual(F(enclosed[i][j].hi), expected[i][j])

    def test_joint_never_exceeds_improved_pure(self):
        for y in (F(1), F(99,100), F(1,2), F(1,10)):
            for tau in (F(1), F(7), F(255)):
                for w in range(17):
                    joint, pure = packet(w,y,tau), packet(w,y,tau,False)
                    self.assertLessEqual(joint[1][0], pure[1][0])
                    self.assertEqual(joint[0],pure[0])
                    self.assertEqual(joint[1][1],pure[1][1])

    def test_two_packet_coefficients_and_normalization(self):
        degree, y, tau = 3, F(3,4), F(5)
        local = [[[comb(16,w)*entry for entry in row] for row in packet(w,y,tau)]
                 for w in range(degree+1)]
        reference = [[[F(0),F(0)],[F(0),F(0)]] for _ in range(degree+1)]
        for w in range(degree+1):
            for left in range(w+1):
                product = multiply(local[left],local[w-left])
                for i in range(2):
                    for j in range(2):
                        reference[w][i][j] += product[i][j]
        actual = np.exp(log_coefficients(degree,float(y),float(tau),32))
        self.assertTrue(np.allclose(actual,np.asarray(reference,dtype=float),rtol=2e-13,atol=1e-15))
        with np.errstate(divide='ignore'):
            regional = np.exp(regional_from_coefficients(np.log(actual),degree,32,log(float(tau))))
        expected = [[sum(reference[h][i][j]*F(comb(degree,h),comb(32,h))*tau**h
                         for h in range(degree+1)) for j in range(2)] for i in range(2)]
        self.assertTrue(np.allclose(regional,np.asarray(expected,dtype=float),rtol=2e-13,atol=1e-15))

    def test_narrow_gap_frozen_witness(self):
        point = (6.30781915507741,1.85171677017025)
        self.assertAlmostEqual(-bound(256,point,20,128,88),4.589434914362993,places=6)
        self.assertAlmostEqual(bound(256,point,20,128,88),
                               bound(256,point,20,128,88,False),places=10)

    def test_invalid_region_geometry(self):
        for active, rows in ((-1,32),(33,32),(2,31),(0,0)):
            with self.assertRaises(ValueError):
                log_coefficients(active,.9,7,rows)

    def test_exact_odd_dimension_shortening(self):
        coefficients = shortening_coefficients(87)
        for degree in range(87):
            self.assertEqual(xor_sum(mul(c,power(j,degree)) for j,c in enumerate(coefficients)),
                             power(87,degree))
        receipt = shortening_cost(87)
        self.assertEqual(receipt['basis_products'],5)
        self.assertEqual(receipt['coefficient_multiplicities'][1],7)
        self.assertEqual(sorted(receipt['coefficient_multiplicities'].values()),[7]+[8]*10)

    def test_generic_six_density_gap_frozen_witness(self):
        point = (-.47159007457111185,-.5857630212846654,2.3631604486331756,.21219579148909212)
        self.assertAlmostEqual(generic_six_density_bound(point,125174),5881.714303298001,places=6)


if __name__ == '__main__':
    unittest.main()
