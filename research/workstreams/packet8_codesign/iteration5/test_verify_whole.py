"""Exact small checks for the independent sequential whole-code replay."""
from fractions import Fraction as F
import unittest
import numpy as np
import scaled_positive as sp
import verify_whole as v


def exact_matrix(matrix):
    return np.array([[F(float(x))*F(2)**matrix.exponent for x in row]
                     for row in matrix.value], dtype=object)


class WholeReplay(unittest.TestCase):
    def test_integer_dyadic_roundtrip_and_strict_union(self):
        value = F((1 << 79)+1593, 1 << 211)
        self.assertEqual(v.decode(v.encode(value)), value)
        endpoints = {str(q):v.encode(F(1, 64)) for q in range(1, 5)}
        decoded, total = v.check_endpoint_union(endpoints, v.encode(F(1, 16)), groups=4, target=4)
        self.assertEqual(total, F(1, 16))
        self.assertEqual(len(decoded), 4)
        with self.assertRaises(ArithmeticError):
            v.check_endpoint_union(endpoints, v.encode(F(1, 32)), groups=4, target=4)
        with self.assertRaises(ArithmeticError):
            v.check_endpoint_union(endpoints, v.encode(F(1, 16)), groups=4, target=5)
        del endpoints['3']
        with self.assertRaises(ValueError):
            v.check_endpoint_union(endpoints, v.encode(F(1, 16)), groups=4, target=4)

    def test_sequential_state_continuity_against_exact_product(self):
        matrix = sp.Scaled(np.array([[.125, .25], [.0625, .5]]))
        exact = exact_matrix(matrix)
        for regions in (1, 3, 8):
            product = np.eye(2, dtype=object)
            for _ in range(regions):
                product = product@exact
            expected = sum(product[0])
            bound = sp.scalar_fraction(v.sequential_moment(matrix, regions))
            self.assertGreaterEqual(bound, expected)
            self.assertLess(float(bound-expected), 1e-12)
        reset = sum(exact[0])**3
        self.assertNotEqual(sum((exact@exact@exact)[0]), reset)

    def test_sequential_prefactor_exact_inequality(self):
        matrix = sp.Scaled(np.array([[.125, .25], [.0625, .5]]))
        exact = exact_matrix(matrix)
        moment = sum((exact@exact@exact)[0])
        for alpha in (F(1), F(2, 5), F(1, 2)):
            bound = v.sequential_bound(matrix, 2, F(3, 4), alpha,
                groups=3, regions=3, cutoff=5, beta=F(5, 4), a=F(2, 3))
            normalized = sp.scalar_fraction(bound)/(3*moment)
            scalar = F(5, 4)**2*F(4, 3)**5*F(2, 3)**6
            self.assertGreaterEqual(normalized**alpha.denominator, scalar**alpha.numerator)

    def test_authentication_rejects_missing_inputs(self):
        with self.assertRaises(ValueError):
            v.authenticate({})
        with self.assertRaises(ArithmeticError):
            v.authenticate({'nonexistent-independent-test-source': '0'*64})


if __name__ == '__main__':
    unittest.main()

