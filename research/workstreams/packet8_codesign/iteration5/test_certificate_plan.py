"""Tiny exact checks for CERTIFICATE_PLAN.md; not a code certificate."""
from fractions import Fraction as F
from math import nextafter, inf
import unittest


def zeros(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def identity(n):
    result = zeros(n)
    for i in range(n):
        result[i][i] = F(1)
    return result


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, column)) for column in zip(*b)] for row in a]


def rebase_exact(family):
    """The plan's rational construction, independent of the float helper."""
    n = len(family[0])
    q = identity(n)
    masses = []
    for j, matrix in enumerate(family[1:], 1):
        mass = sum(matrix[0][2:])
        if mass <= 0:
            raise ValueError('positive selected upper birth mass required')
        masses.append(mass)
        q[j+1] = [F(0), F(0)]+[x/mass for x in matrix[0][2:]]
    rebased = []
    for j, matrix in enumerate(family):
        product = multiply(q, matrix)
        out = zeros(n)
        for i in range(n):
            out[i][:2] = product[i][:2]
        if j:
            out[0][j+1] = masses[j-1]
        rebased.append(out)
    return rebased, q


def upper_root(x, p, d, bits=32):
    """Exact fixed-grid ceiling of x**(p/d), using integer comparisons."""
    x = F(x)
    if x == 0:
        return F(0)
    scale = 1 << bits
    lhs_factor = x.denominator**p
    rhs = x.numerator**p * scale**d
    lo, hi = 0, scale
    while hi**d * lhs_factor < rhs:
        hi *= 2
    while hi-lo > 1:
        mid = (lo+hi)//2
        if mid**d * lhs_factor >= rhs:
            hi = mid
        else:
            lo = mid
    return F(hi, scale)


class CertificatePlan(unittest.TestCase):
    def family(self, duplicate=False):
        family = []
        for j in range(3):
            matrix = zeros(4)
            matrix[0][0] = F(1, 2**j)
            if j:
                matrix[0][2:] = [F(1, 8), F(j if not duplicate else 1, 16)]
            for i in range(1, 4):
                matrix[i][0] = F(j*(i+1), 256)
                matrix[i][1] = F(i+j+1, 16)
            family.append(matrix)
        return family

    def test_exact_rebase_and_all_short_products(self):
        from itertools import product
        for duplicate in (False, True):
            family = self.family(duplicate)
            new, q = rebase_exact(family)
            self.assertEqual([sum(row) for row in q], [1]*4)
            self.assertEqual(q[0], [1, 0, 0, 0])
            for old, rebased in zip(family, new):
                self.assertEqual(multiply(rebased, q), multiply(q, old))
            for sequence in product(range(3), repeat=4):
                old_product = new_product = identity(4)
                for j in sequence:
                    old_product = multiply(old_product, family[j])
                    new_product = multiply(new_product, new[j])
                self.assertEqual(sum(old_product[0]), sum(new_product[0]))

    def test_exact_occupancy_scaling_commutes_with_rebase(self):
        a = F(7, 8)
        old = self.family()
        scaled = [[[x/a**j for x in row] for row in matrix] for j, matrix in enumerate(old)]
        rebased, q = rebase_exact(old)
        scaled_rebased, scaled_q = rebase_exact(scaled)
        self.assertEqual(q, scaled_q)
        self.assertEqual(scaled_rebased,
                         [[[x/a**j for x in row] for row in matrix] for j, matrix in enumerate(rebased)])

    def test_checked_roots_and_positive_tangents(self):
        for p, d in ((1, 2), (2, 5), (3, 10), (7, 20), (1, 1)):
            for t in (F(1, 2), F(3, 4), F(1), F(3, 2), F(2)):
                u = upper_root(t, p, d)
                self.assertGreaterEqual(u**d, t**p)
                self.assertLess((u-F(1, 1 << 32))**d, t**p)
                alpha = F(p, d)
                for x in (F(0), F(1, 64), F(1, 2), F(1), F(7, 4), F(8)):
                    bound = (1-alpha)*u+alpha*u*x/t
                    self.assertGreaterEqual(bound**d, x**p)

    def test_underflow_correction_not_just_relative_gamma(self):
        tau = F(1, 1 << 1074)
        exact_product = tau/F(4)
        rounded = float(tau)*.25
        self.assertEqual(rounded, 0.)
        gamma2 = F(2, (1 << 53)-2)
        upper = (F.from_float(rounded)+tau)/(1-gamma2)
        self.assertGreaterEqual(upper, exact_product)
        self.assertGreater(F.from_float(nextafter(0., inf)), 0)


if __name__ == '__main__':
    unittest.main()
