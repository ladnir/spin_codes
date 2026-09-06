"""Exact mean/variance of the squared norm of sums of fixed-weight rows."""
from fractions import Fraction as F
import itertools
import unittest


class EnergyTests(unittest.TestCase):
    def test_exhaustive_slices(self):
        for n, w, q in [(4, 2, 2), (4, 2, 3), (6, 2, 3)]:
            p = F(w, n)
            rows = [tuple(F(i in support)-p for i in range(n))
                    for support in itertools.combinations(range(n), w)]
            total = total_squared = F(0)
            for sample in itertools.product(rows, repeat=q):
                energy = sum(sum(row[i] for row in sample)**2 for i in range(n))
                total += energy
                total_squared += energy**2
            count = len(rows)**q
            v = p*(1-p)
            self.assertEqual(total/count, q*n*v)
            self.assertEqual(total_squared/count-(total/count)**2,
                             2*q*(q-1)*F(n*n, n-1)*v*v)


if __name__ == '__main__':
    unittest.main()
