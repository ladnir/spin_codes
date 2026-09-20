"""Exact exhaustive checks for the shared-permutation fourth tensor."""
import itertools
import unittest
from fractions import Fraction as F


class FourthTests(unittest.TestCase):
    def test_cyclic_fixed_weight_family(self):
        n = 6
        word = (1, 1, 0, 0, 0, 0)
        family = [word[j:]+word[:j] for j in range(n)]
        p = F(sum(word), n)
        centered = [tuple(F(x)-p for x in u) for u in family]
        cov = [[sum(u[i]*u[j] for u in centered)/len(family)
                for j in range(n)] for i in range(n)]
        v = p*(1-p)
        d = [[cov[i][j]-(v if i == j else -v/(n-1))
              for j in range(n)] for i in range(n)]
        self.assertTrue(all(sum(row) == 0 for row in d))
        tau = sum(d[i][j]**2 for i in range(n) for j in range(n) if i != j)/(n*(n-1))
        self.assertGreater(tau, 0)
        patterns = [(0, 1, 0, 1), (0, 1, 0, 2), (0, 1, 2, 3), (0, 0, 1, 2)]
        expected = [tau, -tau/(n-2), 2*tau/((n-2)*(n-3)), F(0)]
        totals = [F(0)]*len(patterns)
        for perm in itertools.permutations(range(n)):
            for u, z in itertools.product(centered, repeat=2):
                for a, (i, j, k, l) in enumerate(patterns):
                    totals[a] += u[perm[i]]*u[perm[j]]*z[perm[k]]*z[perm[l]]
        count = math_factorial(n)*len(family)**2
        for total, (i, j, k, l), want in zip(totals, patterns, expected):
            product = (v if i == j else -v/(n-1))*(v if k == l else -v/(n-1))
            self.assertEqual(total/count-product, want)
        inner_second = sum(sum(x*y for x, y in zip(u, z))**2
                           for u, z in itertools.product(centered, repeat=2))/len(family)**2
        self.assertEqual(inner_second, F(n*n, n-1)*v*v+n*(n-1)*tau)


def math_factorial(n):
    import math
    return math.factorial(n)


if __name__ == '__main__':
    unittest.main()
