"""Exact toy checks of the convex reduction, including nonmonotone weights."""
import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction as F


def product_weight(counts, h):
    return math.prod(h[k] for k in counts)


def multinomial_mean(n, balls, h):
    total = F(0)
    for draws in itertools.product(range(n), repeat=balls):
        counts = Counter(draws)
        total += product_weight([counts[i] for i in range(n)], h)
    return total/n**balls


class ConvexReductionTests(unittest.TestCase):
    def test_nonmonotone_logconvex_weights(self):
        n = 4
        # U-shaped, so a monotonicity/negative-association shortcut is invalid.
        h = [F(2)**((k-2)**2) for k in range(7)]
        subsets = list(itertools.combinations(range(n), 2))
        actual = F(0)
        for first, second in itertools.product(subsets, repeat=2):
            counts = Counter(first+second)
            actual += product_weight([counts[i] for i in range(n)], h)
        actual /= len(subsets)**2
        self.assertLessEqual(actual, multinomial_mean(n, 4, h))
        means = [multinomial_mean(n, a, h) for a in range(7)]
        self.assertTrue(all(means[a+2]-2*means[a+1]+means[a] >= 0 for a in range(5)))

    def test_endpoint_and_tail_replacement(self):
        # Base endpoint2; the true tail has mass1/10 at3, raised to1/5.
        actual = {0:F(1, 5), 1:F(2, 5), 2:F(3, 10), 3:F(1, 10)}
        mean = sum(a*p for a, p in actual.items())
        upper = {3:F(1, 5)}
        upper[2] = (mean-3*upper[3])/2
        upper[0] = 1-upper[2]-upper[3]
        self.assertTrue(all(p >= 0 for p in upper.values()))
        self.assertEqual(sum(a*p for a, p in upper.items()), mean)
        # All stop-loss functions and equal means characterize convex order.
        for t in range(4):
            self.assertLessEqual(sum(max(a-t, 0)*p for a, p in actual.items()),
                                 sum(max(a-t, 0)*p for a, p in upper.items()))

    def test_conditioned_poisson_coefficient(self):
        n, balls, degree = 3, 4, 4
        h = [F(2)**((k-2)**2) for k in range(degree+1)]
        poly = [F(1)]+[F(0)]*degree
        for _ in range(n):
            poly = [sum(poly[j]*h[k-j]/math.factorial(k-j) for j in range(k+1))
                    for k in range(degree+1)]
        exact = math.factorial(balls)*poly[balls]/n**balls
        self.assertEqual(exact, multinomial_mean(n, balls, h))
        # Positive-coefficient bound is valid even with the series truncated
        # at the target degree; higher terms cannot affect that coefficient.
        tilt = F(3, 2)
        series = sum(h[k]*tilt**k/math.factorial(k) for k in range(degree+1))
        bound = math.factorial(balls)*series**n/(n*tilt)**balls
        self.assertLessEqual(exact, bound)

    def test_optional_binomial_support_replacement(self):
        n, q, regions = 6, 3, 3
        h = [F(2)**((k-2)**2) for k in range(2*q+1)]
        values = [multinomial_mean(regions, a, h) for a in range(2*q+1)]
        # nu has equal mass at0 and2. Preserve its random jumps; comparing
        # r alone is not a license to replace all jumps by their mean.
        def conditional(r):
            return sum(F(math.comb(r, j), 2**r)*values[2*j] for j in range(r+1))
        actual = sum(F(math.comb(q, r)*math.comb(n-q, q-r), math.comb(n, q))*conditional(r)
                     for r in range(q+1))
        upper = sum(F(math.comb(q, r), 2**q)*conditional(r) for r in range(q+1))
        direct = sum(F(math.comb(q, j)*3**(q-j), 4**q)*values[2*j] for j in range(q+1))
        self.assertLessEqual(actual, upper)
        self.assertEqual(upper, direct)


if __name__ == '__main__':
    unittest.main()
