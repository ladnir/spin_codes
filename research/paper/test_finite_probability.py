"""Exact small-instance checks of finite-appendix probability reductions.

These regression checks do not replace the analytic proof or replay certificates.
The transfer matrices are deliberately noncommuting synthetic matrices.
"""
from fractions import Fraction as F
from itertools import product
from math import comb, isqrt
import unittest


ZERO = ((F(0), F(0)), (F(0), F(0)))
IDENTITY = ((F(1), F(0)), (F(0), F(1)))


def add(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(2)) for i in range(2))


def scale(a, x):
    return tuple(tuple(x * v for v in row) for row in a)


def mul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(2))
                       for j in range(2)) for i in range(2))


def power(a, n):
    result = IDENTITY
    for _ in range(n):
        result = mul(result, a)
    return result


def poisson_binomial(probabilities):
    masses = [F(1)]
    for p in probabilities:
        nxt = [F(0)] * (len(masses) + 1)
        for j, mass in enumerate(masses):
            nxt[j] += (1-p)*mass
            nxt[j+1] += p*mass
        masses = nxt
    return masses


class FiniteProbabilityTests(unittest.TestCase):
    def test_uniform_slice_keeps_matrix_order(self):
        width, epochs = 2, 3
        transfers = [((F(1), F(j+1, 7)), (F(j, 9), F(2, 3)))
                     for j in range(width+1)]
        self.assertNotEqual(mul(transfers[0], transfers[1]),
                            mul(transfers[1], transfers[0]))
        coefficients = [IDENTITY]
        for _ in range(epochs):
            nxt = [ZERO] * (len(coefficients) + width)
            for a, previous in enumerate(coefficients):
                for b, transfer in enumerate(transfers):
                    nxt[a+b] = add(nxt[a+b], scale(mul(previous, transfer), comb(width, b)))
            coefficients = nxt
        enumerated = [ZERO] * (width*epochs+1)
        for bits in product((0, 1), repeat=width*epochs):
            value = IDENTITY
            for e in range(epochs):
                value = mul(value, transfers[sum(bits[e*width:(e+1)*width])])
            enumerated[sum(bits)] = add(enumerated[sum(bits)], value)
        for j in range(width*epochs+1):
            self.assertEqual(scale(coefficients[j], F(1, comb(width*epochs, j))),
                             scale(enumerated[j], F(1, comb(width*epochs, j))))

    def test_shuffled_bernoulli_factor_including_endpoints(self):
        grid = (F(0), F(1, 4), F(1, 2), F(1))
        for n in range(1, 6):
            ceil_sqrt = isqrt(n) + (isqrt(n)**2 != n)
            factor = min(n+1, 4*ceil_sqrt)
            for probabilities in product(grid, repeat=n):
                mean = sum(probabilities)/n
                for j, actual in enumerate(poisson_binomial(probabilities)):
                    reference = comb(n, j)*mean**j*(1-mean)**(n-j)
                    self.assertLessEqual(actual, factor*reference)

    def test_sparse_triangle_including_all_one_label(self):
        labels = ((F(1, 4), F(3, 2)), (F(2, 3), F(5, 4)), (F(1), F(1)))
        for q in range(5):
            regions = [((F(1, j+1), F(j+1, 8)), (F(1, 5), F(2, j+3)))
                       for j in range(q+1)]
            envelope = regions[:]
            for r in range(q):
                nxt = []
                for j in range(q-r):
                    candidates = [scale(add(scale(envelope[j], 1-p),
                                            scale(envelope[j+1], p)), rho)
                                  for p, rho in labels]
                    nxt.append(tuple(tuple(max(a[i][k] for a in candidates)
                                           for k in range(2)) for i in range(2)))
                envelope = nxt
            actual_sum = ZERO
            for assignment in product(labels, repeat=q):
                law = poisson_binomial([p for p, _ in assignment])
                rho = F(1)
                for _, factor in assignment:
                    rho *= factor
                region = ZERO
                for j, probability in enumerate(law):
                    region = add(region, scale(regions[j], rho*probability))
                actual_sum = add(actual_sum, power(region, 3))
            bound = scale(power(envelope[0], 3), len(labels)**q)
            for i, j in product(range(2), repeat=2):
                self.assertLessEqual(actual_sum[i][j], bound[i][j])

    def test_exponential_tilt_ratio(self):
        probabilities = (F(0), F(1, 4), F(1, 2), F(3, 4), F(1))
        n = len(probabilities)
        theta = sum(probabilities)/n
        original = poisson_binomial(probabilities)
        for j in range(1, n):
            y = F(j, n)
            a = y*(1-theta)/(theta*(1-y))
            tilted = [p*a/(1-p+p*a) for p in probabilities]
            normalizer = F(1)
            for p in probabilities:
                normalizer *= 1-p+p*a
            normalizer /= (1-theta+theta*a)**n
            reference = comb(n, j)*theta**j*(1-theta)**(n-j)
            at_mean = comb(n, j)*y**j*(1-y)**(n-j)
            self.assertEqual(original[j]/reference,
                             normalizer*poisson_binomial(tilted)[j]/at_mean)


if __name__ == '__main__':
    unittest.main()
