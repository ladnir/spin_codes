import itertools
from fractions import Fraction as F
import unittest
from fixed_slice_additive_moments import moments, additive_mean_variance


class MomentTests(unittest.TestCase):
    def test_exhaustive(self):
        n, w, q, degree = 6, 2, 3, 3
        p = F(w, n)
        rows = [tuple(F(i in support)-p for i in range(n))
                for support in itertools.combinations(range(n), w)]
        table = moments(n, w, q, degree)
        direct = {key:F(0) for key in table}
        coefficients = [F(2), F(-3, 5), F(7, 8), F(-1, 17)]
        total = square = F(0)
        for sample in itertools.product(rows, repeat=q):
            counts = [sum(row[i] for row in sample) for i in range(n)]
            for a, b in direct:
                direct[a, b] += counts[0]**a*counts[1]**b
            value = sum(c*x**j for x in counts for j, c in enumerate(coefficients))
            total += value
            square += value*value
        count = len(rows)**q
        self.assertEqual(table, {key:value/count for key, value in direct.items()})
        mean, variance = additive_mean_variance(n, w, q, coefficients)
        self.assertEqual(mean, total/count)
        self.assertEqual(variance, square/count-mean*mean)


if __name__ == '__main__':
    unittest.main()
