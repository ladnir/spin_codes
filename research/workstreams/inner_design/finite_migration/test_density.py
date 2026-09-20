from fractions import Fraction as F
import itertools
import math
import unittest

from flint import arb, ctx
import density_comparison as density


def pmf(probabilities):
    values = [F(1)]
    for p in probabilities:
        updated = [F(0)] * (len(values) + 1)
        for j, value in enumerate(values):
            updated[j] += value * (1 - p)
            updated[j + 1] += value * p
        values = updated
    return values


class DensityTests(unittest.TestCase):
    def test_exhaustive_point_ratios(self):
        ctx.prec = 160
        ps = (F(1, 4), F(1, 2), F(3, 4))
        for rows in (4, 7):
            for q in range(1, 4):
                for active in itertools.combinations_with_replacement(ps + (F(1),), q):
                    nu = sum(active, F(0)) / q
                    if nu == 1:
                        continue
                    theta = q * nu / rows
                    bound = density.factor(rows, q, q, nu, nu, ps)
                    for j, probability in enumerate(pmf(active)):
                        reference = math.comb(rows, j) * theta**j * (1 - theta)**(rows - j)
                        self.assertLessEqual(density.number(probability / reference), bound)

    def test_rectangle_covers_all_compositions(self):
        ctx.prec = 160
        ps = (F(1, 5), F(2, 5), F(4, 5))
        rows, lo, hi, nlo, nhi = 9, 2, 5, F(1, 3), F(5, 6)
        bound = density.factor(rows, lo, hi, nlo, nhi, ps)
        for q in range(lo, hi + 1):
            for active in itertools.combinations_with_replacement(ps + (F(1),), q):
                nu = sum(active, F(0)) / q
                if not nlo <= nu <= nhi:
                    continue
                theta = q * nu / rows
                for j, probability in enumerate(pmf(active)):
                    reference = math.comb(rows, j) * theta**j * (1 - theta)**(rows - j)
                    self.assertLessEqual(density.number(probability / reference), bound)

    def test_variance_lp_encloses_actual_compositions(self):
        ps = (F(1, 4), F(1, 2), F(3, 4))
        for active in itertools.combinations_with_replacement(ps + (F(1),), 4):
            nu = sum(active, F(0)) / len(active)
            floor = density.variance_floor(ps, nu, F(1, 3), F(5, 2))
            for x in (F(1, 3), F(1, 2), F(1), F(5, 2)):
                actual = sum((density.variance_term(p, x) for p in active), F(0)) / len(active)
                self.assertLessEqual(floor, actual)


if __name__ == '__main__':
    unittest.main()
