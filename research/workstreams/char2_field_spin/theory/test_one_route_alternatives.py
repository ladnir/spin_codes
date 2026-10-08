"""Exact small checks of one-route formulas and rank-one refresh laws."""

from collections import Counter
from fractions import Fraction as F
from itertools import product
import unittest

from batched_transfer import multiply
from one_route_alternatives import one_scan_run_upper, scan_step


class OneRouteAlternatives(unittest.TestCase):
    def test_one_scan_closed_open_run_upper(self):
        n, e, g = 8, 3, F(1, 3)
        states = {(0, 0, 0): F(1)}
        for _ in range(n):
            next_states = Counter()
            for (old, inputs, outputs), count in states.items():
                for source in (0, 1):
                    for new, returns in scan_step(old, source):
                        key = new, inputs + source, outputs + new
                        next_states[key] += count * g ** returns
            states = next_states
        for source_weight in range(1, 2 * e + 1):
            actual = sum(count for (_, inputs, outputs), count in states.items()
                         if inputs == source_weight and 1 <= outputs <= e)
            self.assertLessEqual(actual, one_scan_run_upper(n, source_weight, e, g))

    def test_exact_rank_one_projective_refresh_law(self):
        q, modulus = 4, 0x7
        mul = lambda a, b: multiply(a, b, modulus)
        vectors = list(product(range(q), repeat=2))
        nonzero = [value for value in vectors if any(value)]
        dot = lambda a, b: mul(a[0], b[0]) ^ mul(a[1], b[1])
        pairs = [(u, v) for u in nonzero for v in vectors if dot(u, v) != 1]
        self.assertEqual(len(pairs), (q * q - 1) * (q - 1) * q)
        for source in nonzero:
            ordinary, scaled = Counter(), Counter()
            for u, v in pairs:
                coefficient = dot(v, source)
                target = tuple(x ^ mul(a, coefficient) for x, a in zip(source, u))
                self.assertTrue(any(target))
                ordinary[target] += 1
                for scalar in range(1, q):
                    scaled[tuple(mul(scalar, x) for x in target)] += 1
            for target in nonzero:
                same_line = mul(source[0], target[1]) == mul(source[1], target[0])
                if target == source:
                    expected = F(q, q * q - 1)
                elif same_line:
                    expected = F(1, q * q - 1)
                else:
                    expected = F(q - 1, q * (q * q - 1))
                self.assertEqual(F(ordinary[target], len(pairs)), expected)
                scaled_expected = F(q - 1, q * (q * q - 1))
                if same_line:
                    scaled_expected += F(1, q * (q - 1))
                self.assertEqual(F(scaled[target], len(pairs) * (q - 1)), scaled_expected)

    def test_rank_one_refresh_does_not_give_full_rank_return_equations(self):
        # T-I = u v^T always has rank <=1, although T is invertible.
        q, modulus = 4, 0x7
        mul = lambda a, b: multiply(a, b, modulus)
        for u in product(range(q), repeat=2):
            if not any(u):
                continue
            for v in product(range(q), repeat=2):
                if mul(u[0], v[0]) ^ mul(u[1], v[1]) == 1:
                    continue
                difference = [[mul(a, b) for b in v] for a in u]
                self.assertEqual(mul(difference[0][0], difference[1][1])
                                 ^ mul(difference[0][1], difference[1][0]), 0)


if __name__ == '__main__':
    unittest.main()
