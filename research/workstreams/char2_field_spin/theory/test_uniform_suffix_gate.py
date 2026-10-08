"""Exhaustive small-layout checks for the uniform-suffix obstruction."""

from fractions import Fraction
import itertools
import math
import unittest

from uniform_suffix_gate import obstruction, screen, tail_probability, success_bound_exponent


class UniformSuffixTests(unittest.TestCase):
    def test_tail_against_exact_subset_count(self):
        for n in range(1, 11):
            for b in range(1, min(n, 4)+1):
                for e in range(n+1):
                    subsets = tuple(itertools.combinations(range(n), e))
                    for d in range(1, b+1):
                        count = sum(sum(x < b for x in subset) >= d for subset in subsets)
                        self.assertEqual(tail_probability(n, e, b, d), Fraction(count, len(subsets)))

    def test_pair_covariance_and_union_lower_bound(self):
        for b, k, rows in ((2, 1, 3), (3, 2, 3), (4, 2, 2), (4, 3, 2)):
            n, d = b*rows, b-k+1
            for e in range(n+1):
                both = any_bad = total = 0
                for subset in itertools.combinations(range(n), e):
                    counts = [0]*rows
                    for x in subset:
                        counts[x//b] += 1
                    flags = [c >= d for c in counts]
                    both += flags[0] and flags[1]
                    any_bad += any(flags)
                    total += 1
                p, lower = obstruction(rows, b, k, e)
                self.assertLessEqual(Fraction(both, total), p*p)
                self.assertLessEqual(lower, Fraction(any_bad, total))
                exponent = success_bound_exponent(rows, b, k, e)
                self.assertLessEqual(Fraction(total-any_bad, total), Fraction(n+1, 1 << exponent))

    def test_mds_shortening_witness_and_causal_support(self):
        # GRS[4,2] over GF4: f(z)=z+a has exactly the chosen root a.
        # Every chosen 3-subset is consequently an exact codeword support.
        # Exhaust all permutations and all tap coefficients of a width-one
        # recurrence; whenever a suffix contains a 3-subset, no prefix leaks.
        def mul(a, b):
            out = 0
            while b:
                if b & 1:
                    out ^= a
                b >>= 1
                a <<= 1
                if a & 4:
                    a ^= 7
            return out
        for order in itertools.permutations(range(4)):
            root = order[0]
            source = [point ^ root for point in order]
            self.assertEqual(source[0], 0)
            for coefficients in itertools.product(range(4), repeat=3):
                y = [source[0]]
                for value, a in zip(source[1:], coefficients):
                    y.append(value ^ mul(a, y[-1]))
                self.assertEqual(y[0], 0)
                self.assertLessEqual(sum(x != 0 for x in y), 3)

    def test_no_event_and_full_suffix_endpoints(self):
        self.assertEqual(obstruction(10, 16, 8, 8), (0, 0))
        self.assertEqual(obstruction(10, 16, 8, 160), (1, 1))

    def test_natural_rounding_and_exact_threshold(self):
        result = screen(1 << 16, 16, 12, Fraction(1, 10))
        self.assertEqual(result['k'], 65544)
        self.assertEqual(result['n'], 87392)
        self.assertEqual(result['forbidden_suffix'], 8739)
        self.assertTrue(result['rules_out_target'])
        self.assertEqual(result['forbidden_suffix'], math.ceil(Fraction(1, 10)*result['n'])-1)

    def test_invalid_geometry(self):
        for args in ((0, 0, 1, 1), (10, 11, 2, 1), (10, 4, 3, 4)):
            with self.assertRaises(ValueError):
                tail_probability(*args)


if __name__ == '__main__':
    unittest.main()
