"""Exact finite checks for the accumulator identities in the BA tail proof.

These regression tests supplement, not replace, the analytic proof and the
two outward-rounded verifiers linked from artifact/PAPER_MAP.md.
"""

from fractions import Fraction
from math import comb
import unittest


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def transition_count(n, v, w):
    if v == 0:
        return int(w == 0)
    r = (v + 1) // 2
    return choose(w - 1, r - 1) * choose(n - w, v - r)


class AccumulatorTailTests(unittest.TestCase):
    def test_enumerator_against_all_short_inputs(self):
        for n in range(1, 11):
            counts = [[0] * (n + 1) for _ in range(n + 1)]
            for word in range(1 << n):
                state = weight = 0
                for i in range(n):
                    state ^= (word >> i) & 1
                    weight += state
                counts[word.bit_count()][weight] += 1
            for v in range(n + 1):
                self.assertEqual(
                    counts[v], [transition_count(n, v, w) for w in range(n + 1)]
                )

    def test_both_tail_bounds_with_exact_rationals(self):
        delta = Fraction(13, 125)
        base = Fraction(13, 28)
        for n in range(1, 97):
            cutoff = (delta * n).__floor__()
            low = [
                Fraction(
                    sum(transition_count(n, v, w) for w in range(cutoff + 1)),
                    comb(n, v),
                )
                for v in range(n + 1)
            ]
            for v, probability in enumerate(low):
                self.assertLessEqual(probability, base ** (v // 2), (n, v))
            for v in range(1, n + 1):
                high = Fraction(
                    sum(transition_count(n, v, w) for w in range(n - cutoff, n + 1)),
                    comb(n, v),
                )
                bound = Fraction(v, n - v + 1) * low[v - 1]
                if v < n:
                    bound += Fraction(n - v, v + 1) * low[v + 1]
                self.assertLessEqual(high, bound, (n, v))

    def test_parity_reflection_and_support(self):
        for n in range(1, 65):
            for v in range(1, n + 1):
                for w in range(n + 1):
                    count = transition_count(n, v, w)
                    if count:
                        self.assertGreaterEqual(w, (v + 1) // 2)
                        self.assertLessEqual(w, n - v // 2)
                    if v % 2:
                        # Odd input weight reflects about (n+1)/2, not n/2.
                        self.assertEqual(count, transition_count(n, v, n + 1 - w))
                    else:
                        self.assertEqual(
                            w * count,
                            (n - w) * transition_count(n, v, n - w),
                        )


if __name__ == "__main__":
    unittest.main()
