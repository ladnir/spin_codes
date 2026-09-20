#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest

from flint import arb, ctx

from exact_row_certificate import (
    accumulator_slice_tail_bound_arb,
    verify_exact_row_parameters,
)


def exact_slice_tail(n: int, cutoff: int, u: int) -> arb:
    if u == 0:
        return arb(1)
    numerator = sum(
        math.comb(cutoff, j) * math.comb(n - cutoff, u - j)
        for j in range((u + 1) // 2, min(u, cutoff) + 1)
    )
    return arb(numerator) / math.comb(n, u)


def normalized_krawtchouk(n: int, d: int, j: int) -> float:
    lo = max(0, d - (n - j))
    hi = min(d, j)
    numerator = sum(
        (-1) ** t * math.comb(j, t) * math.comb(n - j, d - t)
        for t in range(lo, hi + 1)
    )
    return numerator / math.comb(n, d)


class ExactRowCertificateTests(unittest.TestCase):
    def setUp(self) -> None:
        ctx.prec = 192

    def test_geometric_slice_bound_dominates_exact_tail(self) -> None:
        for n in range(4, 18):
            cutoff = max(1, n // 4)
            for u in range(n + 1):
                bound = accumulator_slice_tail_bound_arb(n, cutoff, u)
                exact = exact_slice_tail(n, cutoff, u)
                self.assertTrue(bound.contains(exact) or bound > exact)

    def test_krawtchouk_variance_bound_on_small_spheres(self) -> None:
        for n in range(2, 24):
            for d in range(1, n):
                for j in range(n + 1):
                    value = normalized_krawtchouk(n, d, j)
                    variance = d * (j / n) * (1 - j / n) * (n - d) / (n - 1)
                    self.assertLessEqual(value * value, math.exp(-4 * variance) + 1e-14)

    def test_three_checked_parameter_sets_prove_twenty_bits(self) -> None:
        cases = ((20, 31, 7), (25, 35, 16), (30, 39, 48))
        for k_log2, row_weight, exact_limit in cases:
            with self.subTest(k_log2=k_log2):
                k = 1 << k_log2
                result = verify_exact_row_parameters(
                    k=k,
                    n=5 * k,
                    cutoff=k // 4,
                    row_weight=row_weight,
                    exact_limit=exact_limit,
                )
                self.assertTrue(result.success)
                self.assertTrue(result.security_bits > 20)
                self.assertEqual(result.largest_exact_r, 1)


if __name__ == "__main__":
    unittest.main()
