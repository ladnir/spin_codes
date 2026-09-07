from fractions import Fraction as F
import math
import unittest

import numpy as np
from flint import arb, ctx

import scalable_split_grid as fast


class SelectedGridTests(unittest.TestCase):
    def test_exact_count_recurrence(self):
        for d in range(10):
            for qs in (list(range(d, 20)), list(range(d, 20, 3))):
                self.assertEqual(list(fast.counts(20, d, qs)),
                    [math.comb(20, q)*math.comb(q, d)*12**d for q in qs])

    def test_matches_full_grid_and_replays_each_case_once(self):
        ctx.prec = 256
        lo, hi = 2, 8
        spec = dict(rows=8, cutoff=2)
        region = [tuple(arb((i+1)*(j+1))/256 for j in range(16)) for i in range(hi+1)]
        ps = [fast.base.encode(F(1, 2))]*12
        expected = fast.frozen.evaluate(region, ps, spec, lo, hi, -10)
        best, owners = fast.table(lo, hi), np.full((hi-lo+1, hi+1), -1, dtype=np.int64)
        count = fast.evaluate(region, ps, spec, lo, hi, -10, best, owners, 0)
        self.assertEqual(count, sum(q+1 for q in range(lo, hi+1)))
        self.assertEqual(fast.pack(best, lo, hi), expected)
        self.assertEqual(fast.evaluate(region, ps, spec, lo, hi, -10, best, owners, 0, True), count)
        best[0, 0] -= 1
        with self.assertRaises(ValueError):
            fast.evaluate(region, ps, spec, lo, hi, -10, best, owners, 0, True)

    def test_skip_capped_cases(self):
        lo, hi = 2, 3
        best = fast.unpack([[-80]*3, [-80]*4], lo, hi)
        owners = np.zeros_like(best)
        region = [tuple(arb(1)/256 for _ in range(16)) for _ in range(hi+1)]
        ps = [fast.base.encode(F(1, 2))]*12
        self.assertEqual(fast.evaluate(region, ps, dict(rows=4, cutoff=0), lo, hi, -10, best, owners, 1), 0)


if __name__ == '__main__':
    unittest.main()
