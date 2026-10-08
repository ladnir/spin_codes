from fractions import Fraction
from itertools import combinations
from math import comb
import unittest

from one_route_tail_obstruction import obstruction, pair_tail, row_tail


class OneRouteTailTests(unittest.TestCase):
    def test_literal_tiny_windows(self):
        n, b, minimum = 12, 4, 3
        for w in range(1, 11):
            one = pair = either = 0
            for window in combinations(range(n), w):
                counts = [sum(i//b == row for i in window) for row in range(3)]
                events = [c >= minimum for c in counts]
                one += events[0]
                pair += events[0] and events[1]
                either += any(events)
            p = row_tail(n, w, b, minimum)
            p2 = pair_tail(n, w, b, minimum)
            self.assertEqual(p, Fraction(one, comb(n, w)))
            self.assertEqual(p2, Fraction(pair, comb(n, w)))
            self.assertLessEqual(3*p-3*p2, Fraction(either, comb(n, w)))

    def test_k16_grs16_is_not_40_bit(self):
        record = obstruction()
        self.assertEqual(record['window'], 13104)
        self.assertEqual(record['target'], 13108)
        self.assertGreater(record['code_failure_lower'], Fraction(1, 100))
        self.assertLess(record['code_failure_lower'], Fraction(1, 10))

    def test_whole_step_window_is_below_target(self):
        for b in (8, 16, 32, 64):
            for step in (8, 16, 32):
                result = obstruction(block=b, dimension=b//2, step=step)
                self.assertEqual(result['window'] % step, 0)
                self.assertLess(result['window'], result['target'])
                self.assertLessEqual(result['code_failure_lower'], 1)


if __name__ == '__main__':
    unittest.main()
