"""Independent exhaustive/seeded checks of exact midpoint box sums."""

from random import Random
import unittest

import hd_field128_gamma_midpoint_boxes as mode
from hd_shifted_pairs_count import exact_profiles


class BoxTests(unittest.TestCase):
    def test_all_tiny_rectangles_and_supports(self):
        for p in range(1, 8):
            oracle = exact_profiles(p)
            for r in range(2*p+1):
                boxes = mode.MidpointBoxes(p, r)
                lo, hi = boxes.counts.full.start, boxes.counts.full.stop-1
                for al in range(lo, hi+1):
                    for ah in range(al, hi+1):
                        for bl in range(lo, hi+1):
                            for bh in range(bl, hi+1):
                                expected = sum(boxes.counts.midpoint_upper(a, b)
                                    for a in range(al, ah+1) for b in range(bl, bh+1))
                                actual = sum(oracle.get((r, a, b), 0)
                                    for a in range(al, ah+1) for b in range(bl, bh+1))
                                self.assertEqual(boxes.rectangle(al, ah, bl, bh), expected)
                                self.assertGreaterEqual(expected, actual)

    def test_random_larger_boxes_partition_and_complement(self):
        random = Random(20261010)
        for p in (8, 15, 32, 79):
            for r in (1, p//2, p-1, p, p+1, 2*p-1):
                boxes, inverse = mode.MidpointBoxes(p, r), mode.MidpointBoxes(p, 2*p-r)
                lo, hi = boxes.counts.full.start, boxes.counts.full.stop-1
                for _ in range(20):
                    al, ah = sorted(random.choices(range(lo, hi+1), k=2))
                    bl, bh = sorted(random.choices(range(lo, hi+1), k=2))
                    value = boxes.rectangle(al, ah, bl, bh)
                    expected = sum(boxes.counts.midpoint_upper(a, b)
                        for a in range(al, ah+1) for b in range(bl, bh+1))
                    self.assertEqual(value, expected)
                    mid = (al+ah)//2
                    self.assertEqual(value, boxes.rectangle(al, mid, bl, bh)+boxes.rectangle(mid+1, ah, bl, bh))
                    self.assertEqual(value, inverse.rectangle(*(x+p-r for x in (al, ah, bl, bh))))
                    self.assertEqual(value, boxes.rectangle(bl, bh, al, ah))

    def test_clipping_and_guards(self):
        boxes = mode.MidpointBoxes(8, 7)
        self.assertEqual(boxes.rectangle(-20, 20, -20, 20), boxes.rectangle(0, 3, 0, 3))
        self.assertEqual(boxes.rectangle(4, 3, 0, 3), 0)
        self.assertEqual(boxes.rectangle(0, 3, 10, 12), 0)
        with self.assertRaises(ValueError):
            boxes.rectangle(True, 3, 0, 3)
        with self.assertRaises(ValueError):
            mode.MidpointBoxes(0, 0)
        with self.assertRaises(ValueError):
            mode.MidpointBoxes(8, 6, counts=boxes.counts)


if __name__ == '__main__':
    unittest.main()
