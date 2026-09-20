"""Exercise exact integer-simplex coverage, including malformed covers."""
import itertools
import unittest

import type_box_coverage as coverage


class CoverageTests(unittest.TestCase):
    def test_count_against_enumeration(self):
        for dimensions in (2,3,4):
            for total in range(6):
                for lower in itertools.product((0,1),repeat=dimensions):
                    upper = [a+(j%3) for j,a in enumerate(lower)]
                    exact = sum(sum(row) == total for row in itertools.product(
                        *(range(a,b+1) for a,b in zip(lower,upper))))
                    self.assertEqual(coverage.lattice_count(list(lower),upper,total),exact)

    def test_partition(self):
        boxes = [dict(lower=[0,0,0],upper=[1,6,6]),
                 dict(lower=[2,0,0],upper=[4,6,6])]
        self.assertEqual(coverage.check(boxes,6,2,3),25)

    def test_reject_holes_overlap_and_outside(self):
        box = dict(lower=[0,0,0],upper=[4,6,6])
        for boxes in ([dict(lower=[0,0,0],upper=[3,6,6])],
                      [box,box], [dict(lower=[0,0,0],upper=[5,6,6])]):
            with self.assertRaises(ValueError):
                coverage.check(boxes,6,2,3)

    def test_real_overlap_without_integer_overlap_is_allowed(self):
        # These rectangular ranges intersect, but their intersection does
        # not contain a point on the required integer-count simplex.
        boxes = [dict(lower=[0,0],upper=[1,1]),dict(lower=[0,1],upper=[0,2])]
        self.assertEqual(coverage.check(boxes,2,1,2),2)


if __name__ == '__main__':
    unittest.main()
