import itertools
import unittest
import inner_candidate_boxes as boxes


class CoverTests(unittest.TestCase):
    def test_count_by_enumeration(self):
        for total in range(7):
            for lo,hi in [([0,0,0],[3,4,5]),([1,0,1],[2,3,2]),([2,0,0],[1,4,4])]:
                exact=sum(sum(v)==total for v in itertools.product(*(range(a,b+1) for a,b in zip(lo,hi))))
                self.assertEqual(boxes.count_box(lo,hi,total),exact)

    def test_cover_rejects_overlap_and_gaps(self):
        a=dict(lower=[0,0,0],upper=[1,4,4])
        b=dict(lower=[2,0,0],upper=[2,4,4])
        self.assertEqual(boxes.check_cover([a,b],4,2),12)
        with self.assertRaises(AssertionError): boxes.check_cover([a],4,2)
        with self.assertRaises(AssertionError): boxes.check_cover([a,a,b],4,2)


if __name__=='__main__': unittest.main()
