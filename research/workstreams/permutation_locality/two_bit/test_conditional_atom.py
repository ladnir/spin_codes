from collections import Counter
from fractions import Fraction as Q
import unittest

from model import conditional_atom
from test_model import linear


class ConditionalAtomTests(unittest.TestCase):
    def test_exposed_subset_bounds_every_small_shape_and_syndrome(self):
        columns=[1,2,3,5,4,7,6,1,5,3];windows=5
        all_counts={}
        for word in range(1<<(2*windows)):
            weights=[((word>>(2*w))&3).bit_count() for w in range(windows)]
            shape=(weights.count(1),weights.count(2))
            all_counts.setdefault(shape,Counter())[linear(columns,word)]+=1
        census={shape:(None,counts[0],max((c for s,c in counts.items() if s),default=0),sum(counts.values()),{})
                for shape,counts in all_counts.items() if 1<=sum(shape)<=2}
        for shape,counts in all_counts.items():
            if not sum(shape):continue
            cap=conditional_atom(shape,census,windows)
            self.assertGreaterEqual(cap,Q(max(counts.values()),sum(counts.values())))

    def test_census_subset_can_improve_last_packet_bound(self):
        # Any census cap may be used; its validity is independently supplied
        # by the exact feedback census. Check the combinatorial factor here.
        census={(4,4):(None,1,1,1<<18,{})}
        self.assertEqual(conditional_atom((5,4),census),Q(1,1<<18)*Q(64,56))
        self.assertLess(conditional_atom((5,4),census),Q(1,112))


if __name__=='__main__':unittest.main()
