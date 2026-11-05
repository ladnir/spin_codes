from fractions import Fraction
from itertools import combinations, product
import unittest

from hd_bit_syndrome_nonquiet_expurgation import event_bound, nonquiet_array_count, rank_profile
from hd_bit_syndrome_scheduled_expurgation import EXPECTED, event_union


class NonquietExpurgationTests(unittest.TestCase):
    def test_exact_toy_array_enumeration(self):
        # Quiet subspaces require the first rank coordinates to vanish.
        q,a,ranks=4,2,(1,2,0)
        for h in range(4):
            count=0
            for flat in product(range(q),repeat=a*len(ranks)):
                nonquiet=sum(any(flat[a*j+i] for i in range(r)) for j,r in enumerate(ranks))
                count += any(flat) and nonquiet<=h
            self.assertEqual(nonquiet_array_count(ranks,a,h,q=q),count)

    def test_all_columns_and_zero_event_limits(self):
        for ranks in ((1,2,0),(3,3),(0,0)):
            a=max(1,max(ranks))
            self.assertEqual(nonquiet_array_count(ranks,a,len(ranks),q=4),4**(a*len(ranks))-1)
            self.assertEqual(nonquiet_array_count(ranks,a,0,q=4),4**sum(a-r for r in ranks)-1)

    def test_no_worse_than_padding_union(self):
        ranks=(1,2,2,3)
        for h in range(5):
            union=sum(4**(sum(3-r for r in ranks)+sum(ranks[j] for j in chosen))-1
                      for chosen in combinations(range(4),h))
            self.assertLessEqual(nonquiet_array_count(ranks,3,h,q=4),union)

    def test_zero_nonquiet_matches_frozen_expurgation(self):
        hist={1:{0:16},2:{0:120},3:{0:560},**EXPECTED}
        previous,_=event_union(hist,9)
        for a in range(4,10):
            self.assertEqual(event_bound(hist,a,0),previous[a])
            for d in hist[a]:
                self.assertEqual(sum(a-r for r in rank_profile(a,d)),d)

    def test_combined_threshold_and_next_floor_obstruction(self):
        chosen={4:4,5:4,6:3,7:2,8:1,9:0}
        total=sum((event_bound(EXPECTED,a,h) for a,h in chosen.items()),Fraction(0))
        self.assertLess(total,Fraction(1,2**45))
        self.assertGreater(event_bound(EXPECTED,6,4),Fraction(1,2**45))
        self.assertGreater(event_bound(EXPECTED,7,3),Fraction(1,2**45))

    def test_guards(self):
        with self.assertRaises(ValueError):
            nonquiet_array_count((1,),1,2)
        with self.assertRaises(ValueError):
            rank_profile(3,0)
        with self.assertRaises(ValueError):
            event_bound({6:{16:1}},6,3)


if __name__ == '__main__':
    unittest.main()
