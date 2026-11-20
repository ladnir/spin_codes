import unittest
from collections import Counter
import numpy as np
from inner_pair_spectrum import pair_counts,walsh


class PairTests(unittest.TestCase):
    def test_transform_involution(self):
        for bits in range(1,9):
            n=1<<bits;x=np.arange(3*n,dtype=np.int64).reshape(3,n)-n
            np.testing.assert_array_equal(walsh(walsh(x)),n*x)

    def test_exhaustive_pairs(self):
        for generators in ([3,5],[15,51,85],[255,15,51,85]):
            words=[0]
            for g in generators:words += [w^g for w in words]
            expected=Counter((x.bit_count(),y.bit_count(),(x^y).bit_count()) for x in words for y in words)
            actual,marginal=pair_counts(generators)
            self.assertEqual(actual,expected)
            self.assertEqual(marginal,Counter(w.bit_count() for w in words))


if __name__=='__main__':unittest.main()
