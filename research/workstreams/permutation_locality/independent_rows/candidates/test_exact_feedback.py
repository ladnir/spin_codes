"""Check limb inversion against independent unbounded-integer transforms."""
import unittest
from itertools import product
import numpy as np

from exact_feedback import inverse_counts,build


def transform(values):
    return [sum(value*(-1 if (i&a).bit_count()%2 else 1) for i,value in enumerate(values))
            for a in range(len(values))]


class ExactFeedback(unittest.TestCase):
    def check(self,counts):
        coefficients=np.array(transform(counts),dtype=np.int64)
        self.assertTrue(np.array_equal(inverse_counts(coefficients,sum(counts)),counts))

    def test_all_small_distributions(self):
        for counts in product(range(3),repeat=4):
            if sum(counts): self.check(counts)

    def test_limb_reconstruction(self):
        # Each unnormalized inverse needs more than 64 signed bits.
        for location in range(8):
            counts=[0]*8; counts[location]=(1<<62)+17
            self.check(counts)
        self.check([((1<<58)+3)*a for a in (1,2,0,3,1,0,2,4)])
        self.check([((1<<57)-1)*a for a in (7,0,3,2,0,5,1,0)])

    def test_rejects_invalid_inputs(self):
        for values,denominator in (([1.,1.],1),([1,1,1],1),([2,0],1),([1,2],1),([1,0],1)):
            with self.assertRaises(ValueError): inverse_counts(np.array(values),denominator)
        with self.assertRaises(ValueError): inverse_counts(np.array([2,2,2,-2],dtype=np.int64),2)
        for maximum in (0,11,True):
            with self.assertRaises(ValueError): build(maximum)

    def test_production_size_signed_limbs(self):
        size=1<<19; denominator=(1<<61)+123
        states=np.arange(size,dtype=np.uint32)
        signs=1-2*(np.bitwise_count(states&(size-1))&1).astype(np.int64)
        counts=inverse_counts(signs*denominator,denominator)
        self.assertEqual(int(counts[-1]),denominator)
        self.assertEqual(np.count_nonzero(counts),1)


if __name__=='__main__':
    unittest.main()
