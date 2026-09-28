from fractions import Fraction as Q
from itertools import permutations,product
from math import factorial
import unittest

from heterogeneous import density_loss,tilted_reference,capped_density_loss


class HeterogeneousTests(unittest.TestCase):
    def test_capped_loss_is_exact_maximum_over_small_profiles(self):
        for n in range(1,13):
            self.assertEqual(capped_density_loss(n,(n,n,0)),density_loss(n,2))
            for caps in ((n,n,n),(n,n,1),(n//2,n,n//3)):
                values=[]
                for a in range(n+1):
                    for b in range(n-a+1):
                        counts=[a,b,n-a-b]
                        if any(x>c for x,c in zip(counts,caps)):continue
                        value=Q(n**n,factorial(n))
                        for c in counts:value*=Q(factorial(c),c**c)
                        values.append(value)
                self.assertEqual(capped_density_loss(n,caps),max(values))

    def test_balanced_loss_dominates_every_small_profile(self):
        for n in range(1,20):
            bound=density_loss(n)
            for a in range(n+1):
                for b in range(n-a+1):
                    counts=[a,b,n-a-b]
                    actual=Q(n**n,factorial(n))
                    for c in counts:actual*=Q(factorial(c),c**c)
                    self.assertGreaterEqual(bound,actual)

    def test_shuffled_heterogeneous_law_pointwise(self):
        laws=[(Q(1),Q(0),Q(0)),(Q(1,4),Q(1,2),Q(1,4)),(Q(3,5),Q(1,10),Q(3,10))]
        for tilt in ([1,1,1],[1,Q(2,7),Q(3,11)],[1,2,3]):
            reference,factor=tilted_reference(laws,[1,1,1],tilt)
            for pattern in product(range(3),repeat=3):
                actual=Q(0)
                for order in permutations(range(3)):
                    term=Q(1)
                    for i,j in enumerate(order):term*=laws[j][pattern[i]]
                    actual+=term/6
                bound=factor
                for j in pattern:bound*=reference[j]
                self.assertGreaterEqual(bound,actual)


if __name__=='__main__':unittest.main()
