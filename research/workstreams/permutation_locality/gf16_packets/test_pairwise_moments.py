import unittest
from fractions import Fraction as Q
from itertools import product
from math import comb
import pairwise_support  # Establish the existing proof-module search paths.
import pairwise_moments as moments
from joint_support import span


class PairwiseMomentTests(unittest.TestCase):
    def test_recurrence_and_orthogonality(self):
        for n in range(1,9):
            rows = [moments.krawtchouk4(n,u,n) for u in range(n+1)]
            for u in range(n+1):
                for j in range(n+1):
                    exact = sum((-1)**h*3**(j-h)*comb(u,h)*comb(n-u,j-h)
                                for h in range(max(0,j-(n-u)),min(u,j)+1))
                    self.assertEqual(rows[u][j],exact)
            for a,b in product(range(n+1),repeat=2):
                norm = sum(comb(n,u)*3**u*rows[u][a]*rows[u][b] for u in range(n+1))
                self.assertEqual(norm,4**n*3**a*comb(n,a) if a==b else 0)

    def test_moments_and_shell_caps_on_small_codes(self):
        for basis,n in (([3,5],4),([3,5,9],4),([15,51,85],7),([1,2,4],3)):
            words = span(basis)
            dual = [x for x in range(1<<n) if all((x&r).bit_count()%2==0 for r in basis)]
            distance = min((x.bit_count() for x in dual if x),default=n+1)
            counts = [0]*(n+1)
            for x,y in product(words,repeat=2):
                counts[(x|y).bit_count()] += 1
            for d in range(distance):
                actual = Q(sum(c*u**d for u,c in enumerate(counts)),len(words)**2)
                ideal = Q(sum(comb(n,u)*3**u*u**d for u in range(n+1)),4**n)
                self.assertEqual(actual,ideal)
            caps = moments.shell_caps(n,len(basis),distance)
            self.assertTrue(all(c<=v for c,v in zip(counts,caps)))


if __name__ == '__main__':
    unittest.main()
