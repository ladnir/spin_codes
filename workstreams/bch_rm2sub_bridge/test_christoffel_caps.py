import math
import unittest
from fractions import Fraction as F
import christoffel_caps as caps


class KernelTests(unittest.TestCase):
    def test_orthogonality(self):
        for n in (4,8,12):
            for a in range(n//2+1):
                for b in range(n//2+1):
                    total=sum(math.comb(n,w)*caps.kraw(n,a,w)*caps.kraw(n,b,w) for w in range(n+1))
                    self.assertEqual(total,(1<<n)*math.comb(n,a) if a==b else 0)

    def test_reproducing_polynomial_and_even_parity_code(self):
        n=8;degree=3
        for target in range(n+1):
            values=[sum((F(caps.kraw(n,j,target)*caps.kraw(n,j,w),math.comb(n,j)) for j in range(degree+1)),F(0)) for w in range(n+1)]
            full=sum((math.comb(n,w)*v*v for w,v in enumerate(values)),F(0))/(1<<n)
            even=sum((math.comb(n,w)*values[w]**2 for w in range(0,n+1,2)),F(0))/(1<<(n-1))
            self.assertEqual(full,values[target]);self.assertEqual(even,full)
            actual=math.comb(n,target) if target%2==0 else 0
            self.assertLessEqual(actual,caps.shell_bound(n,n-1,degree,target))

    def test_caps_only_tighten(self):
        original=caps.q2.deterministic_caps();new=caps.deterministic_caps()
        self.assertEqual(set(new),set(original))
        self.assertTrue(all(0<=new[w]<=original[w] for w in new))
        self.assertTrue(all(new[w]==new[256-w] for w in range(38,130,2)))


if __name__=='__main__':unittest.main()
