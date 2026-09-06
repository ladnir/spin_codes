"""Exact small checks for the parity argument and zero-state branch."""
import itertools
import unittest
from fractions import Fraction as F


class ParityTests(unittest.TestCase):
    def test_even_convolution_return(self):
        # Arbitrary even-support row laws, including nonuniform mixtures.
        n=4;support=[x for x in range(1<<n) if x.bit_count()%2==0]
        for masses in ((1,)*8,(0,1,2,3,4,5,6,7),(0,0,0,0,0,0,0,1)):
            law={x:F(m,sum(masses)) for x,m in zip(support,masses)}
            for q in (2,4,6):
                current={0:F(1)}
                for _ in range(q):
                    nxt={x:F(0) for x in support}
                    for x,p in current.items():
                        for y,r in law.items():nxt[x^y]+=p*r
                    current=nxt
                fourier=F(0)
                for u in range(1<<n):
                    character=sum((p*(-1)**((u&x).bit_count()%2) for x,p in law.items()),F(0))
                    fourier+=character**q/(1<<n)
                self.assertEqual(current[0],fourier)
                self.assertGreaterEqual(fourier,F(1,1<<(n-1)))

    def test_kernel_region_probability(self):
        # Two epochs, B=parity on two bits. Enumerate all weight-j supports.
        valid=[x for x in range(16) if (x&3).bit_count()%2==0 and (x>>2).bit_count()%2==0]
        counts=[sum(x.bit_count()==j for x in valid) for j in range(5)]
        self.assertEqual(counts,[1,0,2,0,1]) # (1+u^2)^2
        for x in valid:
            state=0;output=0
            for i in range(2):
                epoch=(x>>(2*i))&3
                output|=(epoch^(3*state))<<(2*i)
                state=state^(epoch.bit_count()%2)
            self.assertEqual(output,x)
            self.assertEqual(state,0)


if __name__=='__main__':unittest.main()
