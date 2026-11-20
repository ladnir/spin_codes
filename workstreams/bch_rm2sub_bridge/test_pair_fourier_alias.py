import itertools
import math
import unittest
from collections import Counter
from fractions import Fraction as F


class AliasTests(unittest.TestCase):
    def test_exact_mod_two_alias_and_l1_upper(self):
        kernel=[x for x in range(16) if all((x&g).bit_count()%2==0 for g in [3,5])]
        regional=[a|(b<<4) for a in kernel for b in kernel]
        counts=Counter(tuple(sum(2*((x>>j)&1)+((y>>j)&1)==s for j in range(8)) for s in range(4))
            for x in regional for y in regional)
        p=[F(1,2),F(1,8),F(1,8),F(1,4)]
        masses={m:count*math.prod(x**n for x,n in zip(p,m)) for m,count in counts.items()}
        transform={v:sum((mass*(-1)**sum(a*b for a,b in zip(v,m[1:])) for m,mass in masses.items()),F(0))
            for v in itertools.product(range(2),repeat=3)}
        l1=sum(map(abs,transform.values()))/8
        for m,mass in masses.items():
            alias=sum((value for n,value in masses.items() if all((a-b)%2==0 for a,b in zip(m[1:],n[1:]))),F(0))
            inverse=sum((value*(-1)**sum(a*b for a,b in zip(v,m[1:])) for v,value in transform.items()),F(0))/8
            self.assertEqual(alias,inverse)
            self.assertGreaterEqual(alias,mass)
            self.assertLessEqual(alias,l1)


if __name__=='__main__':unittest.main()
