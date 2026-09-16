from collections import Counter
import itertools
import math
import random
import unittest

import screen_macro as macro
search=macro.search


def accumulate(x,n,transpose=False):
    state=0;out=0
    for i in (range(n-1,-1,-1) if transpose else range(n)):
        state^=(x>>i)&1;out|=state<<i
    return out


def permute(x,p,transpose=False):
    if transpose:return sum(((x>>i)&1)<<j for i,j in enumerate(p))
    return sum(((x>>j)&1)<<i for i,j in enumerate(p))


def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0


def accumulator_count(n,a,b):
    if a==0:return int(b==0)
    return choose(n-b,a//2)*choose(b-1,(a+1)//2-1)


class RedesignTests(unittest.TestCase):
    def test_factored_adjoint_and_round_order(self):
        a0=[1,2,3,4,5];s=3
        # Toy feedback may repeat; correctness identities do not require rank
        # or distinct columns. Actual candidate audits enforce those properties.
        for repeat in (2,4):
            a=a0*repeat;t=len(a);b=[1+(3*i)%7 for i in range(t)]
            rounds=[[(1,2),(3,3)],[(5,5),(2,1)],[(6,1),(7,3)]]
            def forward(words):
                state=0;out=[]
                for x,factors in zip(words,rounds):
                    out.append(x^search.image(a,state))
                    for u,v in factors:state=search.g.tv.update(state,u,v)
                    state^=search.inject(b,x)
                return out
            def transpose(words):
                state=0;out=[]
                for x,factors in zip(words[::-1],rounds[::-1]):
                    out.append(x^search.image(b,state))
                    for u,v in factors[::-1]:state=search.g.tv.update(state,v,u)
                    folded=0
                    for j in range(repeat):folded^=(x>>(5*j))&31
                    self.assertEqual(search.inject(a,x),search.inject(a0,folded))
                    state^=search.inject(a0,folded)
                return out[::-1]
            for i,j in itertools.product(range(3*t),repeat=2):
                x=[0]*3;u=[0]*3;x[i//t]=1<<(i%t);u[j//t]=1<<(j%t)
                lhs=sum((v&w).bit_count() for v,w in zip(forward(x),u))%2
                rhs=sum((v&w).bit_count() for v,w in zip(x,transpose(u)))%2
                self.assertEqual(lhs,rhs)

    def test_accumulator_enumerator(self):
        for n in range(1,11):
            exact=Counter((x.bit_count(),accumulate(x,n).bit_count()) for x in range(1<<n))
            for a,b in itertools.product(range(n+1),repeat=2):
                self.assertEqual(exact[a,b],accumulator_count(n,a,b),(n,a,b))
            for a in range(n+1):self.assertEqual(sum(accumulator_count(n,a,b) for b in range(n+1)),math.comb(n,a))

    def test_accumulator_cascade_adjoint(self):
        n=12;rng=random.Random(4)
        for passes in (2,3):
            permutations=[rng.sample(range(n),n) for _ in range(passes)]
            def forward(x):
                for p in permutations:x=accumulate(permute(x,p),n)
                return x
            def transpose(x):
                for p in permutations[::-1]:x=permute(accumulate(x,n,True),p,True)
                return x
            for i,j in itertools.product(range(n),repeat=2):
                self.assertEqual((forward(1<<i)>>j)&1,(transpose(1<<j)>>i)&1)

    def test_macro_spectrum_wider_than_128(self):
        a=[1,2,3,4,5]*52
        expected=Counter(search.image(a,q).bit_count() for q in range(8))
        self.assertEqual(macro.spectrum(a,3),dict(expected))


if __name__=='__main__':unittest.main()
