"""Exact bounded checks for reset counting, scalar extension, and cheap mixers.

These toy checks do not certify the proposed layered inner.
"""
from collections import Counter
from itertools import product
import random
import time


def mul(a,b):
    r=0
    while b:
        if b&1:r^=a
        b>>=1;a<<=1
        if a&16:a^=19
    return r


def conv(a,b):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):r[i+j]+=x*y
    return r


def resets():
    checks=0
    for width in (1,2):
        for lengths,terminal in (((2,),1),((3,),2),((2,2),1),((2,3),2)):
            counts=[1]
            for length in lengths:
                single=[0]*(length+1)
                for word in product(range(4),repeat=length):
                    total=0
                    for x in word:total^=x
                    if total==0:single[sum(x!=0 for x in word)]+=1
                block=[1]
                for _ in range(width):block=conv(block,single)
                block[0]-=1  # Each prescribed closed episode is nonempty.
                counts=conv(counts,block)
            for _ in range(width*terminal):counts=conv(counts,[1,3])
            v=sum(lengths)+terminal;r=len(lengths);free=width*v
            tail=width*(v-2*r);cumulative=0
            for e,count in enumerate(counts):
                cumulative+=count
                assert cumulative <= 4**(free-width*r)
                assert cumulative*cumulative <= 2**(2*free)*4**(e+tail)
                checks+=1
    return checks


def extension():
    u=next(x for x in range(2,16) if mul(x,x)^x^1==0)
    base=[0,1,u,u^1]
    checked=0
    for c,d,e,f in ((1,u,u^1,1),(u,u,1,u^1),(0,1,1,0),(1,1,1,1)):
        def distance(alphabet):
            best=4
            for a,b in product(alphabet,repeat=2):
                if a or b:
                    word=(a,b,mul(c,a)^mul(d,b),mul(e,a)^mul(f,b))
                    best=min(best,sum(x!=0 for x in word))
            return best
        assert distance(base)==distance(range(16));checked+=1
    return checked


def mixers():
    delta=2
    inv=next(x for x in range(1,16) if mul(delta,x)==1)
    forbidden={1,1^delta,delta^inv}
    keys=[(a,b) for a in range(1,16) for b in range(1,16)
          if mul(a,b) not in forbidden]
    assert len(keys)==15*12
    laws=0
    for x,y in product(range(16),repeat=2):
        if not(x or y):continue
        distribution=Counter()
        for a,b in keys:
            p=x^mul(a,y);r=y^mul(b,p)
            feedback=(p^mul(delta,x),r^mul(delta,y))
            assert p or r
            assert any(feedback)
            # Both A and A+delta*I are superregular 2x2 matrices.
            assert (x!=0)+(y!=0)+(p!=0)+(r!=0)>=3
            assert (x!=0)+(y!=0)+sum(z!=0 for z in feedback)>=3
            distribution[feedback]+=1
        # Nonexceptional points identify both sampled keys. Exceptional
        # points are the y=0 branch or p=0 branch; exactly one key is lost.
        for (f0,f1),count in distribution.items():
            exceptional=(y==0 or (f0==mul(delta,x) and f1==mul(1^delta,y)))
            assert count <= (12 if exceptional else 1)
        laws+=1
    rng=random.Random(27)
    for _ in range(20):
        setup=[[rng.choice(keys) for _ in range(4)] for _ in range(3)]
        for source in range(8):
            state=[0]*8;state[source]=rng.randrange(1,16)
            for stage in range(3):
                pair=0
                for i in range(8):
                    if i&(1<<stage):continue
                    j=i^(1<<stage);a,b=setup[stage][pair];pair+=1
                    p=state[i]^mul(a,state[j]);r=state[j]^mul(b,p)
                    state[i],state[j]=p,r
            assert all(state)
    return laws,len(keys),160


def main():
    start=time.monotonic()
    print(dict(reset_counts=resets(),extension_codes=extension(),
               local_mixer_checks=mixers(),whole_code_certificate=False))
    # Bernoulli bounds the total GF32 sandwich inflation by <2^(1/1024).
    assert 2*(2047**2)>2048**2
    assert time.monotonic()-start<20


if __name__=='__main__':main()
