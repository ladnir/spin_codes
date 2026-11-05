"""Exact turnover-history and conserved-route checks; serial <=30s."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
import math
import time


def multiply(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


def bernoulli_law(probs):
    out=[F(1)]
    for p in probs:out=multiply(out,[1-p,p])
    return out


def run():
    began=time.monotonic();histories=0;coefficient_cases=0
    for w in (1,2,3):
        masks=[sum(1<<i for i in c) for c in combinations(range(4),w)]
        S=[F(math.comb(w,b)*math.comb(4-w,b)) for b in range(min(w,4-w)+1)]
        polynomial=[F(1)]
        for T in range(1,5):
            polynomial=multiply(polynomial,S);observed=Counter()
            for states in product(masks,repeat=T+1):
                B=sum((b&~a).bit_count() for a,b in zip(states,states[1:]))
                area=states[0].bit_count()+sum((a|b).bit_count() for a,b in zip(states,states[1:]))
                common=sum((a&b).bit_count() for a,b in zip(states,states[1:]))
                assert area==(T+1)*w+B and common==T*w-B
                assert area-common==w+2*B
                observed[B]+=1;histories+=1
            assert [observed[b] for b in range(len(polynomial))]==[len(masks)*v for v in polynomial]
            for x in (F(1,4),F(1,2)):
                value=len(masks)*sum((v*x**j for j,v in enumerate(S)),F(0))**T
                for cap in range(len(polynomial)):
                    assert sum(observed[b] for b in range(cap+1))<=value/x**cap
            coefficient_cases+=1
    smoothing=0
    for probs in product((F(1,4),F(1,2),F(3,4)),repeat=4):
        law=bernoulli_law(probs);mean=sum(probs,F(0))/4
        envelope=bernoulli_law([mean]*4)
        for c in range(5):
            values=[F(3,2)**max(0,b-c)*F(4,3)**(2-b) for b in range(5)]
            assert sum(x*y for x,y in zip(law,values))<=sum(x*y for x,y in zip(envelope,values))
            smoothing+=1
    conditioned=0
    for sizes in ((1,2,3),(2,2,2),(0,2,4),(1,1,4)):
        choices=[list(combinations(range(4),a)) for a in sizes]
        probs=[F(a,4) for a in sizes];area=sum(sizes);center=area//4
        probability=math.prod(math.comb(4,a)*p**a*(1-p)**(4-a) for a,p in zip(sizes,probs))
        law=bernoulli_law(probs)
        for c in range(4):
            for r,t in ((F(3,2),F(1)),(F(3,2),F(4,3)),(F(2),F(3,2))):
                total=F(0)
                for regions in product(*choices):
                    counts=[sum(i in region for region in regions) for i in range(4)]
                    total+=r**sum(max(0,b-c) for b in counts)
                actual=total/math.prod(len(x) for x in choices)
                moment=sum((p*r**max(0,b-c)*t**(center-b) for b,p in enumerate(law)),F(0))
                bound=moment**4*t**(area-4*center)/probability
                assert actual<=bound
                assert probability>=F(1,5)**3
                conditioned+=1
    assert time.monotonic()-began<30
    print(dict(histories=histories,coefficient_cases=coefficient_cases,
               convex_smoothing_cases=smoothing,conserved_conditioning_cases=conditioned,
               seconds=time.monotonic()-began))


if __name__=='__main__':run()
