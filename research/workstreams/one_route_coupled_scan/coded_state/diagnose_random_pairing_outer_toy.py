"""Full two-step comparison with a tiny repeated-row MDS outer.

Tests the pin/outer count, NOT a certificate for the RS256 construction.
q4,n4,L2; exact fractions; serial <=30s, no output files.
"""
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import permutations,product
import math
import time
from diagnose_random_pairing import all_matchings,matchings
from diagnose_uniform_pairs import field


@lru_cache(maxsize=256)
def colour_factor(n,s,pins,zeros,q):
    # Exactly s source colours in this batch. Every pin must be coloured;
    # every coloured fixed zero costs one outer dimension in the relaxation.
    free=n-pins-zeros
    if pins>s:return F(0)
    return sum((F(math.comb(zeros,j)*math.comb(free,s-pins-j),q**(pins+j))
                for j in range(max(0,s-pins-free),min(zeros,s-pins)+1)),F(0))/math.comb(n,s)


def transfer(n,s,q):
    out=[[[F(0) for weight in range(n+1)] for nxt in range(n+1)]
         for current in range(n+1)]
    for w in range(n+1):
        active=(1<<w)-1
        for matching in all_matchings(tuple(range(n))):
            live=0
            for i,j in matching:
                if active&((1<<i)|(1<<j)):live|=(1<<i)|(1<<j)
            D=live.bit_count()
            for nxt in range(1<<n):
                v=nxt.bit_count();common=(active&nxt).bit_count()
                pins=(nxt&~live).bit_count();zeros=n-(nxt|live).bit_count()
                colour=colour_factor(n,s,pins,zeros,q)
                for z in range(common+1):
                    weight=(active^nxt).bit_count()+common-z
                    value=F(q)**(v-z-D)*colour*math.comb(common,z)/matchings(n)
                    out[w][v][weight]+=value
    return out


def bound(n,s,q):
    K=transfer(n,s,q);answer=[F(0)]*(2*n+1)
    for middle in range(n+1):
        for final in range(n+1):
            for first,a in enumerate(K[0][middle]):
                for second,b in enumerate(K[middle][final]):answer[first+second]+=q**s*a*b
    return answer


def direct(n,s,q):
    # First route can be fixed by symmetry. Second independent region route
    # is uniform. Each active pair's image is uniform J^2 in this comparison.
    mul=field(q,7);answer=Counter();operations=0
    for values in product(range(1,q),repeat=s):
        message=values+(0,)*(n-s);active=(1<<s)-1
        for permutation in permutations(range(n)):
            X=tuple(message[i] for i in permutation)
            for matching in all_matchings(tuple(range(n))):
                live=[]
                for i,j in matching:
                    if active&((1<<i)|(1<<j)):live.extend((i,j))
                frequency=Counter()
                for image in product(range(q),repeat=len(live)):
                    AQ=[0]*n
                    for i,x in zip(live,image):AQ[i]=x
                    Y=tuple(AQ[i]^X[i]^mul[2][message[i]] for i in range(n))
                    frequency[s+sum(x!=0 for x in Y)]+=1;operations+=1
                factor=F(1,math.factorial(n)*matchings(n)*q**len(live))
                for weight,count in frequency.items():answer[weight]+=count*factor
    assert sum(answer.values())==(q-1)**s
    return answer,operations


def run():
    began=time.monotonic();n=q=4
    for s in (1,2):
        upper=bound(n,s,q);actual,operations=direct(n,s,q)
        for weight in range(2*n+1):assert actual[weight]<=upper[weight]
        print(dict(source_rows=s,enumerated_images=operations,
                   exact_expected_low_weight=str(sum(actual[w] for w in range(4))),
                   pin_bound_low_weight=str(sum(upper[:4])),
                   all_weight_coefficients_pass=True,whole_code_certificate=False),flush=True)
        if time.monotonic()-began>30:raise TimeoutError('outer toy 30s cap')
    print(dict(seconds=time.monotonic()-began),flush=True)


if __name__=='__main__':run()
