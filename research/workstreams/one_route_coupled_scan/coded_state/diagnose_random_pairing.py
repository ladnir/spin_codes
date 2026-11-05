"""Exact random-matching law and inverse-count diagnostics. No code cert.

Independent uniform perfect matchings at setup, accepted SL2 arithmetic
unchanged. Standard library; serial, <=30 seconds; no data files.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
import math
import time


@lru_cache(maxsize=1024)
def matchings(n):
    assert n>=0 and n%2==0
    return math.prod(range(1,n,2))


def crossing_law(n,w):
    """Exact count by number k of active--inactive pairs."""
    assert n%2==0 and 0<=w<=n
    total=matchings(n)
    result={k:F(math.comb(w,k)*math.comb(n-w,k)*math.factorial(k)
                *matchings(w-k)*matchings(n-w-k),total)
            for k in range(w%2,min(w,n-w)+1,2)}
    assert sum(result.values())==1
    return result


def marked_law(n,w,marked):
    """Joint crossing count k and j marked inactive coordinates covered."""
    assert 0<=marked<=n-w
    result={}
    for k,prob in crossing_law(n,w).items():
        for j in range(max(0,k-(n-w-marked)),min(k,marked)+1):
            result[k,j]=prob*F(math.comb(marked,j)*math.comb(n-w-marked,k-j),
                              math.comb(n-w,k))
    assert sum(result.values())==1
    return result


def all_matchings(vertices):
    if not vertices:
        yield ();return
    first=vertices[0]
    for i in range(1,len(vertices)):
        rest=vertices[1:i]+vertices[i+1:]
        for tail in all_matchings(rest):yield ((first,vertices[i]),)+tail


def log2_fraction(x):
    """Display only; certificate arithmetic stays rational."""
    return math.log2(x.numerator)-math.log2(x.denominator)


def test():
    began=time.monotonic();cases=joint=0
    for n in (2,4,6,8,10):
        for w in range(n+1):
            observed=Counter();subsets=Counter()
            for pairs in all_matchings(tuple(range(n))):
                partners=frozenset(b if a<w else a for a,b in pairs if (a<w)!=(b<w))
                observed[len(partners)]+=1;subsets[partners]+=1;cases+=1
            law=crossing_law(n,w)
            assert {k:F(v,matchings(n)) for k,v in observed.items()}==law
            for subset,count in subsets.items():
                assert F(count,matchings(n))==law[len(subset)]/math.comb(n-w,len(subset))
            for marked in range(n-w+1):
                counts=Counter()
                for subset,count in subsets.items():
                    counts[len(subset),sum(w<=x<w+marked for x in subset)]+=count
                assert {key:F(v,matchings(n)) for key,v in counts.items()}==marked_law(n,w,marked)
                joint+=1
            assert sum(k*p for k,p in law.items())==F(w*(n-w),n-1)
            if w%2==0:
                assert law[0]==F(math.comb(n//2,w//2),math.comb(n,w))
    # Exact affine parity-check probability in the uniform-image comparison.
    # One parity row H=(1,1,1,1) over GF2: rank(F)=1 unless F is empty.
    rank_cases=0
    for fixed in range(16):
        for free in range(16):
            positions=[i for i in range(4) if free>>i&1]
            success=0
            for bits in range(1<<len(positions)):
                value=fixed&~free
                for j,i in enumerate(positions):value|=((bits>>j)&1)<<i
                success+=int(value.bit_count()%2==0)
            rank=int(bool(free));prob=F(success,1<<len(positions))
            assert prob in (0,F(1,1<<rank))
            rank_cases+=1
    assert time.monotonic()-began<30
    print(dict(matching_support_cases=cases,marked_profiles=joint,
               affine_rank_checks=rank_cases,seconds=time.monotonic()-began),flush=True)


def report():
    began=time.monotonic();n=512;q=1<<32
    for w in (1,2,16,64,128,256,384,510):
        law=crossing_law(n,w)
        # Full rank on all selected coordinates: favorable RS regime. The
        # moment still cannot replace the outer's actual regional rank.
        moment=sum((p/F(q**k) for k,p in law.items()),F(0))
        minimum=min(law)
        dominant=(law[minimum]/F(q**minimum))/moment
        print(dict(weight=w,expected_live=float(w+sum(k*p for k,p in law.items())),
                   zero_cross_log2=log2_fraction(law[0]) if 0 in law else None,
                   cross_at_most_half_mean_log2=log2_fraction(sum(p for k,p in law.items()
                       if k<=w*(n-w)/(2*(n-1)))) if w>1 else None,
                   inverse_extra_rank_bits=-log2_fraction(moment),
                   minimum_cross_fraction_of_inverse_moment=float(dominant),
                   diagnostic_only=True),flush=True)
        if time.monotonic()-began>30:raise TimeoutError('random pairing report 30s cap')
    # One specified birth set outside H must be reached by live pairs if
    # those input coordinates are constrained to zero by the outer support.
    for marked in (1,16,64):
        law=marked_law(n,256,marked)
        probability=sum(p for (k,j),p in law.items() if j==marked)
        print(dict(active=256,required_births=marked,
                   all_births_hit_log2=log2_fraction(probability)),flush=True)
    print(dict(seconds=time.monotonic()-began,whole_code_certificate=False),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--test',action='store_true')
    args=parser.parse_args()
    if args.test:test()
    else:report()
