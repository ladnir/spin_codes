"""Independent exact matching enumeration and directed series checks."""
from collections import Counter
from fractions import Fraction as F
from math import comb
import time
from layered_matching_moments import moment_bounds,moment_log_upper,factorial_envelope,TILT
from layered_matching_chernoff import matching_envelope
from test_layered_matching_chernoff import matchings
from rational_log_bounds import log2_bounds


def run():
    began=time.monotonic();enumerated=0;factorials=0;tails=0
    for M in (2,4,6,8,10):
        mats=list(matchings(tuple(range(M))))
        for s in range(M+1):
            counts=Counter(sum(a<s and b<s for a,b in mat) for mat in mats)
            for t in (F(1),F(5,4),F(3,2),F(2)):
                exact=sum(F(c,len(mats))*t**z for z,c in counts.items())
                # Independent exact factorial-moment series.
                coefficient=F(1);series=F(1)
                for j in range(s//2):
                    coefficient*=F((s-2*j)*(s-2*j-1),2*(j+1)*(M-2*j-1))
                    series+=coefficient*(t-1)**(j+1)
                assert series==exact;factorials+=1
                lo,hi=moment_bounds(M,s,t,bits=16)
                assert F(lo,1<<16)<=exact<=F(hi,1<<16);enumerated+=1
                if t>1:
                    # Direct rational Markov checks at small field sizes.
                    for R in (1,3,7):
                        lowlog=log2_bounds(t,1<<24)[0]
                        uplog=log2_bounds(F(hi,1<<16),1<<24)[1]
                        k=(R*uplog+5*(1<<24)+lowlog-1)//lowlog
                        assert exact**R/t**k<=F(1,32);tails+=1
    for s in (1,512,2048,2415,4096):
        l96,u96=moment_bounds(8192,s,bits=96)
        l128,u128=moment_bounds(8192,s,bits=128)
        assert max(F(l96,1<<96),F(l128,1<<128))<=min(F(u96,1<<96),F(u128,1<<128))
        # The tighter interval may be checked against the certified log cap.
        assert log2_bounds(F(u128,1<<128),1<<24)[1]<=moment_log_upper(8192,s)
    old=matching_envelope(8192,194);new=factorial_envelope(8192,194)
    assert all(0<=a<=b for a,b in zip(new,old))
    assert all(new[s]==194*(s-4096)+new[8192-s] for s in range(4097,8193))
    print(dict(enumerated_moments=enumerated,exact_series_identities=factorials,
        rational_tail_checks=tails,precision_cross_checks=5,
        improved_cardinalities=sum(a<b for a,b in zip(new,old)),elapsed_seconds=time.monotonic()-began))
    assert time.monotonic()-began<60


if __name__=='__main__':run()
