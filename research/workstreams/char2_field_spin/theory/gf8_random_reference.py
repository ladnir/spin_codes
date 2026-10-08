"""Finite random-linear-code first-moment reference at 10% symbol distance.

The GV/random-linear-code curve is a comparison, not an upper bound for
all structured codes. The finite calculation below is an outward-rounded
sufficient condition for a uniformly random GF256 generator matrix.
"""
from fractions import Fraction as F
import json
import math
from flint import arb,ctx


def reference_rate():
    h=(.1*math.log2(255)-.1*math.log2(.1)-.9*math.log2(.9))/8
    return 1-h


def finite_margin(k,n):
    if not 1<=k<=n:
        raise ValueError('positive K<=N required')
    threshold=(n+9)//10-1
    # For z=1/2295: V_256(N,D)<=z^-D(1+255z)^N.
    # Union over all nonzero messages is bounded by256^K.
    return 8*(n-k)-(n*(arb(10)/9).log()+threshold*arb(2295).log())/arb(2).log()


def smallest_sufficient_length(k=2**20,bits=40,precision=192):
    if bits<1 or precision<80:
        raise ValueError('positive margin and adequate precision required')
    previous=ctx.prec
    try:
        ctx.prec=precision
        lo,hi=k,2*k
        while not finite_margin(k,hi).lower()>=bits:
            hi*=2
        # Strict threshold jumps every10 outputs, so the margin is not
        # monotone at single symbols. Search blocks of10, then all ten
        # positions at the first viable block and its neighbors.
        first,last=lo//10,hi//10+1
        while first<last:
            middle=(first+last)//2
            if finite_margin(k,10*middle+10).lower()>=bits:last=middle
            else:first=middle+1
        candidates=range(max(k,10*(first-1)+1),10*(first+2)+11)
        n=next(value for value in candidates if finite_margin(k,value).lower()>=bits)
        return dict(k=k,n=n,distance='1/10',minimum_distance=(n+9)//10,
                    sufficient_margin_bits=bits,requested_rate=str(F(k,n)),
                    rate_float=k/n,random_linear_asymptotic_rate=reference_rate(),
                    ensemble='uniform GF256 generator matrix',rigorous=True,
                    applies_to_spin=False)
    finally:ctx.prec=previous


if __name__=='__main__':print(json.dumps(smallest_sufficient_length()))
