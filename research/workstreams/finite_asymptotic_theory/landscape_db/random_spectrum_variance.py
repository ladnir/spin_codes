"""Simultaneous shell caps using exact first and second subspace moments.

For binary vectors, any two distinct nonzero vectors are linearly
independent. Their joint inclusion probability in one uniform D-subspace
is M(M-1)/(T(T-1)), where M=2^D-1 and T=2^B-1. Consequently a shell of
size n has the same first two counting moments as a hypergeometric sample.
No independence between shells, rows, or multiple codewords is assumed.
"""
import math
from fractions import Fraction


def shell_moments(block,dimension,weight):
    if not 0<dimension<block or not 1<=weight<=block:
        raise ValueError('invalid subspace or shell')
    total=(1<<block)-1;mass=(1<<dimension)-1;n=math.comb(block,weight)
    mean=Fraction(n*mass,total)
    variance=Fraction(n*mass*(total-mass)*(total-n),total*total*(total-1))
    return mean,variance


def ceil_sqrt_fraction(numerator,denominator):
    if numerator<0 or denominator<=0:raise ValueError('invalid nonnegative rational')
    ceiling=(numerator+denominator-1)//denominator
    root=math.isqrt(ceiling)
    return root+int(root*root<ceiling)


def caps(block,dimension,failure_bits=60):
    if not 0<dimension<block or failure_bits<0:
        raise ValueError('invalid spectrum event parameters')
    total=(1<<block)-1;mass=(1<<dimension)-1
    allocation=block*(1<<failure_bits)
    result={}
    for w in range(1,block+1):
        n=math.comb(block,w)
        mean_numerator=n*mass
        markov=mean_numerator*allocation//total
        mean_ceiling=(mean_numerator+total-1)//total
        variance_numerator=n*mass*(total-mass)*(total-n)
        variance_denominator=total*total*(total-1)
        deviation=ceil_sqrt_fraction(variance_numerator*allocation,variance_denominator)
        # If A_w exceeds this integer cap, its deviation from the mean is
        # at least sqrt(variance/delta). Chebyshev therefore costs <=delta.
        variance_cap=mean_ceiling+deviation-1 if variance_numerator else mean_ceiling
        cap=min(mass,n,markov,variance_cap)
        if cap:result[w]=cap
    return result


def evidence(block,dimension,failure_bits=60):
    return dict(block_bits=block,dimension=dimension,setup_failure_bits=failure_bits,
                counts=caps(block,dimension,failure_bits),
                method='simultaneous minimum of Markov and exact-variance Chebyshev shell caps; one reused full-rank subspace')
