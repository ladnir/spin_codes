"""Exact route expectation and a deterministic state-reset obstruction.

All arithmetic deciding the obstruction is rational/integer. No samples,
generated data, encoder timings, or independence assumption across blocks.
"""
from fractions import Fraction
import time


def occupancy_coefficients(numerators,denominator):
    coeff = [1]
    for p in numerators:
        if not 0<=p<=denominator:
            raise ValueError('Invalid inclusion probability')
        nxt = [0]*(len(coeff)+1)
        for j,a in enumerate(coeff):
            nxt[j] += a*(denominator-p)
            nxt[j+1] += a*p
        coeff = nxt
    total = denominator**len(numerators)
    assert sum(coeff)==total
    return coeff,total


def gate(state,intervals=128):
    # In each pair of regions, select a suffix of the first and a prefix
    # of the second. The selections meet, so there are only 128 intervals.
    numerators = [32]*108+[31]*148
    coeff,denom = occupancy_coefficients(numerators,64)
    mean = Fraction(8192*sum(max(0,j-128)*a for j,a in enumerate(coeff)),denom)
    size = 128*sum(numerators)
    reset_rank = state*intervals
    max_dimension = size//2
    lower = max(Fraction(0),(mean-reset_rank)/(max_dimension-reset_rank))
    return size,mean,reset_rank,lower


def intervals():
    """Half-open batch intervals; each touches exactly two adjacent regions."""
    widths = [32]*108+[31]*148
    return [(2*i*64+64-widths[2*i],(2*i+1)*64+widths[2*i+1]) for i in range(128)]


if __name__=='__main__':
    started = time.monotonic()
    for s in (64,128,139,192,256):
        size,mean,rank,lower = gate(s)
        if s in (64,128):
            assert size<=491*2**21//1000 and lower>Fraction(1,1000)
            assert lower>(Fraction(188,10000) if s==64 else Fraction(269,100000))
        print(f's={s}: support={size}, E[shortened dimension]={float(mean):.9f}, '
              f'reset-rank cap={rank}, rigorous failure lower={float(lower):.9f}',flush=True)
    assert len(intervals())==128
    assert sum(b-a for a,b in intervals())*128==gate(64)[0]
    assert all(b<c for (a,b),(c,d) in zip(intervals(),intervals()[1:]))
    assert intervals()[0][0]>0 and intervals()[-1][1]<16384
    assert time.monotonic()-started<10
    print('Exact threshold decisions PASS; decimal values are display only.')
