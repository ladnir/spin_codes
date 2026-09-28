"""Outward bound for dense messages whose group histograms are typical.

Not a full SPIN certificate: the product of active-group histogram
probabilities must be >=2^(-bq) under independent uniform four-bit columns.
The other histograms and smaller group occupancies remain uncovered.
"""
import argparse
from fractions import Fraction as Q
from math import factorial,comb
from flint import arb,ctx
from joint_support import span


def histogram_probability(histogram):
    assert len(histogram)==5 and sum(histogram)==256 and min(histogram)>=0
    numerator=factorial(256)
    denominator=16**256
    for w,n in enumerate(histogram):
        numerator*=comb(4,w)**n;denominator*=factorial(n)
    return Q(numerator,denominator)


def eligible(histogram,b):
    p=histogram_probability(histogram)
    return p.numerator*(1<<b)>=p.denominator


def eligible_collection(histograms,b):
    from collections import Counter
    histograms=[tuple(row) for row in histograms]
    assert histograms and all(row[0]<256 for row in histograms)
    probability=Q(1)
    for row,count in Counter(histograms).items():
        probability*=histogram_probability(row)**count
    return probability.numerator*(1<<(b*len(histograms)))>=probability.denominator


def contribution(q,b,threshold=209715):
    dimension=1024*q
    assert 2*threshold<dimension
    # Output information set has dimension 1024q; its coordinates are
    # independent unbiased bits. The remaining coordinates can only
    # increase weight. Use the exact rational Chernoff witness z<1.
    z=arb(threshold)/(dimension-threshold)
    return (arb(comb(2048,q))*arb(2)**(512*q+b*q)
            *z**(-threshold)*((1+z)/2)**dimension)


def self_test():
    checks=0
    for basis in ([1,2,4],[3,5],[15,51,85],[0x97,0x4b,0x2d,0x1e]):
        words=span(basis);d=len(basis)
        for z in (Q(0),Q(1,2),Q(3,4),Q(1)):
            actual=sum(z**w.bit_count() for w in words)/len(words)
            assert actual<=((1+z)/2)**d
            checks+=1
    # Explicitly count small histogram classes to check the multinomial
    # and lane multiplicities used by the 256-column formula.
    from itertools import product
    from collections import Counter
    for length in (1,2,3):
        counts=Counter()
        for columns in product(range(16),repeat=length):
            counts[tuple(sum(c.bit_count()==w for c in columns) for w in range(5))]+=1
        for histogram,count in counts.items():
            numerator=factorial(length)
            denominator=1
            for w,n in enumerate(histogram):
                numerator*=comb(4,w)**n;denominator*=factorial(n)
            assert numerator==count*denominator;checks+=1
    assert eligible((16,64,96,64,16),20)
    assert not eligible((128,32,56,32,8),20)
    assert eligible_collection([(16,64,96,64,16)]*3,20)
    assert not eligible_collection([(128,32,56,32,8)]*3,20)
    print('Dense histogram exact information-set/multinomial checks:',checks,'passed',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bits',type=int,default=20)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--budget',type=int,default=60)
    args=parser.parse_args();assert args.bits>=0
    self_test();ctx.prec=args.precision
    total=arb(0);first=None;accepted=None
    for q in range(2048,409,-1):
        term=contribution(q,args.bits)
        candidate=total+term
        if candidate<arb(2)**(-args.budget):
            first=q;accepted=candidate
        total=candidate
    assert first is not None
    print('VERIFIED restricted histogram class: q from',first,'through 2048; b',args.bits,
          'precision',args.precision,'upper',accepted,'margin',-accepted.log()/arb(2).log(),flush=True)
    for q in (first-1,first,2048):
        value=contribution(q,args.bits)
        print('Occupancy',q,'log2 upper',value.log()/arb(2).log(),flush=True)
    p=histogram_probability((16,64,96,64,16))
    probability=arb(p.numerator)/p.denominator
    print('Example histogram (16,64,96,64,16): log2 probability',probability.log()/arb(2).log(),flush=True)
    print('Excluded: atypical histograms and occupancies below the stated threshold. Not full-code distance.',flush=True)


if __name__=='__main__':main()
