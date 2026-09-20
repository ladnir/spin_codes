"""Exact pair-kernel evaluation under explicitly IID coordinate pairs.

This evaluator is NOT a replacement for the construction's shared fixed-type
region permutation. It supplies its complete weight-generating polynomial.
"""
from fractions import Fraction as F
from collections import Counter


def probability(masses,length,dimension,pairs):
    assert len(masses)==4 and all(p>=0 for p in masses) and sum(masses)==1
    p00,p01,p10,p11=masses
    r10=p00+p01-p10-p11;r01=p00-p01+p10-p11;r11=p00-p01-p10+p11
    value=F(0)
    for (a,b,c),count in pairs.items():
        n10=(a+c-b)//2;n01=(b+c-a)//2;n11=(a+b-c)//2
        assert (a+b+c)%2==0 and min(n10,n01,n11)>=0 and n10+n01+n11<=length
        value+=count*r10**n10*r01**n01*r11**n11
    return value/(1<<(2*dimension))


def single(p,spectrum,dimension):
    return sum((count*(1-2*p)**w for w,count in spectrum.items()),F(0))/(1<<dimension)
