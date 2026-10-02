"""Positive four-row comparison measure, not a replacement setup law."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import factorial,log
import sys

import numpy as np
import kernel

sys.path.insert(0,str(kernel.ROOT/'two_bit'))
from row_mixture import envelope
from heterogeneous import density_loss,capped_density_loss
from bch_joint_support import authenticated_caps

G=2048
REGIONS=256
PACKETS=G*REGIONS
EPOCHS=16384
N=4*PACKETS


def group_components(rows):
    result=[]
    for types in combinations_with_replacement(range(len(rows)),4):
        coefficient=Q(factorial(4))
        for count in Counter(types).values():coefficient/=factorial(count)
        probabilities=[Q(1)]
        for i in types:
            c,p=rows[i];coefficient*=c
            next_=[Q(0)]*(len(probabilities)+1)
            for j,v in enumerate(probabilities):
                next_[j]+=v*(1-p);next_[j+1]+=v*p
            probabilities=next_
        if coefficient:
            result.append((''.join(map(str,types)),coefficient,tuple(probabilities),int(any(types))))
    return result


def actual_components(central_bits=130,theta=Q(2,5)):
    return group_components(envelope(authenticated_caps(),1<<central_bits,Q(theta)))


def logq(x):return log(x.numerator)-log(x.denominator)


def log_power(matrix,n=EPOCHS):
    """Binary powering with rescaling; proposals only."""
    v=np.array([1.,0.]);logv=0.;a=matrix.copy();loga=0.
    while n:
        if n&1:
            v=v@a;s=v.max()
            if s<=0:return -np.inf
            v/=s;logv+=loga+np.log(s)
        n>>=1
        if n:
            a=a@a;s=a.max();a/=s;loga=2*loga+np.log(s)
    return logv+np.log(v.sum())
