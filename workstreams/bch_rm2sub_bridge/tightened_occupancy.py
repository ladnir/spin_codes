"""All-weight envelope with both valid termination-mass bounds intersected.

This is a new implementation: retained general_occupancy certificates are
unchanged. The minimum is taken between upper bounds on the same weighted
termination mass, not between differently conditioned probability laws.
"""
import math
from fractions import Fraction as F

from flint import arb
import general_occupancy as original
from general_occupancy import BANDS, NAME, kernel_spectrum, compositions, distribution, multiplicity
import activation_bridge as q1


def upper_min(a,b,number):
    if number is arb:
        from audit_bch_q1_full_arb import rational
        value=min(rational(a.upper()),rational(b.upper()))
        return arb(value.numerator)/value.denominator
    return min(a,b)


def epoch_matrices(t,s,spectrum,kernel,z,number,maximum):
    matrices=original.epoch_matrices(t,s,spectrum,kernel,z,number,maximum)
    den=(1<<s)-1
    result=[]
    for row in matrices:
        row=list(row)
        # For each entering state, alpha*q+B*x vanishes with probability
        # 1/M on B*x!=0. Its weighted mass is at most both that event's
        # probability and the full (unrestricted) emission moment, /M.
        row[3]=upper_min(row[3],row[5]/den,number)
        row[6]=upper_min(row[6],row[8]/den,number)
        result.append(tuple(row))
    return result


def regions(t,s,spectrum,kernel,z,number,maximum,length=8192):
    assert 0<=maximum<=length and length%t==0
    epoch=epoch_matrices(t,s,spectrum,kernel,z,number,maximum)
    weighted=[tuple(v*math.comb(t,j) for v in m) for j,m in enumerate(epoch)]
    current=[tuple(number(int(i==j)) for i in range(3) for j in range(3))]
    for step in range(length//t):
        limit=min(maximum,(step+1)*t);updated=[]
        for degree in range(limit+1):
            terms=[q1.positive_mul(current[degree-j],weighted[j])
                   for j in range(max(0,degree-len(current)+1),min(degree,len(weighted)-1)+1)]
            updated.append(tuple(sum((v[k] for v in terms),number(0)) for k in range(9)))
        current=updated
    return [tuple(v/math.comb(length,j) for v in current[j]) for j in range(maximum+1)]
