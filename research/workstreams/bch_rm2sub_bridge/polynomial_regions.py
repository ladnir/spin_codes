"""Binary polynomial-matrix powering for certified region coefficients.

Arb polynomial products enclose the same positive coefficient formula as
the sequential implementation. Truncation cannot affect lower coefficients.
"""
import math
from flint import arb,arb_poly
import tightened_occupancy as tight


def multiply(left,right,maximum):
    return tuple(sum((left[3*i+k]*right[3*k+j] for k in range(3)),arb_poly()).truncate(maximum+1)
                 for i in range(3) for j in range(3))


def regions(t,s,spectrum,kernel,z,maximum,length=8192):
    assert 0<=maximum<=length and length%t==0
    epoch=tight.epoch_matrices(t,s,spectrum,kernel,z,arb,maximum)
    power=tuple(arb_poly([row[k]*math.comb(t,j) for j,row in enumerate(epoch)]) for k in range(9))
    current=tuple(arb_poly([int(i==j)]) for i in range(3) for j in range(3))
    exponent=length//t
    while exponent:
        if exponent&1:current=multiply(current,power,maximum)
        exponent>>=1
        if exponent:power=multiply(power,power,maximum)
    # Each exact coefficient is nonnegative. Preserve an explicit point
    # upper bound even if interval arithmetic gives a negative lower end.
    return [tuple(max(arb(0),(p[j]/math.comb(length,j)).upper()) for p in current) for j in range(maximum+1)]
