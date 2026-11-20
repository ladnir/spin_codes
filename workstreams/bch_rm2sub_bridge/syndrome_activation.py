"""Fourier density bound for activation from uniform fixed-weight inputs.

The original arbitrary-state activation is kept at weights 0,1,2. For j>=3,
the nonzero syndrome measure is dominated pointwise by rho_j times uniform
nonzero state. This pays its density cost on activation instead of treating
every new state as the worst A-codeword on the next epoch.
"""
import math
from fractions import Fraction as F
from flint import arb,arb_poly
import tightened_occupancy as tight
from general_occupancy import BANDS,NAME,kernel_spectrum,compositions,distribution,multiplicity
from christoffel_caps import kraw
import activation_bridge as q1


def syndrome_density(t,s,spectrum,kernel,j):
    count=math.comb(t,j);den=(1<<s)-1
    absolute=count+sum(c*abs(kraw(t,j,w)) for w,c in spectrum.items())
    fourier=F(den*absolute,(1<<s)*count)
    nonkernel=F(count-kernel.get(j,0),count)
    return min(fourier,den*nonkernel)


def epoch_matrices(t,s,spectrum,kernel,z,number,maximum):
    original=tight.epoch_matrices(t,s,spectrum,kernel,z,number,maximum)
    result=[]
    for j,row in enumerate(original):
        row=list(row)
        if j>=3:
            rho=syndrome_density(t,s,spectrum,kernel,j)
            row[1]=number(0)
            row[2]=(number(rho.numerator)/rho.denominator)*z**j
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


def polynomial_regions(t,s,spectrum,kernel,z,maximum,length=8192):
    import polynomial_regions as poly
    assert 0<=maximum<=length and length%t==0
    epoch=epoch_matrices(t,s,spectrum,kernel,z,arb,maximum)
    power=tuple(arb_poly([row[k]*math.comb(t,j) for j,row in enumerate(epoch)]) for k in range(9))
    current=tuple(arb_poly([int(i==j)]) for i in range(3) for j in range(3))
    exponent=length//t
    while exponent:
        if exponent&1:current=poly.multiply(current,power,maximum)
        exponent>>=1
        if exponent:power=poly.multiply(power,power,maximum)
    return [tuple(max(arb(0),(p[j]/math.comb(length,j)).upper()) for p in current) for j in range(maximum+1)]
