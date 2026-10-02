"""Zero-state activation split by occupancy, with a Fourier density tail.

Coordinates and all nonzero-source rows are those of rank_return.
The chosen cutoff is a proof parameter, not an encoder parameter.
"""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb
import rank_return as base
import birth_refresh_probe as proposal
import scalar_cover as sc

aq,up=base.aq,base.up


def prepare(images,columns,bits,updates=2):
    return proposal.prepare(base.prepare(images,columns,bits,updates))


def actual(updates=2):
    data=proposal.prepare(base.actual(updates))
    print('GF16 BIRTH DENSITY KERNEL',len(data['records']),'character profiles',flush=True)
    return data


def select(data,p,z,original):
    options=proposal.candidates(data,p,z,original)
    name=min(options,key=lambda key:sc.log_power(options[key]))
    return None if name=='original' else int(name.removeprefix('birth_from_j')),options[name]


def floating(data,probabilities,tilt):
    original=base.floating(data,probabilities,tilt)
    return select(data,float(1-probabilities[0]),np.exp(-tilt),original)[1]


def fourier_parts(data,p,z,degree):
    """Full character moment and its occupancy terms through degree."""
    p=Q(p);W=data['windows'];S=1<<data['bits']
    if not 0<=p<=1 or not 0<z<=1 or type(degree) is not int or not 0<=degree<=W:
        raise ValueError('valid activity, output weight, and occupancy degree required')
    values=[((1+z)**(4-r)*(1-z)**r-1)/15 for r in range(5)]
    factors=[aq(1-p)+aq(p)*value for value in values]
    powers=[[value**j for j in range(degree+1)] for value in values]
    weights=[aq(1-p)**(W-j)*aq(p)**j for j in range(degree+1)]
    fulls=[];lows=[]
    for row in data['records']:
        sums=[sum((int(n)*powers[r][i] for r,n in enumerate(row)),arb(0)) for i in range(1,degree+1)]
        coefficients=[arb(1)]
        for j in range(1,degree+1):
            coefficients.append(sum(((-1)**(i-1)*coefficients[j-i]*sums[i-1] for i in range(1,j+1)),arb(0))/j)
        full=arb(1)
        for factor,n in zip(factors,row):full*=factor**int(n)
        fulls.append(full);lows.append(sum((weight*value for weight,value in zip(weights,coefficients)),arb(0)))
    return fulls,lows


def outward_candidate(data,p,z,cutoff,original):
    p=Q(p);W=data['windows'];S=1<<data['bits']
    if not 0<=p<=1 or not 0<z<=1 or type(cutoff) is not int or not 1<=cutoff<=W+1:
        raise ValueError('valid activity, output weight, and occupancy cutoff required')
    fulls,lows=fourier_parts(data,p,z,cutoff-1)
    return _candidate_from_parts(data,p,z,cutoff,original,fulls,lows)


def _candidate_from_parts(data,p,z,cutoff,original,fulls,lows):
    """Finish a candidate from Fourier parts freshly built in this call."""
    W=data['windows'];S=1<<data['bits']
    tails=[full-low for full,low in zip(fulls,lows)]
    # Any constant may be subtracted at a nonzero syndrome. A floating
    # median merely chooses an exact dyadic constant; no ordering claim
    # from floating arithmetic enters the inequality.
    order=sorted(range(len(tails)),key=lambda i:float(tails[i]));seen=0;center=arb(0)
    for i in order:
        seen+=int(data['multiplicities'][i])
        if 2*seen>=S:center=arb(float(tails[i]));break
    cap=up(sum((int(n)*abs(value-center) for n,value in zip(data['multiplicities'],tails)),arb(0))/S)
    packet=((1+z)**4-1)/15
    masses=[comb(W,j)*aq(p)**j*aq(1-p)**(W-j)*packet**j for j in range(W+1)]
    cap=min(cap,up(sum(masses[cutoff:],arb(0))))
    matrix=original*arb(1)
    matrix[0,1]=up(sum(masses[1:cutoff],arb(0)))
    matrix[0,2]=up((S-1)*cap)
    return matrix


def outward_at_z(data,probabilities,z):
    original=base.outward_at_z(data,probabilities,z)
    p=1-Q(probabilities[0]);zf=float(z)
    initial=base.floating(data,probabilities,-np.log(zf))
    cutoff,_=select(data,float(p),zf,initial)
    return original if cutoff is None else outward_candidate(data,p,z,cutoff,original)


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
