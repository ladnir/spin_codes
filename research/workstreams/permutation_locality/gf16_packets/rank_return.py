"""Output-weight bound on a feedback fiber using its number of free bits."""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb
import trimmed_return as base

aq,up=base.aq,base.up


def attach(data):
    W=data['windows'];counts=[0]*(W+1)
    for row,n in zip(data['records'],data['multiplicities']):counts[int(row[0])]+=int(n)
    result=[]
    for histogram in data['histograms']:
        weights=[w for w in range(4,-1,-1) for _ in range(int(histogram[w]))]
        total=sum(weights);subset=np.zeros((W+1,total+1),dtype=object);subset[0,0]=1
        row=np.zeros((W+1,4*W+1),dtype=object)
        for c in range(W+1):
            if counts[c]:row[:,:total+1]+=counts[c]*subset[:,::-1]
            if c<W:
                v=weights[c]
                if v:subset[1:,v:]+=subset[:-1,:-v].copy()
                else:subset[1:]+=subset[:-1].copy()
        for j,coefficients in enumerate(row):
            expected=sum(n*comb(c,j) for c,n in enumerate(counts) if c>=j)
            if any(x<0 for x in coefficients) or sum(coefficients)!=expected:
                raise ArithmeticError('annihilator-weighted support census failed')
        result.append(row)
    floats=np.array(result,dtype=float)/np.array([comb(W,j)*15**j for j in range(W+1)],dtype=float)[None,:,None]
    return dict(data,rank_counts=result,rank_float=floats)


def prepare(images,columns,bits,updates=2):return attach(base.prepare(images,columns,bits,updates))


def actual(updates=2):
    data=attach(base.actual(updates))
    print('GF16 RANK KERNEL exact annihilator-weighted support counts',flush=True)
    return data


def float_bounds(data,p,z):
    W=data['windows'];L=(1<<data['bits'])-1
    factors=(1+z)**(4*np.arange(W+1)-data['bits'])
    powers=z**np.arange(4*W+1)
    values=np.minimum((data['rank_float']@powers)*factors[None,:],data['trimmed_float']@powers)
    maximum=values.max(axis=0)
    mean=np.minimum(data['histogram_multiplicities']@values/L,data['trimmed_uniform_float']@powers)
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    return float(weights@maximum),float(weights@mean)


def fixed_profile_bounds(data,z):
    W=data['windows']
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    powers=[z**w for w in range(4*W+1)];result=[]
    for j in range(W+1):
        factor=(1+z)**(4*j-data['bits'])/(comb(W,j)*15**j)
        D=comb(W,j)*15**j
        values=[min(up(factor*sum((int(n)*power for n,power in zip(rows[j],powers) if n),arb(0))),
                    up(sum((int(n)*power for n,power in zip(trimmed[j],powers) if n),arb(0))/D))
                for rows,trimmed in zip(data['rank_counts'],data['trimmed'])]
        result.append(values)
    return result


def fixed_bounds(data,z):
    W=data['windows'];L=(1<<data['bits'])-1
    powers=[z**w for w in range(4*W+1)];maxima=[];means=[]
    for j,values in enumerate(fixed_profile_bounds(data,z)):
        D=comb(W,j)*15**j;maxima.append(max(values))
        means.append(min(up(sum((int(n)*v for n,v in zip(data['histogram_multiplicities'],values)),arb(0))/L),
                         up(sum((int(n)*power for n,power in zip(data['trimmed_uniform'][j],powers) if n),arb(0))/(L*D))))
    return maxima,means


def bounds(data,p,z):
    p=Q(p);W=data['windows']
    if not 0<=p<=1:raise ValueError('probability required')
    maxima,means=fixed_bounds(data,z)
    weights=[comb(W,j)*aq(p)**j*aq(1-p)**(W-j) for j in range(W+1)]
    return up(sum((w*v for w,v in zip(weights,maxima)),arb(0))),up(sum((w*v for w,v in zip(weights,means)),arb(0)))


def floating(data,probabilities,tilt):
    matrix=base.floating(data,probabilities,tilt)
    maximum,mean=float_bounds(data,float(1-probabilities[0]),np.exp(-tilt))
    alpha=2.**-data['updates'];L=(1<<data['bits'])-1
    matrix[1,0]=min(matrix[1,0],alpha*maximum+matrix[1,2]/L)
    matrix[2,0]=min(matrix[2,0],alpha*mean+matrix[2,2]/L)
    return matrix


def outward_at_z(data,probabilities,z):
    matrix=base.outward_at_z(data,probabilities,z)
    maximum,mean=bounds(data,base.base.base.base.base.activity(probabilities),z)
    alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
    matrix[1,0]=min(matrix[1,0],up(alpha*maximum+matrix[1,2]/L))
    matrix[2,0]=min(matrix[2,0],up(alpha*mean+matrix[2,2]/L))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
