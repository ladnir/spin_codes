"""Bound returns by the lightest outputs permitted by their event count.

All occupancy/output enumerators and truncation budgets use exact integers.
The event may correlate arbitrarily with emitted weight.
"""
from fractions import Fraction as Q
from math import comb
from numbers import Integral
import numpy as np
from flint import arb
import conditioned_return as base

aq,up=base.aq,base.up


def output_counts(histogram):
    if (len(histogram)!=5 or any(not isinstance(n,Integral) or n<0 for n in histogram)
            or not 0<=sum(histogram)<=32):raise ValueError('at most 32 four-bit packet weights required')
    W=int(sum(histogram));values=np.ones((1,1),dtype=object);done=0
    for v,n in enumerate(histogram):
        for _ in range(int(n)):
            new=np.zeros((done+2,4*(done+1)+1),dtype=object)
            rows,cols=values.shape
            new[:rows,v:v+cols]+=values
            for w in range(5):
                count=comb(4,w)-int(w==v)
                if count:new[1:rows+1,w:w+cols]+=count*values
            values=new;done+=1
    for j,row in enumerate(values):
        if any(x<0 for x in row) or sum(row)!=comb(W,j)*15**j:
            raise ArithmeticError('output census has incorrect occupancy mass')
    return values


def lightest(counts,budget):
    """Select budget outcomes in nondecreasing output-weight order."""
    if type(budget) is not int or budget<0:raise ValueError('nonnegative integer budget required')
    result=[];left=budget
    for value in counts:
        count=int(value)
        if count!=value or count<0:raise ValueError('nonnegative integer counts required')
        take=min(count,left);result.append(take);left-=take
    if left:raise ValueError('event budget exceeds the sample space')
    if sum(result)!=budget:raise ArithmeticError('trimming lost event mass')
    return tuple(result)


def attach(data):
    W=data['windows'];L=(1<<data['bits'])-1
    denoms=[comb(W,j)*15**j for j in range(W+1)]
    budgets=[int(Q(cap)*D) for cap,D in zip(data['nonzero_atom_caps'],denoms)]
    all_counts=[output_counts(h) for h in data['histograms']]
    trimmed=[[lightest(row,budgets[j]) for j,row in enumerate(counts)] for counts in all_counts]
    aggregate=sum((int(n)*counts for n,counts in zip(data['histogram_multiplicities'],all_counts)),
                  np.zeros((W+1,4*W+1),dtype=object))
    uniform=[]
    for j,(row,D,zero) in enumerate(zip(aggregate,denoms,data['zero_probabilities'])):
        count=D*(1-zero)
        if count.denominator!=1 or sum(row)!=L*D:raise ArithmeticError('uniform-state event count is not integral')
        uniform.append(lightest(row,int(count)))
    floats=np.array(trimmed,dtype=float)/np.array(denoms,dtype=float)[None,:,None]
    uniform_floats=np.array(uniform,dtype=float)/(L*np.array(denoms,dtype=float)[:,None])
    return dict(data,trimmed=trimmed,trimmed_uniform=uniform,trimmed_denominators=denoms,
                trimmed_float=floats,trimmed_uniform_float=uniform_floats)


def prepare(images,columns,bits,updates=2):return attach(base.prepare(images,columns,bits,updates))


def actual(updates=2):
    data=attach(base.actual(updates))
    print('GF16 TRIMMED KERNEL',len(data['trimmed']),'exact occupancy/output profiles',flush=True)
    return data


def float_bounds(data,p,z):
    W=data['windows'];L=(1<<data['bits'])-1;powers=z**np.arange(4*W+1)
    values=data['trimmed_float']@powers
    maximum=values.max(axis=0)
    mean=np.minimum(data['histogram_multiplicities']@values/L,data['trimmed_uniform_float']@powers)
    weights=np.array([comb(W,j)*p**j*(1-p)**(W-j) for j in range(W+1)])
    return float(weights@maximum),float(weights@mean)


def bounds(data,p,z):
    p=Q(p);W=data['windows'];L=(1<<data['bits'])-1
    if not 0<=p<=1 or not 0<z<=1:raise ValueError('probability and output weight required')
    powers=[z**w for w in range(4*W+1)]
    def evaluate(row):return sum((n*power for n,power in zip(row,powers) if n),arb(0))
    maxima=[];means=[]
    for j,D in enumerate(data['trimmed_denominators']):
        values=[up(evaluate(rows[j])/D) for rows in data['trimmed']]
        maxima.append(max(values))
        means.append(min(up(sum((int(n)*v for n,v in zip(data['histogram_multiplicities'],values)),arb(0))/L),
                         up(evaluate(data['trimmed_uniform'][j])/(L*D))))
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
    maximum,mean=bounds(data,base.base.base.base.activity(probabilities),z)
    alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
    matrix[1,0]=min(matrix[1,0],up(alpha*maximum+matrix[1,2]/L))
    matrix[2,0]=min(matrix[2,0],up(alpha*mean+matrix[2,2]/L))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
