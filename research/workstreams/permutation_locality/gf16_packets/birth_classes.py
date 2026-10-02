"""Outward weighted birth classes, retaining expansion weight after activation.

The extra coordinates are arbitrary nonzero measures restricted to one
expansion-weight class, not conditional uniform-state assertions.
"""
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb,arb_mat
import birth_classes_probe as proposal
import birth_refresh as base
import rank_return
import scalar_cover as sc

aq,up=base.aq,base.up


def prepare(images,columns,bits,updates=2):
    return dict(proposal.prepare(images,columns,bits,updates),class_returns=True)


def actual(updates=2):
    data=proposal.actual(updates);data['class_returns']=True
    return data


def select(data,probabilities,tilt):
    options=proposal.candidates(data,probabilities,tilt,class_returns=True)
    name=min(options,key=lambda key:sc.log_power(options[key]))
    return name,options[name]


def floating(data,probabilities,tilt):return select(data,probabilities,tilt)[1]


def weighted_classes(data,p,z,cutoff=None):
    """Outward class masses for all births, or J below the cutoff."""
    fulls,lows=_birth_parts(data,p,z,cutoff)
    return _class_masses_from_values(data,fulls if cutoff is None else lows)


def _birth_parts(data,p,z,cutoff):
    W=data['windows']
    if cutoff is not None and (type(cutoff) is not int or not 1<=cutoff<=W+1):
        raise ValueError('valid occupancy cutoff required')
    return base.fourier_parts(data,p,z,0 if cutoff is None else cutoff-1)


def _class_masses_from_values(data,values):
    S=1<<data['bits']
    return [max(arb(0),up(sum((int(n)*v for n,v in zip(row,values) if n),arb(0))/S))
            for row in data['birth_class_census']]


def outward_candidate(data,probabilities,z,name,original=None):
    p=1-Q(probabilities[0]);W=data['windows'];L=(1<<data['bits'])-1
    if name=='original':return base.outward_at_z(data,probabilities,z)
    if name=='all_classes':cutoff=None
    elif name.startswith('classes_below_'):cutoff=int(name.removeprefix('classes_below_'))
    else:raise ValueError('unknown birth-class transfer candidate')
    # These two output envelopes use exactly the same Fourier moments.
    # Build them once per outward call, without any saved numeric cache.
    fulls,lows=_birth_parts(data,p,z,cutoff)
    classes=_class_masses_from_values(data,fulls if cutoff is None else lows)
    # The birth-density variant changes only the zero-source M/U entries,
    # which are replaced below. Start from the unchanged rank kernel.
    if original is None:original=rank_return.outward_at_z(data,probabilities,z)
    levels=list(map(int,data['birth_class_levels']));n=3+len(levels)
    rows=[[arb(0)]*n for _ in range(n)]
    for i in range(3):
        for j in range(3):rows[i][j]=original[i,j]
    rows[0][1]=arb(0);rows[0][2]=arb(0);rows[0][3:]=classes
    if cutoff is not None:
        rows[0][2]=base._candidate_from_parts(data,p,z,cutoff,original,fulls,lows)[0,2]
    alpha=arb(2)**-data['updates'];beta=1-alpha;empty=aq(1-p)**W
    factors=[aq(1-16*p/15)*z**w+aq(p/15)*(1+z)**4 for w in range(5)]
    moments=[]
    for histogram in data['histograms']:
        value=arb(1)
        for factor,count in zip(factors,histogram):value*=factor**int(count)
        moments.append(up(value))
    fibers=rank_return.fixed_profile_bounds(data,z)
    weights=[comb(W,j)*aq(p)**j*aq(1-p)**(W-j) for j in range(W+1)]
    for i,level in enumerate(levels,3):
        selected=[k for k,w in enumerate(data['image_histogram_weights']) if w==level]
        total=max(moments[k] for k in selected);quiet=empty*z**level
        lazy_return=up(sum((weight*max(row[k] for k in selected) for weight,row in zip(weights,fibers)),arb(0)))
        rows[i][0]=min(original[1,0],total,up(alpha*lazy_return+beta*total/L))
        rows[i][1]=max(arb(0),up(alpha*(total-quiet)))
        rows[i][2]=up(beta*total);rows[i][i]=up(alpha*quiet)
    return arb_mat(rows)


def outward_at_z(data,probabilities,z):
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    name,_=select(data,probabilities,-np.log(float(z)))
    return outward_candidate(data,probabilities,z,name)


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
