"""Use the expansion weight when bounding a cancelling Fourier coefficient.

For activity p<=15/16, a packet is a positive mixture of a zero packet
and uniform four bits. Rearrangement and convexity bound the character
product without allowing every packet to choose weight zero.
"""
from fractions import Fraction as Q
import numpy as np
from flint import arb
import weighted_return as base

aq,up=base.aq,base.up
arithmetic=base.base.base.base


def pair_profile(records,histogram):
    """Pair a fixed weight histogram with each character histogram."""
    records=np.asarray(records,dtype=np.int64)
    if (records.ndim!=2 or records.shape[1]!=5 or np.any(records<0)
            or not len(records) or not np.all(records.sum(axis=1)==records[0].sum())):
        raise ValueError('five-category character histograms required')
    W=int(records[0].sum())
    if len(histogram)!=5 or any(type(n) is not int or n<0 for n in histogram) or sum(histogram)!=W:
        raise ValueError('matching expansion histogram required')
    segments=[(w,histogram[w]) for w in range(4,-1,-1) if histogram[w]]
    end=np.cumsum(records,axis=1);start=end-records;offset=0;parts=[]
    for weight,count in segments:
        parts.append(np.maximum(0,np.minimum(end,offset+count)-np.maximum(start,offset)))
        offset+=count
    if not np.array_equal(sum(parts),records):raise ArithmeticError('profile pairing lost a packet')
    return [weight for weight,count in segments],np.concatenate(parts,axis=1)


def profiles(records,distance):
    """Canonical weight profile, paired with each character histogram."""
    W=int(np.asarray(records)[0].sum())
    if type(distance) is not int or not 0<=distance<=4*W:raise ValueError('valid expansion distance required')
    full,remainder=divmod(distance,4)
    histogram=[W-full-int(bool(remainder)),0,0,0,full]
    if remainder:histogram[remainder]+=1
    return pair_profile(records,histogram)


def attach(data):
    weights,records=profiles(data['records'],data['distance'])
    return dict(data,profile_weights=weights,profile_records=records,
                profile_tuples=[tuple(map(int,row)) for row in records])


def prepare(images,columns,bits,updates=2):return attach(base.prepare(images,columns,bits,updates))


def actual(updates=2):return attach(base.actual(updates))


def float_lazy(data,p,z):
    a=1-16*p/15;b=p/15;W=data['windows']
    factors=np.array([a*z**w+b*(1+z)**(4-r)*(1-z)**r for w in data['profile_weights'] for r in range(5)])
    powers=np.ones((len(factors),W+1));powers[:,1:]=np.cumprod(np.broadcast_to(factors[:,None],(len(factors),W)),axis=1)
    products=np.prod(powers[np.arange(len(factors))[None,:],data['profile_records']],axis=1)
    return max(0.,float(data['multiplicities']@products)/(1<<data['bits'])-a**W*z**data['distance'])


def lazy_bound(data,p,z):
    p=Q(p)
    if not 0<=p<=Q(15,16) or not 0<z<=1:raise ValueError('positive-mixture range and output weight in (0,1] required')
    a=aq(1-16*p/15);b=aq(p/15);W=data['windows']
    factors=[a*z**w+b*(1+z)**(4-r)*(1-z)**r for w in data['profile_weights'] for r in range(5)]
    products=arithmetic.arb_products(factors,data['profile_tuples'],W)
    bound=up(sum((int(n)*value for n,value in zip(data['multiplicities'],products)),arb(0))/(1<<data['bits'])
             -a**W*z**data['distance'])
    if bound<0:raise ArithmeticError('negative profile cancellation bound')
    return bound


def floating(data,probabilities,tilt):
    matrix=base.floating(data,probabilities,tilt);p=float(1-probabilities[0])
    if p<=15/16:
        value=float_lazy(data,p,np.exp(-tilt));alpha=2.**-data['updates'];L=(1<<data['bits'])-1
        matrix[1,0]=min(matrix[1,0],alpha*value+matrix[1,2]/L)
    return matrix


def outward_at_z(data,probabilities,z):
    matrix=base.outward_at_z(data,probabilities,z);p=base.activity(probabilities)
    if p<=Q(15,16):
        value=lazy_bound(data,p,z);alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
        matrix[1,0]=min(matrix[1,0],up(alpha*value+matrix[1,2]/L))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
