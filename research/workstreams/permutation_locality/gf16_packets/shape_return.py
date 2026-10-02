"""Profile-specific return bounds, with separate exceptional-state bounds.

The rare minimum-weight full/zero packet states use Cauchy--Schwarz on
their two character factors. Signed enumeration is an optional diagnostic.
Other states retain their actual weight histogram.
"""
from fractions import Fraction as Q
import numpy as np
from flint import arb
import profile_return as base

aq,up=base.aq,base.up
arithmetic=base.arithmetic


def signed_records(columns,bits,state,image,exact=False):
    W=len(columns)//4;chars=np.arange(1<<bits,dtype=np.uint32);radix=W+1
    if W>32 or any((image>>(4*i))&15 not in (0,15) for i in range(W)):
        raise ValueError('full/zero packet image of at most 32 packets required')
    powers=np.array([radix**i for i in range(10)],dtype=np.uint64)
    codes=np.zeros(len(chars),dtype=np.uint64)
    for i in range(W):
        r=sum((np.bitwise_count(chars&col)&1 for col in columns[4*i:4*i+4]))
        category=5 if (image>>(4*i))&15 else 0
        codes+=powers[category+r]
    marginals=[]
    for part in (codes%(radix**5),codes//(radix**5)):
        keys,mults=np.unique(part,return_counts=True)
        rows=np.array([[int(code)//radix**i%radix for i in range(5)] for code in keys],dtype=np.int64)
        marginals.append(dict(records=rows,tuples=[tuple(map(int,row)) for row in rows],counts=mults))
    if not exact:return dict(state=state,marginals=marginals)
    unique,inverse=np.unique(codes,return_inverse=True)
    counts=np.zeros(len(unique),dtype=np.int64)
    signs=1-2*(np.bitwise_count(chars&state)&1).astype(np.int64)
    np.add.at(counts,inverse,signs)
    keep=counts!=0;unique=unique[keep];counts=counts[keep]
    if int(counts.sum())!=0:raise ArithmeticError('nonzero character did not cancel constant term')
    records=np.array([[int(code)//radix**i%radix for i in range(10)] for code in unique],dtype=np.int64).reshape(-1,10)
    if not np.all(records.sum(axis=1)==W):raise ArithmeticError('signed profile encoding failed')
    return dict(state=state,records=records,tuples=[tuple(map(int,row)) for row in records],counts=counts,marginals=marginals)


def prepare(images,columns,bits,updates=2,exact_exceptions=False):
    data=base.prepare(images,columns,bits,updates);d=data['distance'];W=data['windows']
    special=[]
    if d%4==0:
        special=[state for state,image in enumerate(images) if state and image.bit_count()==d
                 and all((image>>(4*i))&15 in (0,15) for i in range(W))]
    if len(special)>16:special=[]
    canonical=(W-d//4,0,0,0,d//4)
    shapes=[]
    for histogram,n in zip(data['histograms'],data['histogram_multiplicities']):
        histogram=tuple(map(int,histogram));count=int(n)
        if special and histogram==canonical:count-=len(special)
        if not count:continue
        weights,records=base.pair_profile(data['records'],histogram)
        shapes.append(dict(weights=weights,records=records,tuples=[tuple(map(int,row)) for row in records],
                           multiplicity=count,weight=sum(i*n for i,n in enumerate(histogram))))
    exceptions=[signed_records(columns,bits,state,images[state],exact_exceptions) for state in special]
    if sum(row['multiplicity'] for row in shapes)+len(exceptions)!=(1<<bits)-1:
        raise ArithmeticError('state profiles do not cover every nonzero state')
    return dict(data,shapes=shapes,exceptions=exceptions)


def actual(updates=2):
    images,columns,_=arithmetic.maps();data=prepare(images,columns,19,updates)
    print('GF16 SHAPE KERNEL',len(data['shapes']),'ordinary profiles;',len(data['exceptions']),
          'exceptional states with Cauchy bounds',flush=True)
    return data


def float_products(factors,records,W):
    factors=np.asarray(factors,dtype=float)
    # Positive factors dominate the search. Contract directly, without
    # materializing a records-by-factors gather on every output-tilt probe.
    if np.all(factors>0):return np.exp(np.einsum('ij,j->i',records,np.log(factors),optimize=False))
    powers=np.ones((len(factors),W+1));powers[:,1:]=np.cumprod(np.broadcast_to(factors[:,None],(len(factors),W)),axis=1)
    return np.prod(powers[np.arange(len(factors))[None,:],records],axis=1)


def float_bounds(data,p,z):
    W=data['windows'];S=1<<data['bits'];c=1-16*p/15;b=p/15
    factors=np.array([[c*z**w+b*(1+z)**(4-r)*(1-z)**r for r in range(5)] for w in range(5)])
    values=[];multiplicities=[]
    for row in data['shapes']:
        value=float(data['multiplicities']@float_products(factors[row['weights']].ravel(),row['records'],W))/S-c**W*z**row['weight']
        values.append(max(0.,value));multiplicities.append(row['multiplicity'])
    signed=[c*z**w+b*((-1)**r if w else 1)*(1+z)**(4-r)*(1-z)**r for w in (0,4) for r in range(5)]
    for row in data['exceptions']:
        moments=[float(m['counts']@float_products(np.array(signed[5*i:5*i+5])**2,m['records'],W))/S
                 for i,m in enumerate(row['marginals'])]
        value=np.sqrt(moments[0]*moments[1])
        if 'counts' in row:value=min(value,max(0.,float(row['counts']@float_products(signed,row['records'],W))/S))
        values.append(value);multiplicities.append(1)
    return max(values),float(np.array(multiplicities)@values)/(S-1)


def bounds(data,p,z):
    p=Q(p)
    if not 0<=p<=Q(15,16) or not 0<z<=1:raise ValueError('positive-mixture range required')
    W=data['windows'];S=1<<data['bits'];c=aq(1-16*p/15);b=aq(p/15)
    factors=[[c*z**w+b*(1+z)**(4-r)*(1-z)**r for r in range(5)] for w in range(5)]
    values=[];multiplicities=[]
    for row in data['shapes']:
        product=arithmetic.arb_products([x for w in row['weights'] for x in factors[w]],row['tuples'],W)
        value=up(sum((int(n)*v for n,v in zip(data['multiplicities'],product)),arb(0))/S-c**W*z**row['weight'])
        values.append(value);multiplicities.append(row['multiplicity'])
    signed=[c*z**w+b*((-1)**r if w else 1)*(1+z)**(4-r)*(1-z)**r for w in (0,4) for r in range(5)]
    for row in data['exceptions']:
        moments=[]
        for i,m in enumerate(row['marginals']):
            products=arithmetic.arb_products([x*x for x in signed[5*i:5*i+5]],m['tuples'],W)
            moments.append(up(sum((int(n)*v for n,v in zip(m['counts'],products)),arb(0))/S))
        value=up((moments[0]*moments[1]).sqrt())
        if 'counts' in row:
            product=arithmetic.arb_products(signed,row['tuples'],W)
            value=min(value,up(sum((int(n)*v for n,v in zip(row['counts'],product)),arb(0))/S))
        values.append(value);multiplicities.append(1)
    if any(v<0 for v in values):raise ArithmeticError('negative state-specific cancellation bound')
    return max(values),up(sum((n*v for n,v in zip(multiplicities,values)),arb(0))/(S-1))


def floating(data,probabilities,tilt):
    matrix=base.floating(data,probabilities,tilt);p=float(1-probabilities[0])
    if p<=15/16:
        maximum,mean=float_bounds(data,p,np.exp(-tilt));alpha=2.**-data['updates'];L=(1<<data['bits'])-1
        matrix[1,0]=min(matrix[1,0],alpha*maximum+matrix[1,2]/L)
        matrix[2,0]=min(matrix[2,0],alpha*mean+matrix[2,2]/L)
    return matrix


def outward_at_z(data,probabilities,z):
    matrix=base.outward_at_z(data,probabilities,z);p=base.base.activity(probabilities)
    if p<=Q(15,16):
        maximum,mean=bounds(data,p,z);alpha=arb(2)**-data['updates'];L=(1<<data['bits'])-1
        matrix[1,0]=min(matrix[1,0],up(alpha*maximum+matrix[1,2]/L))
        matrix[2,0]=min(matrix[2,0],up(alpha*mean+matrix[2,2]/L))
    return matrix


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
