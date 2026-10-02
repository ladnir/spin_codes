"""Four-bit iid reference kernel; conditioning costs belong to the caller.

Two coordinates distinguish zero and nonzero state. Inputs are independent
packets, uniform within each specified weight class. This is a comparison
measure, not the actual routed BCH input distribution.
"""
from collections import Counter
from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sys

import numpy as np
from flint import arb, arb_mat

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from group_moment import maps
from group_rank_one_verify import up

DIST=np.array([[(x^y).bit_count() for x in range(16)] for y in range(16)])
WEIGHT=np.array([x.bit_count() for x in range(16)])
SIGN=np.array([[(-1)**((x&((1<<b)-1)).bit_count()%2) for x in range(16)] for b in range(5)])
KRAW=np.array([[sum((-1)**h*comb(r,h)*comb(4-r,b-h)
                   for h in range(max(0,b-4+r),min(b,r)+1))/comb(4,b)
                for b in range(5)] for r in range(5)])


def aq(value):
    value=Q(value)
    return arb(value.numerator)/value.denominator


def prepare(images,columns,bits,updates=2):
    if (type(updates) is not int or not 1<=updates<=32 or not columns or len(columns)%4
            or len(images)!=1<<bits or images[0]!=0 or any(x==0 for x in images[1:])
            or any(not 0<=c<1<<bits for c in columns)):
        raise ValueError('injective expansion, complete packet geometry and update count required')
    windows=len(columns)//4
    if (any(not 0<=x<1<<(4*windows) for x in images)
            or any(images[i]!=images[i^(i&-i)]^images[i&-i] for i in range(1,len(images)))):
        raise ValueError('complete linear expansion images in the output geometry required')
    chars=np.arange(1<<bits,dtype=np.uint32)
    hist=np.zeros((len(chars),5),dtype=np.uint16)
    for w in range(windows):
        weight=sum((np.bitwise_count(chars&c)&1 for c in columns[4*w:4*w+4]))
        for r in range(5):hist[:,r]+=weight==r
    records,mult=np.unique(hist,axis=0,return_counts=True)
    patterns=sorted({tuple(Counter(((word>>(4*w))&15).bit_count()
                                  for w in range(windows))[r] for r in range(5))
                     for word in images[1:]})
    assert np.all(records.sum(axis=1)==windows) and sum(mult)==1<<bits
    return dict(bits=bits,windows=windows,updates=updates,records=records.astype(np.int64),
                multiplicities=mult,histograms=np.array(patterns,dtype=np.int64))


def actual(updates=2):
    images,columns,_=maps()
    data=prepare(images,columns,19,updates)
    print('FOUR-BIT IID KERNEL',len(data['records']),'character histograms;',
          len(data['histograms']),'nonzero output patterns',flush=True)
    return data


def float_products(factors,records,windows):
    powers=np.ones((5,windows+1))
    powers[:,1:]=np.cumprod(np.broadcast_to(factors[:,None],(5,windows)),axis=1)
    return np.prod(powers[np.arange(5)[None,:],records],axis=1)


def arb_products(factors,records,windows):
    powers=[[arb(1)] for _ in factors]
    for x,p in zip(factors,powers):
        for _ in range(windows):p.append(p[-1]*x)
    result=[]
    for h in records:
        value=arb(1)
        for p,n in zip(powers,h):value*=p[n]
        result.append(value)
    return result


def floating(data,probabilities,tilt):
    probabilities=np.asarray(probabilities,dtype=float)
    z=np.exp(-tilt);size=1<<data['bits'];W=data['windows']
    word=probabilities[WEIGHT]/np.array([comb(4,int(w)) for w in WEIGHT])
    weighted=word[None,:]*z**DIST
    totals=weighted.sum(axis=1)
    h=float(np.exp(np.max(data['histograms']@np.log(totals[[0,1,3,7,15]]))))
    factors=KRAW@(probabilities*z**np.arange(5))
    transformed=float_products(factors,data['records'],W)
    total=float(factors[0]**W)
    zero=float(data['multiplicities']@transformed)/size
    normal=KRAW@probabilities
    spectrum=float_products(normal,data['records'],W)
    atom=min(float(data['multiplicities']@np.abs(spectrum)),
             float(data['multiplicities']@np.abs(spectrum-probabilities[0]**W)))/size
    rhos=np.max(np.abs(weighted@SIGN.T)/totals[:,None],axis=0)
    tilted_atom=float(data['multiplicities']@float_products(rhos,data['records'],W))/size
    factors2=(word[None,:]*z**(2*DIST)).sum(axis=1)[[0,1,3,7,15]]
    h2=float(np.exp(np.max(data['histograms']@np.log(factors2))))
    cancel=min(h,np.sqrt(max(0,atom*h2)),h*tilted_atom)
    zero=max(0,min(total,zero));alpha=2.**-data['updates']
    return np.array([[zero,total-zero],[alpha*cancel+(1-alpha)*h/(size-1),h]])


def outward_at_z(data,probabilities,z):
    """Arb enclosure; z can be a ball enclosing exp(-lambda)."""
    probs=list(map(Q,probabilities))
    if len(probs)!=5 or min(probs)<0 or sum(probs)!=1 or not 0<z<=1:
        raise ValueError('five probabilities and output weight in (0,1] required')
    probs=list(map(aq,probs));W=data['windows'];size=1<<data['bits']
    records=[tuple(map(int,h)) for h in data['records']]
    mult=list(map(int,data['multiplicities']))
    word=[probs[x.bit_count()]/comb(4,x.bit_count()) for x in range(16)]
    weighted=[[p*z**(x^y).bit_count() for x,p in enumerate(word)] for y in range(16)]
    totals=[sum(row,arb(0)) for row in weighted]
    def moment(power):
        factors=[sum((p*z**(power*(x^y).bit_count()) for x,p in enumerate(word)),arb(0))
                 for y in (0,1,3,7,15)]
        return max(up(np_product(factors,h)) for h in data['histograms'])
    def fourier(z):
        factors=[sum((probs[b]*z**b*sum((-1)**h*comb(r,h)*comb(4-r,b-h)
                         for h in range(max(0,b-4+r),min(b,r)+1))/comb(4,b)
                      for b in range(5)),arb(0)) for r in range(5)]
        # Arb's real-power operation may return NaN when a ball straddles
        # zero, even for an integral exponent. These are polynomials.
        return arb_products(factors,records,W)
    total=totals[0]**W
    zero=sum((n*f for n,f in zip(mult,fourier(z))),arb(0))/size
    transform=fourier(arb(1));empty=probs[0]**W
    atom=min(up(sum((n*abs(f) for n,f in zip(mult,transform)),arb(0))/size),
             up(sum((n*abs(f-empty) for n,f in zip(mult,transform)),arb(0))/size))
    rhos=[]
    for b in range(5):
        mask=(1<<b)-1
        rhos.append(min(arb(1),max(up(abs(sum((p*((-1)**((x&mask).bit_count()%2))
                          for x,p in enumerate(row)),arb(0)))/total)
                        for row,total in zip(weighted,totals))))
    tilted_atom=up(sum((n*f for n,f in zip(mult,arb_products(rhos,records,W))),arb(0))/size)
    h,h2=moment(1),moment(2)
    cancel=min(h,up((h2*atom).sqrt()),up(h*tilted_atom))
    zu,nu=up(zero),up(total-zero)
    if min(zu,nu)<0:raise ArithmeticError('negative zero/nonzero feedback contribution')
    alpha=arb(2)**-data['updates']
    return arb_mat([[min(up(total),zu),min(up(total),nu)],
                    [up(alpha*cancel+(1-alpha)*h/(size-1)),h]])


def np_product(values,counts):
    value=arb(1)
    for x,n in zip(values,counts):
        n=int(n)
        while n:
            if n&1:value*=x
            n>>=1
            if n:x=x*x
    return value


def outward(data,probabilities,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp())
