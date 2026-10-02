"""IID packet comparison retaining a uniform-density state component.

Coordinates represent zero mass, arbitrary nonzero mass, and a measure
dominated pointwise by a multiple of uniform nonzero mass. No assertion
that the actual conditioned state is uniform is made.
"""
from math import comb
from fractions import Fraction as Q
import numpy as np
from flint import arb,arb_mat
import scalar_cover as sc

base=sc.kernel
aq,up=base.aq,base.up


def prepare(images,columns,bits,updates=2):
    if (type(bits) is not int or not 1<=bits<=24 or type(updates) is not int or not 1<=updates<=32
            or not columns or len(columns)%4 or len(columns)>128 or len(images)!=1<<bits
            or images[0]!=0 or any(x==0 for x in images[1:])
            or any(not 0<=c<1<<bits for c in columns)):
        raise ValueError('injective binary expansion and complete packet geometry required')
    W=len(columns)//4
    if (any(not 0<=x<1<<(4*W) for x in images)
            or any(images[i]!=images[i^(i&-i)]^images[i&-i] for i in range(1,len(images)))):
        raise ValueError('linear expansion in the output geometry required')
    radix=W+1;powers=np.array([radix**j for j in range(5)],dtype=np.uint64)
    chars=np.arange(1<<bits,dtype=np.uint32)
    char_codes=np.zeros(len(chars),dtype=np.uint64)
    image_codes=np.zeros(len(images),dtype=np.uint64)
    halves=[np.array([x&((1<<64)-1) for x in images],dtype=np.uint64),
            np.array([x>>64 for x in images],dtype=np.uint64)]
    for w in range(W):
        r=sum((np.bitwise_count(chars&c)&1 for c in columns[4*w:4*w+4]))
        char_codes+=powers[r]
        r=np.bitwise_count((halves[w//16]>>(4*(w%16)))&15)
        image_codes+=powers[r]
    def summarize(codes):
        unique,mult=np.unique(codes,return_counts=True)
        rows=np.array([[int(code)//radix**j%radix for j in range(5)] for code in unique],dtype=np.int64)
        order=np.lexsort(tuple(rows[:,j] for j in range(4,-1,-1)))
        return rows[order],mult[order]
    records,mult=summarize(char_codes);histograms,hist_mult=summarize(image_codes[1:])
    if not np.all(records.sum(axis=1)==W) or not np.all(histograms.sum(axis=1)==W):
        raise ArithmeticError('histogram encoding failed')
    return dict(bits=bits,windows=W,updates=updates,columns=list(columns),records=records,multiplicities=mult,
                histograms=histograms,histogram_multiplicities=hist_mult)


def actual(updates=2):
    images,columns,_=base.maps()
    data=prepare(images,columns,19,updates)
    print('REFRESH IID KERNEL',len(data['records']),'characters;',len(data['histograms']),'image histograms',flush=True)
    return data


def floating(data,probabilities,tilt,*,feedback_aware=False):
    p=np.asarray(probabilities,dtype=float);z=np.exp(-tilt)
    S=1<<data['bits'];W=data['windows'];alpha=2.**-data['updates'];beta=1-alpha
    word=p[base.WEIGHT]/np.array([comb(4,int(w)) for w in base.WEIGHT])
    weighted=word[None,:]*z**base.DIST;totals=weighted.sum(axis=1)
    moments=np.exp(data['histograms']@np.log(totals[[0,1,3,7,15]]))
    h=float(moments.max());mean=float(data['histogram_multiplicities']@moments)/(S-1)
    factors=base.KRAW@(p*z**np.arange(5))
    transformed=base.float_products(factors,data['records'],W)
    total=float(factors[0]**W);zero=float(data['multiplicities']@transformed)/S
    normal=base.float_products(base.KRAW@p,data['records'],W)
    nonzero=max(0,min(1,1-float(data['multiplicities']@normal)/S)) if feedback_aware else 1.
    atom=min(float(data['multiplicities']@np.abs(normal)),
             float(data['multiplicities']@np.abs(normal-p[0]**W)))/S
    rhos=np.max(np.abs(weighted@base.SIGN.T)/totals[:,None],axis=0)
    tilted=float(data['multiplicities']@base.float_products(rhos,data['records'],W))/S
    factors2=(word[None,:]*z**(2*base.DIST)).sum(axis=1)[[0,1,3,7,15]]
    moments2=np.exp(data['histograms']@np.log(factors2))
    h2=float(moments2.max());mean2=float(data['histogram_multiplicities']@moments2)/(S-1)
    cancel=min(h,np.sqrt(max(0,atom*h2)),h*tilted)
    uniform_cancel=min(mean,nonzero/(S-1),np.sqrt(max(0,mean2*nonzero/(S-1))),mean*tilted)
    refresh_h=min(h,np.sqrt(max(0,h2*nonzero)))
    refresh_mean=min(mean,np.sqrt(max(0,mean2*nonzero)))
    zero=max(0,min(total,zero))
    return np.array([[zero,total-zero,0.],
                     [alpha*cancel+beta*refresh_h/(S-1),alpha*h,beta*h],
                     [alpha*uniform_cancel+beta*refresh_mean/(S-1),alpha*mean,beta*mean]])


def outward_at_z(data,probabilities,z,*,feedback_aware=False):
    p=list(map(Q,probabilities))
    if len(p)!=5 or min(p)<0 or sum(p)!=1 or not 0<z<=1:
        raise ValueError('five probabilities and output weight in (0,1] required')
    p=list(map(aq,p));W=data['windows'];S=1<<data['bits'];alpha=arb(2)**-data['updates'];beta=1-alpha
    records=[tuple(map(int,h)) for h in data['records']];mult=list(map(int,data['multiplicities']))
    word=[p[x.bit_count()]/comb(4,x.bit_count()) for x in range(16)]
    weighted=[[prob*z**(x^y).bit_count() for x,prob in enumerate(word)] for y in range(16)]
    totals=[sum(row,arb(0)) for row in weighted]
    def moments(power):
        factors=[sum((prob*z**(power*(x^y).bit_count()) for x,prob in enumerate(word)),arb(0)) for y in (0,1,3,7,15)]
        values=[base.np_product(factors,h) for h in data['histograms']]
        return max(map(up,values)),up(sum((int(n)*v for n,v in zip(data['histogram_multiplicities'],values)),arb(0))/(S-1))
    def fourier(value):
        factors=[sum((p[b]*value**b*sum((-1)**j*comb(r,j)*comb(4-r,b-j)
                         for j in range(max(0,b-4+r),min(b,r)+1))/comb(4,b)
                      for b in range(5)),arb(0)) for r in range(5)]
        return base.arb_products(factors,records,W)
    total=totals[0]**W
    zero=sum((n*f for n,f in zip(mult,fourier(z))),arb(0))/S
    normal=fourier(arb(1));empty=p[0]**W
    nonzero=min(arb(1),max(arb(0),up(1-sum((n*f for n,f in zip(mult,normal)),arb(0))/S))) if feedback_aware else arb(1)
    atom=min(up(sum((n*abs(f) for n,f in zip(mult,normal)),arb(0))/S),
             up(sum((n*abs(f-empty) for n,f in zip(mult,normal)),arb(0))/S))
    rhos=[]
    for b in range(5):
        mask=(1<<b)-1
        rhos.append(min(arb(1),max(up(abs(sum((prob*((-1)**((x&mask).bit_count()%2))
                          for x,prob in enumerate(row)),arb(0)))/mass)
                        for row,mass in zip(weighted,totals))))
    tilted=up(sum((n*f for n,f in zip(mult,base.arb_products(rhos,records,W))),arb(0))/S)
    h,mean=moments(1);h2,mean2=moments(2)
    cancel=min(h,up((atom*h2).sqrt()),up(h*tilted))
    uniform_cancel=min(mean,up(nonzero/(S-1)),up((mean2*nonzero/(S-1)).sqrt()),up(mean*tilted))
    refresh_h=min(h,up((h2*nonzero).sqrt()))
    refresh_mean=min(mean,up((mean2*nonzero).sqrt()))
    zu,nu=up(zero),up(total-zero)
    if min(zu,nu)<0:raise ArithmeticError('negative zero/nonzero feedback contribution')
    return arb_mat([[min(up(total),zu),min(up(total),nu),0],
                    [up(alpha*cancel+beta*refresh_h/(S-1)),up(alpha*h),up(beta*h)],
                    [up(alpha*uniform_cancel+beta*refresh_mean/(S-1)),up(alpha*mean),up(beta*mean)]])


def outward(data,probabilities,tilt,*,feedback_aware=False):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,probabilities,(-aq(tilt)).exp(),feedback_aware=feedback_aware)
