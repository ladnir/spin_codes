"""Two-state output-moment envelope for an iid two-bit input reference.

This reference is used only after paying explicit profile/placement
conditioning. It must never replace the actual structured input law
without those factors. The default two-update inner has lazy probability
1/4. Optional update counts support isolated parameter diagnostics; changing
them requires a new sparse proof too. Output precedes refresh/feedback.
"""
from fractions import Fraction as Q

import numpy as np
from flint import arb,arb_mat

import model

DISTANCES=np.array([[(x^y).bit_count() for x in range(4)] for y in range(4)])
SIGNS=np.array([[(-1)**((x&c).bit_count()%2) for x in range(4)] for c in (1,3)])


def floating_tilted_atom(data,v,r,tilt):
    """Uniform syndrome cap after tilting by distance from any state word."""
    probabilities=np.array([1-v,v*(1-r)/2,v*(1-r)/2,v*r])
    weighted=probabilities[None,:]*np.exp(-tilt*DISTANCES)
    rhos=np.max(np.abs(weighted@SIGNS.T)/weighted.sum(axis=1)[:,None],axis=0)
    n1,n2=data['records'][:,1],data['records'][:,2]
    return float(data['multiplicities']@(rhos[0]**n1*rhos[1]**n2))/(1<<data['bits'])


def outward_tilted_atom(data,v,r,tilt):
    def aq(x):
        x=Q(x);return arb(x.numerator)/x.denominator
    if not 0<=Q(v)<=1 or not 0<=Q(r)<=1 or Q(tilt)<0:raise ValueError('valid probabilities and nonnegative tilt required')
    v,r=aq(v),aq(r);z=(-aq(tilt)).exp()
    probabilities=[1-v,v*(1-r)/2,v*(1-r)/2,v*r];rhos=[]
    for c in (1,3):
        bounds=[]
        for y in range(4):
            weights=[probabilities[x]*z**(x^y).bit_count() for x in range(4)]
            total=sum(weights,arb(0))
            if not total>0:raise ArithmeticError('tilted packet mass must be positive')
            character=sum((weight*((-1)**((x&c).bit_count()%2)) for x,weight in enumerate(weights)),arb(0))
            bounds.append(min(arb(1),model.up(abs(character)/total)))
        rhos.append(max(bounds))
    value=sum((int(count)*rhos[0]**int(row[1])*rhos[1]**int(row[2])
               for row,count in zip(data['records'],data['multiplicities'])),arb(0))/(1<<data['bits'])
    return min(arb(1),model.up(value))


def prepare(images,columns,state_bits,*,updates=2):
    if type(updates) is not int or not 1<=updates<=32:raise ValueError('integer update count in 1..32 required')
    windows=len(columns)//2
    if len(images)!=1<<state_bits or images[0]!=0:raise ValueError('complete linear expansion images required')
    selected,table,_,multiplicities,_=model.characters(columns,state_bits,1)
    n1=(windows-table[selected.index((0,1))])//2
    n2=(windows-n1-table[selected.index((1,0))]//2)//2
    records=np.column_stack((windows-n1-n2,n1,n2))
    if np.any(records<0) or not np.all(records.sum(axis=1)==windows):raise ArithmeticError('invalid character census')
    return dict(windows=windows,bits=state_bits,updates=updates,records=records,multiplicities=multiplicities,
                histograms=np.array(sorted({model.pattern_histogram(x,windows) for x in images[1:]})))


def _float_output(histograms,v,r,tilt):
    z=np.exp(-tilt)
    factors=np.array([(1-v)+v*(1-r)*z+v*r*z*z,
                      (1-v)*z+v*(1-r)*(1+z*z)/2+v*r*z,
                      (1-v)*z*z+v*(1-r)*z+v*r])
    return float(np.exp(np.max(histograms@np.log(factors))))


def floating(data,v,r,tilt):
    z=np.exp(-tilt);W=data['windows'];size=1<<data['bits']
    records=data['records'];mult=data['multiplicities']
    def fourier(z):
        factors=np.array([(1-v)+v*(1-r)*z+v*r*z*z,
                          (1-v)-v*r*z*z,
                          (1-v)-v*(1-r)*z+v*r*z*z])
        return np.prod(factors[None,:]**records,axis=1)
    total=((1-v)+v*(1-r)*z+v*r*z*z)**W
    zero=float(mult@fourier(z))/size
    spectrum=fourier(1.)
    atom=min(float(mult@np.abs(spectrum)),float(mult@np.abs(spectrum-(1-v)**W)))/size
    moment=_float_output(data['histograms'],v,r,tilt)
    moment2=_float_output(data['histograms'],v,r,2*tilt)
    cancel=min(moment,np.sqrt(max(0.,atom*moment2)),moment*floating_tilted_atom(data,v,r,tilt))
    # Numerical clamps are proposals only; outward() uses enclosures instead.
    zero=max(0.,min(total,zero))
    alpha=2.**-data['updates']
    return np.array([[zero,total-zero],[alpha*cancel+(1-alpha)*moment/(size-1),moment]])


def outward(data,v,r,tilt):
    def aq(x):
        x=Q(x);return arb(x.numerator)/x.denominator
    if not 0<=Q(v)<=1 or not 0<=Q(r)<=1 or Q(tilt)<=0:raise ValueError('reference probabilities and positive tilt required')
    tilted_atom=outward_tilted_atom(data,v,r,tilt)
    v,r,lam=aq(v),aq(r),aq(tilt);z=(-lam).exp();W=data['windows'];size=1<<data['bits']
    records=[tuple(map(int,row)) for row in data['records']]
    mult=list(map(int,data['multiplicities']))
    def fourier(z):
        factors=[(1-v)+v*(1-r)*z+v*r*z*z,(1-v)-v*r*z*z,(1-v)-v*(1-r)*z+v*r*z*z]
        return [factors[0]**n0*factors[1]**n1*factors[2]**n2 for n0,n1,n2 in records]
    total=((1-v)+v*(1-r)*z+v*r*z*z)**W
    zero=sum((count*f for count,f in zip(mult,fourier(z))),arb(0))/size
    spectrum=fourier(arb(1));empty=(1-v)**W
    atom=min(model.up(sum((count*abs(f) for count,f in zip(mult,spectrum)),arb(0))/size),
             model.up(sum((count*abs(f-empty) for count,f in zip(mult,spectrum)),arb(0))/size))
    def moment(z):
        factors=[(1-v)+v*(1-r)*z+v*r*z*z,(1-v)*z+v*(1-r)*(1+z*z)/2+v*r*z,(1-v)*z*z+v*(1-r)*z+v*r]
        return max(model.up(factors[0]**int(h[0])*factors[1]**int(h[1])*factors[2]**int(h[2])) for h in data['histograms'])
    h,h2=moment(z),moment(z*z)
    cancel=min(h,model.up((h2*atom).sqrt()),model.up(h*tilted_atom))
    # Exact boundary kernels may have zero transition mass. Intersect its
    # enclosure with [0,total] instead of rejecting an interval around zero.
    zupper,lupper=model.up(zero),model.up(total-zero)
    if zupper<0 or lupper<0:raise ArithmeticError('negative feedback mass contradicts positivity')
    alpha=arb(2)**-data['updates']
    return arb_mat([[min(model.up(total),zupper),min(model.up(total),lupper)],
                    [model.up(alpha*cancel+(1-alpha)*h/(size-1)),h]])
