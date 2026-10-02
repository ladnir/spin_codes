"""GF16 epoch envelope for a fixed number of active packet positions.

Packet positions are a uniform subset of the epoch; each selected packet
is independently uniform nonzero. Coordinates are those of REFRESH_BOUND.
No iid activity approximation or worst packet-weight choice is used here.
"""
from collections import Counter
from fractions import Fraction as Q
from math import comb
import numpy as np
from flint import arb,arb_mat,ctx
import refresh_kernel as refresh

aq,up=refresh.aq,refresh.up
TERMINAL=np.ones(3)


def prepare(data,exact_feedback=False):
    W=data['windows'];S=1<<data['bits']
    counts=Counter()
    for row,n in zip(data['records'],data['multiplicities']):counts[int(row[0])]+=int(n)
    zero=[];atoms=[]
    for q in range(W+1):
        fourier=[(sum((Q(comb(c,j)*comb(W-c,q-j))*(-Q(1,15))**(q-j)
                       for j in range(max(0,q-W+c),min(c,q)+1)),Q(0))/comb(W,q),n)
                 for c,n in counts.items()]
        z=sum((f*n for f,n in fourier),Q(0))/S
        # A constant character term vanishes at every nonzero syndrome.
        # A weighted median minimizes this exact l1 majorant.
        seen=0;median=Q(0)
        for f,n in sorted(fourier):
            seen+=n
            if 2*seen>=S:median=f;break
        b=min(Q(1),sum((abs(f)*n for f,n in fourier),Q(0))/S,
              sum((abs(f-median)*n for f,n in fourier),Q(0))/S)
        if not 0<=z<=1 or not 0<=b<=1:raise ArithmeticError('invalid Fourier probability bound')
        zero.append(z);atoms.append(b)
    exact=[]
    if exact_feedback:
        from feedback_exact import census
        exact=census(data['columns'],data['bits'],min(W,8))
        for row in exact:
            q=row['occupancy'];z=Q(row['zero'],row['denominator']);b=Q(row['nonzero_peak'],row['denominator'])
            if z!=zero[q] or b>atoms[q]:raise ArithmeticError('exact feedback census disagrees with Fourier bounds')
            atoms[q]=b
    return dict(data,zero_probabilities=zero,nonzero_atom_caps=atoms,exact_feedback=exact)


def actual(updates=2,exact_feedback=False):return prepare(refresh.actual(updates),exact_feedback)


def polynomial(histogram,z):
    """All fixed-occupancy output moments for one expansion histogram."""
    result=[arb(1)]
    whole=(1+z)**4
    for w,n in enumerate(histogram):
        inactive=z**w;active=(whole-inactive)/15
        for _ in range(int(n)):
            new=[arb(0)]*(len(result)+1)
            for q,value in enumerate(result):
                new[q]+=inactive*value;new[q+1]+=active*value
            result=new
    W=sum(map(int,histogram))
    return [value/comb(W,q) for q,value in enumerate(result)]


def moments(data,z):
    rows=[polynomial(h,z) for h in data['histograms']]
    L=(1<<data['bits'])-1
    maximum=[max(up(row[q]) for row in rows) for q in range(data['windows']+1)]
    mean=[up(sum((int(n)*row[q] for n,row in zip(data['histogram_multiplicities'],rows)),arb(0))/L)
          for q in range(data['windows']+1)]
    return maximum,mean


def outward_at_z(data,z):
    if not 0<z<=1:raise ValueError('output weight in (0,1] required')
    L=(1<<data['bits'])-1;alpha=arb(2)**-data['updates'];beta=1-alpha
    h,mean=moments(data,z);h2,mean2=moments(data,z*z)
    packet=((1+z)**4-1)/15
    packet2=((1+z*z)**4-1)/15
    distance=min(sum(w*int(n) for w,n in enumerate(row)) for row in data['histograms'])
    result=[]
    for q,(zero,atom) in enumerate(zip(data['zero_probabilities'],data['nonzero_atom_caps'])):
        zero=aq(zero);atom=aq(atom);nonzero=1-zero
        mass=packet**q
        zz=min(up(mass),up(z**q*zero),up((packet2**q*zero).sqrt()))
        zn=up(mass-z**(4*q)*zero)
        pointwise=z**max(0,distance-4*q)
        cancel=min(h[q],up((h2[q]*atom).sqrt()),up(atom*pointwise))
        uniform_cancel=min(mean[q],up(nonzero*pointwise/L),up((mean2[q]*nonzero/L).sqrt()))
        fresh=min(h[q],up((h2[q]*nonzero).sqrt()))
        fresh_mean=min(mean[q],up((mean2[q]*nonzero).sqrt()))
        matrix=arb_mat([[zz,zn,0],
                        [up(alpha*cancel+beta*fresh/L),up(alpha*h[q]),up(beta*h[q])],
                        [up(alpha*uniform_cancel+beta*fresh_mean/L),up(alpha*mean[q]),up(beta*mean[q])]])
        if any(matrix[i,j]<0 for i in range(3) for j in range(3)):
            raise ArithmeticError('negative occupancy envelope entry')
        result.append(matrix)
    return result


def outward(data,tilt):
    if Q(tilt)<=0:raise ValueError('positive output tilt required')
    return outward_at_z(data,(-aq(tilt)).exp())


def build_operators(args):
    from occupancy_model import placement
    from occupancy_memory import rounded
    data=actual(2,getattr(args,'exact_feedback',False));ctx.prec=args.precision;result={}
    for tilt in args.tilts:
        local=outward(data,tilt)
        exact=placement(local,rounding=rounded,maximum_groups=args.groups)
        arrays=[np.array([[float(matrix[i,j]) for j in range(3)] for i in range(3)]) for matrix in exact]
        result[tilt,'1']=exact,arrays
        print('GF16 fixed-occupancy operators',tilt,'degree',args.groups,flush=True)
    return result
