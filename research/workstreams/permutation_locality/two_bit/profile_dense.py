"""Dense two-bit diagnostic retaining both packet types by conditioning.

Artificial iid packet types are a coefficient witness, not a new encoder.
The scalar envelope is deliberately independent of the feedback dynamics:
it bounds the tilted output for every fixed incoming state. Listed boxes
restrict *every* active pair to the indicated union-support interval.
They are not a complete mixed-support or full-code cover.
"""
import argparse
from fractions import Fraction as Q
from math import comb,log

import numpy as np
from flint import arb,ctx
from scipy.optimize import minimize
from scipy.special import logsumexp

import model
from bch_joint_support import authenticated_caps
from profiles import pair_profiles,conditioned_shells
from probe import aq
import iid_kernel


def float_epoch(histograms,v,r,tilt):
    z=np.exp(-tilt)
    factors=np.array([(1-v)+v*(1-r)*z+v*r*z*z,
                      (1-v)*z+v*(1-r)*(1+z*z)/2+v*r*z,
                      (1-v)*z*z+v*(1-r)*z+v*r])
    return float(np.max(histograms@np.log(factors)))


def outward_epoch(histograms,v,r,tilt):
    v,r=aq(v),aq(r);z=(-aq(tilt)).exp()
    factors=[(1-v)+v*(1-r)*z+v*r*z*z,
             (1-v)*z+v*(1-r)*(1+z*z)/2+v*r*z,
             (1-v)*z*z+v*(1-r)*z+v*r]
    return max(model.up(factors[0]**h[0]*factors[1]**h[1]*factors[2]**h[2]) for h in histograms)


def profile_arrays(profiles,lo,hi):
    selected=[(u,b,c) for (u,b),c in profiles.items() if lo<=u<=hi]
    us=np.array([u for u,_,_ in selected]);bs=np.array([b for _,b,_ in selected])
    # Retain exact counts until taking logarithms; this is only a proposal.
    logs=np.array([log(c.numerator)-log(c.denominator)-log(comb(256,u)*comb(u,b)) for u,b,c in selected])
    return us,bs,logs


def point_proposal(histograms,profiles,q,lo,hi,stateful=None,density=False):
    us,bs,logs=profile_arrays(profiles,lo,hi);choose=log(comb(4096,q))
    def objective(w):
        lam,a,p,r=w
        if not (lam>0 and 0<a<=1 and 0<p<1 and 0<r<1):return np.inf
        if a==1 and q!=4096:return np.inf
        condition=choose+q*log(a)+(4096-q)*np.log1p(-a) if a<1 else 0.
        profile_costs=logs-us*log(p)-(256-us)*np.log1p(-p)-bs*log(r)-(us-bs)*np.log1p(-r)
        outer=np.max(profile_costs) if density else logsumexp(profile_costs)
        inner=(16384*float_epoch(histograms,a*p,r,lam) if stateful is None else
               model.log_power_moment(iid_kernel.floating(stateful,a*p,r,lam),16384,np.ones(2)))
        return inner+lam*model.THRESHOLD+choose+q*outer-256*condition
    best=None
    for p0,r0 in ((.75,1/3),(.9,.5),(.99,.75),(.5,.25)):
        # Optimize probabilities directly; fixed q=4096 uses a=1 exactly.
        if q==4096:
            fit=minimize(lambda x:objective([x[0],1,x[1],x[2]]),[1.5,p0,r0],method='L-BFGS-B',
                         bounds=[(.00001,5),(.0001,.999999),(.0001,.9999)],options={'maxiter':150,'ftol':1e-12})
            w=[fit.x[0],1,fit.x[1],fit.x[2]]
        else:
            fit=minimize(objective,[1.5,q/4096,p0,r0],method='L-BFGS-B',
                         bounds=[(.00001,5),(.00001,.999999),(.0001,.999999),(.0001,.9999)],
                         options={'maxiter':150,'ftol':1e-12})
            w=fit.x
        if best is None or fit.fun<best[0]:best=fit.fun,w
    witness=[Q(round(float(x)*10**9),10**9) for x in best[1]]
    value=objective(list(map(float,witness)))/log(2)
    if not np.isfinite(value):raise ArithmeticError('nonfinite profile proposal')
    return value,witness


def outward_fixed(histograms,profiles,lo,hi,witness,stateful=None,density=False):
    lam,a,p,r=map(Q,witness)
    if not (lam>0 and 0<a<=1 and 0<p<1 and 0<r<1):raise ValueError('invalid frozen witness')
    if density:
        support_masses={u:comb(256,u)*aq(p)**u*aq(1-p)**(256-u) for u in range(lo,hi+1)}
        one=[aq(1-r)**a for a in range(257)];two=[aq(r)**b for b in range(257)]
        outer=max(model.up(aq(count)/(support_masses[u]*comb(u,b)*two[b]*one[u-b]))
                  for (u,b),count in profiles.items() if lo<=u<=hi)
    else:
        shells=conditioned_shells(profiles,256,r)
        outer=model.up(sum((aq(shells[u])/(comb(256,u)*aq(p)**u*aq(1-p)**(256-u)) for u in range(lo,hi+1)),arb(0)))
    if stateful is None:inner=outward_epoch(histograms,a*p,r,lam)**16384
    else:
        matrix=iid_kernel.outward(stateful,a*p,r,lam)**16384
        inner=matrix[0,0]+matrix[0,1]
    base=model.up((aq(lam)*model.THRESHOLD).exp()*inner)
    def bound(q):
        if not 1<=q<=4096 or (a==1 and q!=4096):raise ValueError('occupancy outside witness domain')
        choose=comb(4096,q);condition=choose*aq(a)**q*aq(1-a)**(4096-q)
        if not condition>0:raise ArithmeticError('nonpositive conditioning mass')
        return model.up(base*choose*outer**q/condition**256)
    return bound


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,nargs='+',default=[512,1024,2048,4096])
    parser.add_argument('--interval',type=int,nargs=2,action='append',default=None)
    parser.add_argument('--outward',action='store_true')
    parser.add_argument('--stateful',action='store_true',help='Retain zero-state activation and the two-update lazy/refresh split')
    parser.add_argument('--density-envelope',action='store_true',help='Bound the complete pair input measure, avoiding a sum of profile conditioning charges')
    args=parser.parse_args()
    intervals=args.interval or [(128,128),(160,160),(192,192),(224,224),(256,256),(38,256)]
    if any(not 1<=q<=4096 for q in args.groups) or any(not 38<=lo<=hi<=256 for lo,hi in intervals):
        parser.error('valid occupancy and nonzero-pair support intervals required')
    caps=authenticated_caps();ctx.prec=192
    profiles=pair_profiles(caps);profiles.pop((0,0))
    images,columns,_=model.maps();histograms=sorted({model.pattern_histogram(x) for x in images})
    stateful=iid_kernel.prepare(images,columns,19) if args.stateful else None
    print('PROFILE DENSE: all',len(histograms),'expansion histograms including zero;',len(profiles),'nonzero pair profiles',flush=True)
    floating=np.array(histograms)
    for q in args.groups:
        for lo,hi in intervals:
            value,witness=point_proposal(floating,profiles,q,lo,hi,stateful,args.density_envelope)
            print('PROFILE BOX q/support',q,(lo,hi),'binary64 log2',value,'lambda/a/p/r',list(map(str,witness)),flush=True)
            if args.outward:
                bound=outward_fixed(histograms,profiles,lo,hi,witness,stateful,args.density_envelope)(q)
                print('PROFILE BOX OUTWARD q/support',q,(lo,hi),'log2 upper',model.up(bound.log()/arb(2).log()),flush=True)
    print('Selected support boxes only; no mixed-support cover or full-code certificate.',flush=True)


if __name__=='__main__':main()
