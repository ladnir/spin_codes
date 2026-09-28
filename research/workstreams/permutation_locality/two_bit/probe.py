"""Replay isolated two-bit bounds; no full-code certificate is claimed.

--one-group covers every message on exactly one active row pair.
Otherwise only the listed homogeneous support events are examined.
Binary64 proposes witnesses; --outward replays selected witnesses in Arb.
"""
import argparse
from fractions import Fraction as Q
from math import comb,log

import numpy as np
from flint import arb,arb_mat,ctx
from scipy.optimize import minimize_scalar

import model
from bch_joint_support import authenticated_caps
from support import weighted_cdf_upper,weighted_union_shells
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_binomial_mass
from model import log_power_moment


def aq(value):
    value=Q(value)
    return arb(value.numerator)/value.denominator


def score(region,q,u,count,tilt):
    def objective(p):
        return (log_power_moment(matrix_for_probabilities(region,[p]*q),256,model.terminal())
                +float(tilt)*model.THRESHOLD-q*log_binomial_mass(256,u,p))
    if u==256:p=Q(1)
    else:
        fit=minimize_scalar(lambda z:objective(1/(1+np.exp(-z))),bounds=(-12.,18.),method='bounded')
        p=Q(max(1,min(10**9-1,round(10**9/(1+np.exp(-fit.x))))),10**9)
    return ((objective(float(p))+q*(log(count.numerator)-log(count.denominator))
             +log(comb(model.GROUPS,q)))/log(2),p)


def replay(region,q,u,count,tilt,p):
    p=aq(p);n=region[0].nrows()
    matrix=sum((comb(q,r)*p**r*(1-p)**(q-r)*region[r] for r in range(q+1)),arb_mat(n,n))**256
    moment=sum((matrix[0,j] for j in range(n) if model.terminal()[j]),arb(0))
    mass=comb(256,u)*p**u*(1-p)**(256-u)
    assert mass>0
    return model.up(moment*(aq(tilt)*model.THRESHOLD).exp()*comb(model.GROUPS,q)*(aq(count)/mass)**q)


def one_group_moments(regions):
    """Exact support coefficient extraction, avoiding a support Chernoff loss."""
    n=regions[0].nrows();current=[arb_mat([[int(i==0) for i in range(n)]])]
    for _ in range(256):
        following=[]
        for u in range(len(current)+1):
            value=arb_mat(1,n)
            if u<len(current):value+=current[u]*regions[0]
            if u:value+=current[u-1]*regions[1]
            following.append(model.rounded(value))
        current=following
    return [model.up(sum((v[0,j] for j in range(n) if model.terminal()[j]),arb(0))/comb(256,u))
            for u,v in enumerate(current)]


def fold(cdf,shells,bounds):
    """Two separately valid outer sums; CDF increments are not shell caps."""
    majorant=list(bounds)
    for u in range(255,-1,-1):majorant[u]=max(majorant[u],majorant[u+1])
    prefix=sum((aq(cdf[u]-(cdf[u-1] if u else 0))*majorant[u] for u in range(257)),arb(0))
    shell=sum((aq(shells[u])*bounds[u] for u in range(257)),arb(0))
    return model.up(min(prefix,shell))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--one-group',action='store_true')
    parser.add_argument('--updates',type=int,default=2)
    parser.add_argument('--tilts',nargs='+')
    parser.add_argument('--penalty',default='1')
    parser.add_argument('--cutoff',type=int,default=8)
    parser.add_argument('--groups',type=int,nargs='+',default=[32,64,80])
    parser.add_argument('--supports',type=int,nargs='+',default=[128,160,192,200,224,256])
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--outward',action='store_true')
    parser.add_argument('--output-degree',type=int,default=None)
    args=parser.parse_args()
    if (not 1<=args.cutoff<=10 or not 1<=args.updates<=32 or not 0<Q(args.penalty)<=1 or args.precision<128
            or any(not 1<=q<=4096 for q in args.groups)
            or any(not 38<=u<=256 for u in args.supports)):
        parser.error('valid geometry, precision >=128, census <=10 and positive penalty <=1 required')
    tilts=args.tilts or (['.00016','.00025','.0004','.00064','.001','.0016','.0025'] if args.one_group
                        else ['.032','.048','.064'])
    if any(Q(t)<=0 for t in tilts):parser.error('positive output tilts required')
    ctx.prec=args.precision
    print('ISOLATED TWO-BIT ensemble: K=2^20, BCH[256,128], IMT(128,19),',args.updates,'updates;',
          '4096 pairs, 256 regions, 64 epochs/region, 64 two-bit slots/epoch.',flush=True)
    caps=authenticated_caps()
    # The older BCH verifier sets the global Arb context while authenticating
    # its certificate. Restore the precision requested for this new replay.
    ctx.prec=args.precision
    cdf=weighted_cdf_upper(caps,1<<128,rows=2,full_weight=1/Q(args.penalty))
    shells=weighted_union_shells(caps,rows=2,full_weight=1/Q(args.penalty));shells[0]-=1
    assert cdf[37]==0 and cdf[-1]>0
    data=model.census(1 if args.one_group else args.cutoff)
    assert ctx.prec==args.precision
    print('Outward replay precision:',ctx.prec,'bits',flush=True)
    if args.one_group:
        best=[arb(1) for _ in range(257)]
        for tilt in tilts:
            ops=model.epoch_operators(data,tilt,args.penalty,maximum=1,updates=args.updates)
            moments=one_group_moments(model.region(ops,1))
            factor=model.up((aq(tilt)*model.THRESHOLD).exp())
            best=[min(old,model.up(value*factor)) for old,value in zip(best,moments)]
            bound=model.up(model.GROUPS*fold(cdf,shells,best))
            print('OUTWARD ONE-GROUP COMPLETE SUPPORT COVER through tilt',tilt,
                  'log2 upper',model.up(bound.log()/arb(2).log()),flush=True)
        print('Covers all 4096 pair locations and every nonzero pair message, only q=1.',flush=True)
        return
    counts={u:min(cdf[u],shells[u]) for u in args.supports};winners={};operators={}
    model.prepare_density(data,tilts)
    for tilt in tilts:
        ops=model.epoch_operators(data,tilt,args.penalty,maximum=min(64,max(args.groups)),output_degree=args.output_degree,updates=args.updates)
        operators[tilt]=ops;regions=model.float_region(ops,max(args.groups))
        for q in args.groups:
            for u in args.supports:
                if not counts[u]:continue
                value,p=score(regions,q,u,counts[u],tilt)
                print('BINARY64 SELECTED EVENT q/u',q,u,'tilt',tilt,'log2',value,'p',p,flush=True)
                if (q,u) not in winners or value<winners[q,u][0]:winners[q,u]=value,tilt,p
    if args.outward:
        for tilt in tilts:
            selected={key:record for key,record in winners.items() if record[1]==tilt}
            if not selected:continue
            regions=model.region(operators[tilt],max(q for q,_ in selected))
            for (q,u),(_,_,p) in selected.items():
                bound=replay(regions,q,u,counts[u],tilt,p)
                print('OUTWARD SELECTED EVENT q/u',q,u,'tilt',tilt,'p',p,
                      'log2 upper',model.up(bound.log()/arb(2).log()),flush=True)
    print('Selected homogeneous supports only. No mixed-support cover or full-code certificate.',flush=True)


if __name__=='__main__':main()
