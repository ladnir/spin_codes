"""A positive-coefficient bound that can cover whole occupancy intervals.

An auxiliary Bernoulli marking removes the expensive regional coefficient
array. Its conditioning probability is paid once per region. For a frozen
witness, the resulting log bound is convex in the active-pair count q.
"""
import argparse
from fractions import Fraction as Q
from math import comb,log

import numpy as np
from flint import arb,arb_mat,ctx
from scipy.optimize import minimize,minimize_scalar
from scipy.special import gammaln

import model
from probe import aq
from shell_cover import IntervalFolds
from bch_joint_support import authenticated_caps
from support import weighted_cdf_upper,weighted_union_shells
from model import log_power_moment

DENOMINATOR=10**9


def binomial_mixture(operators,p):
    if len(operators)!=65 or not 0<=p<=1:raise ValueError('complete slot degrees and a probability required')
    return sum((comb(64,j)*p**j*(1-p)**(64-j)*t for j,t in enumerate(operators)),arb_mat(9,9))


def float_mixture(operators,p):
    weights=np.array([comb(64,j)*p**j*(1-p)**(64-j) for j in range(65)])
    return np.einsum('r,rij->ij',weights,operators)


def point_proposal(operators,fold,q,tilt):
    def objective(a,p):
        moment=log_power_moment(float_mixture(operators,a*p),64*256,model.terminal())
        choose=log(comb(model.GROUPS,q))
        if q==model.GROUPS and a==1:condition=0.
        else:condition=choose+q*log(a)+(model.GROUPS-q)*np.log1p(-a)
        return moment+float(tilt)*model.THRESHOLD+q*fold(p)+choose-256*condition
    if q==model.GROUPS:
        fit=minimize_scalar(lambda z:objective(1.,1/(1+np.exp(-z))),bounds=(-10.,18.),method='bounded')
        a,p=1.,1/(1+np.exp(-fit.x))
    else:
        a=q/model.GROUPS;best=None
        for p0 in (.5,.8,.95):
            def function(logits):
                probabilities=1/(1+np.exp(-logits))
                return objective(*probabilities)
            start=np.log(np.array([a,p0])/(1-np.array([a,p0])))
            fit=minimize(function,start,method='L-BFGS-B',bounds=[(-12.,18.)]*2,
                         options={'maxiter':100,'ftol':1e-12})
            if best is None or fit.fun<best.fun:best=fit
        a,p=1/(1+np.exp(-best.x))
    def rational(x,endpoint=False):
        if endpoint and x==1:return Q(1)
        return Q(max(1,min(DENOMINATOR-1,round(x*DENOMINATOR))),DENOMINATOR)
    a,p=rational(a,True),rational(p)
    value=objective(float(a),float(p))/log(2)
    if not np.isfinite(value):raise ArithmeticError('nonfinite dense proposal is not proof evidence')
    return value,a,p


def outward_fixed(operators,folds,tilt,a,p):
    """Build one fixed-witness bound evaluable at every compatible occupancy."""
    a,p=Q(a),Q(p)
    if not 0<a<=1 or not 0<p<1:raise ValueError('valid auxiliary and support probabilities required')
    matrix=binomial_mixture(operators,aq(a*p))**(64*256)
    moment=sum((matrix[0,j] for j in range(9) if model.terminal()[j]),arb(0))
    base=model.up(moment*(aq(tilt)*model.THRESHOLD).exp())
    outer=folds.outward(38,256,aq(p))
    def bound(q):
        if not 1<=q<=model.GROUPS or (a==1 and q!=model.GROUPS):
            raise ValueError('occupancy outside fixed-witness domain')
        choose=comb(model.GROUPS,q)
        probability=arb(choose)*aq(a)**q*aq(1-a)**(model.GROUPS-q)
        if not probability>0:raise ArithmeticError('conditioning mass is not positive')
        return model.up(base*choose*outer**q/probability**256)
    return bound


def interval_bound(bound,lo,hi):
    """Log convexity bounds every integer between these fixed-witness endpoints."""
    if not 1<=lo<=hi<=model.GROUPS:raise ValueError('valid occupancy interval required')
    return model.up((hi-lo+1)*max(bound(lo),bound(hi)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,nargs='+',default=[32,64,128,256,512,1024,2048,4095,4096])
    parser.add_argument('--tilts',nargs='+',default=['.032','.064','.096','.128','.192','.256'])
    parser.add_argument('--cutoff',type=int,default=8)
    parser.add_argument('--output-degree',type=int,default=None)
    parser.add_argument('--precision',type=int,default=192)
    parser.add_argument('--outward',action='store_true')
    parser.add_argument('--mass-columns',action='store_true')
    args=parser.parse_args()
    if args.precision<128 or any(not 1<=q<=4096 for q in args.groups) or any(Q(t)<=0 for t in args.tilts):
        parser.error('bounded occupancy, positive tilt and outward precision >=128 required')
    caps=authenticated_caps();ctx.prec=args.precision
    cdf=weighted_cdf_upper(caps,1<<128,rows=2)
    shells=weighted_union_shells(caps,rows=2);shells[0]-=1
    folds=IntervalFolds(cdf,shells);fold=folds.function(38,256)
    data=model.census(args.cutoff)
    if not args.mass_columns:model.prepare_density(data,args.tilts)
    winners={};operators={}
    for tilt in args.tilts:
        ops=model.epoch_operators(data,tilt,output_degree=args.output_degree,mass_columns=args.mass_columns);operators[tilt]=ops
        floating=np.array([[[float(t[i,j]) for j in range(9)] for i in range(9)] for t in ops])
        for q in args.groups:
            value,a,p=point_proposal(floating,fold,q,tilt)
            print('DENSE COMPLETE-SUPPORT PROPOSAL q',q,'tilt',tilt,'log2',value,'a/p',a,p,flush=True)
            if q not in winners or value<winners[q][0]:winners[q]=value,tilt,a,p
    if args.outward:
        for q,(_,tilt,a,p) in winners.items():
            bound=outward_fixed(operators[tilt],folds,tilt,a,p)(q)
            print('DENSE OUTWARD COMPLETE SUPPORT q',q,'tilt',tilt,'a/p',a,p,
                  'log2 upper',model.up(bound.log()/arb(2).log()),flush=True)
    print('Only listed occupancies examined; interval and full-range closure not presumed.',flush=True)


if __name__=='__main__':main()
