"""Floating diagnostic: retain the packet count in each shuffled region.

LP proposals, floating conditional operators, and these output files are
not certificates. The production encoder and existing verifiers are unchanged.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from scipy.special import gammaln,logsumexp
from flint import ctx
import birth_classes
import occupancy_birth_classes
import scalar_cover as sc
import variance_partition as variance


def log_choose(n,k):
    return gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)


def placement(local,epochs):
    """All conditional matrices, normalized at each step to avoid overflow."""
    local=np.asarray(local);W=len(local)-1;n=local.shape[1]
    if epochs<1 or W<1 or local.shape!=(W+1,n,n) or np.any(local<0):
        raise ValueError('nonnegative square local operators and positive geometry required')
    result=np.eye(n)[None,:,:]
    for step in range(epochs):
        previous=step*W;total=previous+W
        output=np.zeros((total+1,n,n))
        for j in range(W+1):
            indices=np.arange(j,previous+j+1)
            chance=np.exp(log_choose(W,j)+log_choose(previous,indices-j)-log_choose(total,indices))
            output[indices]+=chance[:,None,None]*(result@local[j])
        result=output
    return result


def scaled_placement(local,epochs):
    """Conditional count matrices with a separate log scale per count."""
    local=np.asarray(local);W=len(local)-1;size=local.shape[1]
    if epochs<1 or W<1 or local.shape!=(W+1,size,size) or np.any(local<0):
        raise ValueError('nonnegative square operators and positive geometry required')
    peaks=local.max(axis=(1,2));valid=peaks>0
    local_logs=np.full(W+1,-np.inf);local_logs[valid]=np.log(peaks[valid])
    normalized=np.zeros_like(local);normalized[valid]=local[valid]/peaks[valid,None,None]
    result=np.eye(size)[None,:,:];logs=np.zeros(1)
    for step in range(epochs):
        previous=step*W;total=previous+W
        output=np.zeros((total+1,size,size));out_logs=np.full(total+1,-np.inf)
        for j in range(W+1):
            if not valid[j]:continue
            indices=np.arange(j,previous+j+1)
            product=result@normalized[j];product_peaks=product.max(axis=(1,2))
            good=(product_peaks>0)&np.isfinite(logs)
            if not np.any(good):continue
            ix=indices[good];products=product[good]/product_peaks[good,None,None]
            term_logs=(logs[good]+local_logs[j]+np.log(product_peaks[good])
                       +log_choose(W,j)+log_choose(previous,ix-j)-log_choose(total,ix))
            high=np.maximum(out_logs[ix],term_logs)
            output[ix]=output[ix]*np.exp(out_logs[ix]-high)[:,None,None]+products*np.exp(term_logs-high)[:,None,None]
            peak=output[ix].max(axis=(1,2));output[ix]/=peak[:,None,None]
            out_logs[ix]=high+np.log(peak)
        result,logs=output,out_logs
    return result,logs


def count_ratio(features,active,mean,interval,minimum_groups,groups,cap):
    """Pointwise PB/binomial ratio cap from floating MGF LPs; diagnostic."""
    features=np.asarray(features,dtype=float);active=np.asarray(active,dtype=float)
    if not 0<mean<1 or not 0<=interval[0]<=interval[1] or cap<=0:
        raise ValueError('interior mean, variance interval, and positive cap required')
    variances=features*(1-features);counts=np.arange(groups+1)
    base=log_choose(groups,counts)+counts*np.log(mean)+(groups-counts)*np.log1p(-mean)
    upper=np.full(groups+1,np.log(cap));successful=0
    for tilt in (-16.,-8.,-4.,-2.,-1.,-.5,-.25,0.,.25,.5,1.,2.,4.,8.,16.):
        with np.errstate(divide='ignore'):
            mgf=np.logaddexp(np.log1p(-features),np.log(features)+tilt)
        fit=linprog(-mgf,A_ub=np.array([variances,-variances,-active]),
            b_ub=[interval[1],-interval[0],-minimum_groups/groups],
            A_eq=np.array([np.ones(len(features)),features]),b_eq=[1.,mean],
            bounds=(0,None),method='highs')
        if fit.success:
            successful+=1
            upper=np.minimum(upper,-groups*fit.fun-tilt*counts-base)
    return upper,successful


def weighted_region(region,log_weights):
    shift=float(np.max(log_weights))
    matrix=np.tensordot(np.exp(log_weights-shift),region,axes=1)
    peak=float(matrix.max())
    if not np.isfinite(peak) or peak<=0:raise ArithmeticError('floating regional moment underflowed')
    # A whole-region moment may be 1e-200. Normalize before its first
    # squaring; scaling only after squaring would silently produce zero.
    return matrix/peak,shift+np.log(peak)


def probe(model,p,tilt_scale=Q(1),baseline=None):
    x=model.tilt*p/(1-p+model.tilt*p);cell=(x,x)
    score,witness=model.proposal(cell) if baseline is None else baseline
    if Q(witness['tilt'])!=model.tilt:
        raise ValueError('diagnostic requires a base-tilt witness')
    if 'variance_partition' not in witness:
        _,witness=variance.propose(model,cell,model.propose_with(cell,model.tilt))
    if tilt_scale<=0:raise ValueError('positive output-tilt scale required')
    lam=Q(witness['parameters'][0])*tilt_scale;z=(-sc.kernel.aq(lam)).exp()
    local=occupancy_birth_classes.outward_at_z(model.data,z)
    n=local[0].nrows()
    arrays=np.array([[[float(m[i,j]) for j in range(n)] for i in range(n)] for m in local])
    region=placement(arrays,64)
    weights,_=model.weights(cell,model.tilt);scale=sum(weights);reference=float(weights[1]/scale)
    counts=np.arange(sc.G+1)
    log_binomial=log_choose(sc.G,counts)+counts*np.log(reference)+(sc.G-counts)*np.log1p(-reference)
    basic,shift=weighted_region(region,log_binomial)
    basic_log=sc.log_power(basic,sc.REGIONS)+sc.REGIONS*shift
    W=len(local)-1;indices=np.arange(W+1)
    local_weights=np.exp(log_choose(W,indices)+indices*np.log(reference)+(W-indices)*np.log1p(-reference))
    iid_log=sc.log_power(np.tensordot(local_weights,arrays,axes=1))
    if abs(iid_log-basic_log)>1e-4:raise ArithmeticError('conditional placement disagrees with direct iid product')
    _,_,logs=model.family(model.tilt);features=np.array(list(map(float,model.features)))
    active=np.array(model.active);parts=[];baseline_terms=[];new_terms=[]
    for interval,dual in variance.validate(cell,witness['variance_partition']):
        lo,hi=map(float,interval);eta,mu,gamma=map(float,dual)
        count=(sc.G*logsumexp(logs+eta*features+mu*active+gamma*features*(1-features))
            -sc.G*(eta*float(x)+min(gamma*lo,gamma*hi))-mu*model.q_min)
        cap=variance.factor(model,cell,interval[0],witness['variance_dual'])
        ratio,success=count_ratio(features,active,float(x),(lo,hi),model.q_min,sc.G,float(cap))
        matrix,shift=weighted_region(region,log_binomial+ratio)
        moment=sc.log_power(matrix,sc.REGIONS)+sc.REGIONS*shift
        baseline_terms.append(count+basic_log+sc.REGIONS*sc.logq(cap));new_terms.append(count+moment)
        parts.append(dict(interval=list(map(str,interval)),lp_successes=success,
            original_count_factor=float(cap),regional_gain_bits=float((basic_log+sc.REGIONS*sc.logq(cap)-moment)/np.log(2))))
    offset=sc.PACKETS*sc.logq(scale)+float(lam)*model.threshold
    return dict(activity=str(p),x=str(x),output_tilt=str(lam),tilt_scale=str(tilt_scale),existing_proposal=score,
        uniform_ratio_proposal=float((logsumexp(baseline_terms)+offset)/np.log(2)),
        count_sensitive_proposal=float((logsumexp(new_terms)+offset)/np.log(2)),parts=parts)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distance',default='19/200')
    parser.add_argument('--activities',nargs='+',default=['.15','.2'])
    parser.add_argument('--tilt-scales',nargs='+',default=['1'])
    parser.add_argument('--variance-bins',type=int,default=8)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();ctx.prec=192
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),birth_classes.actual(),
        int(Q(args.distance)*sc.N),97,tilt=Q(3,16),inner=birth_classes,variance_shuffle=True,variance_bins=args.variance_bins)
    rows=[]
    for p in map(Q,args.activities):
        x=model.tilt*p/(1-p+model.tilt*p);baseline=model.proposal((x,x))
        for scale in map(Q,args.tilt_scales):
            row=probe(model,p,scale,baseline);rows.append(row);print(json.dumps(row,allow_nan=False),flush=True)
            if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,distance=args.distance,
                kernel='gf-birth-classes',minimum_groups=97,base_tilt='3/16',variance_bins=args.variance_bins,
                rows=rows),indent=2,allow_nan=False)+'\n')


if __name__=='__main__':main()
