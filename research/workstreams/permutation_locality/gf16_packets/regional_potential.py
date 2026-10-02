"""Diagnostic regional composition of shared-mass future-cost bounds.

Keeps exact conditional packet counts, with floating arithmetic. Positive
potentials propose a regional super-solution, not a distance certificate.
"""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from flint import ctx
import mass_density_potential as potential
import regional_count_probe as placement


def conditional(bound,values,epochs=64):
    """Backward cost bounds conditioned on each exact regional count."""
    v=np.asarray(values,dtype=float)
    if (type(epochs) is not int or epochs<1 or v.shape!=(bound.n,)
            or not np.isfinite(v).all() or (v<=0).any()):raise ValueError('positive vector and step count required')
    W=bound.W;result=(v/v.max())[None,:];logs=np.array([np.log(v.max())])
    for step in range(epochs):
        previous=step*W;total=previous+W
        output=np.zeros((total+1,bound.n));out_logs=np.full(total+1,-np.inf)
        for j in range(W+1):
            indices=np.arange(j,previous+j+1)
            product=bound.floating(j,result);peaks=product.max(axis=1)
            good=(peaks>0)&np.isfinite(logs)
            if not good.any():continue
            ix=indices[good];products=product[good]/peaks[good,None]
            term_logs=(logs[good]+np.log(peaks[good])+placement.log_choose(W,j)
                +placement.log_choose(previous,ix-j)-placement.log_choose(total,ix))
            high=np.maximum(out_logs[ix],term_logs)
            output[ix]=output[ix]*np.exp(out_logs[ix]-high)[:,None]+products*np.exp(term_logs-high)[:,None]
            peak=output[ix].max(axis=1);output[ix]/=peak[:,None];out_logs[ix]=high+np.log(peak)
        result,logs=output,out_logs
    return result,logs


def weighted_image(region,logs,log_weights):
    terms=logs+log_weights;shift=float(np.max(terms))
    image=np.exp(terms-shift)@region;peak=float(image.max())
    if peak<=0 or not np.isfinite(image).all():raise ArithmeticError('nonpositive regional image')
    return image/peak,shift+np.log(peak)


def propose(bound,weighted,offset,start,regions=256,epochs=64,iterations=6):
    """Each variance part may retain its own positive rational potential."""
    if type(iterations) is not int or iterations<1:raise ValueError('positive proposal iteration count required')
    v=np.asarray(start,dtype=float);best=[None]*len(weighted)
    for iteration in range(iterations):
        region,logs=conditional(bound,v,epochs);images=[];scores=[]
        for i,(weights,count) in enumerate(weighted):
            image,shift=weighted_image(region,logs,weights)
            log_rho=shift+np.log(np.max(image/v))
            value=count+regions*log_rho+np.log(v[0]/v.min())
            if best[i] is None or value<best[i]['term']:
                coordinate=int(np.argmax(image/v))
                with np.errstate(divide='ignore'):
                    contributions=weights+logs+np.log(region[:,coordinate])
                mass=np.exp(contributions-logsumexp(contributions));cdf=np.cumsum(mass)
                best[i]=dict(term=float(value),log_rho=float(log_rho),iteration=iteration,
                    potential=[str(Q(float(x))) for x in v],outer_count_log=float(count),
                    bottleneck_coordinate=coordinate,peak_count=int(np.argmax(contributions)),
                    mean_count=float(np.arange(len(mass))@mass),
                    count_quantiles=[int(np.searchsorted(cdf,q)) for q in (.1,.5,.9)])
            images.append(image);scores.append(best[i]['term'])
        score=float((logsumexp(scores)+offset)/np.log(2))
        print('SHARED MASS REGIONAL iteration',iteration,'log2',score,flush=True)
        # Refine the currently dominant variance part; retain independent
        # previous proposals for every other part.
        v=images[int(np.argmax(scores))]
        if (v<=0).any() or not np.isfinite(v).all():raise ArithmeticError('positive regional potential lost')
    return score,best


def main():
    import birth_classes
    import fiber_density
    import scalar_cover as sc
    import regional_count as regional
    import regional_frontier_probe as frontier
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path);parser.add_argument('--activities',nargs='+',default=['.4'])
    parser.add_argument('--precision',type=int,default=256);parser.add_argument('--iterations',type=int,default=6)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();source=json.loads(args.input.read_text())
    if source.get('schema')!='gf16-regional-count-screen-1' or not source.get('diagnostic_only'):
        parser.error('regional floating screen required')
    if args.precision<128:parser.error('precision at least 128 required')
    data=fiber_density.attach(birth_classes.actual(source['updates']));ctx.prec=args.precision
    model=sc.Model(sc.actual_components(129,Q(21,50),Q(1001,1000),True),data,source['threshold'],
        source['minimum_groups'],tilt=Q(3,16),inner=birth_classes,variance_shuffle=True,variance_bins=16,regional_count=True)
    rows=[]
    for p in map(Q,args.activities):
        matching=[row for row in source['rows'] if Q(row['activity'])==p]
        if not matching:parser.error('requested activity missing')
        row=min(matching,key=lambda r:r['proposal']);witness=row['witness'];lam=Q(witness['parameters'][0])
        cell=tuple(map(Q,row['cell']));parts,checked=regional.prepare_witness(model,cell,witness)
        baseline=regional.local_operators(model,witness)
        raw={k:v for k,v in witness.items() if k not in (
            'regional_lazy_density_through','regional_feedback_classes_from','regional_feedback_classes_through',
            'regional_feedback_uniform_classes','regional_feedback_uniform_replace')}
        original=regional.local_operators(model,raw)
        bound=potential.Bound(data,original,(-potential.aq(lam)).exp(),[baseline])
        weighted,_,scale=frontier.weighted_parts(model,cell,witness,parts,checked)
        offset=sc.PACKETS*sc.logq(scale)+float(lam)*model.threshold
        # Confirm the same regional normalization and saved baseline first.
        existing=frontier.evaluate(bound.float_linear[1],weighted,offset)
        if abs(existing-row['proposal'])>1e-5:raise ArithmeticError('saved baseline mismatch')
        initial=potential.iid_proposal(bound,float(p))['potential']
        score,proposals=propose(bound,weighted,offset,initial,iterations=args.iterations)
        result=dict(activity=str(p),cell=row['cell'],witness=checked,baseline=existing,proposal=score,potentials=proposals)
        rows.append(result);print('SHARED MASS REGIONAL result',str(p),'log2',score,flush=True)
        if args.output:args.output.write_text(json.dumps(dict(diagnostic_only=True,partial_only=True,
            scope='selected regional mean intervals, floating potentials',precision=args.precision,
            threshold=source['threshold'],minimum_groups=source['minimum_groups'],updates=source['updates'],
            rows=rows),indent=2)+'\n')
    print('Floating regional diagnostic only; no outward or full-domain certificate.',flush=True)


if __name__=='__main__':main()
