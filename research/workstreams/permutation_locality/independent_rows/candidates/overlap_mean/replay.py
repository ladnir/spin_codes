"""Selected-point overlap refinement using an already built local memo.

This driver does not generate a missing memo. It applies the new bound
after loading the unchanged source-bound baseline family, then searches
fixed convex column witnesses and reconstructs them with outward arithmetic.
The result concerns one homogeneous support vector, not a full occupancy.
"""
import argparse
import hashlib
from fractions import Fraction as Q
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from flint import ctx

from mean import build,zero_refine
from spectral_feedback import build as feedback_build
from mass_density import blend
from mass_density_screen import as_array,score,baseline
from mass_optimize import family,objective
from local_sensitivity import placement_tape
from occupancy_memory import Z,C
import local_family


def refinement_sources():
    digest=hashlib.sha256()
    for path in sorted(Path(__file__).resolve().parent.glob('*.py')):
        digest.update(path.name.encode()+b'\0'+path.read_bytes()+b'\0')
    return digest.hexdigest()


def refined_family(directory,tilt,penalty,density_chord=False):
    sources=local_family.source_digest(); refinement=refinement_sources()
    key=local_family.parameters(tilt,penalty,192,16,10,0,True,10,sources)
    new_key=dict(key,stage='independent-row-overlap-family',overlap_sources=refinement,
                 overlap_minimum=5,overlap_maximum=10,density_chord=bool(density_chord))
    cached=local_family.load(directory,new_key)
    if cached is not None:
        print('OVERLAP FAMILY exact refined memo hit',tilt,'density chord',density_chord,flush=True)
        return cached
    cached=local_family.load(directory,key)
    if cached is None and not density_chord:
        raise ValueError('matching exact10/joint4/density10 local memo is not ready; run mass_optimize first')
    if cached is None:
        print('Baseline memo still pending; computing the independent density refinement first',flush=True)
    overlap=build(10)
    if density_chord:
        from density import build as density_build,refine as density_refine
        from group_moment import maps
        records,feedback=density_build(10,[tilt],overlap=overlap)
    else:
        feedback=feedback_build(10)
    if local_family.source_digest()!=sources or refinement_sources()!=refinement:
        raise RuntimeError('generator sources changed during refinement; no memo or proof replay')
    if cached is None:
        cached=local_family.load(directory,key)
    if cached is None:
        raise ValueError('local memo still unavailable after refinement; no proof replay')
    base,coefficients=cached
    changed=zero_refine(base,overlap,feedback,tilt,penalty)
    if density_chord:
        changed=density_refine(changed,records[tilt],maps()[2],penalty,minimum=5,maximum=10)
    local_family.save(directory,new_key,(changed,coefficients))
    print('OVERLAP FAMILY saved exact refined memo',tilt,'density chord',density_chord,flush=True)
    return changed,coefficients


def optimize(base,coefficients,q,support,count,tilt,iterations):
    maximum=16
    ones={j:Q(1) for j in range(1,maximum+1)}
    zero=blend(base,coefficients,ones,target=Z)
    density=blend(base,coefficients,ones,target=C)
    arrays=tuple(map(as_array,(base,zero,density)))
    fractions=np.zeros((2,32));fractions[:,5:12]=1
    best,p=score(placement_tape(family(*arrays,fractions),q)[-1],q,support,count,tilt)
    print('OVERLAP initial fixed ZC-6-12 score',best,flush=True)
    bounds=[(0.,1.) if j<maximum else (0.,0.) for _ in range(2) for j in range(32)]
    for turn in range(2):
        def evaluate(vector):
            value,gradient=objective(*arrays,vector.reshape(2,32),q,p,baseline.TAIL_TERMINAL)
            return value,gradient.ravel()
        fit=minimize(evaluate,fractions.ravel(),jac=True,bounds=bounds,method='L-BFGS-B',
                     options={'maxiter':iterations,'ftol':1e-12,'gtol':1e-8})
        proposal=np.clip(fit.x.reshape(2,32),0,1)
        candidate,new_p=score(placement_tape(family(*arrays,proposal),q)[-1],q,support,count,tilt)
        if candidate<best: fractions,best,p=proposal,candidate,new_p
        print('OVERLAP optimized proposal',turn+1,'score',best,'iterations',fit.nit,flush=True)
    rational=[{j+1:Q(round(float(f)*10**6),10**6) for j,f in enumerate(row)
               if j<maximum and round(float(f)*10**6)>0} for row in fractions]
    print('OVERLAP fixed zero/density witnesses',[{j:str(f) for j,f in row.items()} for row in rational],flush=True)
    result=blend(base,coefficients,rational[0],target=Z)
    return blend(result,coefficients,rational[1],target=C)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--epoch-cache',required=True)
    parser.add_argument('--groups',type=int,default=80)
    parser.add_argument('--support',type=int,default=200)
    parser.add_argument('--tilt',default='.072')
    parser.add_argument('--penalty',default='.9')
    parser.add_argument('--iterations',type=int,default=30)
    parser.add_argument('--density-chord',action='store_true',help='Also refine translated density with the overlap budget')
    args=parser.parse_args();ctx.prec=192
    if not 1<=args.groups<=2048 or not 38<=args.support<256 or args.iterations<1:
        parser.error('invalid selected-point parameters')
    changed,coefficients=refined_family(args.epoch_cache,args.tilt,args.penalty,args.density_chord)
    caps=baseline.authenticated_caps()
    cdf=baseline.integer_cdf(baseline.weighted_cdf_upper(caps,1<<128,full_weight=1/Q(args.penalty)))
    shells=baseline.weighted_union_shells(caps,full_weight=1/Q(args.penalty));shells[0]-=1
    count=min(Q(cdf[args.support]),shells[args.support])
    selected=optimize(changed,coefficients,args.groups,args.support,count,args.tilt,args.iterations)
    exact=baseline.placement(selected,rounding=baseline.rounded,maximum_groups=args.groups)
    label=args.penalty+':overlap-mean:optimized-mass'
    if args.density_chord:label+=':density-chord'
    args.probe_supports=[args.support];args.target_bits=52
    baseline.shell_points(args,{(args.tilt,label):(exact,as_array(exact))},{label:cdf},{label:shells})
    print('Other supports and occupancies remain; not a full-code certificate.',flush=True)


if __name__=='__main__':
    main()
