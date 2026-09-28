"""Homogeneous-histogram SPIN diagnostic, with categorical conditioning.

All numerical inner arithmetic is binary64. Profiles are selected points,
not a complete histogram cover. No certificate is claimed.
"""
import argparse
from math import log,comb
from pathlib import Path
from hashlib import sha256
import numpy as np
from scipy.optimize import minimize_scalar,minimize
from scipy.special import gammaln,xlogy
from flint import ctx,fmpz_poly
from categorical_model import build,mix,normalization_test
from occupancy_sensitivity import (prepare,local_data,fresh_census,prepare_inputs,zero_census,
    averages,authenticated_caps,dimension_caps,support_caps,improve_caps,float_placement)
from occupancy_memory import TERMINAL
from occupancy_screen import matrix_for_probabilities
from two_group_screen import log_power_moment,log_binomial_mass
from shortening_moments import improve
from odd_column_caps import marginal
from occupancy_allones import full_column_caps


def code_hash():
    here=Path(__file__).parent
    paths=[here/name for name in ('categorical_model.py','occupancy_memory.py','occupancy_fresh_moment.py',
        'fresh_collision.py','zero_moment.py','full_feedback_census.py','group_moment.py',
        'occupancy_model.py','feedback_convolution.py','random_group_verify.py','two_group_verify.py',
        'two_column_moment.py','two_group_moment.py','window_feedback.py','feedback_character_census.py',
        'pair_tail.py','occupancy_window_average.py')]
    paths += [here.parents[2]/'spin/src/kernels/generated/SelectedMaps.h',
              here.parents[1]/'workstreams/inner_design/NO_CONSTANT_MAP.json']
    return sha256(b''.join(p.read_bytes() for p in paths)).hexdigest()


def profile_cap(profile,spectrum,caps,dimensions):
    u=sum(profile);v=profile[0]+profile[2];w=sum(i*n for i,n in enumerate(profile,1));j=profile[3]
    total=sum(row[u] for row in caps)
    odd=sum(row[v] for row in marginal(spectrum,caps,u,dimensions=dimensions))
    weight=int((fmpz_poly(spectrum)**4)[w])
    allone=full_column_caps(spectrum,caps,dimensions,u,total)[j]
    return min(total,odd,weight,allone)


def score(data,groups,profile,theta,tilt,outer,use_maximum=False):
    u=sum(profile)
    ops=np.array([t.max(axis=0) for _,t,_ in data]) if use_maximum else mix(data,theta)
    region=float_placement(ops,groups)
    log_conditional=(0 if use_maximum else gammaln(u+1)-sum(gammaln(n+1) for n in profile)+sum(xlogy(n,p) for n,p in zip(profile,theta)))
    def objective(z):
        p=1/(1+np.exp(-z))
        value=log_power_moment(matrix_for_probabilities(region,[p]*groups),256,TERMINAL)
        return value+tilt*209715-groups*(log_binomial_mass(256,u,p)+log_conditional)
    fit=minimize_scalar(objective,bounds=(-8,16),method='bounded')
    return (fit.fun+outer)/log(2),1/(1+np.exp(-fit.x))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=int,default=64)
    parser.add_argument('--tilt',default='.032')
    parser.add_argument('--cache')
    parser.add_argument('--load-cache')
    parser.add_argument('--profile',type=int,nargs=4,action='append')
    parser.add_argument('--optimize-theta',action='store_true')
    args=parser.parse_args();ctx.prec=192
    profiles=args.profile or [(32,56,32,8),(64,0,64,0),(0,128,0,0),(96,0,32,0),(32,80,16,0)]
    assert all(0<sum(row)<=256 and min(row)>=0 for row in profiles)
    if args.load_cache:
        with np.load(args.load_cache,allow_pickle=False) as saved:
            assert str(saved['code_hash'])==code_hash() and str(saved['tilt'])==args.tilt
            data=[(saved[f's{j}'],saved[f't{j}'],saved[f'c{j}']) for j in range(33)]
    else:
        prepared=prepare(local_data(4));fresh=fresh_census(prepared)
        inputs=prepare_inputs();window=averages(inputs,args.tilt);zeros=zero_census()
        from full_feedback_census import census
        full=census(6)
        data=build(prepared,fresh,window,zeros,full,float(args.tilt))
        # Compare per-shape minima against the existing shape-maximized
        # nine-coordinate pipeline. Mature-tail lifting is intentionally absent.
        from occupancy_sensitivity import fresh_refine,window_refine,multi_refine,collision_refine,zero_refine,transform
        from full_feedback_refinement import refine as full_refine
        old=fresh_refine(prepared,fresh,args.tilt,32)
        old=window_refine(prepared,window,args.tilt,32,base=old)
        old=multi_refine(old,fresh,args.tilt)
        old=collision_refine(old,prepared,args.tilt,1)
        old=zero_refine(old,zeros,args.tilt,1)
        old=transform(old,inputs[3],args.tilt,2)
        old=full_refine(old,full,inputs[3],args.tilt,1,2)
        from categorical_model import asarray
        for j,(_,matrices,_) in enumerate(data):
            assert (matrices<=asarray(old[j])[None,:,:]+1e-13).all(),j
        print('All per-shape pilot coefficients are below the existing nine-coordinate envelope',flush=True)
        if args.cache:
            values={f'{prefix}{j}':value for j,row in enumerate(data) for prefix,value in zip(('s','t','c'),row)}
            np.savez_compressed(args.cache,**values,code_hash=np.array(code_hash()),tilt=np.array(args.tilt))
    normalization_test(data)
    from categorical_model import state_test
    state_test(data,prepare_inputs(),float(args.tilt))
    spectrum=authenticated_caps();dimensions=dimension_caps()
    caps,_=improve(improve_caps(support_caps(spectrum,g=4,dimensions=dimensions),spectrum),dimensions)
    for profile in profiles:
        cap=profile_cap(profile,spectrum,caps,dimensions)
        if not cap:
            print('Profile excluded by exact outer caps',profile,flush=True);continue
        theta=np.array(profile)/sum(profile)
        outer=args.groups*log(cap)+log(comb(2048,args.groups))
        value,p=score(data,args.groups,profile,theta,float(args.tilt),outer)
        print('CATEGORICAL PILOT profile/log2-count/log2-bound/p',profile,log(cap)/log(2),value,p,flush=True)
        baseline,bp=score(data,args.groups,profile,theta,float(args.tilt),outer,use_maximum=True)
        print('MATCHED NINE-COORDINATE MAXIMUM',baseline,'p',bp,'averaging gain bits',baseline-value,flush=True)
        if args.optimize_theta:
            active=np.flatnonzero(np.array(profile)>0)
            def distribution(z):
                values=np.exp(np.append(z,0)-max(0,max(z,default=0)))
                result=np.zeros(4);result[active]=values/values.sum()
                return result
            if len(active)>1:
                initial=np.log(theta[active[:-1]]/theta[active[-1]])
                fit=minimize(lambda z:score(data,args.groups,profile,distribution(z),float(args.tilt),outer)[0],
                             initial,method='Nelder-Mead',options={'maxiter':60,'xatol':1e-4,'fatol':1e-3})
                print('OPTIMIZED THETA',distribution(fit.x),'log2-bound',fit.fun,
                      'converged',fit.success,'evaluations',fit.nfev,flush=True)
    print('Selected homogeneous profiles; binary64 pilot without mature-tail coordinates. No full certificate.',flush=True)


if __name__=='__main__':main()
