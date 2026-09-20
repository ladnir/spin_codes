"""Full-occupancy pure-band diagnostics, not a mixed-type certificate.

At Q=L the auxiliary Bernoulli inputs are iid throughout each epoch. Thus
the pure-band transfer is the epoch mixture raised to 256*L/t, avoiding
degree-L polynomials. All final values here are binary64 diagnostics.
"""
import argparse
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from scipy.stats import binom
from flint import arb,ctx

import certificate_search_core as core
import certificate_search_backend as backend
import tightened_occupancy as tight

base=core.base
BANDS=backend.sparse.BANDS
LN2=math.log(2)


def logmul(a,b):
    return logsumexp(a[:,:,None]+b[None,:,:],axis=1)


def logpower(a,n):
    result=np.full_like(a,-np.inf);np.fill_diagonal(result,0.)
    while n:
        if n&1:result=logmul(result,a)
        n>>=1
        if n:a=logmul(a,a)
    return result


def epoch_logs(t,s,spectrum,kernel,tilt):
    ctx.prec=256
    lam=(arb(tilt)/10).exp()
    epoch=tight.epoch_matrices(t,s,spectrum,kernel,(-lam).exp(),arb,t)
    return np.array([[float(v.log()) if v>0 else -np.inf for v in row] for row in epoch]).reshape(-1,3,3)


def band_profile(logs,band,p,caps,rows,t,cutoff,tilt):
    w=np.array(band)
    if band==(256,):cost=0.;shell=256
    else:
        values=np.array([math.log(caps[int(v)])-math.log(math.comb(256,int(v))) for v in w])
        costs=values-w*math.log(p)-(256-w)*math.log1p(-p)
        index=int(np.argmax(costs));cost=float(costs[index]);shell=int(w[index])
    mixture=logsumexp(binom.logpmf(np.arange(t+1),t,p)[:,None,None]+logs,axis=0)
    count=256*rows//t
    final=logpower(mixture,count)[0]
    moment=float(logsumexp(final));zero_only=count*float(mixture[0,0])
    correction=cutoff*math.exp(tilt/10)
    value=moment+rows*cost+correction
    return dict(band=list(band),p=p,tilt=tilt,dominant_cost_shell=shell,
        row_cost_bits=rows*cost/LN2,moment_bits=moment/LN2,chernoff_cost_bits=correction/LN2,
        margin_bits=-value/LN2,never_activated_fraction=math.exp(min(0.,zero_only-moment)),
        never_activated_log_fraction_bits=(zero_only-moment)/LN2,
        final_state_fractions=np.exp(final-moment).tolist(),
        epoch_log_matrix=mixture.tolist())


def optimize(logs,band,caps,rows,t,cutoff,tilt):
    if band==(256,):return band_profile(logs,band,1.,caps,rows,t,cutoff,tilt)
    # Multiple starts on a fixed coarse mesh avoid assuming a unimodal objective.
    def profile(theta):return band_profile(logs,band,1/(1+math.exp(-theta)),caps,rows,t,cutoff,tilt)
    mesh=np.linspace(-6,14,41)
    values=[profile(float(x)) for x in mesh]
    index=max(range(len(values)),key=lambda i:values[i]['margin_bits'])
    lo,hi=float(mesh[max(0,index-1)]),float(mesh[min(len(mesh)-1,index+1)])
    fit=minimize_scalar(lambda x:-profile(float(x))['margin_bits'],bounds=(lo,hi),method='bounded',
                        options={'xatol':1e-6})
    core.require(fit.success,'Pure-band optimization failed')
    return max((values[index],profile(float(fit.x))),key=lambda v:v['margin_bits'])


def run(name,m,tilts,output):
    core.require(not output.exists(),'Use a fresh diagnostic output')
    spec=core.instance(name,m);t,s,spectrum,kernel=core.inputs.load(name)
    caps=core.inputs.caps_module.caps();started=time.monotonic();results=[]
    for tilt in tilts:
        logs=epoch_logs(t,s,spectrum,kernel,tilt)
        profiles=[optimize(logs,band,caps,spec['rows'],t,spec['cutoff'],tilt) for band in BANDS]
        worst=min(profiles,key=lambda r:r['margin_bits'])
        results.append(dict(tilt=tilt,bands=profiles,worst_band=worst['band'],
            worst_pure_margin_bits=worst['margin_bits'],
            pure_proxy_after_label_count_bits=worst['margin_bits']-spec['rows']*math.log2(len(BANDS))))
        print(name,'tilt',tilt,'worst band',worst['band'],'margin',round(worst['margin_bits'],3),
              'zero fraction',round(worst['never_activated_fraction'],6),flush=True)
    hashes=backend.sources();hashes[Path(__file__).relative_to(base.ROOT).as_posix()]=base.sha(Path(__file__))
    base.write_new(output,dict(status='PURE_BAND_BINARY64_DIAGNOSTIC_ONLY',instance=spec,
        full_coverage=False,results=results,seconds=time.monotonic()-started,source_sha256=hashes,
        scope='Only homogeneous band assignments. Neither their maximum nor a label-count correction certifies mixed assignments.'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--configuration',required=True);p.add_argument('--m',type=int,default=20)
    p.add_argument('--tilts',nargs='+',type=int,default=[-5,0,2,4,6,8])
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.configuration,a.m,a.tilts,a.output)
