"""Refine common log-odds witnesses for the adaptive range bound."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
import bridge as base
import general_occupancy as general
import adaptive_range as adaptive


def terminal_log(region,ps,roots):
    current=region.copy();scale=0.
    for _ in range(len(region)-1):
        a=current[:-1];b=current[1:]
        updated=roots[0]*((1-ps[0])*a+ps[0]*b)
        for p,r in zip(ps[1:],roots[1:]):np.maximum(updated,r*((1-p)*a+p*b),out=updated)
        normalizer=float(updated.max());current=updated/normalizer;scale+=math.log(normalizer)
    return general.log_power(current[0].copy())+256*scale


def screen(lower,upper):
    t,s,spectrum=base.load_map(general.NAME);caps=general.q2.deterministic_caps();kernel=general.kernel_spectrum()
    weights=[np.array(band) for band in general.BANDS]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in general.BANDS]
    mid=[(band[0]+band[-1])/512 for band in general.BANDS]
    best={q:(math.inf,None,None,None) for q in range(lower,upper+1)}
    for tenth in range(-60,-29):
        lam=math.exp(tenth/10)
        epoch=np.array(general.epoch_matrices(t,s,spectrum,kernel,math.exp(-lam),float,upper)).reshape(-1,3,3)
        region=adaptive.normalized_regions(epoch,upper)
        for q in best:
            def objective(shift,ret=False):
                factor=math.exp(shift);ps=[p*factor/(1-p+p*factor) for p in mid]
                gamma=[float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)]
                value=terminal_log(region[:q+1],ps,np.exp(np.array(gamma)/256))+209716*lam
                return (value,ps) if ret else value
            opt=minimize_scalar(objective,bounds=(-0.5,6.),method='bounded',options={'xatol':0.00005})
            value,ps=objective(float(opt.x),True)
            if value<best[q][0]:best[q]=(value,tenth,float(opt.x),ps)
        if tenth%5==0:print('Refined adaptive tilt',tenth,flush=True)
    rows=[]
    for q,(value,tenth,shift,ps) in best.items():
        margin=-(value+math.log(math.comb(8192,q))+q*math.log(5))/math.log(2)
        rows.append(dict(occupation=q,margin_bits_diagnostic=margin,witness_tenth=tenth,
                         shift_diagnostic=shift,p=[base.encode(F.from_float(p)) for p in ps]))
    base.write_new(base.HERE/'generated'/f'adaptive_q{lower}_q{upper}_refined_screen.json',
                   dict(status='REFINED_ADAPTIVE_SCREEN_ONLY',rows=rows,
                        source_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(adaptive.__file__),Path(general.__file__))}))
    print([(r['occupation'],round(r['margin_bits_diagnostic'],4)) for r in rows],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--lower',type=int,default=17);parser.add_argument('--upper',type=int,default=64)
    args=parser.parse_args();screen(args.lower,args.upper)
