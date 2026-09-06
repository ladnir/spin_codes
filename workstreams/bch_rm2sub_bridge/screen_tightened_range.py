"""Exploratory witnesses for the tightened all-weight occupancy envelope."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
import bridge as base
import general_occupancy as general
import tightened_occupancy as tight
import adaptive_range as adaptive
import refine_adaptive_range as refine


def screen(occupancies,tilts,tag):
    t,s,spectrum=base.load_map(tight.NAME);caps=general.q2.deterministic_caps();kernel=tight.kernel_spectrum()
    weights=[np.array(band) for band in tight.BANDS]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in tight.BANDS]
    mid=[(band[0]+band[-1])/512 for band in tight.BANDS]
    best={q:(math.inf,None,None,None) for q in occupancies}
    for tenth in tilts:
        lam=math.exp(tenth/10)
        epoch=np.array(tight.epoch_matrices(t,s,spectrum,kernel,math.exp(-lam),float,max(occupancies))).reshape(-1,3,3)
        region=adaptive.normalized_regions(epoch,max(occupancies))
        for q in best:
            def objective(shift,ret=False):
                factor=math.exp(shift);ps=[p*factor/(1-p+p*factor) for p in mid]
                gamma=[float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)]
                value=refine.terminal_log(region[:q+1],ps,np.exp(np.array(gamma)/256))+209716*lam
                return (value,ps) if ret else value
            opt=minimize_scalar(objective,bounds=(-.5,6.),method='bounded',options={'xatol':.00005})
            value,ps=objective(float(opt.x),True)
            if value<best[q][0]:best[q]=(value,tenth,float(opt.x),ps)
        print('Tightened tilt',tenth,'margins',[(q,round(-(v[0]+math.log(math.comb(8192,q))+q*math.log(5))/math.log(2),2)) for q,v in best.items()],flush=True)
    rows=[]
    for q,(value,tenth,shift,ps) in best.items():
        margin=-(value+math.log(math.comb(8192,q))+q*math.log(5))/math.log(2)
        rows.append(dict(occupation=q,margin_bits_diagnostic=margin,witness_tenth=tenth,
                         shift_diagnostic=shift,p=[base.encode(F.from_float(p)) for p in ps]))
    base.write_new(base.HERE/'generated'/f'tightened_{tag}_screen.json',
                   dict(status='TIGHTENED_SCREEN_ONLY',rows=rows,
                        source_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(tight.__file__),Path(adaptive.__file__),Path(general.__file__),Path(refine.__file__))}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--occupancies',type=int,nargs='+',required=True)
    parser.add_argument('--tilts',type=int,nargs='+',default=list(range(-60,-19,5)))
    parser.add_argument('--tag',required=True)
    args=parser.parse_args();assert all(1<=q<=8192 for q in args.occupancies)
    screen(args.occupancies,args.tilts,args.tag)
