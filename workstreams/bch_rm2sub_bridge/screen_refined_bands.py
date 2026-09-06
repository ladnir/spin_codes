"""Screen finer bands with a floating upper-hull acceleration.

Hull pruning is discovery-only. Certificates must maximize over every band.
"""
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
import occupation_three as q3
import refine_adaptive_range as refine


def active_lines(ps,roots):
    lines=sorted(((r*p,r*(1-p),i) for i,(p,r) in enumerate(zip(ps,roots))))
    hull=[];starts=[]
    for slope,intercept,index in lines:
        if hull and slope==hull[-1][0]:
            if intercept<=hull[-1][1]:continue
            hull.pop();starts.pop()
        start=-math.inf
        while hull:
            start=(hull[-1][1]-intercept)/(slope-hull[-1][0])
            if start>starts[-1]:break
            hull.pop();starts.pop()
        if not hull:start=-math.inf
        hull.append((slope,intercept,index));starts.append(start)
    return [line[2] for k,line in enumerate(hull) if k+1==len(hull) or starts[k+1]>=0]


def screen(qs,tilts,mode,tag):
    bands=q3.BANDS if mode=='13' else tuple((w,) for w in base.WEIGHTS)
    t,s,spectrum=base.load_map(tight.NAME);caps=general.q2.deterministic_caps();kernel=tight.kernel_spectrum()
    weights=[np.array(band) for band in bands]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in bands]
    midpoint=np.array([min((band[0]+band[-1])/512,1-1e-7) for band in bands])
    best={q:(math.inf,None,None) for q in qs}
    for tenth in tilts:
        lam=math.exp(tenth/10)
        epoch=np.array(tight.epoch_matrices(t,s,spectrum,kernel,math.exp(-lam),float,max(qs))).reshape(-1,3,3)
        region=adaptive.normalized_regions(epoch,max(qs))
        for q in qs:
            def objective(shift,ret=False):
                factor=math.exp(shift);ps=midpoint*factor/(1-midpoint+midpoint*factor)
                gamma=np.array([float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)])
                roots=np.exp(gamma/256);keep=active_lines(ps,roots)
                value=refine.terminal_log(region[:q+1],ps[keep],roots[keep])+209716*lam
                return (value,ps.tolist()) if ret else value
            opt=minimize_scalar(objective,bounds=(-.5,6.),method='bounded',options={'xatol':.00005})
            value,ps=objective(float(opt.x),True)
            if value<best[q][0]:best[q]=(value,tenth,ps)
        print('Bands',mode,'tilt',tenth,'best',[(q,round(-(v[0]+math.log(math.comb(8192,q))+q*math.log(len(bands)))/math.log(2),2)) for q,v in best.items()],flush=True)
    rows=[]
    for q,(value,tenth,ps) in best.items():
        rows.append(dict(occupation=q,margin_bits_diagnostic=-(value+math.log(math.comb(8192,q))+q*math.log(len(bands)))/math.log(2),
                         witness_tenth=tenth,p=[base.encode(F.from_float(p)) for p in ps]))
    base.write_new(base.HERE/'generated'/f'refined_bands_{tag}_screen.json',dict(status='REFINED_BANDS_SCREEN_ONLY',bands=bands,rows=rows,
        source_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(tight.__file__),Path(adaptive.__file__),Path(general.__file__),Path(refine.__file__),Path(q3.__file__))}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--occupancies',type=int,nargs='+',required=True)
    p.add_argument('--tilts',type=int,nargs='+',required=True);p.add_argument('--mode',choices=('13','singleton'),default='13')
    p.add_argument('--tag',required=True);a=p.parse_args();screen(a.occupancies,a.tilts,a.mode,a.tag)
