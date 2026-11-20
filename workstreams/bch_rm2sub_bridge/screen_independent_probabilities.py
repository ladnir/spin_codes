"""Choose independent per-band Bernoulli witnesses before adaptive maximization.

Pure-group optimization selects witnesses only. The reported final bound
still uses the adaptive maximum over every band and pays the assignment count.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import binom
import bridge as base
import general_occupancy as general
import tightened_occupancy as tight
import adaptive_range as adaptive
import occupation_three as q3
import refine_adaptive_range as refine
import screen_refined_bands as hull
import christoffel_caps as outer


def search(qs,tilts,tag):
    bands=q3.BANDS
    t,s,spectrum=base.load_map(tight.NAME);caps=outer.deterministic_caps();kernel=tight.kernel_spectrum()
    weights=[np.array(band) for band in bands]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in bands]
    best={q:(math.inf,None,None) for q in qs}
    for tenth in tilts:
        lam=math.exp(tenth/10)
        epoch=np.array(tight.epoch_matrices(t,s,spectrum,kernel,math.exp(-lam),float,max(qs))).reshape(-1,3,3)
        region=adaptive.normalized_regions(epoch,max(qs))
        for q in qs:
            logits=[]
            for w,v in zip(weights,logs):
                def pure(theta):
                    p=1/(1+math.exp(-theta))
                    gamma=float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                    matrix=np.sum(binom.pmf(np.arange(q+1),q,p)[:,None,None]*region[:q+1],axis=0)
                    scale=float(matrix.max())
                    if scale==0:return math.inf
                    try:moment=general.log_power(matrix/scale)
                    except AssertionError:return math.inf
                    return moment+256*math.log(scale)+q*gamma
                opt=minimize_scalar(pure,bounds=(-4,6),method='bounded',options={'xatol':1e-5})
                logits.append(float(opt.x))
            logits=np.array(logits)
            def objective(shift,ret=False):
                ps=1/(1+np.exp(-logits-shift))
                gamma=np.array([float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)])
                roots=np.exp(gamma/256);keep=hull.active_lines(ps,roots)
                value=refine.terminal_log(region[:q+1],ps[keep],roots[keep])+209716*lam
                return (value,ps.tolist()) if ret else value
            opt=minimize_scalar(objective,bounds=(-.25,.25),method='bounded',options={'xatol':1e-5})
            value,ps=min((objective(float(opt.x),True),objective(0,True)),key=lambda pair:pair[0])
            if value<best[q][0]:best[q]=(value,tenth,ps)
        print('Independent-p tilt',tenth,'best',[(q,round(-(v[0]+math.log(math.comb(8192,q))+q*math.log(len(bands)))/math.log(2),2)) for q,v in best.items()],flush=True)
    rows=[]
    for q,(value,tenth,ps) in best.items():
        rows.append(dict(occupation=q,margin_bits_diagnostic=-(value+math.log(math.comb(8192,q))+q*math.log(len(bands)))/math.log(2),
                         witness_tenth=tenth,p=[base.encode(F.from_float(p)) for p in ps]))
    base.write_new(base.HERE/'generated'/f'independent_p_{tag}_screen.json',dict(status='INDEPENDENT_P_SCREEN_ONLY',bands=bands,rows=rows,
        source_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(tight.__file__),Path(adaptive.__file__),Path(general.__file__),Path(refine.__file__),Path(q3.__file__),Path(hull.__file__),Path(outer.__file__))}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--occupancies',type=int,nargs='+',required=True)
    p.add_argument('--tilts',type=int,nargs='+',required=True);p.add_argument('--tag',required=True)
    a=p.parse_args();search(a.occupancies,a.tilts,a.tag)
