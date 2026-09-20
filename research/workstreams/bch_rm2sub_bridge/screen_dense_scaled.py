"""Dense-range screen with Arb region coefficients and binary matrix scales."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import binom
from flint import arb,ctx
import bridge as base
import general_occupancy as general
import tightened_occupancy as tight
import polynomial_regions as poly
import occupation_three as q3
import screen_refined_bands as hull
import christoffel_caps as outer
import scaled_adaptive as scaled


def search(qs,tilts,tag):
    bands=q3.BANDS;ctx.prec=192
    t,s,spectrum=base.load_map(tight.NAME);caps=outer.deterministic_caps();kernel=tight.kernel_spectrum()
    weights=[np.array(band) for band in bands]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in bands]
    best={q:(math.inf,None,None) for q in qs}
    for tenth in tilts:
        lam=(arb(tenth)/10).exp()
        region=poly.regions(t,s,spectrum,kernel,(-lam).exp(),max(qs))
        mantissas,exponents=scaled.initial(region,outward=False)
        with np.errstate(divide='ignore'):logregion=np.log(mantissas)+exponents[:,None,None]*math.log(2)
        print('Dense regions ready, tilt',tenth,flush=True)
        for q in qs:
            ps=[]
            for w,v in zip(weights,logs):
                def pure(theta):
                    p=1/(1+math.exp(-theta))
                    gamma=float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                    terms=binom.logpmf(np.arange(q+1),q,p)[:,None,None]+logregion[:q+1]
                    scale=float(terms.max());matrix=np.exp(terms-scale).sum(axis=0)
                    try:moment=general.log_power(matrix)
                    except AssertionError:return math.inf
                    return moment+256*scale+q*gamma
                opt=minimize_scalar(pure,bounds=(-4,12),method='bounded',options={'xatol':1e-5})
                ps.append(1/(1+math.exp(-float(opt.x))))
            ps=np.array(ps)
            gamma=np.array([float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)])
            roots=np.exp(gamma/256);keep=hull.active_lines(ps,roots)
            value=scaled.terminal_log(mantissas[:q+1],exponents[:q+1],ps[keep],roots[keep])+209716*float(lam)
            if value<best[q][0]:best[q]=(value,tenth,ps.tolist())
            print('Dense Q',q,'tilt',tenth,'margin',-(value+math.log(math.comb(8192,q))+q*math.log(len(bands)))/math.log(2),flush=True)
    rows=[]
    for q,(value,tenth,ps) in best.items():
        rows.append(dict(occupation=q,margin_bits_diagnostic=-(value+math.log(math.comb(8192,q))+q*math.log(len(bands)))/math.log(2),
                         witness_tenth=tenth,p=[base.encode(F.from_float(p)) for p in ps]))
    base.write_new(base.HERE/'generated'/f'dense_scaled_{tag}_screen.json',dict(status='DENSE_SCALED_SCREEN_ONLY',bands=bands,rows=rows,
        source_sha256={p.name:base.sha(p) for p in (Path(__file__),Path(tight.__file__),Path(poly.__file__),Path(general.__file__),Path(q3.__file__),Path(hull.__file__),Path(outer.__file__),Path(scaled.__file__))}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--occupancies',type=int,nargs='+',required=True)
    p.add_argument('--tilts',type=int,nargs='+',required=True);p.add_argument('--tag',required=True)
    a=p.parse_args();search(a.occupancies,a.tilts,a.tag)
