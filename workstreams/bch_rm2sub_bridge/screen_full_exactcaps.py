"""Thirteen-band discovery with current exact caps; records pure-band losses."""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import bridge as base
import occupation_three as groups
import screen_exponential_modes as cap_loader
import screen_refined_bands as hull
import scaled_adaptive as scaled
import region_log_cache as cache
from screen_constant_row_split import moment


def search(qs,tilts,tag):
    bands=groups.BANDS;caps,sources=cap_loader.latest_caps()
    weights=[np.array(b) for b in bands]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in b]) for b in bands]
    best={q:(math.inf,None) for q in qs};trials=[]
    for tilt in tilts:
        region=cache.get(tilt);lam=math.exp(tilt/10)
        for q in qs:
            selected=region[:q+1];ps=[];pure=[]
            offset=209716*lam+math.log(math.comb(8192,q))+q*math.log(len(bands))
            for w,v in zip(weights,logs):
                def objective(theta):
                    p=1/(1+math.exp(-theta))
                    gamma=float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                    return moment(selected,q,p)+q*gamma
                opt=minimize_scalar(objective,bounds=(-6,14),method='bounded',options={'xatol':1e-6})
                ps.append(1/(1+math.exp(-float(opt.x))));pure.append(-(float(opt.fun)+offset)/math.log(2))
            ps=np.array(ps)
            gamma=np.array([float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)])
            roots=np.exp(gamma/256);keep=hull.active_lines(ps,roots)
            exponents=np.floor(selected.max(axis=(1,2))/math.log(2)).astype(np.int64)
            mantissas=np.exp(selected-exponents[:,None,None]*math.log(2))
            value=scaled.terminal_log(mantissas,exponents,ps[keep],roots[keep])+offset
            row=dict(occupation=q,witness_tenth=tilt,p=[base.encode(F.from_float(float(p))) for p in ps],
                margin_bits_diagnostic=-value/math.log(2),pure_margins_with_assignment_cost=pure)
            trials.append(row)
            if value<best[q][0]:best[q]=(value,row)
            print('Full Q',q,'tilt',tilt,'margin',row['margin_bits_diagnostic'],'worst pure',min(pure),
                'band',int(np.argmin(pure)),flush=True)
    local=[Path(__file__),Path(groups.__file__),Path(cache.__file__),Path(cap_loader.__file__),Path(hull.__file__),
        Path(scaled.__file__),base.HERE/'screen_constant_row_split.py']+sources
    base.write_new(base.HERE/'generated'/f'full_exactcaps_{tag}_screen.json',dict(status='FULL_EXACTCAPS_SCREEN_ONLY',
        bands=bands,rows=[row for value,row in best.values()],trials=trials,
        used_caps={str(w):caps[w] for w in base.WEIGHTS},cap_receipts=[str(p.relative_to(base.HERE)) for p in sources],
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in local}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--occupancies',type=int,nargs='+',required=True)
    p.add_argument('--tilts',type=int,nargs='+',required=True);p.add_argument('--tag',required=True);a=p.parse_args()
    search(a.occupancies,a.tilts,a.tag)
