"""Pure-shell diagnostics, not a bound for arbitrary mixed shell sequences."""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
import bridge as base
import region_log_cache as cache
import screen_exponential_modes as caps_source
from screen_constant_row_split import moment


def run(q,tilts,tag):
    caps,sources=caps_source.latest_caps();results=[]
    for tilt in tilts:
        region=cache.get(tilt)[:q+1];lam=math.exp(tilt/10);rows=[]
        offset=209716*lam+math.log(math.comb(8192,q))
        for w in base.WEIGHTS:
            logcap=math.log(caps[w])-math.log(math.comb(256,w))
            def objective(theta):
                p=1/(1+math.exp(-theta))
                return moment(region,q,p)+q*(logcap-w*math.log(p)-(256-w)*math.log1p(-p))
            opt=minimize_scalar(objective,bounds=(-6,14),method='bounded',options={'xatol':1e-5})
            rows.append(dict(weight=w,p=1/(1+math.exp(-float(opt.x))),margin_bits=-(float(opt.fun)+offset)/math.log(2)))
        worst=sorted(rows,key=lambda r:r['margin_bits'])[:6]
        print('Pure singleton Q',q,'tilt',tilt,'NO assignment cost',[(r['weight'],round(r['margin_bits'],2)) for r in worst],flush=True)
        results.append(dict(tilt_tenth=tilt,rows=rows))
    base.write_new(base.HERE/'generated'/f'single_shells_{tag}_screen.json',dict(status='PURE_SHELL_DIAGNOSTIC_ONLY',
        occupation=q,results=results,used_caps={str(w):caps[w] for w in base.WEIGHTS},
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),Path(cache.__file__),
        Path(caps_source.__file__),base.HERE/'screen_constant_row_split.py']+sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--occupancy',type=int,required=True)
    p.add_argument('--tilts',nargs='+',type=int,required=True);p.add_argument('--tag',required=True);a=p.parse_args()
    run(a.occupancy,a.tilts,a.tag)
