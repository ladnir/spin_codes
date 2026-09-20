"""Separate all-one rows exactly from ordinary nonzero rows.

d ordinary rows use the 12 nonconstant weight groups; h all-one rows add
exactly h ones to every region. They require no Bernoulli density cost.
The final adaptive maximum still covers all ordinary group mixtures.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from scipy.stats import binom
import bridge as base
import occupation_three as groups
import screen_exponential_modes as cap_loader
import screen_refined_bands as hull
import scaled_adaptive as scaled
import region_log_cache as cache


def log_power(matrix):
    for _ in range(8):matrix=logsumexp(matrix[:,:,None]+matrix[None,:,:],axis=1)
    return float(logsumexp(matrix[0]))


def moment(region,d,p):
    distribution=binom.logpmf(np.arange(d+1),d,p)
    return log_power(logsumexp(distribution[:,None,None]+region,axis=0))


def search(cases,tilts,tag):
    bands=groups.BANDS[:-1];assert groups.BANDS[-1]==(256,)
    caps,sources=cap_loader.latest_caps();weights=[np.array(band) for band in bands]
    logs=[np.array([math.log(caps[w])-math.log(math.comb(256,w)) for w in band]) for band in bands]
    best={case:(math.inf,None,None) for case in cases};cache_sources=[]
    for tilt in tilts:
        region=cache.get(tilt);lam=math.exp(tilt/10)
        cache_sources.append(base.HERE/'generated'/f'region_logs_t128_s15_tilt_{tilt}.json')
        print('Constant split regions ready, tilt',tilt,flush=True)
        for d,h in cases:
            assert 0<=d and 0<=h and 1<=d+h<=8192
            selected=region[h:h+d+1]
            if d:
                ps=[]
                for w,v in zip(weights,logs):
                    def objective(theta):
                        p=1/(1+math.exp(-theta))
                        gamma=float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p)))
                        return moment(selected,d,p)+d*gamma
                    opt=minimize_scalar(objective,bounds=(-6,6),method='bounded',options={'xatol':1e-5})
                    ps.append(1/(1+math.exp(-float(opt.x))))
                ps=np.array(ps)
                gamma=np.array([float(np.max(v-w*math.log(p)-(256-w)*math.log1p(-p))) for v,w,p in zip(logs,weights,ps)])
                roots=np.exp(gamma/256);keep=hull.active_lines(ps,roots)
                exponents=np.floor(selected.max(axis=(1,2))/math.log(2)).astype(np.int64)
                mantissas=np.exp(selected-exponents[:,None,None]*math.log(2))
                value=scaled.terminal_log(mantissas,exponents,ps[keep],roots[keep])+d*math.log(len(bands))
            else:
                ps=[];value=log_power(selected[0])
            value+=209716*lam+math.log(math.comb(8192,d))+math.log(math.comb(8192-d,h))
            if value<best[d,h][0]:best[d,h]=(value,tilt,list(ps))
            print('Constant split d,h',d,h,'tilt',tilt,'margin',-value/math.log(2),flush=True)
    rows=[dict(ordinary_rows=d,all_one_rows=h,occupation=d+h,margin_bits_diagnostic=-value/math.log(2),
               witness_tenth=tilt,p=[base.encode(F.from_float(float(p))) for p in ps]) for (d,h),(value,tilt,ps) in best.items()]
    base.write_new(base.HERE/'generated'/f'constant_split_{tag}_screen.json',dict(status='CONSTANT_SPLIT_SCREEN_ONLY',rows=rows,bands=bands,
        used_caps={str(w):caps[w] for w in base.WEIGHTS},local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
            [Path(__file__),Path(cache.__file__),Path(groups.__file__),Path(cap_loader.__file__),Path(hull.__file__),Path(scaled.__file__)]+sources+cache_sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cases',nargs='+',required=True,help='ordinary_rows:all_one_rows')
    p.add_argument('--tilts',nargs='+',type=int,required=True);p.add_argument('--tag',required=True);a=p.parse_args()
    search([tuple(map(int,s.split(':'))) for s in a.cases],a.tilts,a.tag)
