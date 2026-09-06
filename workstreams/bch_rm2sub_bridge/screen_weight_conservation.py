"""Discovery only: preserve total input weight across region envelopes.

R_j <= exp(theta*j) E on j in [h,h+d]. Since the total of the
256 region weights equals the sum of the BCH row weights, a cell's
enumerator is bounded by exp(256*h*theta) S_ord(exp(theta))**d.
No Bernoulli comparison or group assignment loss is needed.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
import bridge as base
import region_log_cache as cache
import screen_exponential_modes as cap_loader
from screen_constant_row_split import log_power


def search(cases,tilts,tag):
    caps,sources=cap_loader.latest_caps()
    weights=np.array([w for w in base.WEIGHTS if 0<w<256])
    logs=np.array([math.log(caps[int(w)]) for w in weights])
    best={case:(math.inf,None,None) for case in cases}
    for tilt in tilts:
        region=cache.get(tilt);lam=math.exp(tilt/10)
        for d,h in cases:
            assert 0<=d and 0<=h and 1<=d+h<=8192
            selected=region[h:h+d+1];j=np.arange(h,h+d+1)
            def objective(theta):
                envelope=np.max(selected-theta*j[:,None,None],axis=0)
                spectrum=float(logsumexp(logs+weights*theta))
                # The dimension bound remains valid with either sign of theta.
                spectrum=min(spectrum,math.log((1<<128)-2)+max(weights*theta))
                return log_power(envelope)+256*h*theta+d*spectrum
            grid=np.linspace(-3,3,121)
            values=[objective(theta) for theta in grid]
            i=int(np.argmin(values))
            opt=minimize_scalar(objective,bounds=(grid[max(i-1,0)],grid[min(i+1,120)]),method='bounded')
            value=float(opt.fun)+209716*lam+math.log(math.comb(8192,d))+math.log(math.comb(8192-d,h))
            print('Conserved weight d,h',d,h,'tilt',tilt,'theta',opt.x,'margin',-value/math.log(2),flush=True)
            if value<best[d,h][0]:best[d,h]=(value,tilt,float(opt.x))
    base.write_new(base.HERE/'generated'/f'weight_conservation_{tag}_screen.json',dict(
        status='WEIGHT_CONSERVATION_SCREEN_ONLY',rows=[dict(ordinary_rows=d,all_one_rows=h,
        margin_bits_diagnostic=-v/math.log(2),tilt_tenth=t,theta=theta) for (d,h),(v,t,theta) in best.items()],
        used_caps={str(w):caps[w] for w in base.WEIGHTS},
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),Path(cache.__file__),Path(cap_loader.__file__)]+sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--cases',nargs='+',required=True)
    p.add_argument('--tilts',type=int,nargs='+',required=True);p.add_argument('--tag',required=True);a=p.parse_args()
    search([tuple(map(int,s.split(':'))) for s in a.cases],a.tilts,a.tag)
