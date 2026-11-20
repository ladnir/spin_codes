"""Discovery: shared segment rates avoid activation-entry positive slopes.

Every matrix entry on a segment is covered at the chosen common rate.
Modes are not probabilities. Large off-diagonal coefficients are retained.
Only the final all-entry fit and fixed-weight formula matter, not slopes.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from flint import arb,ctx
from scipy.special import logsumexp
import bridge as base
import region_log_cache as cache
import screen_exponential_modes as previous


def run(tilt,grid,tag):
    ctx.prec=256;region=cache.get(tilt).reshape(-1,9);caps,sources=previous.latest_caps()
    breaks=[0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096,5120,6144,7168,7680,8000,8192]
    logs=np.full((2*grid+1,9),-np.inf);logs[0]=region[0];segments=[]
    for start,end in zip(breaks,breaks[1:]):
        start=max(1,start)
        # Z->Z gives the decaying kernel contribution. For its initial zeros,
        # borrow a later negative slope; coverage is checked independently.
        reference_start=max(32,start);reference_end=max(48,end)
        slope=float(np.median(np.diff(region[reference_start:reference_end+1,0])))
        rate=math.exp(min(math.log(2),max(-100,slope)))
        label=max(1,min(2*grid,round(rate*grid)));lograte=math.log(label/grid)
        intercept=np.max(region[start:end+1]-np.arange(start,end+1)[:,None]*lograte,axis=0)
        logs[label]=np.maximum(logs[label],intercept)
        segments.append(dict(lower=start,upper=end,label=label))
    covered=logsumexp(logs[None,:,:]+np.arange(8193)[:,None,None]*np.log(np.maximum(np.arange(2*grid+1)/grid,1e-300))[None,:,None],axis=1)
    # Rate zero contributes only at j=0, not at positive j.
    assert np.all(covered>=region-1e-8)
    modes=[tuple(arb(float(v)).exp() if math.isfinite(v) else arb(0) for v in row) for row in logs]
    print('Shared-slope modes',[(i,round(float(max(row).log()),2)) for i,row in enumerate(modes) if max(row)>0],flush=True)
    coefficients=previous.mode_coefficients(modes,256);terms=[];costs=[]
    for b,c in enumerate(coefficients):
        if c<=0 or b==0:continue
        terms.append(float(c.log()));x=b/(256*grid)
        costs.append(float(logsumexp([math.log(caps[w])+w*math.log(x) for w in base.WEIGHTS])))
    terms=np.array(terms);costs=np.array(costs);rows=[]
    for q in [1280,1536,1792,2048,3072,4096,6144,8192]:
        value=float(logsumexp(terms+q*costs))+209716*math.exp(tilt/10)+math.log(math.comb(8192,q))
        rows.append(dict(occupation=q,margin_bits_diagnostic=-value/math.log(2)))
    print('Shared-slope margins',rows,flush=True)
    base.write_new(base.HERE/'generated'/f'shared_mode_slopes_{tag}_screen.json',dict(status='SHARED_SLOPE_MODE_SCREEN_ONLY',
        tilt_tenth=tilt,grid=grid,segments=segments,rows=rows,
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in [Path(__file__),Path(cache.__file__),Path(previous.__file__)]+sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tilt',type=int,required=True);p.add_argument('--grid',type=int,default=32)
    p.add_argument('--tag',required=True);a=p.parse_args();run(a.tilt,a.grid,a.tag)
