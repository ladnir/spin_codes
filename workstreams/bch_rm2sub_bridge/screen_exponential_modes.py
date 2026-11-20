"""Screen fixed-row-weight mode envelopes; not an occupancy certificate.

Region coefficients and covering mode coefficients are computed outward in
Arb. The final objective is evaluated with floating logs for discovery only.
Every integer input weight is assigned to a covering exponential segment.
"""
import argparse
import math
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from flint import arb,arb_poly,ctx
import bridge as base
import tightened_occupancy as tight
import polynomial_regions as poly
import christoffel_caps as outer


def latest_caps():
    caps=outer.deterministic_caps();sources=[]
    for path in sorted((base.HERE/'generated').glob('joint_shell_*/cap.json')):
        receipt=base.read(path)
        assert receipt['rational_primal_dual_checks_passed']
        for name,digest in receipt['local_sha256'].items():assert base.sha(base.HERE/name)==digest
        for name,digest in receipt['outer_sha256'].items():assert base.sha(base.BCH/name)==digest
        w=receipt['weight'];value=receipt['cap'];assert isinstance(value,int) and value>0
        caps[w]=caps[256-w]=min(caps[w],value);sources.append(path)
    return caps,sources


def fit(region,grid):
    maximum=len(region)-1
    breaks=sorted(set(min(maximum,j) for j in (0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096,5120,6144,7168,7680,8000,8192)))
    if breaks[-1]!=maximum:breaks.append(maximum)
    modes=[[arb(0)]*9 for _ in range(2*grid+1)]
    modes[0]=list(region[0])
    labels=[]
    for lo,hi in zip(breaks,breaks[1:]):
        lo=max(1,lo)
        if hi<lo:continue
        for k in range(9):
            live=[j for j in range(lo,hi+1) if region[j][k]>0]
            if not live:continue
            left,right=live[0],live[-1]
            slope=float((region[right][k]/region[left][k]).log())/(right-left) if right>left else 0.
            rate=math.exp(min(math.log(2),max(-100,slope)))
            label=min(2*grid,max(1,round(rate*grid)))
            r=arb(label)/grid
            coefficient=max((region[j][k]/r**j).upper() for j in live)
            modes[label][k]=max(modes[label][k],coefficient)
            labels.append(dict(lower=lo,upper=hi,entry=k,label=label))
    # All requested integer weights have a covering segment, not just samples.
    assert {0}|{j for row in labels for j in range(row['lower'],row['upper']+1)}==set(range(maximum+1))
    return [tuple(row) for row in modes],labels


def mode_coefficients(modes,positions):
    maximum=(len(modes)-1)*positions
    current=tuple(arb_poly([int(i==j)]) for i in range(3) for j in range(3))
    power=tuple(arb_poly([row[k] for row in modes]) for k in range(9));exponent=positions
    while exponent:
        if exponent&1:current=poly.multiply(current,power,maximum)
        exponent>>=1
        if exponent:power=poly.multiply(power,power,maximum)
    return [(sum((p[b] for p in current[:3]),arb(0))).upper() for b in range(maximum+1)]


def screen(tilt,grid,qs,tag):
    assert grid>0 and all(1<=q<=8192 for q in qs)
    ctx.prec=256;t,s,spectrum=base.load_map(tight.NAME);kernel=tight.kernel_spectrum()
    caps,sources=latest_caps();lam=(arb(tilt)/10).exp()
    region=poly.regions(t,s,spectrum,kernel,(-lam).exp(),8192)
    print('All region coefficients ready',flush=True)
    modes,labels=fit(region,grid)
    print('Modes fitted; active labels',[i for i,row in enumerate(modes) if any(v>0 for v in row)],flush=True)
    coefficients=mode_coefficients(modes,256)
    print('Mode polynomial powered',flush=True)
    terms=[];log_enumerator=[]
    for b,c in enumerate(coefficients):
        if c<=0 or b==0:continue
        x=F(b,256*grid)
        cost=float(logsumexp([math.log(caps[w])+w*math.log(float(x)) for w in base.WEIGHTS]))
        terms.append(float(c.log()));log_enumerator.append(cost)
    terms=np.array(terms);log_enumerator=np.array(log_enumerator)
    rows=[]
    for q in qs:
        logbound=float(logsumexp(terms+q*log_enumerator))+209716*float(lam)+math.log(math.comb(8192,q))
        rows.append(dict(occupation=q,margin_bits_diagnostic=-logbound/math.log(2)))
    print('Mode margins',rows,flush=True)
    base.write_new(base.HERE/'generated'/f'exponential_modes_{tag}_screen.json',dict(status='MODE_SCREEN_ONLY',tilt=tilt,rate_grid_denominator=grid,
        rows=rows,segments=labels,used_caps={str(w):caps[w] for w in base.WEIGHTS},
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
          [Path(__file__),Path(poly.__file__),Path(tight.__file__),Path(outer.__file__)]+sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tilt',type=int,required=True);p.add_argument('--grid',type=int,default=16)
    p.add_argument('--occupancies',type=int,nargs='+',default=[1024,1536,2048,3072,4096,6144,8192]);p.add_argument('--tag',required=True)
    a=p.parse_args();screen(a.tilt,a.grid,a.occupancies,a.tag)
