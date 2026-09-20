"""Interior-slope exponential modes; explicit coverage of every region weight.

This avoids letting exceptional endpoint drops set a whole segment's rate.
Fitting remains discovery; final floating objectives are not certificates.
"""
import argparse
import math
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from flint import arb,ctx
import bridge as base
import tightened_occupancy as tight
import polynomial_regions as poly
import screen_exponential_modes as previous


def fit(region,grid):
    maximum=len(region)-1
    breaks=sorted(set(min(maximum,j) for j in (0,1,2,3,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024,1536,2048,3072,4096,5120,6144,7168,7680,8000,8192)))
    modes=[[arb(0)]*9 for _ in range(2*grid+1)];modes[0]=list(region[0]);labels=[]
    for start,end in zip(breaks,breaks[1:]):
        start=max(1,start)
        for k in range(9):
            live=[j for j in range(start,end+1) if region[j][k]>0]
            if not live:continue
            slopes=[float((region[j+1][k]/region[j][k]).log()) for j in range(start,end)
                    if region[j][k]>0 and region[j+1][k]>0]
            slope=float(np.median(slopes)) if slopes else -math.log(grid)
            rate=math.exp(min(math.log(2),max(-100,slope)))
            label=min(2*grid,max(1,round(rate*grid)));r=arb(label)/grid
            coefficient=max((region[j][k]/r**j).upper() for j in live)
            modes[label][k]=max(modes[label][k],coefficient)
            labels.append(dict(lower=start,upper=end,entry=k,label=label))
    assert {0}|{j for row in labels for j in range(row['lower'],row['upper']+1)}==set(range(maximum+1))
    return [tuple(row) for row in modes],labels


def screen(tilt,grid,qs,tag):
    ctx.prec=256;t,s,spectrum=base.load_map(tight.NAME);caps,sources=previous.latest_caps();lam=(arb(tilt)/10).exp()
    region=poly.regions(t,s,spectrum,tight.kernel_spectrum(),(-lam).exp(),8192)
    print('Interior fit: region coefficients ready',flush=True)
    modes,segments=fit(region,grid)
    print('Interior mode log maxima',[(i,round(float(max(row).log()),2)) for i,row in enumerate(modes) if max(row)>0],flush=True)
    coefficients=previous.mode_coefficients(modes,256)
    terms=[];costs=[]
    for b,c in enumerate(coefficients):
        if c<=0 or b==0:continue
        x=b/(256*grid)
        terms.append(float(c.log()));costs.append(float(logsumexp([math.log(caps[w])+w*math.log(x) for w in base.WEIGHTS])))
    terms=np.array(terms);costs=np.array(costs);rows=[]
    for q in qs:
        logbound=float(logsumexp(terms+q*costs))+209716*float(lam)+math.log(math.comb(8192,q))
        rows.append(dict(occupation=q,margin_bits_diagnostic=-logbound/math.log(2)))
    print('Interior mode margins',rows,flush=True)
    base.write_new(base.HERE/'generated'/f'exponential_modes_v2_{tag}_screen.json',dict(status='INTERIOR_MODE_SCREEN_ONLY',tilt=tilt,rate_grid_denominator=grid,
        rows=rows,segments=segments,used_caps={str(w):caps[w] for w in base.WEIGHTS},
        local_sha256={str(p.relative_to(base.HERE)):base.sha(p) for p in
            [Path(__file__),Path(previous.__file__),Path(poly.__file__),Path(tight.__file__),Path(previous.outer.__file__)]+sources}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--tilt',type=int,required=True);p.add_argument('--grid',type=int,default=16)
    p.add_argument('--occupancies',type=int,nargs='+',default=[1024,1536,2048,3072,4096,6144,8192]);p.add_argument('--tag',required=True)
    a=p.parse_args();screen(a.tilt,a.grid,a.occupancies,a.tag)
