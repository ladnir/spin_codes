"""Broaden the sparse tilt bank before attributing a failed bound to the map."""
from fractions import Fraction as F
import argparse
from pathlib import Path
import math
from flint import arb,ctx
import refine
import model
import sparse_bch

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--s',type=int,default=15)
p.add_argument('--q',type=int,nargs='+',default=[32,63])
p.add_argument('--output',type=Path,required=True)
p.add_argument('--seed',type=Path)
a=p.parse_args()
if a.output.exists(): raise FileExistsError(a.output)
ctx.prec=256
engine=refine.Engine(64,a.s)
rows=[]
prior={r['q']:r for r in model.base.base.read(a.seed)['rows']} if a.seed else {}
for q in a.q:
    best=None
    tilts=[F(prior[q]['tilt'])+F(j,40) for j in range(-10,11)] if prior else [F(j,2) for j in range(-12,1)]
    for tilt in tilts:
        region=engine.region(tilt,q)
        ps=sparse_bch.choose(engine,region,q)
        value=model.up(engine.adaptive(region,ps,q)*(engine.cutoff*model.number(tilt).exp()).exp())
        margin=float(-value.log()/arb(2).log())
        if best is None or margin>best['margin_bits']:
            best=dict(q=q,tilt=str(tilt),margin_bits=margin,upper=model.base.pack(value),
                probabilities=list(map(model.base.base.encode,ps)))
        print('retune',a.s,q,str(tilt),round(margin,3),'best',round(best['margin_bits'],3),flush=True)
        engine.region.cache_clear();engine.epoch.cache_clear()
    rows.append(best)
model.base.base.write_new(a.output,dict(status='K16_RETUNED_SPARSE_POINTS',
    instance=engine.identity(),rows=rows,full_distance_proved=False,source_sha256=model.base.sources()))
