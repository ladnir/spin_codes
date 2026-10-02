"""Probe a conservative r1-to-r2 bridge for Bernoulli-only dense bounds.

For a nonzero source row, put c>=1 for the unnormalized Fourier cap and
m=2^s-1. The coefficient ratio is
  (c/(4(m+1))+3/(4m))/(c/(2(m+1))+1/(2m)) <= (4m+3)/(4m+2).
Zero rows are unchanged. Positive products of E epochs cost at most f^E.
This does NOT justify transferring a fixed-weight or arbitrary r1 bound.
"""
import argparse
from fractions import Fraction as F
from pathlib import Path
from flint import arb,ctx
import model
import subspace_cover


def run(a):
    if a.output.exists(): raise FileExistsError(a.output)
    ctx.prec=256
    maps=model.base.base.read(a.maps)
    record=next(r['inner'] for r in maps['rows'] if r['s']==12 and r['found'])
    c=subspace_cover.Checker(record)
    seed=model.base.base.read(a.seed)
    model.base.authenticate(seed)
    assert seed['instance']==c.engine.identity()
    m=c.engine.m
    factor=F(4*m+3,4*m+2)
    log_penalty=model.up((c.engine.output_bits//c.engine.t)*model.number(factor).log())
    rows=[]
    for key,node in sorted(seed['leaves'].items(),key=lambda kv:kv[1]['power'],reverse=True)[:a.count]:
        geometry=subspace_cover.budget_dense.geometry.geometry(node)
        # Both expressions use only Bernoulli moments with the unchanged zero
        # activation alternative. Never call bound(), which also tries r1-only
        # fixed-weight moment bounds.
        ordinary=c.fixed_input_bound(*geometry,node['witness'])
        split=c.split_bound(*geometry,node['witness'])
        candidates=[x for x in (ordinary,split) if x is not None]
        if not candidates: raise ValueError('No Bernoulli bound at this box')
        upper=(min(candidates)+log_penalty).exp()
        row=dict(key=key,old_power=node['power'],upper=model.base.pack(upper),
                 margin_bits=float(-upper.log()/arb(2).log()))
        rows.append(row)
        print('dense bridge',key,'old power',node['power'],'new margin',row['margin_bits'],flush=True)
    instance=c.engine.identity()
    instance['inner']['transvection_rounds']=2
    instance['setup']='independent row/region permutations; two independent transvections before feedback per update; zero start; output before update; persistent state; no flush'
    result=dict(status='K16_TARGET49_BERNOULLI_BRIDGE_PROBE',instance=instance,
        factor_per_epoch=model.base.base.encode(factor),penalty_bits=float(log_penalty/arb(2).log()),
        rows=rows,full_distance_proved=False,source_sha256=model.base.sources())
    for path in (Path(__file__),a.maps,a.seed,Path(subspace_cover.__file__)):
        result['source_sha256'][path.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(path)
    model.base.base.write_new(a.output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--count',type=int,default=10)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args())
