"""Recompute full sparse/dense components for the two-round K16 inner."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from types import FunctionType, SimpleNamespace as NS
from flint import ctx
import model
import subspace_cover
import target49_rounds_model as rounds
import sparse_ranges
import budget_dense


def sources(a):
    result=model.base.sources()
    for path in (Path(__file__),Path(rounds.__file__),Path(subspace_cover.__file__),a.maps,a.seed):
        result[path.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(path)
    return result


def sparse(a,record):
    out=a.output.resolve()
    saved=model.base.base.read(out) if a.verify else None
    if saved:
        model.base.authenticate(saved)
        assert saved['status']=='K16_TWO_ROUND_SPARSE_COMPONENT'
    elif out.exists(): raise FileExistsError(out)
    ctx.prec=512 if saved else 256
    e=rounds.Engine(record)
    if saved: assert saved['instance']==e.identity()
    source=model.base.base.read(a.seed)
    model.base.authenticate(source)
    jobs=list(saved['witnesses'] if saved else source['witnesses'])
    best,witnesses={},[]
    for job in jobs:
        tilt=F(job['tilt'])
        ps=list(map(model.base.base.decode,job['probabilities']))
        rows=sparse_ranges.evaluate(e,e.region(tilt,63),ps,tilt,2,63)
        assert [r['occupation'] for r in rows]==list(range(2,64))
        if saved:
            assert all(r['power']<=old['power'] for r,old in zip(rows,job['rows']))
            rows=job['rows']
        witnesses.append(dict(tilt=str(tilt),probabilities=job['probabilities'],rows=rows))
        for row in rows:
            q,power=row['occupation'],row['power']
            best[q]=min(best.get(q,power),power)
        e.region.cache_clear();e.epoch.cache_clear()
        print('two-round sparse',len(witnesses),'/',len(jobs),'below -55',sum(p<=-55 for p in best.values()),flush=True)
    assert set(best)==set(range(2,64))
    total=sum((F(2)**p for p in best.values()),F(0))
    result=dict(status='K16_TWO_ROUND_SPARSE_COMPONENT',instance=e.identity(),
        covered_occupancies=[2,63],witnesses=witnesses,best=best,
        upper=model.base.base.encode(total),margin_bits=math.log2(total.denominator)-math.log2(total.numerator),
        full_distance_proved=False,source_sha256=sources(a))
    print('sparse margin',result['margin_bits'],flush=True)
    if saved:
        assert result['upper']==saved['upper']
        result=dict(status='K16_TWO_ROUND_SPARSE_512_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(out),full_distance_proved=False)
        out=out.with_name(out.stem+'_replay.json')
    model.base.base.write_new(out,result)


def run(a):
    records=model.base.base.read(a.maps)['rows']
    record=next(r['inner'] for r in records if r['s']==12 and r['found'])
    if a.mode=='sparse': return sparse(a,record)
    original=budget_dense.short_dense.mixed_dense.ladder.candidate.sources
    def pinned_sources(): return {**original(),**sources(a)}
    adapter=NS(Checker=lambda exponent:rounds.Checker(record),
        mixed_dense=NS(ladder=NS(candidate=NS(sources=pinned_sources))))
    driver=FunctionType(budget_dense.run.__code__,{**budget_dense.run.__globals__,'short_dense':adapter})
    driver(a.output.resolve(),16,64,a.seed.resolve(),a.nodes,a.seconds,54,a.verify)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--mode',choices=['sparse','dense'],required=True)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--nodes',type=int,default=100)
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
