"""Full-range component bounds for the selected expansion-subspace maps."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import math
from pathlib import Path
from types import SimpleNamespace
from flint import arb,ctx
import model
import subspace
import sparse_ranges
import budget_dense
import fast_combined_dense
import high_band_dense
import short_dense


class Checker(model.Checker):
    def __init__(self,record):
        super().__init__(64,record['s'])
        self.engine=subspace.Engine(record)
        self.constant_comparison=lru_cache(maxsize=2048)(self._constant_comparison)

    _comparison=fast_combined_dense.Checker._comparison
    _constant_comparison=fast_combined_dense.Checker._constant_comparison
    split_bound=high_band_dense.Checker.split_bound

    def bound(self,*args):
        old=short_dense.Checker.bound(self,*args)
        if old < -90*arb(2).log(): return old
        improved=self.split_bound(*args)
        return old if improved is None else min(old,improved)


def sparse(a,record):
    out=a.output.resolve()
    saved=model.base.base.read(out) if a.verify else None
    if saved:
        model.base.authenticate(saved)
        assert saved['status']=='K16_SUBSPACE_SPARSE_COMPONENT'
    elif out.exists(): raise FileExistsError(out)
    ctx.prec=512 if saved else 256
    e=subspace.Engine(record)
    if saved: assert saved['instance']==e.identity()
    source=model.base.base.read(a.seed)
    jobs=saved['witnesses'] if saved else [r for r in source['sparse'] if r['q']<=63]
    best,witnesses={},[]
    for job in jobs:
        tilt=F(job['tilt'])
        ps=list(map(model.base.base.decode,job['probabilities']))
        rows=sparse_ranges.evaluate(e,e.region(tilt,63),ps,tilt,2,63)
        if saved:
            assert [r['occupation'] for r in rows]==list(range(2,64))
            assert all(r['power']<=old['power'] for r,old in zip(rows,job['rows']))
            rows=job['rows']
        witnesses.append(dict(tilt=str(tilt),probabilities=job['probabilities'],rows=rows))
        for row in rows:
            q,power=row['occupation'],row['power']
            best[q]=min(best.get(q,power),power)
        e.region.cache_clear();e.epoch.cache_clear()
        print('sparse cover',len(witnesses),'covered below -50',sum(p<=-50 for p in best.values()),'/62',flush=True)
        if not saved and len(witnesses)==len(jobs) and max(best.values())>-55 and len(jobs)<30:
            anchor=next(q for q in range(2,64) if best[q]>-55)
            proposal=model.study.core.grid.ladder.screen.sparse(e,anchor)
            print('retuned sparse anchor',anchor,proposal['margin_bits'],flush=True)
            if proposal['margin_bits']>55:
                jobs.append(proposal)
    assert set(best)==set(range(2,64))
    total=sum((F(2)**p for p in best.values()),F(0))
    result=dict(status='K16_SUBSPACE_SPARSE_COMPONENT',instance=e.identity(),
        covered_occupancies=[2,63],witnesses=witnesses,best=best,
        upper=model.base.base.encode(total),margin_bits=math.log2(total.denominator)-math.log2(total.numerator),
        full_distance_proved=False,source_sha256=model.base.sources())
    print('sparse margin',result['margin_bits'],flush=True)
    if saved:
        assert result['upper']==saved['upper']
        result=dict(status='K16_SUBSPACE_SPARSE_512_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(out),full_distance_proved=False)
        out=out.with_name(out.stem+'_replay.json')
    model.base.base.write_new(out,result)


def run(a):
    records=model.base.base.read(a.maps)['rows']
    record=next(r['inner'] for r in records if r['s']==a.s and r['found'])
    if a.mode=='sparse': return sparse(a,record)
    original=budget_dense.short_dense
    budget_dense.short_dense=SimpleNamespace(Checker=lambda exponent:Checker(record),
        mixed_dense=original.mixed_dense)
    try:
        budget_dense.run(a.output.resolve(),16,64,a.seed.resolve() if a.seed else None,
            a.nodes,a.seconds,48,a.verify)
    finally:
        budget_dense.short_dense=original


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--maps',type=Path,required=True)
    p.add_argument('--s',type=int,default=12)
    p.add_argument('--mode',choices=['sparse','dense'],required=True)
    p.add_argument('--seed',type=Path)
    p.add_argument('--nodes',type=int,default=100)
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
