"""Recompute old search witnesses for new maps; never transfer old bounds."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
import time
from flint import arb,ctx
import model
import sparse_ranges
import dense_ladder
import certify_dense

PRIOR = model.HERE.parent/'finite_migration'


def run(a):
    out = a.output.resolve()
    saved = model.base.base.read(out) if a.verify else None
    if saved:
        model.base.authenticate(saved)
        assert saved['parameters'] == [a.t,a.s,16] and saved['mode'] == a.mode
    elif out.exists():
        raise FileExistsError(out)
    ctx.prec = 512 if saved else 256
    engine = model.Engine(a.t,a.s)
    result = dict(status='K16_OCCUPANCY_COMPONENT',parameters=[a.t,a.s,16],mode=a.mode,
        instance=engine.identity(),full_distance_proved=False,precision_bits=ctx.prec)
    started = time.monotonic()
    if a.mode == 'sparse':
        seed_path = PRIOR/'SPARSE_M16_v2.json'
        seed = model.base.base.read(seed_path)
        jobs = saved['witnesses'] if saved else seed['witnesses']
        best, witnesses = {}, []
        for job in jobs:
            lo,hi = job['interval']
            tilt = F(job['tilt'])
            ps = list(map(model.base.base.decode,job['probabilities']))
            rows = sparse_ranges.evaluate(engine,engine.region(tilt,hi),ps,tilt,lo,hi)
            if saved:
                assert [r['occupation'] for r in rows] == [r['occupation'] for r in job['rows']]
                assert all(r['power'] <= old['power'] for r,old in zip(rows,job['rows']))
                rows = job['rows']
            witnesses.append(dict(interval=[lo,hi],tilt=str(tilt),probabilities=job['probabilities'],rows=rows))
            for row in rows:
                q,power = row['occupation'],row['power']
                best[q] = min(best.get(q,power),power)
            engine.region.cache_clear();engine.epoch.cache_clear()
            print('sparse',a.t,a.s,lo,hi,'worst power',max(best.values()),flush=True)
        assert set(best) == set(range(2,64))
        result.update(witnesses=witnesses,best=best,covered_occupancies=[2,63])
        total = sum((F(2)**p for p in best.values()),F(0))
    else:
        seed_path = PRIOR/'DENSE_M16_fast_combined_v1.json'
        seed = model.base.base.read(seed_path)
        leaves,splits = (saved['leaves'],saved['splits']) if saved else (seed['leaves'],seed['splits'])
        geometry = dense_ladder.geometry
        geometry.check_partition(leaves,splits,64,512)
        checker = model.Checker(a.t,a.s)
        result.update(leaves={},splits=splits,requested_occupancies=[64,512])
        for i,(key,node) in enumerate(leaves.items()):
            value = checker.bound(*geometry.geometry(node),node['witness'])
            power = certify_dense.power_bound(value)
            if saved:
                assert power <= node['power'],(key,power,node['power'])
                power = node['power']
            result['leaves'][key] = dict(node,power=power)
            if i%50 == 0:
                print('dense',a.t,a.s,i+1,'/',len(leaves),'power',power,flush=True)
            if not saved and time.monotonic()-started > a.seconds:
                break
        complete = len(result['leaves']) == len(leaves)
        result['complete_partition_evaluated'] = complete
        if complete:
            result['covered_occupancies'] = [64,512]
        total = sum((F(2)**n['power'] for n in result['leaves'].values()),F(0)) if complete else None
    result['source_sha256'] = model.base.sources()
    result['source_sha256'][seed_path.relative_to(model.base.ROOT).as_posix()] = model.base.base.sha(seed_path)
    if total is not None:
        result['upper'] = model.base.base.encode(total)
        result['margin_bits'] = math.log2(total.denominator)-math.log2(total.numerator)
        print('component margin',result['margin_bits'],flush=True)
    if saved:
        assert result.get('upper') == saved.get('upper')
        result = dict(status='K16_COMPONENT_512_BIT_REPLAY_PASSED',
            producer_sha256=model.base.base.sha(out),full_distance_proved=False)
        out = out.with_name(out.stem+'_replay.json')
    out.parent.mkdir(parents=True,exist_ok=True)
    model.base.base.write_new(out,result)


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--t',type=int,required=True);p.add_argument('--s',type=int,required=True)
    p.add_argument('--mode',choices=['sparse','dense'],required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seconds',type=float,default=180)
    p.add_argument('--verify',action='store_true')
    run(p.parse_args())
