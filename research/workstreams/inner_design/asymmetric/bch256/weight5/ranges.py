"""Candidate-specific contiguous sparse ranges, retaining every attempted witness."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
from flint import ctx
import candidate
import sparse_bch as discovery
import sparse_ranges as transfer
model = candidate.model


def run(output, exponent, last, name, verify):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        exponent,last,name = (saved[k] for k in ('exponent','last','candidate'))
    else:
        assert not output.exists()
    ctx.prec = 512 if verify else 256
    engine = candidate.Engine(exponent,name)
    assert 2 <= last <= engine.length
    best,witnesses = {},[]
    if saved:
        assert saved['instance'] == engine.identity()
        for job in saved['witnesses']:
            tilt = F(job['tilt'])
            lo,hi = job['interval']
            ps = list(map(model.base.decode,job['probabilities']))
            rows = transfer.evaluate(engine,engine.region(tilt,hi),ps,tilt,lo,hi)
            assert len(rows) == len(job['rows'])
            for actual,old in zip(rows,job['rows']):
                assert actual['occupation'] == old['occupation'] and actual['power'] <= old['power']
                if old['power'] <= -60:
                    q = old['occupation']
                    best[q] = min(best.get(q,old['power']),old['power'])
            engine.region.cache_clear()
            engine.epoch.cache_clear()
            print('replay sparse',exponent,lo,hi,flush=True)
    else:
        while len(best) < last-1:
            anchor = next(q for q in range(2,last+1) if q not in best)
            hi = min(last,anchor+63,max(anchor+7,int(anchor*1.3)))
            prediction = (-76+6*math.log2(anchor)-10*(exponent-20)*math.log(2))/10
            for delta in (0.,-.5,.5,-1.,1.):
                tilt = F.from_float(prediction+delta)
                region = engine.region(tilt,hi)
                ps = discovery.choose(engine,region[:anchor+1],anchor)
                rows = transfer.evaluate(engine,region,ps,tilt,anchor,hi)
                witnesses.append(dict(interval=[anchor,hi],tilt=str(tilt),
                    probabilities=list(map(model.base.encode,ps)),rows=rows))
                for row in rows:
                    if row['power'] <= -60:
                        q = row['occupation']
                        best[q] = min(best.get(q,row['power']),row['power'])
                print('sparse range',name,exponent,anchor,hi,'covered',len(best),'anchor power',rows[0]['power'],flush=True)
                if anchor in best:
                    break
            engine.region.cache_clear()
            engine.epoch.cache_clear()
            if anchor not in best:
                break
    total = sum((F(2)**v for v in best.values()),F(0))
    complete = len(best) == last-1
    if saved:
        assert best == {int(q):power for q,power in saved['best'].items()}
        assert total == model.base.decode(saved['upper'])
        model.base.write_new(output.with_name(output.stem+'_replay.json'),dict(
            status='WEIGHT5_SPARSE_RANGE_REPLAY_PASSED',producer_sha256=model.base.sha(output),full_distance_proved=False))
    else:
        model.base.write_new(output,dict(status='OUTWARD_CANDIDATE_SPARSE_RANGE',
            candidate=name,exponent=exponent,last=last,instance=engine.identity(),
            witnesses=witnesses,best=best,upper=model.base.encode(total),complete_sparse_range=complete,
            source_sha256=candidate.sources(),full_distance_proved=False))
    print('sparse done',exponent,'complete',complete,'covered',len(best),flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,choices=(16,18,20),default=20)
    p.add_argument('--last',type=int,default=511)
    p.add_argument('--candidate',choices=('weight5_seed0','weight5_seed1','balanced'),default='weight5_seed0')
    p.add_argument('--verify',action='store_true')
    a = p.parse_args()
    run(a.output.resolve(),a.m,a.last,a.candidate,a.verify)
