"""Reuse sparse witnesses across sizes, recomputing every bound at the new size.

Transport changes discovery only. A previous certificate's probability bound
is never substituted for a new-size bound. Frozen frontier_sparse replays
the resulting certificate without needing this discovery driver.
"""
import argparse
from functools import lru_cache
from fractions import Fraction as F
import math
from pathlib import Path
from flint import arb,ctx
import frontier_sparse as core


def run(output,m,inputs):
    base=core.base;output=output.resolve();assert not output.exists()
    core.restore_or_verify(base.ROOT);ctx.prec=256
    rows,cutoff=core.parameters(m)
    t,s,spectrum,kernel=core.maps.load('t64_s20');caps=core.caps_module.caps()
    @lru_cache(maxsize=6)
    def region(tilt,hi):
        return core.poly.regions(t,s,spectrum,kernel,(-(arb(tilt)/10).exp()).exp(),hi,length=rows)
    best={};shards=[];source_paths=[];intervals=[]
    def evaluate(witness,lo,hi):
        values=core.evaluate(region(witness['tilt'],hi),witness,lo,hi,rows,cutoff,caps)
        shards.append(dict(occupancy_range=[lo,hi],witness=witness,rows=values))
        for row in values:
            if row['upper_power']<=-70:
                q=row['occupation'];best[q]=min(best.get(q,0),row['upper_power'])
        print('new-size certificate',lo,hi,'anchor margin',round(values[0]['margin_bits_diagnostic'],3),flush=True)
    for source in inputs:
        source=source.resolve();old=base.read(source)
        assert old['status']=='FRONTIER_SPARSE_OUTWARD_CERTIFICATE'
        replay=source.with_name(source.stem+'_replay.json')
        assert base.read(replay)['producer_sha256']==base.sha(source)
        for name,digest in old['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
        source_paths.extend((source,replay));intervals.append(old['occupancy_range'])
        shift=round(-10*(m-old['message_exponent'])*math.log(2))
        for shard in old['shards']:
            lo,hi=shard['occupancy_range']
            witness=dict(shard['witness']);witness['tilt']+=shift
            # The old diagnostic margin is not a new-size result.
            witness.pop('margin_bits_diagnostic',None)
            evaluate(witness,lo,hi)
    intervals.sort();first,last=intervals[0][0],intervals[-1][1]
    assert 2<=first<=last<=min(rows,8192)
    assert all(a[1]+1==b[0] for a,b in zip(intervals,intervals[1:]))
    while len(best)<last-first+1:
        anchor=next(q for q in range(first,last+1) if q not in best)
        hi=min(last,max(anchor,int(anchor*1.12)))
        prediction=round((-76+6*math.log2(anchor)-10*(m-20)*math.log(2))/5)*5
        trials=[core.choose(region(tilt,hi),anchor,tilt,rows,cutoff,caps) for tilt in range(prediction-5,prediction+6)]
        witness=max(trials,key=lambda w:w['margin_bits_diagnostic'])
        evaluate(witness,anchor,hi)
        assert anchor in best,f'Cannot close transported witness at Q{anchor}'
    assert sorted(best)==list(range(first,last+1))
    hashes=core.sources()
    for p in [Path(__file__),*source_paths]:hashes[p.relative_to(base.ROOT).as_posix()]=base.sha(p)
    base.write_new(output,dict(status='FRONTIER_SPARSE_OUTWARD_CERTIFICATE',message_exponent=m,
        parameters=dict(outer_rows=rows,cutoff=cutoff,step_bits=t,state_bits=s),occupancy_range=[first,last],
        upper_powers=[best[q] for q in sorted(best)],range_upper=base.encode(sum((core.as_fraction(v) for v in best.values()),F(0))),
        shards=shards,source_sha256=hashes,full_distance_proved=False))
    print('Complete new-size sparse certificate',m,first,last,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,required=True);p.add_argument('--inputs',type=Path,nargs='+',required=True)
    a=p.parse_args();run(a.output,a.m,a.inputs)
