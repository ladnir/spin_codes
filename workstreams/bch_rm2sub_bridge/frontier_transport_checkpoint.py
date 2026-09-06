"""Checkpointed sparse discovery; frozen frontier_sparse performs the replay.

Reuse auxiliary witnesses, never old probability bounds. Stop searching an
anchor once its recomputed bound meets the 70-bit per-occupancy policy.
"""
import argparse
from functools import lru_cache
from fractions import Fraction as F
import math
from pathlib import Path
from flint import arb,ctx
import frontier_sparse as core


def run(output,m,first,last,source):
    base=core.base;output=output.resolve();source=source.resolve();assert not output.exists()
    core.restore_or_verify(base.ROOT);ctx.prec=256
    rows,cutoff=core.parameters(m);assert 2<=first<=last<=min(rows,8192)
    t,s,spectrum,kernel=core.maps.load('t64_s20');caps=core.caps_module.caps()
    old=base.read(source);replay=source.with_name(source.stem+'_replay.json')
    assert old['status']=='FRONTIER_SPARSE_OUTWARD_CERTIFICATE'
    assert base.read(replay)['producer_sha256']==base.sha(source)
    for name,digest in old['source_sha256'].items():assert base.sha(base.ROOT/name)==digest
    hashes=core.sources()
    for p in (Path(__file__),source,replay):hashes[p.relative_to(base.ROOT).as_posix()]=base.sha(p)
    checkpoint=output.with_name(output.stem+'_checkpoints');checkpoint.mkdir(exist_ok=True)
    @lru_cache(maxsize=3)
    def region(tilt,hi):
        return core.poly.regions(t,s,spectrum,kernel,(-(arb(tilt)/10).exp()).exp(),hi,length=rows)
    best={};shards=[];attempts=[]
    def accept(shard):
        shards.append(shard)
        for row in shard['rows']:
            if row['upper_power']<=-70:
                q=row['occupation'];best[q]=min(best.get(q,0),row['upper_power'])
    # These are discovery checkpoints, not certificates or independent replays.
    for path in sorted(checkpoint.glob('attempt_*.json')):
        data=base.read(path)
        assert data['source_sha256']==hashes and data['message_exponent']==m and data['target_range']==[first,last]
        attempts.append(data['attempt'])
        if data['shard'] is not None:accept(data['shard'])
    def evaluate(witness,lo,hi):
        attempt=dict(occupancy_range=[lo,hi],witness=witness)
        if attempt in attempts:return
        values=core.evaluate(region(witness['tilt'],hi),witness,lo,hi,rows,cutoff,caps)
        shard=dict(occupancy_range=[lo,hi],witness=witness,rows=values)
        useful=any(r['upper_power']<=-70 and r['upper_power']<best.get(r['occupation'],0) for r in values)
        base.write_new(checkpoint/f'attempt_{len(attempts):04d}.json',dict(
            status='DISCOVERY_CHECKPOINT_NOT_A_CERTIFICATE',message_exponent=m,target_range=[first,last],
            attempt=attempt,shard=shard if useful else None,source_sha256=hashes))
        attempts.append(attempt)
        if useful:accept(shard)
        print('recomputed',lo,hi,'margin',round(values[0]['margin_bits_diagnostic'],3),'covered',len(best),flush=True)
    shift=round(-10*(m-old['message_exponent'])*math.log(2))
    for shard in old['shards']:
        a,b=shard['occupancy_range'];lo=max(first,a);hi=min(last,b)
        if lo>hi or all(q in best for q in range(lo,hi+1)):continue
        witness=dict(shard['witness']);witness['tilt']+=shift;witness.pop('margin_bits_diagnostic',None)
        evaluate(witness,lo,hi)
    while len(best)<last-first+1:
        anchor=next(q for q in range(first,last+1) if q not in best)
        hi=min(last,max(anchor,int(anchor*1.12)))
        prediction=round((-76+6*math.log2(anchor)-10*(m-20)*math.log(2))/5)*5
        for offset in (0,-1,1,-2,2,-3,3,-4,4,-5,5,-7,7,-10,10):
            tilt=prediction+offset
            witness=core.choose(region(tilt,hi),anchor,tilt,rows,cutoff,caps)
            if witness['margin_bits_diagnostic']<71:continue
            evaluate(witness,anchor,hi)
            if anchor in best:break
        assert anchor in best,f'Unclosed anchor Q{anchor}; discovery checkpoints retained'
    assert sorted(best)==list(range(first,last+1))
    base.write_new(output,dict(status='FRONTIER_SPARSE_OUTWARD_CERTIFICATE',message_exponent=m,
        parameters=dict(outer_rows=rows,cutoff=cutoff,step_bits=t,state_bits=s),occupancy_range=[first,last],
        upper_powers=[best[q] for q in sorted(best)],range_upper=base.encode(sum((core.as_fraction(v) for v in best.values()),F(0))),
        shards=shards,source_sha256=hashes,full_distance_proved=False))
    print('Complete checkpointed sparse certificate',m,first,last,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--m',type=int,default=28);p.add_argument('--first',type=int,default=8);p.add_argument('--last',type=int,default=8191)
    p.add_argument('--source',type=Path,required=True);a=p.parse_args();run(a.output,a.m,a.first,a.last,a.source)
