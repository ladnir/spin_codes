"""Authenticate component replays and sum a full K16 distance bound exactly."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import model
import budget_dense
import verify_progress


def receipt(path,status):
    record=model.base.base.read(path)
    model.base.authenticate(record)
    replay_path=path.with_name(path.stem+'_replay.json')
    replay=model.base.base.read(replay_path)
    assert replay['status']==status
    assert replay['producer_sha256']==model.base.base.sha(path)
    return record,{p.resolve().relative_to(model.base.ROOT).as_posix():model.base.base.sha(p)
                   for p in (path,replay_path)}


def run(a):
    # Audit the deterministic BCH shell bounds, not an assumed spectrum.
    sources=verify_progress.outer_caps()
    q1,paths=receipt(a.q1,'K16_SUBSPACE_POINTS_512_REPLAY_PASSED');sources.update(paths)
    sparse,paths=receipt(a.sparse,'K16_SUBSPACE_SPARSE_512_REPLAY_PASSED');sources.update(paths)
    dense,paths=receipt(a.dense,'IMT_DENSE_TOTAL_BUDGET_512_BIT_REPLAY_PASSED');sources.update(paths)
    instance=q1['instance']
    assert instance==sparse['instance']==dense['instance']
    assert (instance['message_bits'],instance['output_bits'],instance['cutoff'])==(65536,131072,13107)
    assert (instance['inner']['t'],instance['inner']['s'])==(64,12)
    assert instance['distance']=='1/10'
    weights={int(w):model.base.exact(model.base.unpack(v)) for w,v in q1['q1']['coefficients'].items()}
    assert set(weights)==set(model.base.base.WEIGHTS)
    first,_,_=model.base.base.bch_bound(weights)
    assert first==model.base.base.decode(q1['q1']['upper'])
    best={}
    for witness in sparse['witnesses']:
        assert [r['occupation'] for r in witness['rows']]==list(range(2,64))
        for row in witness['rows']:
            q,power=row['occupation'],row['power']
            assert type(power) is int and -200<=power<=10**7
            best[q]=min(best.get(q,power),power)
    assert best=={int(q):v for q,v in sparse['best'].items()}
    assert sparse['covered_occupancies']==[2,63]
    second=sum((F(2)**p for p in best.values()),F(0))
    assert second==model.base.base.decode(sparse['upper'])
    third=budget_dense.checked_union(dense,64,512)
    total=first+second+third
    assert total < F(2)**-40
    margin=lambda v:math.log2(v.denominator)-math.log2(v.numerator)
    result=dict(status='VERIFIED_K16_SUBSPACE_FULL_DISTANCE_BOUND',instance=instance,
        target_bits=40,covered_occupancies=[1,512],union_upper=model.base.base.encode(total),
        margin_bits=margin(total),component_margin_bits=list(map(margin,(first,second,third))),
        component_upper=list(map(model.base.base.encode,(first,second,third))),
        full_distance_proved=True,asymptotic_claim=False,implementation_bound=False,
        source_sha256={**model.base.sources(),**sources})
    model.base.base.write_new(a.output,result)
    print(result['status'],result['margin_bits'],result['component_margin_bits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('q1','sparse','dense','output'): p.add_argument('--'+name,type=Path,required=True)
    run(p.parse_args())
