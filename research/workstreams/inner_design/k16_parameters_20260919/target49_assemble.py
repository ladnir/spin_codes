"""Authenticate and exactly sum all three two-round K16 components."""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path
import model
import assemble
import budget_dense
import verify_progress


def run(a):
    sources=verify_progress.outer_caps()
    components=[]
    for path,status in ((a.q1,'K16_TWO_ROUND_SUBSPACE_Q1_512_LINEAR_REPLAY_PASSED'),
                        (a.sparse,'K16_TWO_ROUND_SPARSE_512_REPLAY_PASSED'),
                        (a.dense,'IMT_DENSE_TOTAL_BUDGET_512_BIT_REPLAY_PASSED')):
        record,paths=assemble.receipt(path,status)
        components.append(record);sources.update(paths)
    q1,sparse,dense=components
    assert q1['status']=='K16_TWO_ROUND_SUBSPACE_Q1_OUTWARD'
    assert sparse['status']=='K16_TWO_ROUND_SPARSE_COMPONENT'
    instance=q1['instance']
    assert instance==sparse['instance']==dense['instance']
    assert (instance['message_bits'],instance['output_bits'],instance['cutoff'])==(65536,131072,13107)
    assert tuple(instance['inner'][k] for k in ('t','s','transvection_rounds'))==(64,12,2)
    assert instance['distance']=='1/10'
    weights={int(w):model.base.base.decode(v) for w,v in q1['coefficients'].items()}
    assert set(weights)==set(model.base.base.WEIGHTS)
    first,_,_=model.base.base.bch_bound(weights)
    assert first==model.base.base.decode(q1['upper'])
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
    assert total<F(2)**-49
    margin=lambda v:math.log2(v.denominator)-math.log2(v.numerator)
    for path in (Path(__file__),Path(assemble.__file__)):
        sources[path.resolve().relative_to(model.base.ROOT).as_posix()]=model.base.base.sha(path)
    result=dict(status='VERIFIED_K16_TWO_ROUND_FULL_DISTANCE_BOUND',instance=instance,
        target_bits=49,covered_occupancies=[1,512],union_upper=model.base.base.encode(total),
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
