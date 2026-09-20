"""Validate compact receipts and sum the covered K30 occupancies exactly.

This checks provenance and coverage; use the producer --verify commands for
the numerical replays. It never infers missing occupancies from samples.
"""
import argparse
from fractions import Fraction as F
import math
from pathlib import Path

import bridge as base
from certify_k30_q1 import ROWS,CUTOFF
from migrate_legacy_workspace import restore_or_verify


def run(output):
    assert not output.exists()
    restore_or_verify(base.ROOT)
    filenames=('k30_q1_v1.json','k30_sparse_v1.json','k30_kernel_q8_q32_v1.json')
    covered={};inputs={}
    for filename in filenames:
        path=base.HERE/'generated'/filename;data=base.read(path)
        replay_path=path.with_name(path.stem+'_replay.json');replay=base.read(replay_path)
        assert replay['producer_sha256']==base.sha(path)
        assert data['full_distance_proved'] is False and replay['full_distance_proved'] is False
        for name,digest in data['source_sha256'].items():
            assert base.sha(base.ROOT/name)==digest
        p=data['parameters']
        for key,value in dict(message_bits=1<<30,outer_rows=ROWS,cutoff=CUTOFF,step_bits=64,state_bits=20).items():
            assert p[key]==value
        if filename=='k30_q1_v1.json':
            assert replay['status']=='K30_Q1_512_BIT_POLYNOMIAL_REPLAY_PASSED'
            assert replay['coefficients_checked']==92
            bounds={1:base.decode(data['Q1_upper'])}
            co={int(w):base.decode(v) for w,v in data['coefficient_upper'].items()}
            assert base.bch_bound(co)==tuple(base.decode(data[k]) for k in ('Q1_upper','factor','rest'))
        else:
            bounds={row['occupation']:base.decode(row['upper']) for row in data['rows']}
            assert len(bounds)==len(data['rows'])
            expected=list(range(2,8)) if filename=='k30_sparse_v1.json' else list(range(8,33))
            assert sorted(bounds)==expected
            assert sum(bounds.values(),F(0))==base.decode(data['range_upper'])
            assert replay['status']==('K30_Q2_THROUGH_Q7_512_BIT_REPLAY_PASSED'
                if filename=='k30_sparse_v1.json' else 'K30_KERNEL_RANGE_512_BIT_REPLAY_PASSED')
        assert not set(covered).intersection(bounds)
        assert all(v>0 for v in bounds.values())
        covered.update(bounds)
        for source in (path,replay_path):inputs[source.relative_to(base.ROOT).as_posix()]=base.sha(source)
    assert sorted(covered)==list(range(1,33))
    total=sum(covered.values(),F(0));remaining=F(1,1<<40)-total
    assert remaining>0
    margin=lambda v:math.log2(v.denominator)-math.log2(v.numerator)
    payload=dict(status='K30_PARTIAL_Q1_THROUGH_Q32_RECEIPTS_VALIDATED',
        full_distance_proved=False,covered_occupancy_range=[1,32],missing_occupancy_range=[33,ROWS],
        partial_union_upper=base.encode(total),partial_margin_bits=margin(total),
        remaining_failure_budget=base.encode(remaining),remaining_budget_margin_bits=margin(remaining),
        remaining_budget_fraction_of_target=float(remaining*(1<<40)),
        input_sha256=inputs,audit_source_sha256=base.sha(Path(__file__)))
    base.write_new(output,payload)
    print('Certified range: Q1..Q32; NOT full distance',flush=True)
    print('Partial margin',payload['partial_margin_bits'],flush=True)
    print('Remaining budget fraction',payload['remaining_budget_fraction_of_target'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output)
