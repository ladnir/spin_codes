"""Authenticate numerical replays and reconstruct the exact K20 union bound."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path
from flint import ctx
import candidate
import certify_dense
import verify_progress
model = candidate.model


def receipt(name,status):
    path = candidate.HERE/name
    record = model.base.read(path)
    replay_path = path.with_name(path.stem+'_replay.json')
    replay = model.base.read(replay_path)
    model.authenticate(record)
    assert replay['status'] == status and replay['producer_sha256'] == model.base.sha(path)
    return record,{p.relative_to(model.ROOT).as_posix():model.base.sha(p) for p in (path,replay_path)}


def sparse_union(bank,last):
    assert bank['last'] == last and bank['complete_sparse_range']
    best = {}
    for job in bank['witnesses']:
        lo,hi = job['interval']
        assert 2 <= lo <= hi <= last
        assert [r['occupation'] for r in job['rows']] == list(range(lo,hi+1))
        for row in job['rows']:
            q,power = row['occupation'],row['power']
            assert type(power) is int
            if power <= -60:
                best[q] = min(best.get(q,power),power)
    assert best == {int(q):power for q,power in bank['best'].items()}
    assert set(best) == set(range(2,last+1))
    total = sum((F(2)**power for power in best.values()),F(0))
    assert total == model.base.decode(bank['upper'])
    return total


def dense_union(record,minimum,rows):
    assert record['minimum'] == minimum and record['unresolved'] == 0
    coverage = certify_dense.geometry.check_partition(record['leaves'],record['splits'],minimum,rows)
    assert coverage == record['coverage']
    powers = [node['power'] for node in record['leaves'].values()]
    assert all(type(p) is int and p <= -80 for p in powers)
    total = sum((F(2)**p for p in powers),F(0))
    assert total == model.base.decode(record['upper'])
    return total


def verify():
    ctx.prec = 256
    dependencies = verify_progress.outer_caps()
    q1,paths = receipt('Q1_v2.json','CANDIDATE_SCREEN_512_BIT_REPLAY_PASSED')
    dependencies.update(paths)
    sparse,paths = receipt('RANGE_M20.json','WEIGHT5_SPARSE_RANGE_REPLAY_PASSED')
    dependencies.update(paths)
    dense,paths = receipt('COVER_M20.json','CANDIDATE_DENSE_COVER_REPLAY_PASSED')
    dependencies.update(paths)
    assert all(r['candidate'] == 'weight5_seed0' for r in (q1,sparse,dense))
    engine = candidate.Engine(20)
    selected = [r for r in q1['results'] if r['instance']['message_bits'] == 1 << 20]
    assert len(selected) == 1
    first = selected[0]
    assert all(r['instance'] == engine.identity() for r in (first,sparse,dense))
    coefficients = {int(w):model.exact(model.unpack(v)) for w,v in first['q1']['coefficients'].items()}
    assert set(coefficients) == set(model.base.WEIGHTS)
    first_upper,_,_ = model.base.bch_bound(coefficients)
    assert first_upper == model.base.decode(first['q1']['upper'])
    sparse_upper = sparse_union(sparse,511)
    dense_upper = dense_union(dense,512,engine.length)
    total = first_upper+sparse_upper+dense_upper
    assert 0 < total < F(1,1 << 40)
    def margin(value):
        return math.log2(value.denominator)-math.log2(value.numerator)
    return dict(status='VERIFIED_FULL_BCH256_WEIGHT5_FINITE_DISTANCE_BOUND',
        instance=engine.identity(),target_bits=40,covered_occupancies=[1,engine.length],
        q1_upper=model.base.encode(first_upper),sparse_upper=model.base.encode(sparse_upper),
        dense_upper=model.base.encode(dense_upper),union_upper=model.base.encode(total),
        margin_bits=margin(total),component_margin_bits=[margin(v) for v in (first_upper,sparse_upper,dense_upper)],
        dense_leaves=len(dense['leaves']),full_distance_proved=True,
        implementation_bound=False,asymptotic_claim=False,
        scope='Fixed BCH [256,128,d>=38], K=2^20, t=128,s=19; weight5_seed0 feedback; setup probability only',
        source_sha256={**candidate.sources(),**dependencies})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path)
    a = p.parse_args()
    result = verify()
    if a.output:
        model.base.write_new(a.output.resolve(),result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('source_sha256','instance')},indent=2))
