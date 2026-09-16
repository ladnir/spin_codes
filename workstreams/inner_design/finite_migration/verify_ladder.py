"""Authenticate all IMT ladder components and form an exact complete union."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

import ladder
import verify_full as frozen_union
import verify_progress
from flint import ctx

model = ladder.model


def receipt(path, status):
    record = model.base.read(path)
    model.authenticate(record)
    replay_path = path.with_name(path.stem + '_replay.json')
    replay = model.base.read(replay_path)
    assert replay['status'] == status
    assert replay['producer_sha256'] == model.base.sha(path)
    return record, {p.relative_to(model.ROOT).as_posix(): model.base.sha(p)
                    for p in (path, replay_path)}


def verify(exponent, q1_path, sparse_path, dense_path):
    ctx.prec = 256
    dependencies = verify_progress.outer_caps()
    q1, paths = receipt(q1_path, 'IMT_LADDER_Q1_512_BIT_LINEAR_REPLAY_PASSED')
    dependencies.update(paths)
    sparse, paths = receipt(sparse_path, 'IMT_LADDER_SPARSE_512_BIT_REPLAY_PASSED')
    dependencies.update(paths)
    dense, paths = receipt(dense_path, 'IMT_LADDER_DENSE_512_BIT_REPLAY_PASSED')
    dependencies.update(paths)
    assert q1['status'] == 'IMT_LADDER_Q1_ONLY'
    assert sparse['status'] == 'IMT_LADDER_SPARSE_RANGE'
    assert dense['status'] == 'IMT_LADDER_DENSE_COVER'
    engine = ladder.Engine(exponent)
    selected = [row for row in q1['results']
                if row['instance']['message_bits'] == 1 << exponent]
    assert len(selected) == 1
    first = selected[0]
    assert all(row['instance'] == engine.identity() for row in (first, sparse, dense))
    assert sparse['exponent'] == dense['exponent'] == exponent
    coefficients = {int(w): model.exact(model.unpack(v))
                    for w, v in first['q1']['coefficients'].items()}
    assert set(coefficients) == set(model.base.WEIGHTS)
    first_upper, _, _ = model.base.bch_bound(coefficients)
    assert first_upper == model.base.decode(first['q1']['upper'])
    last = sparse['last']
    assert 2 <= last < engine.length
    sparse_upper = frozen_union.sparse_union(sparse, last)
    dense_upper = frozen_union.dense_union(dense, last + 1, engine.length)
    total = first_upper + sparse_upper + dense_upper
    assert 0 < total < F(1, 1 << 40), 'Full union does not clear 40 bits'

    def margin(value):
        return math.log2(value.denominator) - math.log2(value.numerator)

    return dict(status='VERIFIED_IMT_LADDER_FULL_FINITE_DISTANCE_BOUND',
        instance=engine.identity(), target_bits=40,
        covered_occupancies=[1, engine.length],
        union_upper=model.base.encode(total), margin_bits=margin(total),
        component_upper=list(map(model.base.encode, (first_upper, sparse_upper, dense_upper))),
        component_margin_bits=list(map(margin, (first_upper, sparse_upper, dense_upper))),
        full_distance_proved=True, implementation_bound=False, asymptotic_claim=False,
        source_sha256={**ladder.candidate.sources(), **dependencies})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--m', type=int, choices=ladder.EXPONENTS, required=True)
    p.add_argument('--q1', type=Path, required=True)
    p.add_argument('--sparse', type=Path, required=True)
    p.add_argument('--dense', type=Path, required=True)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = verify(a.m, a.q1.resolve(), a.sparse.resolve(), a.dense.resolve())
    if a.output:
        model.base.write_new(a.output.resolve(), result)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('source_sha256', 'instance', 'component_upper', 'union_upper')}, indent=2))
