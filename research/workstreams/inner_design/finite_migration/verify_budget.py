"""Check complete coverage and the actual IMT union, independent of leaf targets."""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

from flint import ctx
import budget_dense
import verify_ladder as prior

model = prior.model


def full_union(parts):
    assert len(parts) == 3 and all(isinstance(x, F) and x > 0 for x in parts)
    total = sum(parts, F(0))
    assert total < F(2)**-40, 'Full union does not clear 40 bits'
    return total


def verify(exponent, q1_path, sparse_path, dense_path):
    ctx.prec = 256
    dependencies = prior.verify_progress.outer_caps()
    q1, paths = prior.receipt(q1_path, 'IMT_LADDER_Q1_512_BIT_LINEAR_REPLAY_PASSED')
    dependencies.update(paths)
    sparse, paths = prior.receipt(sparse_path, 'IMT_LADDER_SPARSE_512_BIT_REPLAY_PASSED')
    dependencies.update(paths)
    dense, paths = prior.receipt(dense_path, 'IMT_DENSE_TOTAL_BUDGET_512_BIT_REPLAY_PASSED')
    dependencies.update(paths)
    assert q1['status'] == 'IMT_LADDER_Q1_ONLY'
    assert sparse['status'] == 'IMT_LADDER_SPARSE_RANGE'
    engine = prior.ladder.Engine(exponent)
    selected = [row for row in q1['results'] if row['instance']['message_bits'] == 1 << exponent]
    assert len(selected) == 1
    first = selected[0]
    assert all(row['instance'] == engine.identity() for row in (first, sparse, dense))
    assert sparse['exponent'] == dense['exponent'] == exponent
    coefficients = {int(w): model.exact(model.unpack(v)) for w, v in first['q1']['coefficients'].items()}
    assert set(coefficients) == set(model.base.WEIGHTS)
    first_upper, _, _ = model.base.bch_bound(coefficients)
    assert first_upper == model.base.decode(first['q1']['upper'])
    last = sparse['last']
    assert 2 <= last < engine.length
    sparse_upper = prior.frozen_union.sparse_union(sparse, last)
    dense_upper = budget_dense.checked_union(dense, last + 1, engine.length)
    replay = model.base.read(dense_path.with_name(dense_path.stem + '_replay.json'))
    assert dense['budget_met'] and replay['budget_met']
    assert replay['leaves_checked'] == len(dense['leaves'])
    parts = (first_upper, sparse_upper, dense_upper)
    total = full_union(parts)
    margin = lambda value: math.log2(value.denominator) - math.log2(value.numerator)
    return dict(status='VERIFIED_IMT_TOTAL_BUDGET_FULL_FINITE_DISTANCE_BOUND',
        instance=engine.identity(), target_bits=40, covered_occupancies=[1, engine.length],
        union_upper=model.base.encode(total), margin_bits=margin(total),
        component_upper=list(map(model.base.encode, parts)), component_margin_bits=list(map(margin, parts)),
        full_distance_proved=True, implementation_bound=False, asymptotic_claim=False,
        source_sha256={**prior.ladder.candidate.sources(), **dependencies})


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--m', type=int, choices=(16, 18), required=True)
    p.add_argument('--q1', type=Path, required=True)
    p.add_argument('--sparse', type=Path, required=True)
    p.add_argument('--dense', type=Path, required=True)
    p.add_argument('--output', type=Path)
    a = p.parse_args()
    result = verify(a.m, a.q1.resolve(), a.sparse.resolve(), a.dense.resolve())
    if a.output:
        model.base.write_new(a.output.resolve(), result)
    print(json.dumps({k: v for k, v in result.items() if k not in ('source_sha256', 'instance', 'component_upper', 'union_upper')}, indent=2))
