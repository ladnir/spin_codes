"""Refine a complete IMT cover against an exact total budget, not a leaf margin."""
import argparse
from fractions import Fraction as F
import heapq
import math
from pathlib import Path
import time

from flint import ctx
import dense_ladder as geometry_driver
import short_dense

model = geometry_driver.model
geometry = geometry_driver.geometry
STATUS = 'IMT_DENSE_TOTAL_BUDGET_COVER'


def total_upper(leaves):
    assert leaves
    powers = [node['power'] for node in leaves.values()]
    assert all(type(p) is int and -200 <= p <= 10**7 for p in powers)
    return sum((F(2)**p for p in powers), F(0))


def checked_union(record, minimum, rows):
    assert record['status'] == STATUS and record['minimum'] == minimum
    assert record['instance']['outer_rows'] == rows
    coverage = geometry.check_partition(record['leaves'], record['splits'], minimum, rows)
    assert coverage == record['coverage']
    total = total_upper(record['leaves'])
    assert total == model.base.decode(record['upper'])
    assert type(record['budget_bits']) is int and 40 <= record['budget_bits'] <= 100
    assert record['budget_met'] is (total < F(2)**-record['budget_bits'])
    return total


def run(output, exponent, minimum, seed_path, nodes, seconds, budget_bits, verify=False):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == STATUS
        exponent, minimum, budget_bits = saved['exponent'], saved['minimum'], saved['budget_bits']
    elif output.exists():
        raise FileExistsError(output)
    assert exponent in (16, 18) and 40 <= budget_bits <= 100
    assert nodes >= 0 and seconds > 0
    ctx.prec = 512 if verify else 256
    checker = short_dense.Checker(exponent)
    target = F(2)**-budget_bits

    def bound(node):
        return geometry_driver.certify_dense.power_bound(checker.bound(
            *geometry.geometry(node), node['witness']))

    if saved:
        assert saved['instance'] == checker.engine.identity()
        leaves, splits = saved['leaves'], saved['splits']
        checked_union(saved, minimum, checker.rows)
    elif seed_path:
        seed = model.base.read(seed_path)
        model.authenticate(seed)
        leaves, splits = geometry_driver.reshape(seed, minimum, checker.rows)
    else:
        leaves = {'': geometry.box(minimum, checker.rows, F(0), F(1),
            checker.witness(minimum, checker.rows, F(0), F(1)))}
        splits = {}
    geometry.check_partition(leaves, splits, minimum, checker.rows)
    for i, (key, node) in enumerate(leaves.items()):
        power = bound(node)
        if saved:
            assert power <= node['power'], (key, power, node['power'])
        else:
            node['power'] = power
        if i % 40 == 0:
            print('budget evaluation', ctx.prec, exponent, i + 1, '/', len(leaves), flush=True)
    total = total_upper(leaves)
    used = 0
    if not saved:
        heap = [(-node['power'], key) for key, node in leaves.items()]
        heapq.heapify(heap)
        started = time.monotonic()
        while total >= target and used < nodes and time.monotonic() - started < seconds:
            _, key = heapq.heappop(heap)
            node = leaves[key]
            before = F(2)**node['power']
            proposal = dict(node, witness=checker.witness(*geometry.geometry(node)))
            proposal['power'] = bound(proposal)
            if proposal['power'] < node['power']:
                node = leaves[key] = proposal
            total += F(2)**node['power'] - before
            used += 1
            if total < target:
                break
            axis = geometry_driver.splitting.split_axis(node, checker)
            children = geometry.children(node, axis)
            total -= F(2)**node['power']
            del leaves[key]
            splits[key] = axis
            for i, child in enumerate(children):
                child['power'] = bound(child)
                child_key = key + str(i)
                leaves[child_key] = child
                total += F(2)**child['power']
                heapq.heappush(heap, (-child['power'], child_key))
            if used % 10 == 0:
                print('budget refinement', exponent, used, 'leaves', len(leaves),
                      'largest power', -heap[0][0], flush=True)
    assert total == total_upper(leaves)
    coverage = geometry.check_partition(leaves, splits, minimum, checker.rows)
    met = total < target
    if saved:
        assert total == model.base.decode(saved['upper']) and met == saved['budget_met']
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_DENSE_TOTAL_BUDGET_512_BIT_REPLAY_PASSED',
            producer_sha256=model.base.sha(output), leaves_checked=len(leaves),
            budget_met=met, full_distance_proved=False))
    else:
        sources = short_dense.mixed_dense.ladder.candidate.sources()
        if seed_path:
            sources[seed_path.relative_to(model.ROOT).as_posix()] = model.base.sha(seed_path)
        result = dict(status=STATUS, exponent=exponent, minimum=minimum,
            instance=checker.engine.identity(), leaves=leaves, splits=splits,
            coverage=coverage, upper=model.base.encode(total), budget_bits=budget_bits,
            budget_met=met, refinements=used, full_distance_proved=False, source_sha256=sources)
        assert checked_union(result, minimum, checker.rows) == total
        model.base.write_new(output, result)
    margin = math.log2(total.denominator) - math.log2(total.numerator)
    print('budget done', exponent, 'met', met, 'margin', margin, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16, 18), default=16)
    p.add_argument('--minimum', type=int, default=64)
    p.add_argument('--seed', type=Path)
    p.add_argument('--nodes', type=int, default=160)
    p.add_argument('--seconds', type=float, default=300)
    p.add_argument('--budget-bits', type=int, default=42)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.output.resolve(), a.m, a.minimum, a.seed.resolve() if a.seed else None,
        a.nodes, a.seconds, a.budget_bits, a.verify)
