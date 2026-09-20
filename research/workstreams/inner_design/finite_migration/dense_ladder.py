"""Fresh IMT dense bounds using a prior split tree only as search geometry."""
import argparse
from fractions import Fraction as F
import heapq
from pathlib import Path
import time

import ladder
import certify_dense
import search_ladder_dense_fixed as splitting
from flint import ctx

model = ladder.model
geometry = certify_dense.geometry


def reshape(seed, minimum, rows):
    """Apply the same binary subdivision choices to a new root rectangle."""
    geometry.check_partition(seed['leaves'], seed['splits'], seed['minimum'],
                             seed['instance']['outer_rows'])
    leaves = {}
    pending = [('', geometry.box(minimum, rows, F(0), F(1), None))]
    while pending:
        key, node = pending.pop()
        if key in seed['leaves']:
            node['witness'] = dict(seed['leaves'][key]['witness'])
            leaves[key] = node
        else:
            pending.extend((key + str(i), child) for i, child in
                           enumerate(geometry.children(node, seed['splits'][key])))
    splits = dict(seed['splits'])
    geometry.check_partition(leaves, splits, minimum, rows)
    return leaves, splits


def run(output, exponent, minimum, seed_path, nodes, seconds, verify):
    saved = model.base.read(output) if verify else None
    if saved:
        model.authenticate(saved)
        assert saved['status'] == 'IMT_LADDER_DENSE_COVER'
        exponent, minimum = saved['exponent'], saved['minimum']
    elif output.exists():
        raise FileExistsError('Use a fresh output path')
    ctx.prec = 512 if verify else 256
    checker = ladder.Checker(exponent)

    def bound(node):
        return certify_dense.power_bound(checker.bound(
            *geometry.geometry(node), node['witness']))

    if saved:
        assert saved['instance'] == checker.engine.identity()
        leaves, splits = saved['leaves'], saved['splits']
    elif seed_path:
        prior = model.base.read(seed_path)
        model.authenticate(prior)
        leaves, splits = reshape(prior, minimum, checker.rows)
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
        if i % 50 == 0:
            print('dense evaluation', ctx.prec, exponent, i + 1, '/', len(leaves), flush=True)
    used = 0
    if not saved:
        heap = [(-node['power'], key) for key, node in leaves.items() if node['power'] > -80]
        heapq.heapify(heap)
        started = time.monotonic()
        while heap and used < nodes and time.monotonic() - started < seconds:
            _, key = heapq.heappop(heap)
            node = leaves[key]
            proposal = dict(node, witness=checker.witness(*geometry.geometry(node)))
            proposal['power'] = bound(proposal)
            if proposal['power'] < node['power']:
                node = leaves[key] = proposal
            used += 1
            if node['power'] > -80:
                axis = splitting.split_axis(node, checker)
                splits[key] = axis
                del leaves[key]
                for i, child in enumerate(geometry.children(node, axis)):
                    child['power'] = bound(child)
                    child_key = key + str(i)
                    leaves[child_key] = child
                    if child['power'] > -80:
                        heapq.heappush(heap, (-child['power'], child_key))
            if used % 10 == 0:
                print('dense refinement', exponent, used, 'weak', len(heap), flush=True)
    coverage = geometry.check_partition(leaves, splits, minimum, checker.rows)
    unresolved = sum(node['power'] > -80 for node in leaves.values())
    total = sum((F(2)**node['power'] for node in leaves.values()), F(0)) if not unresolved else None
    if saved:
        assert coverage == saved['coverage'] and unresolved == saved['unresolved']
        if total is not None:
            assert total == model.base.decode(saved['upper'])
        model.base.write_new(output.with_name(output.stem + '_replay.json'), dict(
            status='IMT_LADDER_DENSE_512_BIT_REPLAY_PASSED',
            producer_sha256=model.base.sha(output),
            complete_dense_certificate=not unresolved, full_distance_proved=False))
    else:
        sources = ladder.candidate.sources()
        if seed_path:
            sources[seed_path.relative_to(model.ROOT).as_posix()] = model.base.sha(seed_path)
        result = dict(status='IMT_LADDER_DENSE_COVER', exponent=exponent,
            minimum=minimum, instance=checker.engine.identity(), leaves=leaves,
            splits=splits, coverage=coverage, unresolved=unresolved,
            refinements=used, full_distance_proved=False, source_sha256=sources)
        if total is not None:
            result['upper'] = model.base.encode(total)
        model.base.write_new(output, result)
    print('dense done', exponent, 'unresolved', unresolved, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=ladder.EXPONENTS, default=22)
    p.add_argument('--minimum', type=int, default=512)
    p.add_argument('--seed', type=Path)
    p.add_argument('--nodes', type=int, default=60)
    p.add_argument('--seconds', type=float, default=120)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    if not 2 <= a.minimum <= 1 << (a.m - 7) or a.nodes < 0 or a.seconds <= 0:
        p.error('Invalid range or search budget')
    run(a.output.resolve(), a.m, a.minimum, a.seed.resolve() if a.seed else None,
        a.nodes, a.seconds, a.verify)
