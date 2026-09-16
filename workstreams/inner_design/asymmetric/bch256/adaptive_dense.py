"""Bounded dense-cover search; every checkpoint retains unresolved rectangles."""
import argparse
from fractions import Fraction as F
import heapq
from pathlib import Path
import time
from flint import ctx
import certify_dense as fixed_cover
import dense_bch as initial
import search_ladder_dense_fixed as splitting

model = fixed_cover.model
geometry = fixed_cover.geometry


def run(output, exponent, minimum, nodes, seconds, seed=None, verify=False):
    prior = model.base.read(output if verify else seed) if verify or seed else None
    if not verify:
        assert not output.exists()
    if prior:
        model.authenticate(prior)
        assert prior['message_exponent'] == exponent and prior['minimum'] == minimum
    ctx.prec = 512 if verify else 256
    checker = fixed_cover.search.Checker(exponent)
    def bound(node):
        return fixed_cover.power_bound(checker.bound(*geometry.geometry(node), node['witness']))
    def tune(node):
        candidates = [node]
        lo, hi, a, b = geometry.geometry(node)
        for witness in (checker.witness(lo, hi, a, b), initial.Checker.witness(checker, lo, hi, a, b)):
            candidate = dict(node, witness=witness)
            candidate['power'] = bound(candidate)
            candidates.append(candidate)
        return min(candidates, key=lambda r:r['power'])
    if prior:
        leaves, splits = prior['leaves'], prior['splits']
        assert prior['instance'] == checker.engine.identity()
    else:
        root = geometry.box(minimum, checker.rows, F(0), F(1), checker.witness(minimum, checker.rows, F(0), F(1)))
        root['power'] = bound(root)
        leaves, splits = {'':root}, {}
    geometry.check_partition(leaves, splits, minimum, checker.rows)
    if verify:
        for i, (key, node) in enumerate(leaves.items()):
            assert bound(node) <= node['power'], key
            if i%100 == 0:
                print('dense adaptive replay', i+1, '/', len(leaves), flush=True)
    else:
        heap = [(-node['power'], key) for key, node in leaves.items() if node['power'] > -80]
        heapq.heapify(heap)
        start = time.monotonic()
        used = 0
        while heap and used < nodes and time.monotonic()-start < seconds:
            _, key = heapq.heappop(heap)
            node = leaves[key] = tune(leaves[key])
            used += 1
            if node['power'] > -80:
                lo, hi, a, b = geometry.geometry(node)
                if lo == hi and b-a <= F(1, 1 << 18):
                    continue
                axis = splitting.split_axis(node, checker)
                children = geometry.children(node, axis)
                del leaves[key]
                splits[key] = axis
                for i, child in enumerate(children):
                    child['power'] = bound(child)
                    leaves[key+str(i)] = child
                    if child['power'] > -80:
                        heapq.heappush(heap, (-child['power'], key+str(i)))
            if used%10 == 0:
                print('adaptive', exponent, used, 'leaves', len(leaves), 'weak', len(heap), flush=True)
    coverage = geometry.check_partition(leaves, splits, minimum, checker.rows)
    unresolved = sum(n['power'] > -80 for n in leaves.values())
    total = sum((F(2)**n['power'] for n in leaves.values()), F(0)) if not unresolved else None
    if verify:
        assert prior['unresolved'] == unresolved
        if total is not None:
            assert total == model.base.decode(prior['upper'])
        model.base.write_new(output.with_name(output.stem+'_replay.json'),
                             dict(status='ASYMMETRIC_ADAPTIVE_DENSE_512_BIT_REPLAY_PASSED',
                                  producer_sha256=model.base.sha(output), complete_dense_certificate=not unresolved,
                                  full_distance_proved=False))
    else:
        result = dict(status='ASYMMETRIC_ADAPTIVE_DENSE_CERTIFICATE' if not unresolved else 'ASYMMETRIC_ADAPTIVE_DENSE_UNRESOLVED',
                      message_exponent=exponent, minimum=minimum, instance=checker.engine.identity(),
                      leaves=leaves, splits=splits, coverage=coverage, unresolved=unresolved,
                      refinements_this_run=used, source_sha256=model.sources(), full_distance_proved=False)
        if total is not None:
            result.update(upper=model.base.encode(total), margin_bits=float(-model.number(total).log()/model.number(2).log()))
        model.base.write_new(output, result)
    print('adaptive complete', exponent, len(leaves), 'unresolved', unresolved, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--m', type=int, choices=(16,18,20), default=20)
    p.add_argument('--minimum', type=int, default=512)
    p.add_argument('--nodes', type=int, default=200)
    p.add_argument('--seconds', type=float, default=180)
    p.add_argument('--seed', type=Path)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    assert 1 <= a.minimum <= 1 << (a.m-7) and a.nodes >= 0 and a.seconds > 0
    run(a.output.resolve(), a.m, a.minimum, a.nodes, a.seconds, a.seed, a.verify)
