"""Refine any authenticated BCH-256 seed cover without changing old producers."""
import argparse
import heapq
import math
from pathlib import Path
import time

import numpy as np

import bch256_dense_types as search

base, core, direct = search.base, search.core, search.direct


def run(seed_path, directory, refinements, seconds, verify=False):
    seed_path, directory = seed_path.resolve(), directory.resolve()
    seed = base.read(seed_path)
    core.authenticate(seed, base.ROOT)
    spec = core.instance(seed['configuration'], seed['message_exponent'])
    core.require(spec['configuration'] == 't128_s19', 'Wrong selected map')
    raw = seed['result']
    minimum = raw['occupation_min']
    core.require(raw['occupation_max'] == spec['rows'] and seconds > 0 and refinements >= 0, 'Wrong geometry/budget')
    directory.mkdir(parents=True, exist_ok=True)
    paths = sorted(directory.glob('cover_*.json'))
    saved = base.read(paths[-1]) if paths else None
    if saved:
        core.authenticate(saved, base.ROOT)
        core.require(saved['seed_sha256'] == base.sha(seed_path) and saved['instance'] == spec, 'Wrong resume seed')
    core.require(not verify or saved is not None, 'No cover to replay')
    model = search.Refiner(spec, raw['bands'], directory)
    boxes = saved['leaves'] if saved else raw['selected_boxes']
    direct.check_cover(boxes, spec['rows'], minimum)
    leaves = {}
    for i, box in enumerate(boxes):
        value = model.direct(box['lower'], box['upper'], box['witness'])
        if saved:
            core.require(abs(value-box['own_log_bound']) < 2e-6, 'Witness replay failed')
        leaves[i] = dict(box, own_log_bound=value)
    def total():
        return float(np.logaddexp.reduce([b['own_log_bound'] for b in leaves.values()]))
    if verify:
        core.require(abs(total()-saved['log_upper']) < 2e-6, 'Union replay failed')
        print('Exact type coverage and direct binary64 witnesses replayed.', flush=True)
        return
    heap = [(-b['own_log_bound'], i) for i, b in leaves.items()]
    heapq.heapify(heap)
    serial = len(leaves)
    used = 0
    started = time.monotonic()
    print('initial complete margin', -total()/search.LN2, flush=True)
    while heap and used < refinements and time.monotonic()-started < seconds and total() > -60*search.LN2:
        _, i = heapq.heappop(heap)
        box = model.refine_box(leaves[i], -80*search.LN2)
        leaves[i] = box
        used += 1
        if box['own_log_bound'] > -80*search.LN2:
            corners = search.initial.engine.typed.vertices(box['lower'], box['upper'], spec['rows'])
            prepared = []
            for lo, hi in search.splitting.split_box(box['lower'], box['upper'], corners):
                lo, hi = list(map(int, lo)), list(map(int, hi))
                if direct.lattice_count(lo, hi, spec['rows']):
                    prepared.append(dict(lower=lo, upper=hi, witness=box['witness'],
                        own_log_bound=model.direct(lo, hi, box['witness'])))
            if prepared:
                core.require(sum(direct.lattice_count(c['lower'], c['upper'], spec['rows']) for c in prepared)
                    == direct.lattice_count(box['lower'], box['upper'], spec['rows']), 'Split volume changed')
                del leaves[i]
                for child in prepared:
                    leaves[serial] = child
                    heapq.heappush(heap, (-child['own_log_bound'], serial))
                    serial += 1
        if used % 25 == 0:
            print('refinements', used, 'leaves', len(leaves), 'margin', -total()/search.LN2, flush=True)
    coverage = direct.check_cover(list(leaves.values()), spec['rows'], minimum)
    result = dict(status='COMPLETE_DENSE_TYPE_BINARY64_DIAGNOSTIC', instance=spec, minimum=minimum,
        seed_sha256=base.sha(seed_path), bands=raw['bands'], leaves=list(leaves.values()),
        log_upper=total(), margin_bits=-total()/search.LN2, coverage=coverage,
        refinements_this_run=used, seconds_refining=time.monotonic()-started,
        source_sha256=search.sources(), outward_certificate=False, full_distance_proved=False)
    path = directory/f'cover_{len(paths):04d}.json'
    base.write_new(path, result)
    print('saved', path, 'margin', result['margin_bits'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed', type=Path, required=True)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--refinements', type=int, default=500)
    p.add_argument('--seconds', type=float, default=180)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.seed, a.directory, a.refinements, a.seconds, a.verify)
