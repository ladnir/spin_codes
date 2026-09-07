"""Local, resumable fixed-type search using the activation-density refinement.

All retained bounds use direct four-state evaluation. The table only proposes
witnesses. This is a binary64 search, not an outward certificate.
"""
import argparse
import hashlib
import heapq
import json
import math
from pathlib import Path
import sys
import time

import numpy as np

import certificate_search_core as core
import inner_candidate_boxes as initial
import joint_dense_witness_v2 as joint
import syndrome_density_v1 as density
import entropy_type_split_v1 as splitting
import refine_bch_dense_v1 as direct

base = core.base
LN2 = math.log(2)


def sources():
    result = core.inputs.source_hashes()
    result.update(core.outer_dependencies())
    for module in list(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path:
            path = Path(path).resolve()
            if path.is_relative_to(base.ROOT) and path.suffix == '.py':
                result[path.relative_to(base.ROOT).as_posix()] = base.sha(path)
    return result


class Refiner(joint.JointCoverRefiner):
    def __init__(self, spec, bands, directory):
        t, s, ac, kernel = core.inputs.load(spec['configuration'])
        counts = core.inputs.caps_module.caps()
        # Initialize the direct engine, avoiding the historical constructor's
        # unused table and writes into the shared landscape directory.
        direct.Refiner.__init__(self, counts, 256, t, s, ac,
                                [kernel.get(j, 0) for j in range(t+1)], spec['rows'], bands)
        self.warm = []
        self.epochs = density.Epochs(t, s, ac, [kernel.get(j, 0) for j in range(t+1)])
        self.xgrid = np.arange(-9., 1.001, .04)
        self.ygrid = np.arange(-14., 8.001, .04)
        self.size = 256*self.length
        identity = hashlib.sha256(json.dumps([spec, sources()], sort_keys=True).encode()).hexdigest()
        path = directory/f'table_{identity}.json'
        if path.exists():
            saved = base.read(path)
            core.require(saved['identity'] == identity, 'Wrong proposal table')
            self.table = np.array(saved['values'])
        else:
            theta = 1/(1+np.exp(-self.ygrid))
            values = []
            for i, x in enumerate(self.xgrid):
                epoch = self.epochs.at(math.exp(x))
                values.append(density.base.terminal_logs(
                    density.base.epoch_mixture(epoch, t, theta), self.size//t)/self.size)
                if i % 50 == 0:
                    print('proposal table', i+1, '/', len(self.xgrid), flush=True)
            self.table = np.array(values)
            base.write_new(path, dict(identity=identity, values=self.table.tolist()))


def seed(spec, minimum, directory):
    path = directory/'seed.json'
    if path.exists():
        saved = base.read(path)
        core.authenticate(saved, base.ROOT)
        core.require(saved['instance'] == spec and saved['minimum'] == minimum, 'Wrong seed instance')
        return saved['result']
    t, s, ac, kernel = core.inputs.load(spec['configuration'])
    counts = core.inputs.caps_module.caps()
    bands = [[w for w in counts if w < 256], [256]]
    model = initial.engine.DenseModel(counts, 256, t, s, ac,
        [kernel.get(j, 0) for j in range(t+1)], spec['rows'],
        [-4., -3., -2., -1., 0., .5, 1.], bands=bands)
    result = initial.plain(model.search(minimum=minimum, maximum_nodes=63, target_bits=60))
    direct.check_cover(result['selected_boxes'], spec['rows'], minimum)
    base.write_new(path, dict(instance=spec, minimum=minimum, result=result, source_sha256=sources()))
    return result


def run(directory, m, minimum, refinements, seconds, verify=False):
    directory = directory.resolve()
    directory.mkdir(parents=True, exist_ok=True)
    spec = core.instance('t128_s19', m)
    core.require(2 <= minimum <= spec['rows'] and refinements >= 0 and seconds > 0, 'Invalid budget or range')
    raw = seed(spec, minimum, directory)
    paths = sorted(directory.glob('cover_*.json'))
    saved = base.read(paths[-1]) if paths else None
    if saved:
        core.authenticate(saved, base.ROOT)
        core.require(saved['instance'] == spec and saved['minimum'] == minimum, 'Wrong resumed instance')
    core.require(not verify or saved is not None, 'No cover to replay')
    model = Refiner(spec, raw['bands'], directory)
    source_boxes = saved['leaves'] if saved else raw['selected_boxes']
    direct.check_cover(source_boxes, spec['rows'], minimum)
    leaves = {}
    for index, box in enumerate(source_boxes):
        value = model.direct(box['lower'], box['upper'], box['witness'])
        if saved:
            core.require(abs(value-box['own_log_bound']) < 2e-6, 'Direct witness replay failed')
        leaves[index] = dict(box, own_log_bound=value)
    def total():
        return float(np.logaddexp.reduce([b['own_log_bound'] for b in leaves.values()]))
    if verify:
        core.require(abs(total()-saved['log_upper']) < 2e-6, 'Union replay failed')
        print('Complete type cover and direct witnesses replayed; binary64 only.', flush=True)
        return
    heap = [(-box['own_log_bound'], i) for i, box in leaves.items()]
    heapq.heapify(heap)
    serial = len(leaves)
    started = time.monotonic()
    used = 0
    print('initial complete dense margin', -total()/LN2, flush=True)
    while heap and used < refinements and time.monotonic()-started < seconds and total() > -60*LN2:
        _, index = heapq.heappop(heap)
        box = model.refine_box(leaves[index], -80*LN2)
        leaves[index] = box
        used += 1
        if box['own_log_bound'] > -80*LN2:
            corners = initial.engine.typed.vertices(box['lower'], box['upper'], spec['rows'])
            children = splitting.split_box(box['lower'], box['upper'], corners)
            prepared = []
            for lo, hi in children:
                lo, hi = list(map(int, lo)), list(map(int, hi))
                if direct.lattice_count(lo, hi, spec['rows']):
                    prepared.append(dict(lower=lo, upper=hi, witness=box['witness'],
                        own_log_bound=model.direct(lo, hi, box['witness'])))
            if prepared:
                core.require(sum(direct.lattice_count(c['lower'], c['upper'], spec['rows']) for c in prepared)
                    == direct.lattice_count(box['lower'], box['upper'], spec['rows']), 'Split volume changed')
                del leaves[index]
                for child in prepared:
                    leaves[serial] = child
                    heapq.heappush(heap, (-child['own_log_bound'], serial))
                    serial += 1
        if used % 10 == 0:
            print('refinements', used, 'leaves', len(leaves), 'margin', -total()/LN2, flush=True)
    coverage = direct.check_cover(list(leaves.values()), spec['rows'], minimum)
    result = dict(status='COMPLETE_DENSE_TYPE_BINARY64_DIAGNOSTIC', instance=spec,
        minimum=minimum, bands=raw['bands'], leaves=list(leaves.values()), log_upper=total(),
        margin_bits=-total()/LN2, coverage=coverage, refinements_this_run=used,
        seconds_refining=time.monotonic()-started, source_sha256=sources(),
        outward_certificate=False, full_distance_proved=False)
    path = directory/f'cover_{len(paths):04d}.json'
    base.write_new(path, result)
    print('saved', path, 'margin', result['margin_bits'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--directory', type=Path, required=True)
    p.add_argument('--m', type=int, default=16)
    p.add_argument('--minimum', type=int, default=150)
    p.add_argument('--refinements', type=int, default=200)
    p.add_argument('--seconds', type=float, default=120)
    p.add_argument('--verify', action='store_true')
    a = p.parse_args()
    run(a.directory, a.m, a.minimum, a.refinements, a.seconds, a.verify)
