"""Refine weak count boxes while retaining authenticated IMT witnesses."""
import argparse
import copy
import heapq
import json
import math
from pathlib import Path

import numpy as np
import parameter_activation_bank as bank
import type_box_coverage


def split_box(box, total):
    """Bisect a longest varying coordinate; discard no feasible integer type."""
    lo,hi = box['lower'],box['upper']
    corners = bank.base.typed.vertices(lo,hi,total)
    widths = np.max(corners,axis=0)-np.min(corners,axis=0)
    axis = int(np.argmax(widths))
    if widths[axis] == 0:
        return []
    lower,upper = int(min(corners[:,axis])),int(max(corners[:,axis]))
    cut = (lower+upper)//2
    children = []
    for left in (True,False):
        a,b = list(lo),list(hi)
        if left:
            b[axis] = cut
        else:
            a[axis] = cut+1
        if type_box_coverage.lattice_count(a,b,total):
            children.append(dict(lower=a,upper=b))
    assert len(children) == 2
    assert sum(type_box_coverage.lattice_count(c['lower'],c['upper'],total) for c in children) == \
           type_box_coverage.lattice_count(lo,hi,total)
    return children


def run(seed_path,replay_path,output,splits):
    if output.exists():
        raise FileExistsError(output)
    model,base = bank.model,bank.base
    saved,replay = [model.base.read(p) for p in (seed_path,replay_path)]
    model.authenticate(saved); model.authenticate(replay)
    assert saved['status'] == 'BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER'
    assert replay['status'] == 'VERIFIED_BINARY64_IMT_DENSE_COVER'
    assert replay['producer_sha256'] == model.base.sha(seed_path)
    assert saved['geometry'] == replay['geometry'] == [128,64,20,20]
    with bank.previous.maps.use():
        checker = bank.previous.Dense(64,20,128,20,[])
        activation = bank.activation.Dense(64,20,128,20,[])
        assert saved['inner'] == json.loads(json.dumps(checker.record))
        dense = copy.deepcopy(saved['dense'])
        type_box_coverage.check(dense['selected_boxes'],checker.length,65,3)
        cache = {}
        def moment(w):
            theta = float(np.array(w['proposal'])@checker.ps)
            key = (w.get('family'),theta,w['tilt'])
            if key not in cache:
                if key[0] == 'syndrome_activation':
                    if len(activation.epoch_cache) > 4:
                        activation.epoch_cache.clear()
                    cache[key] = activation.moment(theta,w['tilt'])
                elif key[0] is None:
                    cache[key] = checker.moment(theta,w['tilt'])
                else:
                    cache[key] = checker.replay_moment(theta,w)
            return cache[key]
        witnesses = []
        for entry in saved['bank']:
            theta,tilt = entry['input_probability'],entry['tilt']
            # Recompute each bank moment, not its saved display value.
            term = activation.moment(theta,tilt)
            assert abs(term-entry['moment_log']) < 1e-7
            witnesses.append((theta,tilt,term))

        def evaluate(child,parent_witness):
            proposal = np.array(parent_witness['proposal'])
            corners = base.typed.vertices(child['lower'],child['upper'],checker.length)
            value = float(max(base.typed.point_logs(corners,checker.length,128,checker.log_gammas,
                proposal,moment(parent_witness),checker.cutoff,math.exp(parent_witness['tilt']))))
            value += base.typed.lattice_log_count(child['lower'],child['upper'])
            selected = parent_witness
            for theta,tilt,term in witnesses:
                candidate,w = bank.proposal_bound(checker,child,theta,tilt,term)
                if candidate < value:
                    value,selected = candidate,w
            child.update(own_log_bound=value,witness=copy.deepcopy(selected))
            return child

        boxes = {i:b for i,b in enumerate(dense['selected_boxes'])}
        queue = [(-b['own_log_bound'],i) for i,b in boxes.items()]
        heapq.heapify(queue)
        next_id = len(boxes)
        unsplittable = []
        used = 0
        while used < splits and queue:
            _,i = heapq.heappop(queue)
            parent = boxes[i]
            if parent['own_log_bound'] < -60*math.log(2):
                break
            children = split_box(parent,checker.length)
            if not children:
                unsplittable.append(i)
                continue
            del boxes[i]
            for child in children:
                child = evaluate(child,parent['witness'])
                boxes[next_id] = child
                heapq.heappush(queue,(-child['own_log_bound'],next_id))
                next_id += 1
            used += 1
            if used % 64 == 0:
                union = float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes.values()]))
                print('activation splits',used,'dense bits',-union/math.log(2),
                      'unsplittable weak types',len(unsplittable),flush=True)
        dense['selected_boxes'] = list(boxes.values())
        count = type_box_coverage.check(dense['selected_boxes'],checker.length,65,3)
        dense['log_union_upper'] = float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes.values()]))
        sources = base.grid.ladder.candidate.sources()
        _,outer_path = base.grid.outer(128)
        for path in (seed_path,replay_path,outer_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER',
            geometry=saved['geometry'],inner=checker.record,dense=dense,bank=saved['bank'],
            refinement_splits=used,weak_singletons=[boxes[i] for i in unsplittable],
            integer_types_checked=str(count),margin_bits=-dense['log_union_upper']/math.log(2),
            full_distance_proved=False,source_sha256=sources))
        print('refined activation partition',len(boxes),'boxes, bits',
              -dense['log_union_upper']/math.log(2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('seed','replay','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--splits',type=int,default=512)
    args = parser.parse_args()
    run(args.seed.resolve(),args.replay.resolve(),args.output.resolve(),args.splits)
