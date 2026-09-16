"""Refine and replay the dense suffix after a verified sparse-range bridge."""
import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
import parameter_range_bridge as bridge
import parameter_activation_partition as partition
import type_box_coverage

bank = bridge.bank
STATUS = 'BINARY64_IMT_SCOPED_ACTIVATION_COVER'


def check(record,moments):
    assert record['geometry'] == [128,64,20,20]
    assert record['inner'] == json.loads(json.dumps(moments.original.record))
    dense = record['dense']
    assert dense['occupation_max'] == 16384 and 65 <= dense['occupation_min'] <= 16384
    count = type_box_coverage.check(dense['selected_boxes'],16384,dense['occupation_min'],3)
    assert str(count) == dense['integer_types_checked']
    return count


def run(seed_path,replay_path,output,rounds):
    model,base = bank.model,bank.base
    saved,replay = [model.base.read(p) for p in (seed_path,replay_path)]
    model.authenticate(saved); model.authenticate(replay)
    statuses = {'BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND':
                    'VERIFIED_BINARY64_IMT_RANGE_BRIDGED_FULL_PARAMETER_BOUND',
                STATUS:'VERIFIED_'+STATUS}
    assert saved['status'] in statuses and replay['status'] == statuses[saved['status']]
    assert replay['producer_sha256'] == model.base.sha(seed_path)
    with bank.previous.maps.use():
        moments = bridge.Moments()
        count = check(saved,moments)
        checker = moments.activation
        dense = copy.deepcopy(saved['dense'])
        boxes = dense['selected_boxes']
        # Freshly reconstruct every starting value, including bounds after clipping.
        for box in boxes:
            value = moments.bound(box)
            assert abs(value-box['own_log_bound']) < 1e-7
        entries,progress,attempted = [],[],set()
        activation_cache = {}
        def moment(theta,tilt):
            key = (theta,tilt)
            if key not in activation_cache:
                if len(checker.epoch_cache)>4:
                    checker.epoch_cache.clear()
                activation_cache[key] = checker.moment(theta,tilt)
                if len(activation_cache)%100==0:
                    print('scoped moments',len(activation_cache),flush=True)
            return activation_cache[key]

        for iteration in range(rounds):
            eligible = [i for i,b in enumerate(boxes) if i not in attempted and b['own_log_bound']>-60*math.log(2)]
            if not eligible:
                break
            index = max(eligible,key=lambda i:boxes[i]['own_log_bound'])
            box = boxes[index]
            before = box['own_log_bound']
            witness = box['witness']
            theta = float(np.array(witness['proposal'])@checker.ps)
            logit = math.log(theta/(1-theta))
            trials = []
            print('scoped target',iteration+1,box['lower'],box['upper'],-before/math.log(2),flush=True)
            def objective(x):
                logit,tilt = map(float,x)
                if not -18<=logit<=18 or not -8<=tilt<=1:
                    return 1e6
                theta = float(expit(logit))
                value,_ = bank.proposal_bound(checker,box,theta,tilt,moment(theta,tilt))
                trials.append((value,theta,tilt))
                return value/(128*16384)
            for tilt in (witness['tilt'],-.60):
                minimize(objective,[logit,tilt],method='Nelder-Mead',
                    options=dict(maxiter=110,xatol=1e-6,fatol=1e-10))
            _,theta,tilt = min(trials)
            term = moment(theta,tilt)
            changed = 0
            for current in boxes:
                value,w = bank.proposal_bound(checker,current,theta,tilt,term)
                if value<current['own_log_bound']:
                    current.update(own_log_bound=value,witness=w)
                    changed+=1
            entries.append(dict(theta=theta,tilt=tilt,moment_log=term,changed=changed))
            after = box['own_log_bound']
            progress.append(dict(lower=box['lower'],upper=box['upper'],
                before_margin_bits=-before/math.log(2),after_margin_bits=-after/math.log(2)))
            if before-after<math.log(2):
                attempted.add(index)
            # A still-weak nonsingleton needs a smaller common-witness domain.
            if after>-60*math.log(2):
                children = partition.split_box(box,16384)
                if children:
                    for child in children:
                        child['witness'] = copy.deepcopy(box['witness'])
                        child['own_log_bound'] = moments.bound(child)
                    boxes[index]=children[0]
                    boxes.append(children[1])
                    attempted.discard(index)
            log_union = float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes]))
            print('scoped step',iteration+1,'dense bits',-log_union/math.log(2),
                  'weak boxes',sum(b['own_log_bound']>-60*math.log(2) for b in boxes),flush=True)
        count = type_box_coverage.check(boxes,16384,dense['occupation_min'],3)
        dense.update(integer_types_checked=str(count),
                     log_union_upper=float(np.logaddexp.reduce([b['own_log_bound'] for b in boxes])))
        dense['margin_bits'] = -dense['log_union_upper']/math.log(2)
        sources = base.grid.ladder.candidate.sources()
        for path in (seed_path,replay_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()] = model.base.sha(path)
        result = dict(status=STATUS,geometry=saved['geometry'],inner=saved['inner'],dense=dense,
            entries=entries,progress=progress,full_distance_proved=False,source_sha256=sources)
        model.base.write_new(output,result)


def verify(path,output):
    model,base = bank.model,bank.base
    saved = model.base.read(path)
    model.authenticate(saved)
    assert saved['status'] == STATUS
    with bank.previous.maps.use():
        moments = bridge.Moments()
        count = check(saved,moments)
        values = [moments.bound(box) for box in saved['dense']['selected_boxes']]
        errors = [abs(v-b['own_log_bound']) for v,b in zip(values,saved['dense']['selected_boxes'],strict=True)]
        total = float(np.logaddexp.reduce(values))
        errors += [abs(total-saved['dense']['log_union_upper']),
                   abs(-total/math.log(2)-saved['dense']['margin_bits'])]
        assert max(errors)<1e-7,max(errors)
        sources = base.grid.ladder.candidate.sources()
        sources[path.relative_to(base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='VERIFIED_'+STATUS,producer_sha256=model.base.sha(path),
            geometry=saved['geometry'],covered_occupancies=[saved['dense']['occupation_min'],16384],
            integer_types_checked=str(count),margin_bits=-total/math.log(2),
            maximum_log_error=max(errors),full_distance_proved=False,source_sha256=sources))
        print('verified scoped dense bits',-total/math.log(2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path)
    p.add_argument('--replay',type=Path)
    p.add_argument('--verify',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rounds',type=int,default=4)
    a=p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    if a.verify:
        if a.seed or a.replay: p.error('--verify does not take seed/replay')
        verify(a.verify.resolve(),a.output.resolve())
    else:
        if not a.seed or not a.replay: p.error('search requires seed and replay')
        run(a.seed.resolve(),a.replay.resolve(),a.output.resolve(),a.rounds)
