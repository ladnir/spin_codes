"""Extend a replayed activation bank at its weakest remaining count boxes."""
import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
import parameter_activation_bank as bank
import type_box_coverage


def run(seed_path,replay_path,output,rounds):
    if output.exists():
        raise FileExistsError(output)
    base,model = bank.base,bank.model
    saved,replay = [model.base.read(p) for p in (seed_path,replay_path)]
    model.authenticate(saved); model.authenticate(replay)
    assert saved['status'] == 'BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER'
    assert replay['status'] == 'VERIFIED_BINARY64_IMT_DENSE_COVER'
    assert replay['producer_sha256'] == model.base.sha(seed_path)
    assert saved['geometry'] == replay['geometry'] == [128,64,20,20]
    with bank.previous.maps.use():
        checker = bank.activation.Dense(64,20,128,20,[])
        assert saved['inner'] == json.loads(json.dumps(checker.record))
        dense = copy.deepcopy(saved['dense'])
        boxes = dense['selected_boxes']
        count = type_box_coverage.check(boxes,checker.length,65,3)
        entries = copy.deepcopy(saved['bank'])
        attempted,cache,progress = set(),{},[]

        def moment(theta,tilt):
            key = (float(theta),float(tilt))
            if key not in cache:
                if len(checker.epoch_cache) > 4:
                    checker.epoch_cache.clear()
                cache[key] = checker.moment(*key)
                if len(cache) % 100 == 0:
                    print('frontier moments',len(cache),flush=True)
            return cache[key]

        for iteration in range(rounds):
            eligible = [i for i,b in enumerate(boxes)
                        if i not in attempted and b['own_log_bound'] > -60*math.log(2)]
            if not eligible:
                break
            index = max(eligible,key=lambda i:boxes[i]['own_log_bound'])
            box = boxes[index]
            before = box['own_log_bound']
            w = box['witness']
            theta = float(np.array(w['proposal'])@checker.ps)
            logit = math.log(theta/(1-theta))
            print('frontier box',iteration+1,box['lower'],box['upper'],
                  'bits',-before/math.log(2),flush=True)
            candidates = []
            def objective(x):
                logit,tilt = map(float,x)
                if not -18 <= logit <= 18 or not -8 <= tilt <= 1:
                    return 1e6
                theta = float(expit(logit))
                value,_ = bank.proposal_bound(checker,box,theta,tilt,moment(theta,tilt))
                candidates.append((value,theta,tilt))
                return value/(checker.block*checker.length)
            answers = []
            for tilt in (w['tilt'],-.60):
                result = minimize(objective,[logit,tilt],method='Nelder-Mead',
                    options=dict(maxiter=110,xatol=1e-6,fatol=1e-10))
                answers.append(dict(success=bool(result.success),evaluations=result.nfev))
            value,theta,tilt = min(candidates)
            term = moment(theta,tilt)
            changed = 0
            for current in boxes:
                value,witness = bank.proposal_bound(checker,current,theta,tilt,term)
                if value < current['own_log_bound']:
                    current.update(own_log_bound=value,best_log_bound=value,witness=witness)
                    changed += 1
            entries.append(dict(input_probability=theta,tilt=tilt,moment_log=term,boxes_improved=changed))
            if before-box['own_log_bound'] < 1e-5:
                attempted.add(index)
            progress.append(dict(lower=box['lower'],upper=box['upper'],
                before_margin_bits=-before/math.log(2),after_margin_bits=-box['own_log_bound']/math.log(2),
                optimizers=answers))
            union = float(np.logaddexp.reduce([r['own_log_bound'] for r in boxes]))
            print('frontier bank',len(entries),'changed',changed,'dense bits',-union/math.log(2),
                  'weak boxes',sum(r['own_log_bound'] > -60*math.log(2) for r in boxes),flush=True)
        dense['log_union_upper'] = float(np.logaddexp.reduce([r['own_log_bound'] for r in boxes]))
        sources = base.grid.ladder.candidate.sources()
        _,outer_path = base.grid.outer(128)
        for path in (seed_path,replay_path,outer_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER',
            geometry=saved['geometry'],inner=checker.record,dense=dense,bank=entries,progress=progress,
            integer_types_checked=str(count),margin_bits=-dense['log_union_upper']/math.log(2),
            full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('seed','replay','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--rounds',type=int,default=6)
    args = parser.parse_args()
    run(args.seed.resolve(),args.replay.resolve(),args.output.resolve(),args.rounds)
