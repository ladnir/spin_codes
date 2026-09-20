"""Reuse full-syndrome activation moments over an existing dense partition.

For each fixed input probability and tilt, optimize the positive proposal
subject to the same input probability. Only whole valid moment bounds are
compared. All evaluations are nearest binary64, not outward certificates.
"""
import argparse
import copy
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit
import parameter_syndrome_activation as activation
import parameter_joint_mixed_dense as previous
import type_box_coverage

base = previous.base
model = base.grid.ladder.model


def proposal_bound(checker, box, theta, tilt, moment):
    corners = base.typed.vertices(box['lower'], box['upper'], checker.length)
    penalty = base.typed.lattice_log_count(box['lower'], box['upper'])
    # ps=(0,1/2,1). Let u be the all-one proposal probability.
    assert np.array_equal(checker.ps, [0., .5, 1.])
    low, high = max(0., 2*theta-1), theta
    epsilon = min(2.**-40, (high-low)*1e-6)
    low, high = low+epsilon, high-epsilon
    assert low < high

    def evaluate(u):
        proposal = np.array([1-2*theta+u, 2*(theta-u), u])
        assert np.all(proposal > 0)
        value = float(max(base.typed.point_logs(corners, checker.length, checker.block,
            checker.log_gammas, proposal, moment, checker.cutoff, math.exp(tilt)))) + penalty
        return value, proposal

    result = minimize_scalar(lambda u: evaluate(u)[0], bounds=(low, high), method='bounded',
                             options=dict(xatol=1e-12, maxiter=100))
    value, proposal = min((evaluate(u) for u in (low, high, result.x)), key=lambda row: row[0])
    return value, dict(family='syndrome_activation', probability_bank=0,
                      proposal=proposal.tolist(), tilt=float(tilt))


def run(cover_path, replay_path, point_path, output, rounds):
    if output.exists():
        raise FileExistsError(output)
    saved, replay, point = [model.base.read(p) for p in (cover_path, replay_path, point_path)]
    for record in (saved, replay, point):
        model.authenticate(record)
    assert saved['status'] == 'BINARY64_IMT_JOINT_MIXED_DENSE_COVER'
    assert replay['status'] == 'VERIFIED_BINARY64_IMT_DENSE_COVER'
    assert replay['producer_sha256'] == model.base.sha(cover_path)
    assert point['status'] == 'BINARY64_IMT_SYNDROME_ACTIVATION_REFINED_POINT'
    assert saved['geometry'] == point['geometry'] == replay['geometry'] == [128,64,20,20]
    with previous.maps.use():
        checker = activation.Dense(64,20,128,20,[])
        assert saved['inner'] == point['inner'] == json.loads(json.dumps(checker.record))
        dense = copy.deepcopy(saved['dense'])
        boxes = dense['selected_boxes']
        count = type_box_coverage.check(boxes,checker.length,dense['occupation_min'],3)
        bank = []
        moment_cache = {}

        def moment(theta, tilt):
            key = (float(theta),float(tilt))
            if key not in moment_cache:
                if len(checker.epoch_cache) > 4:
                    checker.epoch_cache.clear()
                moment_cache[key] = checker.moment(*key)
            return moment_cache[key]

        def apply(theta, tilt):
            term = moment(theta,tilt)
            changed = 0
            for box in boxes:
                value,witness = proposal_bound(checker,box,theta,tilt,term)
                if value < box['own_log_bound']:
                    box.update(own_log_bound=value,best_log_bound=value,witness=witness)
                    changed += 1
            bank.append(dict(input_probability=float(theta),tilt=float(tilt),moment_log=term,
                             boxes_improved=changed))
            union = float(np.logaddexp.reduce([r['own_log_bound'] for r in boxes]))
            print('activation bank',len(bank),'changed',changed,'dense bits',-union/math.log(2),
                  'weak boxes',sum(r['own_log_bound'] > -60*math.log(2) for r in boxes),flush=True)

        first = point['selected']
        apply(float(np.array(first['proposal'])@checker.ps),first['tilt'])
        progress = []
        for iteration in range(rounds):
            worst = max(boxes,key=lambda row:row['own_log_bound'])
            if worst['own_log_bound'] < -60*math.log(2):
                break
            before = worst['own_log_bound']
            w = worst['witness']
            theta = float(np.array(w['proposal'])@checker.ps)
            start = [math.log(theta/(1-theta)),w['tilt']]
            def objective(x):
                logit,tilt = map(float,x)
                if not -18 <= logit <= 18 or not -8 <= tilt <= 1:
                    return 1e6
                theta = float(expit(logit))
                value,_ = proposal_bound(checker,worst,theta,tilt,moment(theta,tilt))
                return value/(checker.block*checker.length)
            answer = minimize(objective,start,method='Nelder-Mead',
                options=dict(maxiter=150,xatol=1e-6,fatol=1e-10))
            theta,tilt = float(expit(answer.x[0])),float(answer.x[1])
            if not -18 <= answer.x[0] <= 18 or not -8 <= tilt <= 1:
                raise ValueError('Optimizer left permitted domain')
            apply(theta,tilt)
            progress.append(dict(lower=worst['lower'],upper=worst['upper'],
                before_margin_bits=-before/math.log(2),
                after_margin_bits=-worst['own_log_bound']/math.log(2),success=bool(answer.success)))
            if before-worst['own_log_bound'] < 1e-5:
                print('worst box not improved; a different witness or partition is needed',flush=True)
                break
        dense['log_union_upper'] = float(np.logaddexp.reduce([r['own_log_bound'] for r in boxes]))
        sources = base.grid.ladder.candidate.sources()
        _,outer_path = base.grid.outer(128)
        for path in (cover_path,replay_path,point_path,outer_path):
            sources[path.relative_to(base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='BINARY64_IMT_ACTIVATION_BANK_DENSE_COVER',
            geometry=saved['geometry'],inner=checker.record,dense=dense,bank=bank,progress=progress,
            integer_types_checked=str(count),margin_bits=-dense['log_union_upper']/math.log(2),
            full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('cover','replay','point','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--rounds',type=int,default=8)
    args = parser.parse_args()
    run(args.cover.resolve(),args.replay.resolve(),args.point.resolve(),args.output.resolve(),args.rounds)
