"""Jointly tune a dense point's probability proposal and output tilt."""
import argparse
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp
import parameter_activated_dense as source


def run(output,seed):
    assert not output.exists()
    model = source.base.grid.ladder.model
    saved = model.base.read(seed)
    model.authenticate(saved)
    assert saved['geometry'] == [128,64,20,20]
    rows = []
    with source.maps.use():
        checker = source.Dense(64,20,128,20,[])
        boxes = sorted(saved['dense']['selected_boxes'],key=lambda x:x['own_log_bound'],reverse=True)[:4]
        for box in boxes:
            point = source.base.typed.vertices(box['lower'],box['upper'],checker.length)[0]
            initial = np.array(box['witness']['proposal'])
            def objective(coordinates,detail=False):
                tilt = float(coordinates[-1])
                if not -12 <= tilt <= 2 or max(abs(coordinates[:-1])) > 40:
                    return 1e6
                logits = np.r_[coordinates[:-1],0.]
                proposal = np.exp(logits-logsumexp(logits))
                if len(checker.epoch_cache) > 4:
                    checker.epoch_cache.clear()
                moment = checker.moment(float(proposal@checker.ps),tilt)
                bound = float(source.base.typed.point_logs(point[None],checker.length,128,checker.log_gammas,
                    proposal,moment,checker.cutoff,math.exp(tilt))[0])
                return (bound,dict(tilt=tilt,probability_bank=0,proposal=proposal.tolist())) if detail else bound/(128*checker.length)
            start = np.r_[np.log(initial[:-1]/initial[-1]),box['witness']['tilt']]
            answer = minimize(objective,start,method='Nelder-Mead',
                options=dict(maxiter=450,xatol=1e-7,fatol=1e-11))
            bound,witness = objective(answer.x,True)
            rows.append(dict(point=point.tolist(),margin_bits=-bound/math.log(2),witness=witness,
                evaluations=answer.nfev,optimizer_success=bool(answer.success)))
            print('joint point',rows[-1],flush=True)
        sources = source.base.grid.ladder.candidate.sources()
        _,path = source.base.grid.outer(128)
        sources[path.relative_to(source.base.grid.ROOT).as_posix()] = source.base.grid.maps.sha(path)
        sources[seed.relative_to(source.base.grid.ROOT).as_posix()] = model.base.sha(seed)
        model.base.write_new(output,dict(status='BINARY64_IMT_JOINT_POINT_RETUNE',geometry=saved['geometry'],
            inner=checker.record,points=rows,full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=Path,required=True)
    a = p.parse_args()
    run(a.output.resolve(),a.seed.resolve())
