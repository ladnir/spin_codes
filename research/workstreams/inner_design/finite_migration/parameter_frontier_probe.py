"""Resolve the dense pilot's worst zero/all-ordinary face before more covers.

This is a nearest-binary64 point diagnostic, never a distance certificate.
The all-one proposal remains positive even when its count is zero.
"""
import argparse
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
import parameter_activated_dense as source


def run(seed, output):
    assert not output.exists()
    model = source.base.grid.ladder.model
    saved = model.base.read(seed)
    model.authenticate(saved)
    assert saved['geometry'] == [128,64,20,20]
    with source.maps.use():
        checker = source.Dense(64,20,128,20,[])
        box = max(saved['dense']['selected_boxes'], key=lambda r:r['own_log_bound'])
        vertices = source.base.typed.vertices(box['lower'],box['upper'],checker.length)
        assert len(vertices) == 1
        point = vertices[0]
        assert point[2] == 0
        epsilon = 2.**-40
        samples = []

        def evaluate(probability, tilt):
            if not 0 < probability < 1 or not -12 <= tilt <= 2:
                return math.inf, None
            proposal = np.array([(1-probability)*(1-epsilon), probability*(1-epsilon),epsilon])
            theta = float(proposal@checker.ps)
            if len(checker.epoch_cache) > 4:
                checker.epoch_cache.clear()
            moment = checker.moment(theta,tilt)
            value = float(source.base.typed.point_logs(point[None],checker.length,128,
                checker.log_gammas,proposal,moment,checker.cutoff,math.exp(tilt))[0])
            row = dict(tilt=float(tilt),proposal=proposal.tolist(),input_probability=theta,
                margin_bits=-value/math.log(2),moment_log=moment,
                cutoff_log=checker.cutoff*math.exp(tilt),
                likelihood_log=value-moment-checker.cutoff*math.exp(tilt))
            return value,row

        for tilt in np.arange(-5.,1.501,.125):
            values = []
            for probability in np.linspace(.05,.95,37):
                value,row = evaluate(float(probability),float(tilt))
                samples.append((value,row)); values.append(value)
            print('frontier tilt',float(tilt),'best bits',-min(values)/math.log(2),flush=True)
        # Start in separated output-tilt basins, not just variants of one minimum.
        seeds = []
        for value,row in sorted(samples,key=lambda r:r[0]):
            if all(abs(row['tilt']-prior['tilt']) >= .375 for prior in seeds):
                seeds.append(row)
            if len(seeds) == 6:
                break
        refined = []
        for row in seeds:
            probability = row['proposal'][1]/(1-epsilon)
            def objective(x):
                value,_ = evaluate(float(expit(x[0])),float(x[1]))
                return value/(128*checker.length)
            answer = minimize(objective,[math.log(probability/(1-probability)),row['tilt']],
                method='Nelder-Mead',options=dict(maxiter=350,xatol=1e-7,fatol=1e-11))
            _,result = evaluate(float(expit(answer.x[0])),float(answer.x[1]))
            result.update(success=bool(answer.success),evaluations=answer.nfev,initial_tilt=row['tilt'])
            refined.append(result)
            print('refined face',result,flush=True)
        sources = source.base.grid.ladder.candidate.sources()
        _,outer_path = source.base.grid.outer(128)
        for path in (seed,outer_path):
            sources[path.relative_to(source.base.grid.ROOT).as_posix()] = model.base.sha(path)
        model.base.write_new(output,dict(status='BINARY64_IMT_DENSE_FACE_MULTISTART',
            geometry=saved['geometry'],inner=checker.record,point=point.tolist(),
            coarse_best=sorted(samples,key=lambda r:r[0])[0][1],refined=refined,
            grid_evaluations=len(samples),full_distance_proved=False,source_sha256=sources))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--seed',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    run(a.seed.resolve(),a.output.resolve())
