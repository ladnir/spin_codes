"""Jointly refine the full-syndrome activation point; binary64 only."""
import argparse
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
import parameter_syndrome_activation as source


def run(seed, output, iterations):
    if output.exists():
        raise FileExistsError(output)
    model = source.prior.base.grid.ladder.model
    saved = model.base.read(seed)
    model.authenticate(saved)
    assert saved['status'] == 'BINARY64_IMT_DENSE_FACE_MULTISTART'
    assert saved['geometry'] == [128, 64, 20, 20]
    with source.prior.maps.use():
        checker = source.Dense(64, 20, 128, 20, [])
        point = np.array(saved['point'])
        assert point[2] == 0
        epsilon = 2. ** -40
        best = max(saved['refined'], key=lambda r: r['margin_bits'])
        probability = best['proposal'][1] / (1 - epsilon)
        start = np.array([math.log(probability / (1 - probability)), best['tilt']])
        evaluated = {}
        started = time.monotonic()

        def evaluate(x):
            logits, tilt = map(float, x)
            if not -20 <= logits <= 20 or not -8 <= tilt <= 1:
                return math.inf
            key = (logits, tilt)
            if key not in evaluated:
                p = float(expit(logits))
                proposal = np.array([(1-p)*(1-epsilon), p*(1-epsilon), epsilon])
                if len(checker.epoch_cache) > 4:
                    checker.epoch_cache.clear()
                moment = checker.moment(float(proposal @ checker.ps), tilt)
                value = float(source.prior.base.typed.point_logs(point[None], checker.length,
                    128, checker.log_gammas, proposal, moment, checker.cutoff, math.exp(tilt))[0])
                evaluated[key] = dict(tilt=tilt, proposal=proposal.tolist(), moment_log=moment,
                    margin_bits=-value/math.log(2), log_bound=value)
                if len(evaluated) % 20 == 0:
                    print('activation evaluations', len(evaluated), 'best bits',
                          max(row['margin_bits'] for row in evaluated.values()),
                          'seconds', round(time.monotonic()-started, 1), flush=True)
            return evaluated[key]['log_bound'] / (128 * checker.length)

        initial_value = evaluate(start)
        results = []
        for offset in (np.array([0., 0.]), np.array([-.25, .20])):
            answer = minimize(evaluate, start+offset, method='Nelder-Mead',
                options=dict(maxiter=iterations, xatol=1e-7, fatol=1e-11))
            results.append(dict(success=bool(answer.success), evaluations=answer.nfev,
                                coordinates=answer.x.tolist(), message=str(answer.message)))
        selected = max(evaluated.values(), key=lambda r: r['margin_bits'])
        sources = source.prior.base.grid.ladder.candidate.sources()
        _, outer_path = source.prior.base.grid.outer(128)
        for path in (seed, outer_path):
            sources[path.relative_to(source.prior.base.grid.ROOT).as_posix()] = model.base.sha(path)
        record = dict(status='BINARY64_IMT_SYNDROME_ACTIVATION_REFINED_POINT',
            geometry=saved['geometry'], inner=checker.record, point=point.tolist(),
            initial_margin_bits=-initial_value*(128*checker.length)/math.log(2),
            selected=selected, optimizer_results=results, evaluations=len(evaluated),
            elapsed_seconds=time.monotonic()-started, full_distance_proved=False,
            source_sha256=sources)
        model.base.write_new(output, record)
        print('refined activation', selected, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seed', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--iterations', type=int, default=180)
    args = parser.parse_args()
    run(args.seed.resolve(), args.output.resolve(), args.iterations)
