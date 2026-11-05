"""Bounded fresh marker calibration followed by complete FLOAT support sums.

Sources propose z/alpha only. No old bound, marker, or construction guarantee
is imported. A geometric interpolation proposes marks, never a probability
bound. All outputs remain numerical screens, not directed certificates.
"""

from fractions import Fraction as F
from hashlib import sha256
import argparse
import json
from math import log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_midpoint_counts as mode


provider = mode.near.provider


def source_pins():
    pins = mode.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run(weight, z, alpha, *, max_seconds=180., source=None):
    if type(weight) is not int or not 2 <= weight < 4096:
        raise ValueError('Fresh occupancy initialization supports weights2..4095')
    pins, started = source_pins(), monotonic()
    budget = provider.old.Budget(max_seconds)
    coeff = provider.AliasCoefficients(z, alpha, state=30, _budget=budget)
    half = weight//2
    targets = sorted({0, 1, 2, 4, 8, 16, *(weight*j//16 for j in range(1, 8)),
                      half-16, half-8, half-4, half-2, half-1, half})
    targets = [f for f in targets if 0 <= f <= half]
    chosen, fits, cache = [], [], {}

    def numerators(point):
        budget.check()
        point = tuple(map(float, point))
        if point not in cache:
            cache[point] = np.asarray(provider.alias.power(coeff.pair_logs(point, _budget=budget), 64)).ravel()
        return cache[point]

    for f in targets:
        singles, empty = weight-2*f, 4096-weight+f
        initial = (log2(singles/empty) if singles else -64., log2(f/empty) if f else -64.)
        for entry in range(4):
            if monotonic()-started > .75*max_seconds:
                raise TimeoutError('Marker calibration exhausted its reserved budget')
            best = min([initial, *chosen], key=lambda point: float(numerators(point)[entry]-singles*point[0]-f*point[1]))
            best_value = float(numerators(best)[entry]-singles*best[0]-f*best[1])
            calls = 0

            def objective(raw):
                nonlocal best, best_value, calls
                calls += 1
                budget.check()
                if monotonic()-started > .75*max_seconds:
                    raise TimeoutError('Marker calibration exhausted its reserved budget')
                point = (float(raw[0]), -64.) if not f else (-64., float(raw[0])) if not singles else tuple(map(float, raw))
                value = float(numerators(point)[entry]-singles*point[0]-f*point[1])
                if value < best_value:
                    best, best_value = point, value
                return value

            start = [best[0]] if not f else [best[1]] if not singles else best
            fitted = minimize(objective, start, method='Powell', bounds=[(-64., 8.)]*len(start),
                              options=dict(maxfev=128, xtol=1e-6, ftol=1e-8))
            if best not in chosen:
                chosen.append(best)
            fits.append(dict(single=singles, full=f, entry=entry, calls=calls,
                             marker_logs_exact=[str(F(v)) for v in best],
                             objective=best_value, optimizer_success=bool(fitted.success)))
    if not 1 <= len(chosen) <= 128:
        raise ArithmeticError('Fresh marker count exceeds the bank interface')
    bank = provider.AliasBankScreen(z, alpha, chosen, 0., state=30, distance=F(2, 5),
              grid=128, floor_policy='positive-fallback', _coefficients=coeff, _budget=budget)
    regions = bank.regional(weight)
    old = mode.near.support_sum(4096, weight, regions, max_singletons=64, budget=budget)
    refined = mode.support_sum(4096, weight, regions, budget=budget)
    prefix = bank.constant+bank.alpha*bank.density*weight
    old['log2_upper_screen'] = prefix+old['log2_support_sum_screen']
    refined['log2_upper_screen'] = prefix+refined['log2_support_sum_screen']
    if pins != source_pins():
        raise RuntimeError('Calibration or proof sources changed')
    return dict(schema='gamma30-fresh-midpoint-marker-calibration-float-v1',
                active_groups=weight, state_dimension=30, target_distance='2/5',
                source=source, source_role='initial z/alpha proposal only; no inherited markers or bounds',
                witness=bank.witness(), source_sha256=pins, fits=fits,
                marker_targets=targets, marker_count=len(chosen), distinct_points=len(cache),
                max_seconds=max_seconds, elapsed_seconds=monotonic()-started,
                baseline_same_regions=old, midpoint_same_regions=refined,
                log2_outer_prefix=prefix, regional_alias_diagnostics=bank.regional_diagnostics,
                regions=regions.tolist(), regional_logs_sha256=sha256(regions.tobytes()).hexdigest(),
                profile_pairs_evaluated=(weight//2+1)**2, all_profiles_retained=True,
                common_potential_used=False, construction=bank.metadata()['construction'],
                directed=False, whole_code_certificate=False, old_numerical_bounds_reused=False,
                scope='complete FLOAT sum at one weight; not a directed or whole-code certificate')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--source-low', type=Path, required=True)
    parser.add_argument('--source-high', type=Path)
    parser.add_argument('--high-fraction', default='1/2')
    parser.add_argument('--max-seconds', type=float, default=180.)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Output must be unused')
    fraction = F(args.high_fraction)
    if not 0 <= fraction <= 1:
        parser.error('Interpolation fraction must be in[0,1]')
    sources, rows = [], []
    for path in (args.source_low, args.source_high or args.source_low):
        raw = path.read_bytes()
        data = json.loads(raw)
        rows.append(data.get('best', data)['witness'])
        sources.append(dict(path=str(path.resolve()), sha256=sha256(raw).hexdigest()))
    marks = [2.**((1-float(fraction))*log2(float(F(rows[0][key])))
                    + float(fraction)*log2(float(F(rows[1][key])))) for key in ('z', 'alpha')]
    provenance = dict(sources=sources, high_fraction=str(fraction),
                      rule='geometric interpolation of z and alpha only; no interpolated probability claim')
    result = run(args.r, *marks, max_seconds=args.max_seconds, source=provenance)
    if any(sha256(Path(item['path']).read_bytes()).hexdigest() != item['sha256'] for item in sources):
        raise RuntimeError('Witness source changed')
    provider.search.write_new(args.output, result)
    print(json.dumps({key: result[key] for key in ('active_groups', 'marker_count', 'elapsed_seconds')}))
    print(json.dumps(dict(baseline=result['baseline_same_regions']['log2_upper_screen'],
                         midpoint=result['midpoint_same_regions']['log2_upper_screen'],
                         largest=result['midpoint_same_regions']['largest_profile'])))


if __name__ == '__main__':
    main()
