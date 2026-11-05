"""Bounded FLOAT z/alpha/common-potential proposals for a frozen state bank.

Every mark pair builds one coefficient bank and one FFT bank. Potential
fitting reuses that bank. Exact source markers stay fixed, including their
selected-weight calibration. Immutable baseline checkpoints survive a
later timeout. Numerical bounds from the source are never reused.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_full_rs_state_alias_screen as provider


def source_pins():
    pins = provider.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def prepare_proposal(artifact, provenance, *, r, state, distance, old_source=False,
                     z=None, alpha=None):
    if type(r) is not int or not 1 <= r <= provider.old.GROUPS:
        raise ValueError('Selected weight must lie in1..8192')
    calibrated = {key: artifact[key] for key in ('center', 'selected_weight') if key in artifact}
    if any(type(value) is not int or value != r for value in calibrated.values()):
        raise ValueError('Source marker calibration must match the requested selected weight')
    proposal = provider.proposal_from_artifact(artifact, provenance, state=state,
                                              distance=F(distance), old_source=old_source, z=z, alpha=alpha)
    proposal.update(selected_weight=r, marker_bank_policy='preserve every source exact marker; no recalibration or thinning',
                    source_weight_calibration=calibrated,
                    marker_calibrations=artifact.get('marker_calibrations', []),
                    marker_count=len(proposal['witness']['markers']))
    return proposal


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def run_search(artifact, provenance, *, r, state=32, distance=F(9, 20), old_source=False,
               z=None, alpha=None, log_ratio=None, grid=None, max_evaluations=24,
               potential_evaluations=40, max_seconds=600., max_table_bytes=384 << 20,
               checkpoint=None):
    proposal = prepare_proposal(artifact, provenance, r=r, state=state, distance=distance,
                                old_source=old_source, z=z, alpha=alpha)
    if (type(max_evaluations) is not int or not 1 <= max_evaluations <= 1000
            or type(potential_evaluations) is not int or not 1 <= potential_evaluations <= 256
            or not isfinite(max_seconds) or not 0 < max_seconds <= 3600
            or type(max_table_bytes) is not int or not 1 <= max_table_bytes <= 512 << 20):
        raise ValueError('Invalid mark, potential, wall-time, or table-memory budget')
    w = proposal['witness']
    z0, alpha0 = float(F(w['z'])), float(F(w['alpha']))
    odds = provider.alias.checked_markers(tuple(tuple(provider.log_rational(x) for x in pair)
                                                for pair in w['markers']))
    ratio0 = provider.log_rational(proposal['state_potential_ratio']) if log_ratio is None else float(log_ratio)
    grid = w['grid'] if grid is None else grid
    if (not isfinite(ratio0) or abs(ratio0) > 4096 or type(grid) is not int or grid not in (64, 128, 256)
            or len(odds)*4*grid*grid*9 > max_table_bytes):
        raise ValueError('Invalid potential, grid, or marker-bank table memory requirement')
    if not 2.**-1022 <= z0 <= 2.**-1e-8 or not 2.**-32 <= alpha0 <= 1:
        raise ValueError('Initial marks must satisfy1e-8<=-log2(z)<=1022and2**-32<=alpha<=1')
    pins, started = source_pins(), monotonic()
    budget = provider.old.Budget(max_seconds)
    history, failures, cache = [], [], {}
    best, resource_limit, mark_evaluations, potential_calls = None, False, 0, 0
    outcome = None

    def save(row, bank, index, phase):
        nonlocal best
        witness = bank.witness()
        # Frozen evaluate returns the fitted ratio but does not update bank.log_ratio.
        ratio = row['log2_potential_ratio']
        witness.update(log2_potential_ratio=ratio, log2_potential_ratio_exact=str(F(ratio)))
        record = dict(row, evaluation_index=index, phase=phase, z=bank.z, alpha=bank.alpha,
                      witness=witness, selected_weight=r, field_order=provider.old.FIELD_ORDER,
                      target_distance=str(F(distance)), source=provenance, source_sha256=pins,
                      old_numerical_bounds_reused=False,
                      marker_bank_policy=proposal['marker_bank_policy'],
                      source_weight_calibration=proposal['source_weight_calibration'],
                      scope='selected-weight FLOAT proposal only; no directed or whole-code certificate')
        if best is None or record['log2_upper_screen'] < best['log2_upper_screen']:
            best = record
        history.append({key: record[key] for key in ('evaluation_index', 'phase', 'z', 'alpha',
                                                   'log2_potential_ratio', 'log2_upper_screen',
                                                   'potential_evaluations', 'regional_alias_diagnostics')})
        if source_pins() != pins:
            raise RuntimeError('A pinned search helper changed during execution')
        if checkpoint:
            checkpoint(record)
        return record

    def objective(point):
        nonlocal mark_evaluations, potential_calls
        budget.check()
        ell, exponent = map(float, point)
        if not 1e-8 <= ell <= 1022 or not -32 <= exponent <= 0:
            return 1e300
        mark, fractional = 2.**-ell, 2.**exponent
        key = mark, fractional
        if key in cache:
            return cache[key]
        index = mark_evaluations
        if index >= max_evaluations:
            raise TimeoutError('Mark evaluation budget exhausted')
        mark_evaluations += 1
        ratio_seed = ratio0 if best is None else best['log2_potential_ratio']
        try:
            coefficients = provider.AliasCoefficients(mark, fractional, state=state, _budget=budget)
            bank = provider.AliasBankScreen(mark, fractional, odds, ratio_seed,
                                            state=state, distance=F(distance), grid=grid,
                                            floor_policy='positive-fallback', max_table_bytes=max_table_bytes,
                                            _coefficients=coefficients, _budget=budget)
            row = bank.evaluate(r, method='runs', optimize_potential=False, _budget=budget)
            current = save(row, bank, index, 'baseline')
            cache[key] = current['log2_upper_screen']
            # Only the common potential changes; all coefficients and FFT tables are reused.
            fitted = bank.evaluate(r, method='runs', optimize_potential=True,
                                   max_evaluations=potential_evaluations, _budget=budget)
            potential_calls += fitted['potential_evaluations']
            candidate = save(fitted, bank, index, 'potential-fit')
            cache[key] = min(current['log2_upper_screen'], candidate['log2_upper_screen'])
            return cache[key]
        except (TimeoutError, ArithmeticError) as error:
            failures.append(dict(evaluation_index=index, z=mark, alpha=fractional,
                                 error_type=type(error).__name__, message=str(error),
                                 completed_baseline_retained=key in cache))
            if isinstance(error, TimeoutError):
                raise
            cache.setdefault(key, 1e300)
            return cache[key]

    initial = np.array((-log2(z0), log2(alpha0)))
    simplex = np.array((initial, (min(1022., initial[0]+max(1., initial[0]/10)), initial[1]),
                        (initial[0], initial[1]+1 if initial[1] <= -1 else initial[1]-1)))
    if simplex[1, 0] == initial[0]:
        simplex[1, 0] = initial[0]*.9
    try:
        fitted = minimize(objective, initial, method='Nelder-Mead',
                          options=dict(maxfev=max_evaluations, initial_simplex=simplex,
                                       xatol=1e-7, fatol=1e-5))
        outcome = dict(success=bool(fitted.success), message=str(fitted.message),
                       function_calls=int(fitted.nfev))
    except TimeoutError as error:
        resource_limit = True
        outcome = dict(success=False, message=str(error))
    if source_pins() != pins:
        raise RuntimeError('A pinned search helper changed during execution')
    return dict(schema='field128-full-rs-state-alias-mark-search-float-v1',
                construction=provider.local.construction(state), state_dimension=state, inner_width=64,
                selected_weight=r, target_distance=str(F(distance)), field_order=provider.old.FIELD_ORDER,
                field_bits=128, initial_proposal=proposal, best=best, proposal_history=history,
                failures=failures, completed_evaluation_phases=len(history), mark_evaluations=mark_evaluations,
                completed_potential_function_calls=potential_calls,
                max_evaluations=max_evaluations, potential_max_evaluations_per_mark=potential_evaluations,
                max_seconds=max_seconds, max_table_bytes=max_table_bytes, grid=grid,
                resource_limit_reached=resource_limit or not outcome['success'], optimizer=outcome,
                elapsed_seconds=monotonic()-started, source=provenance, source_sha256=pins,
                initial_z=z0, initial_alpha=alpha0, initial_log2_potential_ratio=ratio0,
                actual_marker_count=len(odds), marker_bank_policy=proposal['marker_bank_policy'],
                marker_bank_exact_source=proposal['witness']['markers'],
                floor_policy='positive-fallback', old_numerical_bounds_reused=False,
                directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                scope='only the reported selected weight; proposal search, not a distance certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--old-source', action='store_true')
    parser.add_argument('--state', type=int, choices=provider.local.STATES, default=32)
    parser.add_argument('--distance', default='9/20')
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--z')
    parser.add_argument('--alpha')
    parser.add_argument('--log-ratio', type=float)
    parser.add_argument('--grid', type=int)
    parser.add_argument('--max-evaluations', type=int, default=24)
    parser.add_argument('--potential-evaluations', type=int, default=40)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--max-table-mib', type=int, default=384)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    artifact, provenance = provider.cover.read_json(args.source)
    # Resolve the source calibration before creating any output.
    prepare_proposal(artifact, provenance, r=args.r, state=args.state, distance=F(args.distance),
                     old_source=args.old_source, z=args.z, alpha=args.alpha)
    args.output_dir.mkdir(exist_ok=False)

    def checkpoint(row):
        write_new(args.output_dir/f"evaluation_{row['evaluation_index']:03d}_{row['phase']}.json", row)
        print(json.dumps({key: row[key] for key in ('evaluation_index', 'phase', 'z', 'alpha',
                                                   'log2_potential_ratio', 'log2_upper_screen',
                                                   'regional_alias_diagnostics')}, allow_nan=False), flush=True)

    result = run_search(artifact, provenance, r=args.r, state=args.state, distance=F(args.distance),
                        old_source=args.old_source, z=args.z, alpha=args.alpha, log_ratio=args.log_ratio,
                        grid=args.grid, max_evaluations=args.max_evaluations,
                        potential_evaluations=args.potential_evaluations, max_seconds=args.max_seconds,
                        max_table_bytes=args.max_table_mib << 20, checkpoint=checkpoint)
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps(dict(output_dir=str(args.output_dir.resolve()), mark_evaluations=result['mark_evaluations'],
                          best_log2_upper_screen=None if result['best'] is None else result['best']['log2_upper_screen'],
                          resource_limit_reached=result['resource_limit_reached']), allow_nan=False))


if __name__ == '__main__':
    main()
