"""Bounded FLOAT mark search for all-pair, chronological profile bounds.

The frozen marker bank is retained. Only z and alpha vary; no common state
potential is fitted or used in the objective. Each completed selected-weight
evaluation is checkpointed before the next bank is built. Saved numerical
bounds are discarded, and the source marks are preserved at the initial point.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_gamma_all_pair_screen as provider
import hd_field128_full_rs_gamma_alias_search as saved


def source_pins():
    pins = provider.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def run_search(artifact, provenance, *, r, state=24, distance=F(2, 5), z=None, alpha=None,
               marker_weight=None, grid=None, max_evaluations=24, max_seconds=600.,
               max_table_bytes=384 << 20, checkpoint=None):
    proposal = saved.prepare_proposal(artifact, provenance, r=r, state=state, distance=distance,
                                      z=z, alpha=alpha, marker_weight=marker_weight)
    if (type(max_evaluations) is not int or not 1 <= max_evaluations <= 1000
            or not isfinite(max_seconds) or not 0 < max_seconds <= 3600
            or type(max_table_bytes) is not int or not 1 <= max_table_bytes <= 512 << 20):
        raise ValueError('Invalid mark, wall-time, or table-memory budget')
    w = proposal['witness']
    z0, alpha0 = float(F(w['z'])), float(F(w['alpha']))
    odds = tuple(tuple(float(F(x)) for x in pair) for pair in w['marker_log_odds_exact'])
    ratio0 = float(F(w['log2_potential_ratio_exact']))
    grid = w['grid'] if grid is None else grid
    if (type(grid) is not int or grid not in (64, 128, 256)
            or len(odds)*4*grid*grid*9 > max_table_bytes):
        raise ValueError('Invalid grid or marker-bank table memory requirement')
    if not 2.**-1022 <= z0 <= 2.**-1e-8 or not 2.**-32 <= alpha0 <= 1:
        raise ValueError('Initial marks require 1e-8 <= -log2(z) <= 1022 and 2**-32 <= alpha <= 1')
    pins, started = source_pins(), monotonic()
    budget = provider.old.Budget(max_seconds)
    history, failures, cache = [], [], {}
    best, resource_limit, mark_evaluations = None, False, 0

    def save(row, bank, index):
        nonlocal best
        if not isfinite(row['log2_upper_screen']):
            raise ArithmeticError('The refined selected-weight objective must be finite')
        witness = bank.witness()
        witness.update(log2_potential_ratio=ratio0, log2_potential_ratio_exact=str(F(ratio0)))
        record = dict(row, evaluation_index=index, phase='exact-profiles', z=bank.z, alpha=bank.alpha,
                      witness=witness, selected_weight=r, field_order=provider.FIELD_ORDER,
                      target_distance=str(F(distance)), source=provenance, source_sha256=pins,
                      old_numerical_bounds_reused=False, potential_evaluations=0,
                      common_potential_used=False,
                      potential_role='source witness retained for compatibility; unused by the objective',
                      marker_bank_policy=proposal['marker_bank_policy'],
                      marker_representation=proposal['marker_representation'],
                      marker_literals=proposal['marker_literals'],
                      source_weight_calibration=proposal['source_weight_calibration'],
                      marker_calibrations=proposal['marker_calibrations'],
                      scope='one selected weight; all-pair exact-profile FLOAT proposal; no directed or whole-code certificate')
        if source_pins() != pins:
            raise RuntimeError('A pinned refined-search helper changed during execution')
        if best is None or record['log2_upper_screen'] < best['log2_upper_screen']:
            best = record
        history.append({key: record[key] for key in ('evaluation_index', 'phase', 'z', 'alpha',
                                                   'log2_upper_screen', 'largest_profile',
                                                   'regional_alias_diagnostics')})
        if checkpoint:
            checkpoint(record)
        return record

    def objective(point):
        nonlocal mark_evaluations
        budget.check()
        ell, exponent = map(float, point)
        if not 1e-8 <= ell <= 1022 or not -32 <= exponent <= 0:
            return 1e300
        # Do not perturb either original binary64 mark by a log/exp round trip.
        mark = z0 if ell == initial[0] else 2.**-ell
        fractional = alpha0 if exponent == initial[1] else 2.**exponent
        key = mark, fractional
        if key in cache:
            return cache[key]
        index = mark_evaluations
        if index >= max_evaluations:
            raise TimeoutError('Mark evaluation budget exhausted')
        mark_evaluations += 1
        try:
            coefficients = provider.AliasCoefficients(mark, fractional, state=state, _budget=budget)
            bank = provider.AliasBankScreen(mark, fractional, odds, ratio0,
                                            state=state, distance=F(distance), grid=grid,
                                            floor_policy='positive-fallback', max_table_bytes=max_table_bytes,
                                            _coefficients=coefficients, _budget=budget)
            row = provider.selected_screen(bank, r, support_method='exact-profiles',
                                            optimize_potential=False, _budget=budget)
            current = save(row, bank, index)
            cache[key] = current['log2_upper_screen']
            return cache[key]
        except (TimeoutError, ArithmeticError) as error:
            failures.append(dict(evaluation_index=index, z=mark, alpha=fractional,
                                 error_type=type(error).__name__, message=str(error)))
            if isinstance(error, TimeoutError):
                raise
            cache[key] = 1e300
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
        raise RuntimeError('A pinned refined-search helper changed during execution')
    return dict(schema='field128-gamma-all-pair-exact-profile-mark-search-float-v1',
                construction=provider.gamma.construction(state), state_dimension=state, inner_width=64,
                gamma_extension_degree=16, feedback='gamma A transpose',
                selected_weight=r, target_distance=str(F(distance)), field_order=provider.FIELD_ORDER,
                field_bits=128, initial_proposal=proposal, best=best, proposal_history=history,
                failures=failures, mark_evaluations=mark_evaluations, completed_mark_evaluations=len(history),
                max_evaluations=max_evaluations, max_seconds=max_seconds, max_table_bytes=max_table_bytes,
                grid=grid, resource_limit_reached=resource_limit or not outcome['success'], optimizer=outcome,
                elapsed_seconds=monotonic()-started, source=provenance, source_sha256=pins,
                initial_z=z0, initial_alpha=alpha0, retained_log2_potential_ratio=ratio0,
                common_potential_used=False, completed_potential_function_calls=0,
                support_method='exact-profiles', objective='all-pair caps; sharpened profile upper counts; chronological transfer',
                actual_marker_count=len(odds), marker_bank_policy=proposal['marker_bank_policy'],
                marker_representation=proposal['marker_representation'],
                marker_bank_exact_source=proposal['marker_literals'],
                source_weight_calibration=proposal['source_weight_calibration'],
                marker_calibrations=proposal['marker_calibrations'],
                floor_policy='positive-fallback', old_numerical_bounds_reused=False,
                directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                scope='only the reported selected weight; refined-objective proposal search, not a distance certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--state', type=int, choices=provider.gamma.STATES, default=24)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--marker-weight', type=int, help='Explicit source calibration if absent; cannot override a mismatch')
    parser.add_argument('--z')
    parser.add_argument('--alpha')
    parser.add_argument('--grid', type=int)
    parser.add_argument('--max-evaluations', type=int, default=24)
    parser.add_argument('--max-seconds', type=float, default=600.)
    parser.add_argument('--max-table-mib', type=int, default=384)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error('Output directory must be unused')
    raw = args.source.read_bytes()
    artifact = json.loads(raw)
    provenance = dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest())
    proposal = saved.prepare_proposal(artifact, provenance, r=args.r, state=args.state, distance=F(args.distance),
                                       z=args.z, alpha=args.alpha, marker_weight=args.marker_weight)
    args.output_dir.mkdir(exist_ok=False)
    saved.write_new(args.output_dir/'proposal.json', proposal)

    def checkpoint(row):
        if args.source.read_bytes() != raw:
            raise RuntimeError('The pinned witness source changed during execution')
        saved.write_new(args.output_dir/f"evaluation_{row['evaluation_index']:03d}_exact_profiles.json", row)
        print(json.dumps({key: row[key] for key in ('evaluation_index', 'z', 'alpha', 'log2_upper_screen',
                                                   'largest_profile', 'regional_alias_diagnostics')}, allow_nan=False), flush=True)

    result = run_search(artifact, provenance, r=args.r, state=args.state, distance=F(args.distance),
                        z=args.z, alpha=args.alpha, marker_weight=args.marker_weight,
                        grid=args.grid, max_evaluations=args.max_evaluations, max_seconds=args.max_seconds,
                        max_table_bytes=args.max_table_mib << 20, checkpoint=checkpoint)
    if args.source.read_bytes() != raw:
        raise RuntimeError('The pinned witness source changed during execution')
    saved.write_new(args.output_dir/'summary.json', result)
    print(json.dumps(dict(output_dir=str(args.output_dir.resolve()), mark_evaluations=result['mark_evaluations'],
                          completed_mark_evaluations=result['completed_mark_evaluations'],
                          best_log2_upper_screen=None if result['best'] is None else result['best']['log2_upper_screen'],
                          resource_limit_reached=result['resource_limit_reached']), allow_nan=False))


if __name__ == '__main__':
    main()
