"""Bounded FLOAT mark/potential search for the frozen gamma-feedback bank.

Saved FLOAT rows and search summaries supply witness literals only. Their
selected-weight calibration must match the requested weight. Each mark pair
builds one coefficient/FFT bank; potential fitting reuses that bank. Completed
baseline and fitted phases are checkpointed before further numerical work.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import numpy as np
from scipy.optimize import minimize

import hd_field128_full_rs_gamma_screen as provider


def source_pins():
    pins = provider.source_pins()
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, F)):
        raise ValueError('An exact rational string or integer is required')
    return F(value)


def log_positive(value):
    value = exact(value)
    if value <= 0:
        raise ValueError('A positive exact marker or potential is required')
    return log2(value.numerator)-log2(value.denominator)


def prepare_proposal(artifact, provenance, *, r, state, distance, z=None, alpha=None,
                     marker_weight=None):
    name, distance = provider.construction(state), F(distance)
    if type(r) is not int or not 1 <= r <= provider.old.GROUPS or not 0 < distance < F(1, 2):
        raise ValueError('Expected weight in 1..8192 and distance strictly between zero and one half')
    if not isinstance(artifact, dict):
        raise ValueError('Expected a saved FLOAT row or search summary')
    row = artifact.get('best', artifact)
    if not isinstance(row, dict) or not isinstance(row.get('witness'), dict):
        raise ValueError('The saved row or summary must contain a completed witness')
    witness = row['witness']
    calibration = {}
    for label, container in (('artifact', artifact), ('row', row), ('witness', witness)):
        for key in ('active_groups', 'selected_weight', 'center'):
            if key in container:
                calibration[f'{label}.{key}'] = container[key]
        if 'field_order' in container and container['field_order'] != provider.FIELD_ORDER:
            raise ValueError('Source witness must use GF(2^128)')
    if marker_weight is not None:
        calibration['explicit_marker_weight'] = marker_weight
    if not calibration or any(type(value) is not int or value != r for value in calibration.values()):
        raise ValueError('Explicit source marker calibration must match the selected weight')

    # Prefer the exact saved logarithms of the actual bank over any older proposal.
    if 'marker_log_odds_exact' in witness:
        representation, literals = 'log2-exact', witness['marker_log_odds_exact']
        convert = lambda value: float(exact(value))
    elif 'markers' in witness:
        representation, literals = 'positive-exact', witness['markers']
        convert = log_positive
    else:
        raise ValueError('Source witness needs exact marker_log_odds_exact or exact positive markers')
    if (not isinstance(literals, (list, tuple)) or not 1 <= len(literals) <= 128
            or any(not isinstance(pair, (list, tuple)) or len(pair) != 2 for pair in literals)):
        raise ValueError('Expected 1..128 exact marker pairs')
    marker_literals = [[str(value) for value in pair] for pair in literals]
    raw_odds = tuple(tuple(convert(value) for value in pair) for pair in literals)
    odds = provider.frozen_alias.alias.checked_markers(raw_odds)
    if len(odds) != len(raw_odds):
        raise ValueError('Source markers must remain distinct in the frozen FLOAT bank')

    # Older evaluate() calls may leave a stale potential in bank.witness().
    # The selected row's fitted potential therefore takes precedence.
    ratio = None
    for label, container in (('row', row), ('witness', witness)):
        if 'log2_potential_ratio_exact' in container:
            ratio = float(exact(container['log2_potential_ratio_exact']))
            ratio_source = f'{label}.log2_potential_ratio_exact'
            break
        if 'log2_potential_ratio' in container:
            ratio = float(container['log2_potential_ratio'])
            ratio_source = f'{label}.log2_potential_ratio'
            break
    if ratio is None:
        for label, container in (('row', row), ('artifact', artifact)):
            if 'state_potential_ratio' in container:
                ratio = log_positive(container['state_potential_ratio'])
                ratio_source = f'{label}.state_potential_ratio'
                break
    if ratio is None or not isfinite(ratio) or abs(ratio) > 4096:
        raise ValueError('Source witness needs a finite common potential with log magnitude <=4096')
    mark, exponent = exact(witness['z'] if z is None else z), exact(witness['alpha'] if alpha is None else alpha)
    if not 0 < mark < 1 or not 0 < exponent <= 1:
        raise ValueError('Expected admissible exact z and alpha')
    return dict(schema='field128-full-rs-gamma-alias-witness-only-v1', construction=name,
                state_dimension=state, inner_width=64, gamma_extension_degree=16,
                field_order=provider.FIELD_ORDER, field_bits=128, selected_weight=r,
                target_distance=str(distance), source=provenance,
                source_construction=row.get('construction', artifact.get('construction')),
                source_role='exact witness literals only; source bounds are discarded',
                source_weight_calibration=calibration,
                marker_calibrations=row.get('marker_calibrations', artifact.get('marker_calibrations', [])),
                marker_bank_policy='preserve every exact source marker; no recalibration or thinning',
                marker_representation=representation, marker_literals=marker_literals, marker_count=len(odds),
                potential_source=ratio_source, old_numerical_bounds_reused=False,
                directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                witness=dict(z=str(mark), alpha=str(exponent), grid=witness['grid'],
                             marker_log_odds_exact=[[str(F(x)) for x in pair] for pair in odds],
                             log2_potential_ratio_exact=str(F(ratio))))


def write_new(path, value):
    raw = (json.dumps(value, indent=2, allow_nan=False)+'\n').encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(raw)
    return dict(path=str(Path(path).resolve()), sha256=sha256(raw).hexdigest())


def run_search(artifact, provenance, *, r, state=24, distance=F(2, 5), z=None, alpha=None,
               marker_weight=None, log_ratio=None, grid=None, max_evaluations=24,
               potential_evaluations=40, max_seconds=600., max_table_bytes=384 << 20,
               checkpoint=None):
    proposal = prepare_proposal(artifact, provenance, r=r, state=state, distance=distance,
                                z=z, alpha=alpha, marker_weight=marker_weight)
    if (type(max_evaluations) is not int or not 1 <= max_evaluations <= 1000
            or type(potential_evaluations) is not int or not 1 <= potential_evaluations <= 256
            or not isfinite(max_seconds) or not 0 < max_seconds <= 3600
            or type(max_table_bytes) is not int or not 1 <= max_table_bytes <= 512 << 20):
        raise ValueError('Invalid mark, potential, wall-time, or table-memory budget')
    w = proposal['witness']
    z0, alpha0 = float(F(w['z'])), float(F(w['alpha']))
    odds = tuple(tuple(float(F(x)) for x in pair) for pair in w['marker_log_odds_exact'])
    ratio0 = float(F(w['log2_potential_ratio_exact'])) if log_ratio is None else float(log_ratio)
    grid = w['grid'] if grid is None else grid
    if (not isfinite(ratio0) or abs(ratio0) > 4096 or type(grid) is not int or grid not in (64, 128, 256)
            or len(odds)*4*grid*grid*9 > max_table_bytes):
        raise ValueError('Invalid potential, grid, or marker-bank table memory requirement')
    if not 2.**-1022 <= z0 <= 2.**-1e-8 or not 2.**-32 <= alpha0 <= 1:
        raise ValueError('Initial marks require 1e-8 <= -log2(z) <= 1022 and 2**-32 <= alpha <= 1')
    pins, started = source_pins(), monotonic()
    budget = provider.old.Budget(max_seconds)
    history, failures, cache = [], [], {}
    best, resource_limit, mark_evaluations, potential_calls = None, False, 0, 0

    def save(row, bank, index, phase):
        nonlocal best
        witness = bank.witness()
        ratio = row['log2_potential_ratio']
        witness.update(log2_potential_ratio=ratio, log2_potential_ratio_exact=str(F(ratio)))
        record = dict(row, evaluation_index=index, phase=phase, z=bank.z, alpha=bank.alpha,
                      witness=witness, selected_weight=r, field_order=provider.FIELD_ORDER,
                      target_distance=str(F(distance)), source=provenance, source_sha256=pins,
                      old_numerical_bounds_reused=False,
                      marker_bank_policy=proposal['marker_bank_policy'],
                      marker_representation=proposal['marker_representation'],
                      marker_literals=proposal['marker_literals'],
                      source_weight_calibration=proposal['source_weight_calibration'],
                      marker_calibrations=proposal['marker_calibrations'],
                      scope='gamma selected-weight FLOAT proposal only; no directed or whole-code certificate')
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
        # Preserve the source binary64 marks exactly at the initial coordinate.
        mark = z0 if ell == initial[0] else 2.**-ell
        fractional = alpha0 if exponent == initial[1] else 2.**exponent
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
            # Only the common potential changes; reuse coefficients and FFT tables.
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
    return dict(schema='field128-full-rs-gamma-alias-mark-search-float-v1',
                construction=provider.construction(state), state_dimension=state, inner_width=64,
                gamma_extension_degree=16, feedback='gamma A transpose',
                selected_weight=r, target_distance=str(F(distance)), field_order=provider.FIELD_ORDER,
                field_bits=128, initial_proposal=proposal, best=best, proposal_history=history,
                failures=failures, completed_evaluation_phases=len(history), mark_evaluations=mark_evaluations,
                completed_potential_function_calls=potential_calls,
                max_evaluations=max_evaluations, potential_max_evaluations_per_mark=potential_evaluations,
                max_seconds=max_seconds, max_table_bytes=max_table_bytes, grid=grid,
                resource_limit_reached=resource_limit or not outcome['success'], optimizer=outcome,
                elapsed_seconds=monotonic()-started, source=provenance, source_sha256=pins,
                initial_z=z0, initial_alpha=alpha0, initial_log2_potential_ratio=ratio0,
                actual_marker_count=len(odds), marker_bank_policy=proposal['marker_bank_policy'],
                marker_representation=proposal['marker_representation'],
                marker_bank_exact_source=proposal['marker_literals'],
                source_weight_calibration=proposal['source_weight_calibration'],
                marker_calibrations=proposal['marker_calibrations'],
                floor_policy='positive-fallback', old_numerical_bounds_reused=False,
                directed=False, whole_code_certificate=False, arithmetic='FLOAT',
                scope='only the reported selected weight; gamma proposal search, not a distance certificate')


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--state', type=int, choices=provider.STATES, default=24)
    parser.add_argument('--distance', default='2/5')
    parser.add_argument('--r', type=int, required=True)
    parser.add_argument('--marker-weight', type=int, help='Explicit source calibration if absent; cannot override a mismatch')
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
    raw = args.source.read_bytes()
    artifact = json.loads(raw)
    provenance = dict(path=str(args.source.resolve()), sha256=sha256(raw).hexdigest())
    prepare_proposal(artifact, provenance, r=args.r, state=args.state, distance=F(args.distance),
                     z=args.z, alpha=args.alpha, marker_weight=args.marker_weight)
    args.output_dir.mkdir(exist_ok=False)

    def checkpoint(row):
        if args.source.read_bytes() != raw:
            raise RuntimeError('The pinned witness source changed during execution')
        write_new(args.output_dir/f"evaluation_{row['evaluation_index']:03d}_{row['phase']}.json", row)
        print(json.dumps({key: row[key] for key in ('evaluation_index', 'phase', 'z', 'alpha',
                                                   'log2_potential_ratio', 'log2_upper_screen',
                                                   'regional_alias_diagnostics')}, allow_nan=False), flush=True)

    result = run_search(artifact, provenance, r=args.r, state=args.state, distance=F(args.distance),
                        z=args.z, alpha=args.alpha, marker_weight=args.marker_weight,
                        log_ratio=args.log_ratio, grid=args.grid, max_evaluations=args.max_evaluations,
                        potential_evaluations=args.potential_evaluations, max_seconds=args.max_seconds,
                        max_table_bytes=args.max_table_mib << 20, checkpoint=checkpoint)
    if args.source.read_bytes() != raw:
        raise RuntimeError('The pinned witness source changed during execution')
    write_new(args.output_dir/'summary.json', result)
    print(json.dumps(dict(output_dir=str(args.output_dir.resolve()), mark_evaluations=result['mark_evaluations'],
                          best_log2_upper_screen=None if result['best'] is None else result['best']['log2_upper_screen'],
                          resource_limit_reached=result['resource_limit_reached']), allow_nan=False))


if __name__ == '__main__':
    main()
