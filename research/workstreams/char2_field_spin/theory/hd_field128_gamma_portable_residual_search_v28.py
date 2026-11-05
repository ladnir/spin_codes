"""FLOAT third-witness search against two fixed complete profile bounds.

Only completed class terms are compared: count times chronological transfer
times that witness's outer prefix. Two baseline banks are rebuilt once.
Candidate markers stay fixed; the search varies only z and alpha. No
coefficient or proof bound is imported, and no directed claim is produced.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from math import isfinite, log2
from pathlib import Path
from time import monotonic

import hd_field128_gamma_portable_profile_min_v28 as profiles

portable, HERE = profiles.portable, Path(__file__).resolve().parent
NAME = Path(__file__).name
PROFILE_SHA256 = '029cecfc78f697b6b660b3d241165324ee5f09eae6de9e2512e0ee92e84679bf'
MAX_BANKS, SECONDS = 32, 90.


def candidate_witness(proposals):
    if len(proposals) != 2:
        raise ValueError('Exactly two fixed baseline proposals are required')
    first = proposals[0]['witness']
    union = list(dict.fromkeys(tuple(p) for item in proposals for p in item['witness']['markers']))
    use_union = len(union) <= 128
    markers = [list(p) for p in union] if use_union else deepcopy(first['markers'])
    z, alpha = F(first['z']), F(first['alpha'])
    if (F(float(z)) != z or F(float(alpha)) != alpha
            or not 2.**-1022 <= float(z) <= 2.**-1e-8 or not 2.**-32 <= float(alpha) <= 1):
        raise ValueError('Initial marks must be exact binary64 values in the frozen logarithmic chart')
    return dict(z=first['z'], alpha=first['alpha'], markers=markers, grid=128, state_dimension=30), dict(
        policy='deduplicated union' if use_union else 'first baseline; union exceeds 128',
        deduplicated_union_count=len(union), fixed_candidate_marker_count=len(markers))


def count_array(counts, full, budget, *, chunk=32):
    import numpy as np
    logs = np.empty((len(full), len(full)))
    for start in range(0, len(full), chunk):
        budget.check()
        logs[start:start+chunk] = counts.log_counts(full[start:start+chunk, None], full[None], budget=budget)
    if np.any(np.isnan(logs)) or np.any(np.isposinf(logs)):
        raise ArithmeticError('Invalid common count logarithms')
    return logs


def class_terms(regions, prefix, logs, transfer, budget, *, repetitions=256, chunk=32):
    """Complete ordered class terms, never cross-witness matrix minima."""
    import numpy as np
    regions, logs = np.asarray(regions), np.asarray(logs)
    n = len(regions)
    if (regions.shape != (n, 2, 2) or logs.shape != (n, n) or not n
            or not isfinite(prefix) or np.any(np.isnan(regions)) or np.any(np.isposinf(regions))
            or np.any(np.isnan(logs)) or np.any(np.isposinf(logs))
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2
            or type(chunk) is not int or not 1 <= chunk <= 128):
        raise ValueError('Invalid complete class domain or transfer parameters')
    result = np.empty_like(logs, dtype=float)
    for start in range(0, n, chunk):
        budget.check()
        product = transfer.multiply(regions[start:start+chunk, None], regions[None])
        repeated = transfer.power(product, repetitions//2)
        result[start:start+chunk] = logs[start:start+chunk]+np.logaddexp2(
            repeated[..., 0, 0], repeated[..., 0, 1])+prefix
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ArithmeticError('Invalid completed class term')
    return result


def log_union(terms):
    import numpy as np
    value = float(np.logaddexp2.reduce(terms.ravel()))
    if not isfinite(value):
        raise ArithmeticError('Zero or nonfinite FLOAT union is not evidence')
    return value


def residual_score(baseline, candidate, logs, full):
    import numpy as np
    if candidate.shape != baseline.shape or logs.shape != baseline.shape:
        raise ValueError('All witnesses must use the same complete profile partition')
    combined = np.minimum(baseline, candidate)
    position = np.unravel_index(np.argmax(combined), combined.shape)
    score, single = log_union(combined), log_union(candidate)
    return dict(combined_log2_upper_screen=score, candidate_single_log2_upper_screen=single,
        largest_profile=dict(even_full=int(full[position[0]]), odd_full=int(full[position[1]]),
            log2_term=float(combined[position]), log2_count_upper=float(logs[position]),
            winning_source='candidate' if candidate[position] < baseline[position] else 'baseline'),
        candidate_winning_profile_count=int(((candidate < baseline) & np.isfinite(logs)).sum()),
        profile_pairs_evaluated=combined.size, all_profiles_retained=True)


def search_candidates(initial, evaluate, budget, minimize, *, max_banks=MAX_BANKS):
    """Bound unique candidate banks; only fully evaluated rows enter history."""
    import numpy as np
    if type(max_banks) is not int or not 1 <= max_banks <= MAX_BANKS:
        raise ValueError('Expected 1..32 candidate banks')
    z0, alpha0 = float(F(initial['z'])), float(F(initial['alpha']))
    start = np.array((-log2(z0), log2(alpha0)))
    simplex = np.array((start, (min(1022., start[0]+.5), start[1]),
        (start[0], start[1]+.25 if start[1] <= -.25 else start[1]-.25)))
    if simplex[1, 0] == start[0]:
        simplex[1, 0] = start[0]-.5
    history, failures, cache, best, attempted = [], [], {}, None, 0

    def objective(point):
        nonlocal attempted, best
        budget.check()
        ell, exponent = map(float, point)
        if not 1e-8 <= ell <= 1022 or not -32 <= exponent <= 0:
            return 1e300
        z = z0 if ell == start[0] else 2.**-ell
        alpha = alpha0 if exponent == start[1] else 2.**exponent
        key = z, alpha
        if key in cache:
            return cache[key]
        if attempted >= max_banks:
            raise TimeoutError('Candidate-bank limit reached')
        number = attempted
        attempted += 1
        witness = deepcopy(initial)
        witness.update(z=str(F(z)), alpha=str(F(alpha)))
        try:
            row = evaluate(witness)
            score = row['combined_log2_upper_screen']
            if not isfinite(score):
                raise ArithmeticError('Nonfinite residual objective')
            compact = dict(evaluation_index=number, z=witness['z'], alpha=witness['alpha'],
                combined_log2_upper_screen=score,
                candidate_single_log2_upper_screen=row['candidate_single_log2_upper_screen'],
                largest_profile=row['largest_profile'], candidate_winning_profile_count=row['candidate_winning_profile_count'])
            history.append(compact)
            if best is None or score < best['combined_log2_upper_screen']:
                best = dict(row, evaluation_index=number, witness=witness)
            cache[key] = score
            return score
        except (TimeoutError, MemoryError, ArithmeticError) as error:
            failures.append(dict(evaluation_index=number, z=witness['z'], alpha=witness['alpha'],
                error_type=type(error).__name__, message=str(error), complete_profile_sum=False))
            if isinstance(error, (TimeoutError, MemoryError)):
                raise
            cache[key] = 1e300
            return cache[key]
    try:
        fitted = minimize(objective, start, method='Nelder-Mead', options=dict(maxfev=max_banks,
            initial_simplex=simplex, xatol=1e-7, fatol=1e-5))
        outcome = dict(success=bool(fitted.success), message=str(fitted.message), function_calls=int(fitted.nfev))
    except (TimeoutError, MemoryError) as error:
        outcome = dict(success=False, message=str(error), resource_stop=type(error).__name__)
    return dict(best=best, history=history, failures=failures, attempted_candidate_banks=attempted,
                complete_candidate_banks=len(history), optimizer=outcome)


def run(specifications, weight, output_dir):
    profiles.validate_target(weight)
    if len(specifications) != 2:
        raise ValueError('Supply exactly two baseline proposals')
    with portable.remote_limits():
        if sha256((HERE/profiles.NAME).read_bytes()).hexdigest() != PROFILE_SHA256:
            raise RuntimeError('Reviewed profile-min helper changed')
        loaded = profiles.helpers()
        validator, normalizer, neighbors = loaded
        from scipy.optimize import minimize
        import numpy as np
        proposals = profiles.load_proposals(specifications, weight)
        initial, marker_policy = candidate_witness(proposals)
        pins = profiles.source_pins(loaded)
        pins[profiles.NAME], pins[NAME] = PROFILE_SHA256, sha256((HERE/NAME).read_bytes()).hexdigest()
        pins = dict(sorted(pins.items()))
        def check():
            for proposal in proposals:
                portable.check_pins(pins, proposal['source'])
        check()
        portable.reserve_guard(HERE, portable.MAX_OUTPUT_BYTES)
        destination = HERE/validator.relative_name(output_dir)
        destination.mkdir(exist_ok=False)
        outputs = portable.Outputs(destination, normalizer.normalize_new_metadata)
        provider, mode = neighbors.provider, neighbors.mode
        transfer = provider.profiles.transfer
        budget, started = provider.old.Budget(SECONDS), monotonic()
        rows, singles, baseline, result, proposal_pin = [], [], None, None, None
        stage, status = 'common-counts', 'incomplete'
        try:
            outputs.write('driver.json', dict(schema='gamma-portable-residual-search-driver-v28',
                sources=[p['source'] for p in proposals], source_sha256=pins, active_groups=weight,
                maximum_candidate_banks=MAX_BANKS, baseline_banks=2, cooperative_seconds=SECONDS,
                hard_wall_seconds=portable.WALL_SECONDS, memory_bytes=portable.MEMORY_BYTES,
                minimum_free_bytes=portable.RESERVE, maximum_output_bytes=portable.MAX_OUTPUT_BYTES,
                maximum_simultaneous_banks=1, trial_files_written=False, marker_policy=marker_policy,
                directed=False, whole_code_certificate=False))
            try:
                full = transfer.tail.full_counts(profiles.PAIRS, weight)
                counts = mode.MidpointCounts(profiles.PAIRS, weight, max_singletons=64, budget=budget)
                logs = count_array(counts, full, budget)
                stage = 'baseline-banks'
                for proposal in proposals:
                    budget.check()
                    matrix, prefix, metadata = profiles.build_regions(proposal, weight, provider, budget)
                    terms = class_terms(matrix, prefix, logs, transfer, budget)
                    singles.append(log_union(terms)); rows.append(metadata)
                    if baseline is None:
                        baseline = terms
                    else:
                        np.minimum(baseline, terms, out=baseline)
                    del matrix, terms
                baseline_score = log_union(baseline)
                stage = 'candidate-search'
                def evaluate(witness):
                    proposal = dict(source=proposals[0]['source'], witness=witness, source_marker_key='markers')
                    matrix, prefix, metadata = profiles.build_regions(proposal, weight, provider, budget)
                    terms = class_terms(matrix, prefix, logs, transfer, budget)
                    row = residual_score(baseline, terms, logs, full)
                    if (row['profile_pairs_evaluated'] != (weight//2+1)**2
                            or row['combined_log2_upper_screen'] > baseline_score+1e-9):
                        raise ArithmeticError('Incomplete or non-dominating residual objective')
                    row.update(log2_outer_prefix=prefix, regional_logs_sha256=metadata['regional_logs_sha256'],
                        coefficient_table_sha256=metadata['coefficient_table_sha256'],
                        regional_alias_diagnostics=metadata['regional_alias_diagnostics'])
                    return row
                result = search_candidates(initial, evaluate, budget, minimize)
                if result['best'] is not None:
                    status, stage = 'complete-best', None
            except (TimeoutError, MemoryError) as error:
                status = 'timeout' if isinstance(error, TimeoutError) else 'memory-error'
            check(); outputs.check()
            if result is not None and result['best'] is not None:
                proposal_pin = outputs.write('best_proposal.json', dict(
                    schema='gamma-residual-selected-positive-literal-proposal-v28',
                    construction=profiles.CONSTRUCTION, active_groups=weight, state_dimension=30,
                    field_order=profiles.FIELD_ORDER, target_distance='2/5',
                    sources=[p['source'] for p in proposals], source_sha256={},
                    witness=result['best']['witness'], selected_evaluation_index=result['best']['evaluation_index'],
                    source_role='New exact candidate literals only; no bound or calibration history inherited',
                    directed=False, selected_weight_certificate=False, whole_code_certificate=False))
            summary = dict(schema='gamma-portable-residual-search-FLOAT-v28', status=status,
                active_groups=weight, construction=profiles.CONSTRUCTION, state_dimension=30,
                field_order=profiles.FIELD_ORDER, target_distance='2/5', source_sha256=pins,
                proposals=proposals, initial_candidate_witness=initial, marker_policy=marker_policy,
                completed_baseline_banks=rows, baseline_single_log2_upper_screens=singles,
                baseline_combined_log2_upper_screen=log_union(baseline) if len(rows) == 2 else None,
                search=result, best_proposal=proposal_pin, pending_stage=stage,
                complete_selected_weight_screen=status == 'complete-best',
                count_provenance=counts.diagnostics() if len(rows) == 2 else None,
                elapsed_seconds=monotonic()-started, old_numerical_bounds_reused=False,
                original_weight_used_for_all_histograms_and_prefixes=True,
                common_count_cap_outside_alpha=True, entrywise_cross_witness_minimum=False,
                historical_origin_verified=False, no_weight_interpolation=True,
                directed=False, selected_weight_certificate=False, whole_code_certificate=False,
                scope='Best completed three-witness FLOAT sum; not a convergence or directed certificate')
            summary_pin = outputs.write('summary.json', summary)
            check(); outputs.check()
            return dict(status=status, summary=summary_pin, best_proposal=proposal_pin,
                combined_log2_upper_screen=None if result is None or result['best'] is None else result['best']['combined_log2_upper_screen'],
                attempted_candidate_banks=None if result is None else result['attempted_candidate_banks'],
                whole_code_certificate=False, directed=False)
        finally:
            check(); outputs.check()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--proposal', nargs=2, action='append', metavar=('RELATIVE_JSON', 'SHA256'), required=True)
    parser.add_argument('--target-weight', type=int, required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    result = run(args.proposal, args.target_weight, args.output_dir)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 0 if result['status'] == 'complete-best' else 2


if __name__ == '__main__':
    raise SystemExit(main())
