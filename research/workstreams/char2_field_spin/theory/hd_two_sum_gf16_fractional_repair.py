"""Finite saved-FLOAT-bank repair of the fractional live-state envelope.

No regional provider is called. For each of25 frozen initial witnesses,
only the live->live slope is varied over the existing at-most16 secants.
The best complete-band repair then receives at most one four-coordinate
round. All marks, exponents and type-specific floor markers stay frozen.
Saved cap values remain FLOAT proposals, never certificate inputs.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_fractional_state as source
import hd_two_sum_gf16_state_envelope as state


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = DIRECTORY/'hd_two_sum_gf16_fractional_state_r18_32_float.json'
ARTIFACT_SHA256 = 'd1f4a8aa1d57b8c63f8bcd27a9af9dedac5274fc0b9496d3f5e35558cf07d544'
SOURCE_SHA256 = 'bbeca202d77780ea325bbc70e959e73e929beec7c0341c2d80d85023c96273d8'
STATE_SHA256 = '5107b945ab18415191c48a3434a354c101d8fc5aeda66ec83a3fd2a5c3498af7'


def read_pinned(path, expected):
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != expected:
        raise ValueError(f'Frozen source/artifact hash mismatch: {Path(path).name}')
    return raw


def load_frozen():
    read_pinned(source.__file__, SOURCE_SHA256)
    read_pinned(state.__file__, STATE_SHA256)
    data = json.loads(read_pinned(ARTIFACT, ARTIFACT_SHA256))
    if (data.get('source_sha256') != SOURCE_SHA256 or data.get('construction') != source.CONSTRUCTION
            or data.get('active_group_band') != [18, 32] or data.get('directed') is not False
            or data.get('complete_requested_band') is not True or data.get('raw_regional_banks') != 5
            or data.get('all_cohort_splits') is not True or data.get('all_constituent_row_sets') is not True
            or data.get('all_source_supports') is not True
            or data.get('packet_power_before_geometry_average') is not True
            or data.get('route_normalizer_powered') is not False
            or len(data.get('banks', [])) != 5 or len(data.get('witnesses', [])) != 26):
        raise ValueError('Pinned completed FLOAT all-cohort source scope required')
    banks = []
    for index, (bank, (mark, alpha)) in enumerate(zip(data['banks'], source.CONTROLS)):
        if F(bank['mark']) != mark or F(bank['alpha']) != alpha or bank.get('zero_entries_are_none') is not True:
            raise ValueError('Frozen fixed-control identity mismatch')
        curves = np.asarray([[(-np.inf if value is None else value) for value in row]
                             for row in bank['curves']], dtype=float)
        if curves.shape != (4, 33) or np.any(np.isnan(curves)) or np.any(np.isposinf(curves)):
            raise ValueError('Complete finite-or-zero FLOAT cap bank required')
        banks.append(curves)
        for offset, witness in enumerate(data['witnesses'][5*index:5*index+5]):
            if (witness.get('bank_index') != index or F(witness['mark']) != mark
                    or F(witness['alpha']) != alpha or witness.get('coordinate_rounds') != 0
                    or witness.get('floor_template') != data['floor_template_names'][offset]
                    or len(witness.get('log2_slopes', [])) != 4
                    or not np.all(np.isfinite(witness['log2_slopes']))):
                raise ValueError('Exactly the25 unrefined frozen witnesses are required')
            source.floor_array(list(map(F, witness['log2_floor_markers'])))
    return data, banks


def repair_live_slope(curves, weights, baseline, check_time):
    """Re-evaluate one baseline and vary ONLY coordinate3 over its secants."""
    _, candidates = state.slope_candidates(curves)
    if len(candidates[3]) > 16 or 0. not in candidates[3]:
        raise ArithmeticError('Bounded live-state candidate bank must include zero')
    initial = np.asarray(baseline['log2_slopes'])
    floor_logs = list(map(F, baseline['log2_floor_markers']))
    memo = {}
    def objective(slopes):
        check_time()
        key = tuple(float(value) for value in slopes)
        if key not in memo:
            memo[key] = source.evaluate(curves, weights, F(baseline['mark']),
                F(baseline['alpha']), floor_logs, np.asarray(key))
        return memo[key]
    original = objective(initial)
    if abs(original['failure_log2_screen']-baseline['failure_log2_screen']) > 1e-7:
        raise ArithmeticError('Fresh baseline evaluation differs from the frozen FLOAT artifact')
    best = original
    for candidate in candidates[3]:
        slopes = initial.copy()
        slopes[3] = candidate
        row = objective(slopes)
        if row['failure_log2_screen'] < best['failure_log2_screen']:
            best = row
    if best['log2_slopes'][:3] != initial[:3].tolist():
        raise ArithmeticError('Live-state repair modified another coordinate')
    return dict(best, floor_template=baseline['floor_template'],
        floor_anchor=baseline['floor_anchor'], floor_scale=baseline['floor_scale'],
        bank_index=baseline['bank_index'], repaired_coordinate=3,
        live_slope_candidates=list(candidates[3]), live_slope_evaluations=len(memo),
        original_failure_log2_screen=original['failure_log2_screen'],
        original_log2_slopes=initial.tolist(),
        live_slope_trials=[dict(log2_slopes=row['log2_slopes'],
            failure_log2_screen=row['failure_log2_screen']) for row in memo.values()])


def run(*, max_seconds=30, refine=True, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30 or type(refine) is not bool):
        raise ValueError('Explicit at-most30s budget and Boolean refinement flag required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Saved-bank live-slope repair exhausted its total budget')
    result = dict(source_artifact=ARTIFACT.name, source_artifact_sha256=ARTIFACT_SHA256,
        frozen_provider_sha256=SOURCE_SHA256, frozen_state_helper_sha256=STATE_SHA256,
        construction=source.CONSTRUCTION, field_order=256, forbidden_weight=source.LIMIT,
        active_group_band=[18, 32], all_cohort_splits=True, all_constituent_row_sets=True,
        all_source_supports=True, pair_swaps=True, source_counts='C(8192,r)',
        packet_power_before_geometry_average=True, route_normalizer_powered=False,
        saved_caps_are_float_proposals=True, regional_banks_rebuilt=0,
        dense_weights_recomputed_from_alpha=True, frozen_marks_alpha_floors=True,
        repaired_baseline_count=0, initial_repaired_coordinate=3,
        maximum_live_candidates_per_baseline=16, maximum_extra_coordinate_rounds=int(refine),
        maximum_extra_refined_witnesses=int(refine), mark_alpha_floor_fits=0,
        directed=False, proposal_requires_fresh_directed_replay=True,
        whole_code_certificate=False, complete_requested_band=False,
        compute_budget_seconds=max_seconds, witnesses=[])
    try:
        check_time()
        data, banks = load_frozen()
        weights = [source.fractional_dense_weights(alpha) for _, alpha in source.CONTROLS]
        check_time()
        for index, baseline in enumerate(data['witnesses'][:25]):
            bank_index = baseline['bank_index']
            repaired = repair_live_slope(banks[bank_index], weights[bank_index], baseline, check_time)
            repaired['frozen_baseline_index'] = index
            result['witnesses'].append(repaired)
            result['repaired_baseline_count'] += 1
            if progress:
                progress(dict(stage='repair-live-slope', completed=index+1, total=25,
                    mark=repaired['mark'], alpha=repaired['alpha'], floor_template=repaired['floor_template'],
                    failure_log2_screen=repaired['failure_log2_screen'], evaluations=repaired['live_slope_evaluations'],
                    elapsed_seconds=perf_counter()-start))
        check_time()
        best = min(range(25), key=lambda i: result['witnesses'][i]['failure_log2_screen'])
        result['best_repaired_baseline_index'] = best
        if refine:
            baseline = result['witnesses'][best]
            bank_index = baseline['bank_index']
            refined = source.refine_one_round(banks[bank_index], weights[bank_index], baseline, check_time)
            refined.update(bank_index=bank_index, repaired_baseline_index=best,
                           frozen_baseline_index=baseline['frozen_baseline_index'])
            result['witnesses'].append(refined)
            if progress:
                progress(dict(stage='one-global-coordinate-round', baseline_index=best,
                    failure_log2_screen=refined['failure_log2_screen'], evaluations=refined['slope_evaluations'],
                    elapsed_seconds=perf_counter()-start))
        check_time()
        selected = []
        for offset, weight in enumerate(range(18, 33)):
            index = min(range(len(result['witnesses'])),
                        key=lambda i: result['witnesses'][i]['entries'][offset]['failure_log2_screen'])
            selected.append(dict(result['witnesses'][index]['entries'][offset], witness_index=index))
        result.update(status='completed', complete_requested_band=True, entries=selected,
            failure_log2_screen=float(state.logsum([row['failure_log2_screen'] for row in selected])),
            witness_selection_per_whole_weight=True, extra_witness_factor=False)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-start
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--no-refine', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen artifacts must not be overwritten')
    result = run(max_seconds=args.max_seconds, refine=not args.no_refine,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('witnesses', 'entries')}, allow_nan=False))
