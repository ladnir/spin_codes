"""Exact frozen-witness A/B for paired-NZ versus packetwise low-arity caps.

Rebuild only three paired-NZ regional banks. Preserve all45 original z,
alpha, type-floor, total-arity-marker and numeric-slope witnesses. Recompute
intercepts over every feasible occupancy. No fitting or refinement occurs.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_low_arity as low
import hd_two_sum_gf16_paired_nz as paired


DIRECTORY = Path(__file__).resolve().parent
REFERENCE = DIRECTORY/'hd_two_sum_gf16_low_arity_anchors_float.json'
REFERENCE_SHA256 = '206f34cca8021de88ce8806a05c6736523261bf769ff76f147542e7cb6fdc5d5'
LOW_SHA256 = '943dd058318f0f4c5d5a371fa3a84210646e71ca8c01327399d9a60ce7b9064a'
PAIRED_SHA256 = 'e6938a76e750924d03329e091a21582d5686b606e41517d43da6118c3995e048'


def pinned(path, digest):
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != digest:
        raise ValueError(f'Frozen A/B source mismatch: {Path(path).name}')
    return raw


def load_reference():
    pinned(low.__file__, LOW_SHA256)
    pinned(paired.__file__, PAIRED_SHA256)
    low.checked_dependencies()
    data = json.loads(pinned(REFERENCE, REFERENCE_SHA256))
    if (data.get('source_sha256') != LOW_SHA256 or data.get('directed') is not False
            or data.get('complete_requested_anchors') is not True or len(data.get('witnesses', [])) != 45
            or data.get('maximum_extra_refined_anchors') != 0 or data.get('regional_banks_rebuilt') != 3):
        raise ValueError('The frozen completed unrefined low-side comparison is required')
    for index, row in enumerate(data['witnesses']):
        anchor = low.ANCHORS[index//15]
        low.checked_low_scope(row['active_groups'], row['high_arity_cutoff'], F(row['mark']), F(row['alpha']))
        if ((row['active_groups'], row['high_arity_cutoff'], F(row['mark'])) != anchor
                or row['log2_total_arity_marker'] != low.LOG_T[index % 5]
                or row['floor_template_index'] != (index % 15)//5):
            raise ValueError('Frozen45-control order and exact scope must agree')
    return data


def run(*, max_seconds=30, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30):
        raise ValueError('Explicit at-most30s budget required')
    start = perf_counter()
    def check_time():
        if perf_counter()-start >= max_seconds:
            raise TimeoutError('Paired-NZ fixed-witness A/B exhausted its budget')
    result = dict(construction=low.fractional.CONSTRUCTION, field_order=256,
        forbidden_weight=low.fractional.LIMIT, anchors=[dict(active_groups=r, high_arity_cutoff=a0,
            low_total_active_constituents_max=a0-1, mark=str(z), alpha=str(low.ALPHA)) for r, a0, z in low.ANCHORS],
        reference_artifact=REFERENCE.name, reference_artifact_sha256=REFERENCE_SHA256,
        frozen_low_helper_sha256=LOW_SHA256, paired_provider_sha256=PAIRED_SHA256,
        provider=paired.provider_metadata(), all_source_supports=True,
        all_cohorts_within_low_total_arity=True, all_constituent_row_sets_within_low_total_arity=True,
        source_counts='C(8192,r)', total_arity_marker_outside_alpha=True,
        identical_numeric_slopes_to_reference=True, identical_marks_alpha_type_floors_and_t=True,
        fitted_parameters=0, control_count=45, maximum_regional_occupancy=32,
        regional_banks_rebuilt=0, continuous_weight_band_covered=False,
        complete_requested_anchors=False, high_side_included=False, directed=False,
        requires_fresh_directed_replay=True, whole_code_certificate=False,
        compute_budget_seconds=max_seconds, banks=[], witnesses=[])
    best_indices = []
    try:
        check_time()
        reference = load_reference()
        weights = low.fractional.fractional_dense_weights(low.ALPHA)
        for anchor_index, (r, a0, mark) in enumerate(low.ANCHORS):
            check_time()
            curves_all = low.state.regional_curves(paired.paired_nz_regions(mark, low.ALPHA, 32), 32)
            old_curves = np.array([[(-np.inf if x is None else x) for x in row]
                                   for row in reference['banks'][anchor_index]['curves']])
            if not np.all(curves_all <= old_curves+1e-8):
                raise ArithmeticError('Paired cap exceeds the frozen packetwise FLOAT cap')
            curves = curves_all[:, :r+1]
            initial, candidates, points = low.peak_tail_slopes(curves)
            result['banks'].append(dict(active_groups=r, mark=str(mark), alpha=str(low.ALPHA),
                curves=[[None if not np.isfinite(x) else float(x) for x in row] for row in curves_all],
                zero_entries_are_none=True, envelope_maximum=r,
                diagnostic_peak_tail_log2_slopes=initial.tolist(),
                diagnostic_slope_candidates=[list(row) for row in candidates], diagnostic_slope_points=points,
                diagnostic_slopes_not_used_in_control=True, no_worse_than_packetwise_float=True))
            result['regional_banks_rebuilt'] += 1
            cache, indices = {}, []
            for original_index in range(15*anchor_index, 15*(anchor_index+1)):
                check_time()
                original = reference['witnesses'][original_index]
                floors = list(map(F, original['log2_floor_markers']))
                slopes = np.array(original['log2_slopes'])
                key = tuple(floors), tuple(slopes)
                if key not in cache:
                    cache[key] = low.prepared_envelope(curves, weights, floors, slopes)
                row = low.evaluate_prepared(cache[key], r, a0, mark, original['log2_total_arity_marker'])
                if row['failure_log2_screen'] > original['failure_log2_screen']+1e-7:
                    raise ArithmeticError('Identical-witness paired bound exceeds its packetwise reference')
                row.update(anchor_index=anchor_index, floor_template_index=original['floor_template_index'],
                    floor_template=original['floor_template'], reference_witness_index=original_index,
                    reference_failure_log2_screen=original['failure_log2_screen'],
                    change_from_packetwise=row['failure_log2_screen']-original['failure_log2_screen'])
                indices.append(len(result['witnesses']))
                result['witnesses'].append(row)
            best_index = min(indices, key=lambda i: result['witnesses'][i]['failure_log2_screen'])
            best_indices.append(best_index)
            if progress:
                row = result['witnesses'][best_index]
                progress(dict(stage='paired-bank-fixed-controls', active_groups=r,
                    failure_log2_screen=row['failure_log2_screen'], change_from_packetwise=row['change_from_packetwise'],
                    completed_anchors=anchor_index+1, elapsed_seconds=perf_counter()-start))
        check_time()
        rows = [dict(result['witnesses'][i], witness_index=i) for i in best_indices]
        result.update(status='completed', complete_requested_anchors=True, selected_anchor_bounds=rows,
            selected_witness_indices=best_indices,
            failure_log2_screen=float(low.state.logsum([row['failure_log2_screen'] for row in rows])),
            deterministic_witness_selection_per_anchor=True, extra_witness_factor=False)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-start
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Frozen artifacts must not be overwritten')
    result = run(max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('banks', 'witnesses', 'selected_anchor_bounds')}, allow_nan=False))
