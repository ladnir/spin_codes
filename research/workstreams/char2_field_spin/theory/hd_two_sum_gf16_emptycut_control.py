"""Saved-bank structural control: separate empty R00, retain two states.

Frozen paired-NZ curves and original numeric slopes, no regional builds
or slope fitting. Eight fixed type-weight templates at each of three
anchors. All results are FLOAT proposals, never certificates.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_gf16_state_emptycut as cut
import hd_two_sum_gf16_low_arity as low
import hd_two_sum_input_weight as joint
import hd_two_sum_input_weight_control as joint_control

DIRECTORY = Path(__file__).resolve().parent
PAIRED_NAME = 'hd_two_sum_gf16_paired_low_anchors_float.json'
PAIRED_SHA = '10dd04490a95f124767925c44954c0b375ebcc6f3d25506e52f898452c29761b'


def type_templates(r, mark):
    dense = low.fractional.fractional_dense_weights(low.ALPHA)
    for template in low.floor_templates(r, mark):
        floors = template['log2_floor_markers']
        yield dict(name='old-'+template['name'], family='nonempty-column-floor',
                   markers=[str(x) for x in floors]), cut.fractional_cardinality_weights(low.ALPHA, dense, floors)
    for name, markers in joint_control.marker_templates(r, mark):
        cache = {float(v): joint.type_weights(low.ALPHA, float(v)) for v in set(markers)}
        yield dict(name='joint-'+name, family='pre-H-total-support-floor',
                   markers=[str(F(float(v))) for v in markers]), {a: cache[float(markers[a-1])][a] for a in range(1, 17)}


def run(max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30:
        raise ValueError('An explicit at-most30s budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Empty-separated state control exhausted its budget')
    raw = (DIRECTORY/PAIRED_NAME).read_bytes()
    if sha256(raw).hexdigest() != PAIRED_SHA:
        raise ValueError('Frozen paired-NZ bank changed')
    source = json.loads(raw)
    low.checked_dependencies()
    result = dict(source_bank=PAIRED_NAME, source_bank_sha256=PAIRED_SHA,
        regional_banks_rebuilt=0, regional_amplitudes='independent-uniform-nonzero',
        good_pair_composed_before_alpha=True, fixed_numeric_slopes=True,
        fitted_parameters=0, empty_00_separated=True, state_endpoints_retained=True,
        type_weights_fixed_across_columns=True, controls_per_anchor=8,
        exact_arity_restriction_used=False, all_source_supports_per_anchor=True,
        all_cohorts_per_anchor=True, continuous_weight_band_covered=False,
        directed=False, whole_code_certificate=False, compute_budget_seconds=max_seconds,
        witnesses=[], complete_requested_controls=False)
    try:
        for index, (r, a0, mark) in enumerate(low.ANCHORS):
            bank = source['banks'][index]
            if bank['active_groups'] != r or bank['mark'] != str(mark) or bank['alpha'] != str(low.ALPHA):
                raise ValueError('Frozen paired-bank anchor mismatch')
            curves = np.asarray([[(-np.inf if x is None else x) for x in row] for row in bank['curves']])[:, :r+1]
            slopes = np.asarray(source['witnesses'][15*index]['log2_slopes'])
            for metadata, weights in type_templates(r, mark):
                check_time()
                prepared = cut.prepared_envelope(curves, slopes, weights)
                row = cut.evaluate_prepared(prepared, r, mark, low.ALPHA, check_time=check_time)
                terms = (prepared['path_a']+r/16*low.state.logsum(16*prepared['typed'], axis=1)
                         +np.log2([x[3] for x in prepared['classes']]))
                candidates = [i for i, entry in enumerate(prepared['classes']) if entry[:2] == (0, 0)]
                dominant = max(candidates, key=lambda i: terms[i])
                row.update(template=metadata,
                    dominant_00_path_transition_counts=list(prepared['classes'][dominant][2]),
                    dominant_00_path_allowed_counts=prepared['allowed_counts'][dominant].tolist(),
                    dominant_type=int(np.argmax(prepared['typed'][dominant]))+1,
                    dominant_00_column_log_term=float(terms[dominant]), log2_slopes=slopes.tolist())
                result['witnesses'].append(row)
                if progress:
                    progress(dict(active_groups=r, template=metadata['name'],
                        failure_log2_screen=row['failure_log2_screen'],
                        dominant_allowed_counts=row['dominant_00_path_allowed_counts'],
                        dominant_type=row['dominant_type']))
        check_time()
        result.update(status='completed', complete_requested_controls=True,
            selected=[min((row for row in result['witnesses'] if row['active_groups'] == r),
                          key=lambda row: row['failure_log2_screen']) for r, _, _ in low.ANCHORS])
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=30)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen artifact')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['dependency_pins'] = {Path(module.__file__).name: sha256(Path(module.__file__).read_bytes()).hexdigest()
                               for module in (cut, low, joint, joint_control)}
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key not in ('witnesses', 'selected')}))
