"""Bounded saved-bank FLOAT control of the9a pre-H support constraint.

Three existing low-arity anchors, fixed regional envelopes, five fixed
marker templates each. No regional provider is called. The exact total-A
restriction is also computed; neither result is a directed certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_input_weight as joint
import hd_two_sum_gf16_low_arity as low
import hd_two_sum_gf16_exact_arity as exact
import hd_two_sum_gf16_state_envelope as state

DIRECTORY = Path(__file__).resolve().parent
BANK_NAME = 'hd_two_sum_gf16_low_arity_anchors_float.json'
BANK_SHA = '206f34cca8021de88ce8806a05c6736523261bf769ff76f147542e7cb6fdc5d5'


def marker_templates(r, mark):
    prefix = float(low.ALPHA)*(-low.fractional.LIMIT*log2(float(mark))-log2(255))
    q = (prefix+log2(comb(8192, r))+60)/r+4
    arities = np.arange(1, 17)
    empty = np.array([log2(comb(16, a))-8*a*float(low.ALPHA)*log2(255) for a in arities])
    threshold = np.maximum(0., (q+empty)/(9*arities*float(low.ALPHA)))
    return [('zero', np.zeros(16)), ('uniform1', np.ones(16)), ('uniform2', np.full(16, 2.)),
            ('empty-floor1', threshold), ('empty-floor2', 2*threshold)]


def prepared_joint(curves, slopes, log_markers):
    intercepts, slack = state.affine_envelope(curves, slopes)
    classes = state.path_classes(16)
    path_counts = np.asarray([row[2] for row in classes])
    marks = np.asarray([np.repeat(slopes, row[2]) for row in classes])
    elementary = state.elementary_logs(marks)
    markers = np.asarray(log_markers)
    if markers.shape != (16,) or np.any(~np.isfinite(markers)) or np.any(markers < 0):
        raise ValueError('One nonnegative pre-H marker log per fixed arity required')
    # Distinct arities can choose distinct fixed markers, but never adapt
    # them to paths or columns. Cache common markers within this witness.
    weights = {float(value): joint.type_weights(low.ALPHA, float(value)) for value in set(markers)}
    typed = np.stack([state.logsum(elementary+weights[float(markers[a-1])][a], axis=-1)
                      for a in range(1, 17)], axis=-1)
    return dict(classes=classes, path_a=path_counts@intercepts, typed=typed,
                slopes=np.asarray(slopes), intercepts=intercepts, slack=slack)


def unmarked_bound(prepared, r, mark):
    kernel = state.kernel_from_types(prepared['classes'], prepared['path_a'], prepared['typed'], r)
    power = state.matrix_power_log(kernel, 16)
    prefix = float(low.ALPHA)*(-low.fractional.LIMIT*log2(float(mark))-log2(255))
    return float(prefix+log2(comb(8192, r))+state.logsum(power[0]))


def run(max_seconds=60, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 60:
        raise ValueError('An explicit at-most60s control budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Joint-weight control budget exceeded')
    path = DIRECTORY/BANK_NAME
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != BANK_SHA:
        raise ValueError('Pinned packetwise fractional bank changed')
    data = json.loads(raw)
    low.checked_dependencies()
    result = dict(bank=path.name, bank_sha256=BANK_SHA, regional_banks_rebuilt=0,
        paired_packet_operator=False, pre_H_support_floor='9*a', pre_H_marker_before_alpha=True,
        marker_fixed_per_arity_across_all_paths_and_columns=True,
        exact_input_weight_counts=True, rowset_holder_exponent=16,
        all_source_supports_in_each_selected_low_class=True, directed=False,
        high_side_included=False, whole_code_certificate=False, continuous_weight_band_covered=False,
        compute_budget_seconds=max_seconds, witnesses=[], complete_requested_controls=False)
    try:
        for anchor_index, (r, a0, mark) in enumerate(low.ANCHORS):
            bank = data['banks'][anchor_index]
            if bank['active_groups'] != r or bank['mark'] != str(mark) or bank['alpha'] != str(low.ALPHA):
                raise ValueError('Unexpected frozen anchor metadata')
            curves = np.array([[(-np.inf if x is None else x) for x in row] for row in bank['curves']])[:, :r+1]
            slopes = np.array(bank['initial_log2_slopes'])
            for name, markers in marker_templates(r, mark):
                check_time()
                prepared = prepared_joint(curves, slopes, markers)
                unrestricted = unmarked_bound(prepared, r, mark)
                restricted = exact.evaluate_prepared(prepared, r, a0, mark, low.ALPHA, check_time=check_time)
                row = dict(active_groups=r, maximum_total_arity=a0-1, mark=str(mark), alpha=str(low.ALPHA),
                    marker_template=name, log2_pre_H_markers=[str(F(float(x))) for x in markers],
                    log2_slopes=slopes.tolist(), log2_intercepts=prepared['intercepts'].tolist(),
                    unrestricted_failure_log2_screen=unrestricted,
                    exact_low_failure_log2_screen=restricted['failure_log2_screen'])
                result['witnesses'].append(row)
                if progress:
                    progress({key: row[key] for key in ('active_groups', 'marker_template',
                        'unrestricted_failure_log2_screen', 'exact_low_failure_log2_screen')})
        check_time()
        result.update(status='completed', complete_requested_controls=True)
    except TimeoutError as error:
        result.update(status='compute-budget-exceeded', reason=str(error))
    result['elapsed_seconds'] = perf_counter()-started
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-seconds', type=float, default=60)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite research artifacts')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['dependency_pins'] = {Path(module.__file__).name: sha256(Path(module.__file__).read_bytes()).hexdigest()
                               for module in (joint, low, exact, state)}
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'witnesses'}))
