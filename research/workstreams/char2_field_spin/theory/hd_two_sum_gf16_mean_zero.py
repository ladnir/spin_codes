"""Two saved-bank controls: quiet/full outer law, zero regional boundaries.

This deliberately does NOT run the failed common-C full-state comparator.
Both controls use the SAME pinned R1_00 moment law, one alpha1 witness,
fixed-A Holder, and complete word floor. No provider or fit is called.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import log2
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_mean_column import (common_matrix, raw_holder_weights,
    whole_word_screen)
from hd_two_sum_late_alpha_source_f1 import load_source_count
from hd_xor_multigroup import LIMIT


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_two_sum_late_alpha_live_gate_float.json'
ARTIFACT_SHA = '33da1eb9517c7d5665d22d2c8b93f5622199d4e1b70dd2a8e3aa52e28d4e04cf'
SOURCE_PINS = {
    'hd_two_sum_late_alpha_live.py': 'e28bd5473fe60be32a984e67cbd0b64f52078a30aa4f6ff4b5b304497fd4fb59',
    'hd_two_sum_gf16_mean_column.py': 'fa6377002fa52fbfcb980ddd4b72d30df16b6c12c7c6db24eed54554e61a1593',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
    'hd_two_sum_late_alpha_quiet.py': '8a753036a9be44759206b5beb3baa044f28d80b265a9304da81bdde6e8f11ae3',
}


def checked_inputs(data, source_count):
    required = dict(active_groups=18, arity=8, source_full_pair_max=1,
        mark='691/700', alpha='1', regional_alpha=1, all_sources=False,
        all_sources_in_declared_stratum=True, all_fixed_arity_rowset_vectors=True,
        rowsets_resampled=False, source_count=str(source_count),
        initial_state=0, terminal_vector=[1, 1], regional_matrix_count=256,
        matrix_entry_order=['00', '01', '10', '11'], order_zero_moment=1,
        aggregate_input_weight_floor=1296, raw_regional_banks=1,
        new_moment_candidates=0, parameter_fits=0, directed=False,
        whole_code_certificate=False, complete_requested_control=True,
        regional_full_pair_restriction='f<=min(1,floor(m/2))',
        source_sha256=SOURCE_PINS['hd_two_sum_late_alpha_live.py'])
    if (type(source_count) is not int or source_count <= 0
            or any(data.get(key) != value for key, value in required.items())):
        raise ValueError('The exact saved r18/a8/F1 raw regional bank is required')
    envelope = data.get('two_atom_envelope', {})
    if envelope.get('secant_endpoints') != [5, 18] or envelope.get('moment_order_zero') != 1:
        raise ValueError('Only the frozen00 moment law is allowed')
    log_p, log_b = (envelope.get('log2_probability'), envelope.get('log2_nonzero_atom'))
    curves = np.asarray([[(-np.inf if value is None else value) for value in row]
                         for row in data.get('regional_log2_curves', [])], dtype=float)
    if curves.shape != (4, 19):
        raise ValueError('Every raw matrix entry/occupancy must be saved')
    matrix, _ = common_matrix(curves.T.reshape(19, 2, 2), log_p, log_b)
    if matrix[0, 0] != 0. or curves[0, 0] != 0.:
        raise ArithmeticError('The frozen moment law must cover00 with C00=1 exactly')
    # Keep only the all-zero boundary path; do not inherit the other C entries.
    zero_path = np.array([[0., -np.inf], [-np.inf, -np.inf]])
    return zero_path, log_p, log_b


def read_inputs():
    raw = (DIRECTORY/ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != ARTIFACT_SHA:
        raise ValueError('Pinned raw regional artifact changed')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Pinned dependency changed: '+name)
    data = json.loads(raw)
    source_count = load_source_count()
    return data, source_count, checked_inputs(data, source_count)


def run(*, max_seconds=30, progress=None):
    if (type(max_seconds) not in (int, float) or not np.isfinite(max_seconds)
            or not 0 < max_seconds <= 30):
        raise ValueError('Explicit finite positive budget at most30 seconds required')
    data, source_count, (matrix, log_p, log_b) = read_inputs()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Saved zero-boundary outer-law control exhausted its budget')
    result = dict(active_groups=18, arity=8, source_full_pair_max=1,
        mark='691/700', alpha=1, field_order=256, forbidden_weight=LIMIT,
        source_count=str(source_count), columns=16, constituent_dimension=8,
        source_count_once=True, rowset_choices_inside_holder=True,
        row_sets_fixed_across_columns=True, rowsets_resampled=False,
        extra_rowset_factor=False, rs_normalizer_exponent=1152,
        all_sources=False, all_sources_in_declared_stratum=True,
        all_regional_state_paths_included=False, initial_state=0,
        scope='Fmax1; homogeneousa8; all regional boundary states zero only',
        aggregate_surrogate_floor=1296, surrogate='min(output_weight,8)',
        individual_row_floors=False, component_of_positive_majorant=True,
        complete_homogeneous_class=False, whole_code_certificate=False,
        failure_probability_lower_bound=False, directed=False,
        raw_regional_banks_built=0, parameter_fits=0, requested_controls=2,
        artifact=ARTIFACT, artifact_sha256=ARTIFACT_SHA,
        reviewed_source_sha256=SOURCE_PINS, compute_budget_seconds=max_seconds,
        two_atom_envelope=data['two_atom_envelope'], controls=[],
        complete_requested_controls=False)
    try:
        for quiet in (True, False):
            check_time()
            weights = raw_holder_weights(quiet=quiet)
            row = whole_word_screen(weights, matrix, 18, log_p, log_b,
                log_source_count=log2(source_count), mark=float(F('691/700')),
                limit=LIMIT, check_time=check_time)
            row.update(outer_scope='quiet-only' if quiet else 'all-outer-moments',
                nonquiet_outer_moments_included=not quiet,
                all_regional_state_paths_included=False,
                component_of_positive_majorant=True)
            result['controls'].append(row)
            if progress:
                progress(dict(outer_scope=row['outer_scope'], log2_screen=row['log2_screen'],
                              elapsed_seconds=perf_counter()-started))
        if result['controls'][1]['log2_screen']+1e-8 < result['controls'][0]['log2_screen']:
            raise ArithmeticError('Adding nonquiet positive branches reduced the comparison')
        check_time()
        result['complete_requested_controls'] = True
    except TimeoutError as error:
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', default=30, type=float)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    result = run(max_seconds=args.max_seconds,
                 progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    with args.output.open('x', encoding='utf-8', newline='\n') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=result['complete_requested_controls'],
        output_sha256=sha256(args.output.read_bytes()).hexdigest())))
