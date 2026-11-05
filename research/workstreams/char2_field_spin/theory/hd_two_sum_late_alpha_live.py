"""One fresh F1 R1 bank: all quiet regional states via a shared matrix.

Keep the frozen scalar moment envelope M_m, with M_0=1. Define the four
entries D_ij=max_m R1_F1(m)_ij/M_m over EVERY m0..18. Then the chronological
product of256 actual regional matrices is at most prod_m M_m timesD^256.
The earlier fixed-rowset quiet polynomial remains valid. No probability
row normalization, common potential, or matrix commutation is assumed.

This only adds live regional states to the quiet outer component. The
nonquiet outer moment remainder is not included. FLOAT only.
"""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from hd_two_sum_gf16_factorized import forced_nonzero_regions
from hd_two_sum_late_alpha_moment import two_atom_envelope
from hd_two_sum_late_alpha_quiet import GROUPS, ARITY, COLUMNS, MARK
from hd_two_sum_late_alpha_source_f1 import restrict_regions
from hd_two_sum_late_alpha_word_order import read_source


DIRECTORY = Path(__file__).resolve().parent
WORD_ARTIFACT = 'hd_two_sum_late_alpha_word_order_gate_float.json'
WORD_SHA = 'a48e647831b07086fd0e6f7398759f111b718680f1af3c8f5943909b1b40670c'
SOURCE_PINS = {
    'hd_two_sum_late_alpha_word_order.py': '1e549ed2bc7f260c8a45029edbed1727acb50f59b66b114e3401f7d8de13b95e',
    'hd_two_sum_gf16_factorized.py': 'd34f054a6a1a07cec5ae97c75acdf475763a5c5e5ccab85d8be05cdd281e5da2',
}
ENTRY_ORDER = ('00', '01', '10', '11')


def log_matrix_product(left, right):
    left, right = np.asarray(left), np.asarray(right)
    if (left.shape != (2, 2) or right.shape != (2, 2)
            or np.any(np.isnan(left)) or np.any(np.isnan(right))
            or np.any(np.isposinf(left)) or np.any(np.isposinf(right))):
        raise ValueError('Two finite-or-zero log matrices required')
    return np.logaddexp2(left[:, :1]+right[:1, :], left[:, 1:]+right[1:, :])


def log_matrix_power(matrix, exponent):
    if type(exponent) is not int or exponent < 0:
        raise ValueError('Nonnegative integer matrix exponent required')
    answer = np.array([[0., -np.inf], [-np.inf, 0.]])
    base = np.asarray(matrix)
    # Validate even the exponent-zero case.
    log_matrix_product(answer, base)
    while exponent:
        if exponent & 1:
            answer = log_matrix_product(answer, base)
        exponent >>= 1
        if exponent:
            base = log_matrix_product(base, base)
    return answer


def shared_matrix(curves, moments):
    """Return logD, including every occupancy and the exact order-zero moment."""
    curves, moments = np.asarray(curves), np.asarray(moments)
    if (curves.ndim != 2 or curves.shape[0] != 4
            or moments.shape != (curves.shape[1],) or not len(moments)
            or moments[0] != 0. or np.any(~np.isfinite(moments))
            or np.any(np.isnan(curves)) or np.any(np.isposinf(curves))):
        raise ValueError('Complete four log curves and finite moments with M0=1 required')
    values = np.max(curves-moments[None, :], axis=1)
    slack = values[:, None]+moments[None, :]-curves
    # Exact-zero entries may create inf-inf only when an entire curve is zero.
    for row in range(4):
        if values[row] == -np.inf:
            slack[row] = np.inf
    if np.any(slack < -1e-10):
        raise ArithmeticError('Shared matrix did not dominate every occupancy')
    return values.reshape(2, 2), [np.flatnonzero(row < 1e-9).tolist() for row in slack]


def load_inputs():
    raw = (DIRECTORY/WORD_ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != WORD_SHA:
        raise ValueError('Frozen whole-word control changed')
    for name, expected in SOURCE_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != expected:
            raise ValueError('Frozen dependency changed: '+name)
    word = json.loads(raw)
    required = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), endpoint_alpha='1', source_count_outside_alpha=True,
        rowset_sum_inside_alpha=True, whole_word_tail_before_alpha=True,
        complete_requested_control=True, directed=False, whole_code_certificate=False)
    if any(word.get(key) != value for key, value in required.items()):
        raise ValueError('Only the frozen direct-first-moment F1 quiet control is allowed')
    source, curve, count = read_source()
    if word.get('source_count') != str(count):
        raise ValueError('Source-class count mismatch')
    return source, word, curve, count


def run(*, max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit positive budget at most30 seconds required')
    source, word, old_curve, count = load_inputs()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Shared-matrix quiet control exhausted its budget')
    result = dict(active_groups=GROUPS, arity=ARITY, source_full_pair_max=1,
        mark=str(MARK), alpha='1', regional_alpha=1,
        scope='Fmax1 sources; quiet outer moments; all regional boundary-state paths',
        all_sources=False, all_sources_in_declared_stratum=True,
        all_fixed_arity_rowset_vectors=True, rowsets_resampled=False,
        source_count=str(count), rowset_count=source['rowset_count'],
        source_count_once=True, projective_normalizations=1,
        nonquiet_outer_moments_included=False, all_regional_state_paths_included=True,
        initial_state=0, terminal_vector=[1, 1], regional_matrix_count=COLUMNS*COLUMNS,
        matrix_entry_order=list(ENTRY_ORDER), matrix_rows_are_probability_laws=False,
        same_scalar_moment_at_every_state=True, order_zero_moment=1,
        common_potential_used=False, matrix_commutation_assumed=False,
        aggregate_input_weight_floor=source['aggregate_input_weight_floor'],
        individual_row_floors=False, all_rows_seen_constraint=False,
        whole_word_artifact=WORD_ARTIFACT, whole_word_artifact_sha256=WORD_SHA,
        reviewed_source_sha256=SOURCE_PINS, raw_regional_banks=0,
        new_moment_candidates=0, parameter_fits=0, directed=False,
        whole_code_certificate=False, full_quiet_event_bounded=False,
        failure_probability_lower_bound=False, compute_budget_seconds=max_seconds,
        complete_requested_control=False)
    try:
        check_time()
        bank = forced_nonzero_regions(float(MARK), GROUPS)
        result['raw_regional_banks'] = 1
        check_time()
        restricted = restrict_regions(bank)
        curves = np.asarray([restricted[m, 0].ravel() for m in range(GROUPS+1)]).T
        if (not np.array_equal(np.isfinite(curves[0]), np.isfinite(old_curve))
                or np.max(np.abs(curves[0, np.isfinite(old_curve)]-old_curve[np.isfinite(old_curve)])) > 1e-9):
            raise ArithmeticError('Fresh F1 zero-state curve did not reproduce')
        moments, envelope = two_atom_envelope(old_curve, 5, 18)
        matrix, contacts = shared_matrix(curves, moments)
        powered = log_matrix_power(matrix, COLUMNS*COLUMNS)
        loss = float(np.logaddexp2.reduce(powered[0]))
        baseline = word['comparison']['direct_first_moment_alpha1_log2']
        serialize = lambda values: [[None if value == -np.inf else float(value) for value in row] for row in values]
        result.update(regional_full_pair_restriction='f<=min(1,floor(m/2))',
            regional_log2_curves=serialize(curves), zero_entries_are_none=True,
            two_atom_envelope=envelope, shared_matrix_log2=serialize(matrix),
            shared_matrix_contacts=dict(zip(ENTRY_ORDER, contacts)),
            powered_matrix_log2=serialize(powered),
            log2_state_path_factor=loss, old_allzero_majorant_term_log2=baseline,
            all_quiet_state_paths_majorant_term_log2=baseline+loss,
            complete_requested_control=True)
        check_time()
        if progress:
            progress(dict(shared_matrix_log2=result['shared_matrix_log2'],
                contacts=result['shared_matrix_contacts'], log2_state_path_factor=loss,
                all_quiet_state_paths_majorant_term_log2=baseline+loss))
    except TimeoutError as error:
        result['complete_requested_control'] = False
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--max-seconds', default=30, type=float)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(complete=answer['complete_requested_control'], elapsed_seconds=answer['elapsed_seconds'])))
