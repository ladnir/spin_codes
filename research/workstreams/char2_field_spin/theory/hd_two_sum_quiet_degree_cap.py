"""ONE saved-bank quiet control: cap each degree by its outer-only mass.

For each quiet total input degreeK, both the unweighted outer comparison
and the marked inner comparison upper-bound the bad-word first moment.
Take their minimum AFTER each complete degree sum, then sum degrees.
No negative subtraction, provider, fitting, or nonquiet scope is used.
"""
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_matrix_moment as frozen
from hd_two_sum_closed_kernel import zero_two_moments


DIRECTORY = Path(__file__).resolve().parent
ARTIFACT = 'hd_two_sum_matrix_moment_quiet_float.json'
ARTIFACT_SHA = '9263a38cdf5311a8f5ddf10490c8f9472eaecee46bf8629a86f94f65ec53e823'
HELPER_SHA = 'bc192ebdafc9abec59f7977e660c5d4741737eb62cf479c920fdbb83a060ce67'


def degree_cap(unweighted, marked, log_marker_factor, floor):
    """Both input arrays have summed all words/state endpoints at degreeK."""
    unweighted, marked = np.asarray(unweighted), np.asarray(marked)
    if (unweighted.ndim != 1 or marked.shape != unweighted.shape
            or np.any(np.isnan(unweighted)) or np.any(np.isnan(marked))
            or np.any(np.isposinf(unweighted)) or np.any(np.isposinf(marked))
            or not np.isfinite(log_marker_factor) or type(floor) is not int
            or not 0 <= floor < len(marked)):
        raise ValueError('Two complete positive degree comparisons and valid floor required')
    capped = np.minimum(unweighted, marked+log_marker_factor)
    return capped, dict(marked_tail_log2=float(np.logaddexp2.reduce(marked[floor:]))+log_marker_factor,
        unweighted_tail_log2=float(np.logaddexp2.reduce(unweighted[floor:])),
        capped_tail_log2=float(np.logaddexp2.reduce(capped[floor:])),
        peak_capped_degree=floor+int(np.argmax(capped[floor:])),
        capped_degrees=[int(k) for k in range(floor, len(marked))
                        if unweighted[k] < marked[k]+log_marker_factor])


def outer_word(groups, columns, arity, *, size=16, check_time=lambda: None):
    """Unweighted quiet count for fixed A, summed over A ONCE per group.

    Once output positions are unmarked, every fixed A has the same column
    polynomial sum_b C(arity,b) Z2(b) X^b. No Holder relaxation is needed.
    Its columns-th power is multiplied by C(size,arity), not by that
    multiplicity in each column. RS normalization remains outside.
    """
    if (any(type(x) is not int for x in (groups, columns, arity, size))
            or not 1 <= arity <= size <= 16 or groups < 1 or columns < 1):
        raise ValueError('Positive fixed-row-set parameters required')
    counts = [comb(arity, w)*zero_two_moments(w) for w in range(arity+1)]
    group = np.array([log2(count) if count else -np.inf for count in counts])
    return (frozen.log_power(group, groups*columns, check_time)
            +groups*log2(comb(size, arity)))


def run(*, max_seconds=30):
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 30:
        raise ValueError('Explicit budget at most30 seconds required')
    raw = (DIRECTORY/ARTIFACT).read_bytes()
    if sha256(raw).hexdigest() != ARTIFACT_SHA or sha256(Path(frozen.__file__).read_bytes()).hexdigest() != HELPER_SHA:
        raise ValueError('Frozen quiet comparator changed')
    source = json.loads(raw)
    if source.get('scope') != 'quiet-only' or source.get('complete_requested_controls') is not True:
        raise ValueError('Completed quiet-only comparison required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Quiet degree-cap control exhausted its budget')
    result = dict(active_groups=18, arity=8, source_full_pair_max=1, mark='691/700', alpha=1,
        scope='quiet outer moments; all regional state paths; homogeneous arity8/F1',
        nonquiet_outer_moments_included=False, marked_fixed_rowset_holder_exponent=16,
        outer_rowsets_counted_once_per_group=True, outer_holder_used=False,
        marked_extra_rowset_factor=False, source_count_once=True, projective_normalizations=1,
        rs_normalizer_exponent=1152, minimum_input_degree=1296,
        minimum_after_complete_degree_sums=True, signed_subtractions=False,
        raw_regional_banks_built=0, parameter_fits=0, directed=False,
        whole_code_certificate=False, failure_probability_lower_bound=False,
        reference_artifact=ARTIFACT, reference_sha256=ARTIFACT_SHA,
        frozen_helper_sha256=HELPER_SHA, compute_budget_seconds=max_seconds,
        complete_requested_control=False)
    try:
        data, curves = frozen.load_saved()
        old = data['two_atom_envelope']
        measure = frozen.entry_measure(curves, old['log2_probability'], old['log2_nonzero_atom'])
        weights = frozen.raw_holder_weights(quiet=True)
        matrix, _ = frozen.column_matrix(weights, frozen.GROUPS, measure,
            arity=frozen.ARITY, check_time=check_time)
        word = frozen.matrix_poly_power(matrix, frozen.COLUMNS, check_time)
        marked = np.logaddexp2.reduce(word[0], axis=0)
        outer = outer_word(frozen.GROUPS, frozen.COLUMNS, frozen.ARITY,
                           check_time=check_time)
        marker = -frozen.LIMIT*log2(float(frozen.MARK))
        capped, diagnosis = degree_cap(outer, marked, marker, frozen.FLOOR)
        common_prefix = (log2(int(data['source_count']))-log2(255)
                         -frozen.DIMENSION*frozen.ARITY*frozen.GROUPS*log2(255))
        old_bound = next(row['log2_screen'] for row in source['controls'] if row['mode'] == 'polynomial')
        if abs(common_prefix+diagnosis['marked_tail_log2']-old_bound) > 1e-8:
            raise ArithmeticError('Frozen quiet bound did not reproduce')
        top = sorted(range(frozen.FLOOR, len(capped)), key=lambda k: capped[k], reverse=True)[:12]
        result.update(diagnosis=diagnosis, common_prefix_log2=common_prefix,
            old_marked_screen=old_bound, outer_only_screen=common_prefix+diagnosis['unweighted_tail_log2'],
            capped_screen=common_prefix+diagnosis['capped_tail_log2'],
            peak_degrees=[dict(degree=int(k), outer_only_log2=float(common_prefix+outer[k]),
                marked_log2=float(common_prefix+marker+marked[k]),
                capped_log2=float(common_prefix+capped[k])) for k in top],
            source_count=data['source_count'], complete_requested_control=True)
        check_time()
    except TimeoutError as error:
        result['complete_requested_control'] = False
        result['incomplete_reason'] = str(error)
    result['elapsed_seconds'] = perf_counter()-started
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    return result


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-seconds', type=float, default=30)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(max_seconds=args.max_seconds)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(answer))
