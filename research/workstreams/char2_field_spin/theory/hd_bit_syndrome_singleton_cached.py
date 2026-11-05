"""Cached-measure singleton gate retaining simultaneous regional activity.

Every nonzero singleton maps to all16 outer output positions. Consequently
K=1+255*Y*product(t_i) exactly, independently of the labelled source row.
For r groups, column coefficient k is C(r,k)*255^k*E(k)^16, where E(k)
is the common matrix measure's kth moment. Ordered16-column multiplication
retains both inner states and the total nonzero-column degree. Degree>=9r
is a necessary relaxation of each group's RS[16,8] minimum distance.

No uniform-cover hypothesis, provider, per-row-set DP, or construction
change is used. Cached regional parameters and all arithmetic are FLOAT.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

DIRECTORY = Path(__file__).resolve().parent
INPUT = 'hd_mds_exact_empty_mixed789_9927_pinned_float.json'
INPUT_SHA = 'c89ef4973de29007fe293817d2e3b9f582246a45d19ca19c7fcce2554e7d33b8'
COVER = 'hd_uniform_cover_conditional_9927_float.json'
COVER_SHA = 'd56a6ab7542a9e85fd4a83a71a3a386aa121788bbeeb1500ff5e6a895947e984'


def log_mm(left, right):
    return np.logaddexp2(left[:, 0, None]+right[None, 0, :],
                        left[:, 1, None]+right[None, 1, :])


def column_coefficients(moment, groups, regions=16, q=256):
    if (type(groups) is not int or groups < 1 or type(regions) is not int or regions < 1
            or type(q) is not int or q < 2):
        raise ValueError('Positive group/region counts and a field size required')
    result = np.full((2, 2, groups+1), -np.inf)
    for k in range(groups+1):
        matrix = np.asarray(moment(k), dtype=float)
        if matrix.shape != (2, 2) or np.isnan(matrix).any() or np.isposinf(matrix).any():
            raise ValueError('A finite-or-zero2x2 moment is required')
        word = np.full((2, 2), -np.inf)
        word[0, 0] = word[1, 1] = 0.
        for _ in range(regions):
            word = log_mm(word, matrix)
        result[:, :, k] = word+log2(comb(groups, k))+k*log2(q-1)
    return result


def ordered_word_coefficients(column, columns=16):
    column = np.asarray(column, dtype=float)
    if (column.ndim != 3 or column.shape[:2] != (2, 2) or column.shape[2] < 1
            or type(columns) is not int or columns < 1
            or np.isnan(column).any() or np.isposinf(column).any()):
        raise ValueError('Positive column count and2x2 coefficient array required')
    row = np.array([[0.], [-np.inf]])
    degree = column.shape[2]-1
    for _ in range(columns):
        following = np.full((2, row.shape[1]+degree), -np.inf)
        for i in range(2):
            for j in range(2):
                for k in range(degree+1):
                    following[j, k:k+row.shape[1]] = np.logaddexp2(
                        following[j, k:k+row.shape[1]], row[i]+column[i, j, k])
        row = following
    return np.logaddexp2(row[0], row[1])


def marker_minimum(coefficients, floor, maximum=16., rounds=48):
    values = np.asarray(coefficients, dtype=float)
    if not 0 <= floor < len(values) or maximum <= 0:
        raise ValueError('An in-scope floor and positive marker bound required')
    degrees = np.arange(len(values))
    def evaluate(y):
        marked = values+degrees*y
        mass = np.logaddexp2.reduce(marked)
        mean = float(np.sum(degrees*np.exp2(marked-mass)))
        return float(mass-floor*y), mean-floor
    lo, hi = 0., float(maximum)
    if evaluate(lo)[1] >= 0:
        best = lo
    elif evaluate(hi)[1] <= 0:
        best = hi
    else:
        for _ in range(rounds):
            middle = (lo+hi)/2
            if evaluate(middle)[1] < 0:
                lo = middle
            else:
                hi = middle
        best = (lo+hi)/2
    return dict(log2_marker=best, marked_tail_log2=evaluate(best)[0],
                derivative=evaluate(best)[1], interval=[0., maximum], bisection_rounds=rounds)


def run(*, marker_maximum=16.):
    start = perf_counter()
    raw = (DIRECTORY/INPUT).read_bytes()
    cover_raw = (DIRECTORY/COVER).read_bytes()
    if sha256(raw).hexdigest() != INPUT_SHA or sha256(cover_raw).hexdigest() != COVER_SHA:
        raise ValueError('Frozen cached inputs changed')
    saved, cover = json.loads(raw), json.loads(cover_raw)
    if (saved['groups'] != 18 or saved['mark'] != '9927/10000'
            or saved['source_class'] != 'Fmax<=1' or saved['input_marker'] != 1):
        raise ValueError('Wrong cached scope')
    measure, groups = saved['state_measure'], saved['groups']
    def moment(k):
        return np.array([[0. if k == 0 else measure['log2_probability']+k*measure['log2_atom00'],
                          measure['log2_scale01']+k*measure['log2_atom01']],
                         [measure['log2_scale10']+k*measure['log2_atom10'], measure['log2_scale11']]])
    coefficients = ordered_word_coefficients(column_coefficients(moment, groups))
    floor = 9*groups
    tail = float(np.logaddexp2.reduce(coefficients[floor:]))
    optimized = marker_minimum(coefficients, floor, maximum=marker_maximum)
    # One group has16 possible labelled singleton rows, fixed across columns.
    prefix = (log2(int(saved['source_count']))-734003*log2(F(saved['mark']))-log2(255)
              +groups*log2(16)-8*groups*log2(255))
    old = next(row for row in cover['arities'] if row['arity'] == 1)
    return dict(construction=saved['construction'], groups=groups, arity=1,
        source_class='Fmax<=1', mark=saved['mark'], common_measure=saved['state_measure'],
        fixed_source_row_choices=16, row_choices_counted_once=True,
        simultaneous_all16_regional_activity=True, all_inner_state_paths=True,
        aggregate_active_column_floor=floor, individual_row_floors_relaxed=True,
        exact_positive_tail_log2=prefix+tail,
        optimized_scalar_marker_log2=prefix+optimized['marked_tail_log2'],
        unmarked_log2=prefix+float(np.logaddexp2.reduce(coefficients)),
        marker_control=optimized, prefix_log2=prefix,
        hypothetical_uniform_cover_log2=old['conditional_log2'],
        word_log2_coefficients=[None if value == -np.inf else float(value) for value in coefficients],
        input_receipt=INPUT, input_sha256=INPUT_SHA, conditional_cover_receipt_sha256=COVER_SHA,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
        elapsed_seconds=perf_counter()-start, regional_providers_run=0, per_rowset_contractions_run=0,
        regional_moment_fit=False, scalar_activity_fit=True, uniform_cover_hypothesis_used=False,
        float_only=True, directed=False, whole_code_certificate=False)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--marker-maximum', type=float, default=16.)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen control')
    answer = run(marker_maximum=args.marker_maximum)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(answer, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: answer[key] for key in ('exact_positive_tail_log2',
        'optimized_scalar_marker_log2', 'unmarked_log2', 'marker_control', 'elapsed_seconds')}))
