"""Single-anchor FLOAT inner gate for pre-H Toeplitz2 pairs.

Only r18,F1,homogeneous arity8,alpha1/5,z691/700. Incidence double
counting recovers the fixed-T sum over all row sets before column Holder.
All16 matchings and all binary state paths are retained chronologically.
Saved regional entries only; ordinary affine envelopes (NOT empty-cut).
"""
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import comb, log2
from pathlib import Path
from time import perf_counter

import numpy as np

import hd_two_sum_pre_toeplitz_fractional as masks
from hd_two_sum_input_weight import marked_count
from hd_two_sum_toeplitz_pairs import xor_schedule


DIRECTORY = Path(__file__).resolve().parent
R, ARITY, GROUP, COLUMNS, FULL_MAX = 18, 8, 16, 16, 1
ALPHA, MARK, LIMIT = F(1, 5), F(691, 700), 734003
MASK_HELPER_SHA = 'a1e2b46ce992b5ab433b2c6bbe57e9407df338b792306be548cd72539aa2ab7b'
OLD_CONSTRUCTION = 'two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half'
CONSTRUCTION = 'pre-toeplitz2-two-sum-gf16-g16-rs16x8-h4-swapped-pairs-k20-rate-half'


def orbit_size(weight, full, group=16):
    if (type(group) is not int or group < 2 or group % 2
            or any(type(x) is not int for x in (weight, full))
            or not 0 <= weight <= group or not max(0, weight-group//2) <= full <= weight//2):
        raise ValueError('A nonempty matching orbit is required')
    single = weight-2*full
    return comb(group//2, full)*comb(group//2-full, single)*2**single


def full_pairs(target, matching):
    return sum(bool(target >> left & 1) and bool(target >> right & 1) for left, right in matching)


def incidence_weights(arity=8, group=16, alpha=ALPHA, columns=16, dimension=8, check_time=None):
    """FLOAT fixed-target Holder weights from exact representative C values.

    Sum_A C(A,T)^p = sum_u n_A(u) sum_{T' in orbit(T)} C(A_u,T')^p / n_T.
    Normalization255^(-(columns-dimension)*arity*alpha) is INSIDE the root.
    """
    alpha = F(alpha)
    if (group not in (8, 16) or type(arity) is not int or not 0 <= arity <= group
            or type(columns) is not int or type(dimension) is not int
            or not 0 < dimension < columns or not 0 < alpha <= 1):
        raise ValueError('Fixed row arity, group, columns, dimension, and exponent required')
    tables = [np.full((group+1, group//2+1), -np.inf) for _ in range(2)]
    total_inputs, counts = 0, []
    for internal in range(max(0, arity-group//2), arity//2+1):
        if check_time:
            check_time()
        multiplicity = orbit_size(arity, internal, group)
        total_inputs += multiplicity
        rows = masks.exact_type_masses(arity, internal, group=group, check_time=check_time)
        matching = masks.matching_for_internal(arity, internal, group)
        for target, target_multiplicity, new, old in rows:
            weight, full = target.bit_count(), full_pairs(target, matching)
            factor = log2(multiplicity)+log2(target_multiplicity)-log2(orbit_size(weight, full, group))
            for table, value in zip(tables, (new, old)):
                value = F(value)
                if value:
                    term = factor+float(columns*alpha)*(log2(value.numerator)-log2(value.denominator))
                    table[weight, full] = np.logaddexp2(table[weight, full], term)
        counts.append(dict(internal_pairs=internal, input_sets=multiplicity, target_types=len(rows)))
    if total_inputs != comb(group, arity):
        raise ArithmeticError('Input orbits must partition the fixed arity')
    normalization = float(alpha)*(columns-dimension)*arity*log2(255)
    new, old = ((table-normalization)/columns for table in tables)
    # Independent old-law recovery: fixed T, explicit intersection sum over A.
    for weight in range(group+1):
        reference = -np.inf
        for overlap in range(max(0, arity+weight-group), min(arity, weight)+1):
            count = marked_count(arity, overlap, weight, F(1), group)
            if count:
                choices = comb(weight, overlap)*comb(group-weight, arity-overlap)
                term = log2(choices)+float(columns*alpha)*log2(count)
                reference = np.logaddexp2(reference, term)
        reference = (reference-normalization)/columns
        for full in range(max(0, weight-group//2), weight//2+1):
            if not np.isclose(old[weight, full], reference, atol=1e-11, rtol=0, equal_nan=False):
                raise ArithmeticError('Incidence formula failed the old fixed-target identity')
    if not np.isclose(new[0, 0], (log2(comb(group, arity))-normalization)/columns, atol=1e-12, rtol=0):
        raise ArithmeticError('The all-empty fixed-target moment must count all A once')
    return new, old, counts


def mask_weight_table(weights, matching):
    """Expand matching-dependent per-mask weights, without any multiplicity."""
    group = 2*len(matching)
    matching = masks.checked_matching(matching, group)
    if np.shape(weights) != (group+1, group//2+1):
        raise ValueError('Complete matching-orbit weight table required')
    index = np.arange(1 << group, dtype=np.uint32)
    cardinality = np.zeros(len(index), dtype=np.uint8)
    pairs = np.zeros(len(index), dtype=np.uint8)
    for left, right in matching:
        a, b = (index >> left) & 1, (index >> right) & 1
        cardinality += (a+b).astype(np.uint8)
        pairs += (a & b).astype(np.uint8)
    result = np.asarray(weights)[cardinality, pairs]
    if np.any(np.isnan(result)) or np.any(np.isposinf(result)):
        raise ValueError('Selected mask weights must be finite or exact zero')
    return result


def path_transform(log_mask_weights, log_b):
    """G(sigma)=sum_T weight(T)*prod_i B[sigma_i,sigma_i+1]^1[i in T].

    Sequential tensor contraction keeps <=2^(g+1) doubles. Bit i of T is
    coordinate i; C-order path axes remain sigma0,...,sigmag. No sorting.
    """
    weights, slopes = np.asarray(log_mask_weights, float), np.asarray(log_b, float)
    size = len(weights)
    group = size.bit_length()-1
    if (weights.shape != (size,) or not 1 <= group <= 16 or size != 1 << group
            or slopes.shape != (2, 2) or np.any(~np.isfinite(slopes))
            or np.any(np.isnan(weights)) or np.any(np.isposinf(weights))):
        raise ValueError('A full binary mask table and four finite slopes required')
    # Transpose bit axes so coordinate0 is consumed FIRST, despite integer
    # mask indexing putting its bit on the least-significant tensor axis.
    values = weights.reshape((2,)*group).transpose(tuple(reversed(range(group)))).copy().reshape(-1)
    values = np.stack((values, values))
    for index in range(group):
        suffix = 1 << (group-index-1)
        values = values.reshape(1 << index, 2, 2, suffix)
        values = np.logaddexp2(values[:, :, 0, None, :],
                              values[:, :, 1, None, :]+slopes[None, :, :, None])
    return values.reshape((2,)*(group+1))


@lru_cache(None)
def path_geometry(group):
    """Flattened chronological path endpoints and four transition counts."""
    if type(group) is not int or not 1 <= group <= 16:
        raise ValueError('Path length1..16 required')
    paths = np.arange(1 << (group+1), dtype=np.uint32)
    counts = np.zeros((len(paths), 4), dtype=np.uint8)
    for i in range(group):
        transition = ((paths >> (group-i-1)) & 3).astype(np.uint8)
        for entry in range(4):
            counts[:, entry] += transition == entry
    return (paths >> group) & 1, paths & 1, counts


def kernel_from_transform(nonempty, log_empty, intercepts, groups, log_activity=0., *, columns=16, minimum=9):
    """No additional Holder step: independent groups give G(path)^r."""
    nonempty, intercepts = np.asarray(nonempty, float), np.asarray(intercepts, float)
    group = nonempty.ndim-1
    if (nonempty.shape != (2,)*(group+1) or not 1 <= group <= 16 or intercepts.shape != (4,)
            or np.any(~np.isfinite(intercepts)) or type(groups) is not int or groups < 1
            or not np.isfinite(log_empty) or not np.isfinite(log_activity) or log_activity < 0
            or not 1 <= minimum <= columns):
        raise ValueError('Complete transformed paths, fixed envelope, and nonnegative activity mark required')
    initial, final, counts = path_geometry(group)
    log_g = np.logaddexp2(log_empty, nonempty.ravel()+log_activity)-minimum*log_activity/columns
    terms = counts @ intercepts+groups*log_g
    matrix = np.empty((2, 2))
    for left in range(2):
        for right in range(2):
            matrix[left, right] = np.logaddexp2.reduce(terms[(initial == left) & (final == right)])
    return matrix


def log_matmul(left, right):
    return np.logaddexp2(left[:, 0, None]+right[None, 0, :],
                        left[:, 1, None]+right[None, 1, :])


def load_frozen_anchor():
    if sha256(Path(masks.__file__).read_bytes()).hexdigest() != MASK_HELPER_SHA:
        raise ValueError('Reviewed exact mask helper changed')
    raw = (DIRECTORY/masks.FROZEN_BANK).read_bytes()
    if sha256(raw).hexdigest() != masks.FROZEN_SHA:
        raise ValueError('Frozen regional bank changed')
    source = json.loads(raw)
    if (source.get('construction') != OLD_CONSTRUCTION or source.get('field_order') != 256
            or source.get('source_size') != 8192 or source.get('forbidden_weight') != LIMIT
            or source.get('directed') is not False or source.get('complete_requested_anchors') is not True):
        raise ValueError('Wrong saved regional construction or scope')
    banks = [bank for bank in source['banks'] if bank['active_groups'] == R]
    if len(banks) != 1 or (banks[0]['mark'], banks[0]['alpha']) != (str(MARK), str(ALPHA)):
        raise ValueError('Exactly the frozen r18 anchor required')
    bank = banks[0]
    counts = {entry['full_max']: entry['count'] for entry in bank['source_strata']}
    if (len(counts) != len(bank['source_strata']) or set(counts) != set(range(R//2+1))
            or sum(counts.values()) != comb(8192, R) or any(type(n) is not int or n <= 0 for n in counts.values())):
        raise ValueError('Complete exact source partition required')
    regions = {(row['occupancy'], row['full_pairs']): np.asarray(
        [-np.inf if x is None else x for x in row['log2_matrix']], float) for row in bank['regions']}
    required = {(m, f) for m in range(33) for f in range(m//2+1)}
    if len(regions) != len(bank['regions']) or set(regions) != required:
        raise ValueError('Complete unique saved regional entries required')
    curves = np.stack([np.max([regions[m, f] for f in range(min(FULL_MAX, m//2)+1)], axis=0)
                       for m in range(R+1)], axis=1)
    slopes = np.asarray(bank['log2_slopes'], float)
    if (curves.shape != (4, R+1) or np.any(np.isnan(curves)) or np.any(np.isposinf(curves))
            or slopes.shape != (4,) or np.any(~np.isfinite(slopes))):
        raise ValueError('Finite-or-zero saved curves and four finite slopes required')
    intercepts = np.max(curves-slopes[:, None]*np.arange(R+1), axis=1)
    if np.any(~np.isfinite(intercepts)):
        raise ValueError('Every entry needs a finite affine upper envelope')
    slack = intercepts[:, None]+slopes[:, None]*np.arange(R+1)-curves
    templates = [row for row in bank['type_templates'] if row['name'].startswith('old-')]
    if len(templates) != 3 or any(row['family'] != 'nonempty-column-floor' or len(row['markers']) != 16
                                  for row in templates):
        raise ValueError('The three frozen activity templates are required')
    markers = [(row['name'], F(row['markers'][ARITY-1])) for row in templates]
    if any(value < 0 for _, value in markers):
        raise ValueError('Activity markers must be at least1')
    return dict(count=counts[FULL_MAX], curves=curves, slopes=slopes, intercepts=intercepts,
                minimum_slack=float(np.min(slack)), markers=markers)


def run(max_seconds=30, progress=None):
    if type(max_seconds) not in (int, float) or not np.isfinite(max_seconds) or not 0 < max_seconds <= 30:
        raise ValueError('An explicit at-most30s budget is required')
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('The fixed-schedule gate exhausted its budget')
    result = dict(construction=CONSTRUCTION, old_control_construction=OLD_CONSTRUCTION,
        active_groups=R, source_full_pair_max=FULL_MAX, source_size=8192,
        homogeneous_constituent_arity=ARITY, all_cohorts=False, all_source_supports=False,
        all_source_supports_in_declared_stratum=True, all_constituent_sets_of_declared_arity=True,
        row_sets_fixed_across_columns=True, incidence_double_counting=True,
        holder_columns=COLUMNS, extra_rowset_factor=False, source_count_outside_alpha=True,
        input_marker='1', mark=str(MARK), alpha=str(ALPHA), forbidden_weight=LIMIT,
        columns=COLUMNS, regions_per_column=GROUP, initial_state=0, terminal_vector=[1, 1],
        regional_envelope='ordinary-four-entry-affine', refined_empty00_retained=False,
        exact_offdiagonal_zero_retained=False, regional_banks_built=0, witness_fits=0,
        fixed_XOR_directions=list(range(1, 16))+[1], schedule_resampled=False,
        whole_word_RS_normalization='255^(-64*alpha*r), inside W then products',
        floor='at least9 nonempty columns per active group; frozen Chernoff activity marker',
        no_endpoint_potential=True, directed=False, whole_code_certificate=False,
        fresh_directed_replay_required=True, compute_budget_seconds=max_seconds,
        complete_requested_controls=False, completed_matching_columns=0, results=[])
    try:
        anchor = load_frozen_anchor()
        result.update(source_supports=anchor['count'], source_bank=masks.FROZEN_BANK,
            source_bank_sha256=masks.FROZEN_SHA, log2_slopes=anchor['slopes'].tolist(),
            log2_intercepts=anchor['intercepts'].tolist(), envelope_minimum_slack=anchor['minimum_slack'])
        new, old, orbit_counts = incidence_weights(check_time=check_time)
        result['input_orbits'] = orbit_counts
        result['log2_matching_weights'] = {name: [[None if np.isneginf(x) else float(x) for x in row]
                                                 for row in table] for name, table in (('new', new), ('old', old))}
        markers = anchor['markers']
        accumulators = {(mode, name): np.array([[0., -np.inf], [-np.inf, 0.]])
                        for mode in ('new', 'old') for name, _ in markers}
        kernel_cache = {}
        schedule = xor_schedule()
        for column, matching in enumerate(schedule):
            check_time()
            for mode, weights in (('new', new), ('old', old)):
                # Old diagonal weights depend only on cardinality, independently
                # verified above; one cached transform covers every matching.
                key = (mode, matching if mode == 'new' else ())
                if key not in kernel_cache:
                    mask_table = mask_weight_table(weights, matching)
                    empty = mask_table[0]
                    mask_table[0] = -np.inf
                    nonempty = path_transform(mask_table, anchor['slopes'].reshape(2, 2))
                    kernel_cache[key] = {name: kernel_from_transform(nonempty, empty, anchor['intercepts'],
                        R, float(value)) for name, value in markers}
                for name, _ in markers:
                    accumulators[mode, name] = log_matmul(accumulators[mode, name], kernel_cache[key][name])
            result['completed_matching_columns'] = column+1
            if progress:
                progress(dict(completed_matching_columns=column+1, elapsed_seconds=perf_counter()-started))
        check_time()
        prefix = log2(anchor['count'])-float(ALPHA)*(LIMIT*log2(float(MARK))+log2(255))
        for name, log_activity in markers:
            values = {mode: float(prefix+np.logaddexp2.reduce(accumulators[mode, name][0])) for mode in ('new', 'old')}
            result['results'].append(dict(template=name, log2_activity=str(log_activity),
                new_failure_log2_screen=values['new'], old_failure_log2_screen=values['old'],
                improvement_bits=values['old']-values['new']))
        result.update(status='completed', complete_requested_controls=True,
                      selected=min(result['results'], key=lambda row: row['new_failure_log2_screen']))
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
        parser.error('Output must be a new artifact')
    result = run(args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    result['dependency_pins'] = {name: sha256((DIRECTORY/name).read_bytes()).hexdigest() for name in (
        'hd_two_sum_pre_toeplitz_fractional.py', 'hd_two_sum_pre_toeplitz_pairs.py',
        'hd_two_sum_toeplitz_pairs.py', 'hd_two_sum_input_weight.py', 'hd_two_sum_closed_kernel.py')}
    with args.output.open('x', encoding='utf-8', newline='\n') as destination:
        destination.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output), sha256=sha256(args.output.read_bytes()).hexdigest(),
                         status=result['status'], results=result['results'], elapsed_seconds=result['elapsed_seconds'])))
