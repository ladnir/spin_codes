"""Fresh directed replay of the frozen affine g8 sparse r17 family.

All81 ordered source classes and all18 one/two-row cohort splits are
included. Exact positive hole mixtures preserve regional matrix order.
The original row sets remain fixed through32 columns. No saved FLOAT
estimate, old XOR support kernel, or old certificate enters the arithmetic.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

from hd_xor_group_geometry import matching_joint_counts
from hd_xor_pair_route import arb_pair_regions, feasible_full_pairs


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_two_sum_affine_sparse_r17_all81_inherited_float.json'
PROPOSAL_SHA256 = '234901cce79ea03cc21f1fe8da1b528d35cd7eef15d4052304c157dbffd7bffc'
HELPER_PINS = {
    'hd_two_sum_affine_outer.py': '52d9c3ad412e79d4ac424270734490d3831fc4bbab02e238fedf7ce24fa3d2af',
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
    'hd_shifted_swap_sparse.py': '9b8875ec1433a2bc79dbd1e42d3b90815528a70ce0234069cfd03d61b030c58c',
    'hd_xor_multigroup.py': 'e64a98c2db2a04d34ee61a7814825a85fba6724d281e0963215c6e59442d1bd0',
    'hd_xor_group_geometry.py': '880b51f316341ed4d77d5020694627805cfd708cdf6278a95e142b8ddf5aff17',
}
MARK, LOG_POTENTIAL = F('996248217/1000000000'), F('-127/2')
LIMIT, GROUPS, LENGTH = 734003, 17, 8192
EXPECTED_SCOPE = dict(
    construction='two-sum-affine-gf256-h4-alternating-swapped-pairs-k20-rate-half',
    labels=list(range(8, 16)), field_order=256, field_polynomial='0x11b',
    forbidden_weight=LIMIT, active_groups=GROUPS, total_columns=32,
    source_classes=81, source_supports=comb(LENGTH, GROUPS), cohort_count=18,
    cohort_order=['a1', 'a2', 'a3-through-a8'], cells=1458,
    mark=str(MARK), log2_potential=str(LOG_POTENTIAL), witness_banks=1,
    pair_swaps=True, fixed_source_geometry=True, fixed_rows_across_columns=True,
    local_mask_mixture_exact=True, ordered_regional_matrices=True,
    regional_matching_count_envelopes_used=True,
    column_occupancy_resource_maximum_used=False, pair_multiplicity_strata_used=False,
    potential_quotient_after_hole_mixture=True, one_common_potential=True,
    projective_normalizations=1, endpoint_potential_factors=1,
    raw_regional_banks=1, distinct_regional_full_count_envelopes=9,
    source_conditioned_cap_banks=81, positive_coefficient_power=32,
    alpha=1, fitted_new_witnesses=0, saved_float_estimates_used=False,
    directed=True, complete_source_coverage=True, complete_sparse_cohort_coverage=True,
    complete_dense_cohort_coverage=False, old_certificates_transferred=False,
    whole_code_certificate=False)


def rational(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def ceiling(value):
    return F(int((value.upper()*10**12).upper().ceil().fmpz()), 10**12)


def require_fields(data, expected):
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError(f'Wrong frozen scope field: {key}')


def checked_proposal():
    for name, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed helper changed: {name}')
    raw = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen all81 proposal changed')
    proposal = json.loads(raw)
    keys = ('construction', 'labels', 'field_order', 'field_polynomial', 'forbidden_weight',
        'active_groups', 'total_columns', 'source_classes', 'cohort_count', 'cohort_order',
        'mark', 'log2_potential', 'witness_banks', 'pair_swaps', 'fixed_source_geometry',
        'fixed_rows_across_columns', 'local_mask_mixture_exact', 'ordered_regional_matrices',
        'regional_matching_count_envelopes_used', 'column_occupancy_resource_maximum_used',
        'pair_multiplicity_strata_used', 'potential_quotient_after_hole_mixture',
        'projective_normalizations', 'endpoint_potential_factors', 'raw_regional_banks',
        'fitted_new_witnesses', 'complete_source_coverage', 'complete_sparse_cohort_coverage',
        'old_certificates_transferred', 'whole_code_certificate')
    require_fields(proposal, {key: EXPECTED_SCOPE[key] for key in keys})
    require_fields(proposal, dict(directed=False, budget_exhausted=False, selected_cells=1458,
        source_count_partition=comb(LENGTH, GROUPS), completed_source_supports=comb(LENGTH, GROUPS)))
    counts = matching_joint_counts(LENGTH, GROUPS)
    if (set(counts) != {(a, b) for a in range(9) for b in range(9)}
            or any(type(count) is not int or count <= 0 for count in counts.values())
            or sum(counts.values()) != comb(LENGTH, GROUPS)):
        raise ArithmeticError('Exact source counts must partition C(8192,17)')
    seen = set()
    for source in proposal['sources']:
        key = tuple(source['source_matching_counts'])
        if key not in counts or key in seen:
            raise ValueError('Duplicate or unknown source class')
        seen.add(key)
        require_fields(source, dict(source_supports=counts[key]))
        if len(source['entries']) != 18:
            raise ValueError('Incomplete sparse cohort partition')
        for n1, row in enumerate(source['entries']):
            require_fields(row, cell_scope(key, counts[key], n1, include_source=False))
    if seen != counts.keys():
        raise ValueError('Incomplete ordered source partition')
    return counts


def cell_scope(source, source_count, n1, *, include_source=True):
    n2 = GROUPS-n1
    result = dict(cohort_counts=[n1, n2, 0], source_supports=source_count,
        source_colorings=comb(GROUPS, n1), singleton_row_choices=8**n1,
        pair_row_choices=28**n2, sparse_support_floor=17*(n1+2*n2),
        sparse_normalizer_exponent=16*(n1+2*n2))
    if include_source:
        result['source_matching_counts'] = list(source)
    return result


def regional_envelopes(raw, maximum):
    """Reuse each fixed-source full-count envelope; upper endpoints only."""
    return tuple(tuple(tuple(max(raw[m, f][entry].upper()
        for f in feasible_full_pairs(m, maximum, full)) for entry in range(4))
        for m in range(maximum+1)) for full in range(maximum//2+1))


def hole_vectors_arb(regions, potential, total, *, q=256, check_time=None):
    """Return all Q_(total-d,d)v, sharing the positive suffix recurrence.

    V_i[k]=sum_h C(k,h) R_i(total-h)V_(i+1)[k-h]. Initially
    V_L[0]=v and V_L[k>0]=0. At the end mix C(d,k)(q-1-L)^(d-k)
    V_0[k]/(q-1)^d. Each V[k] is reused for all d>=k, with no max.
    """
    if (type(total) is not int or total < 0 or not 1 <= len(regions) < q-1
            or any(len(table) <= total for table in regions)
            or not potential.is_finite() or not potential > 0):
        raise ValueError('Complete ordered tables and a positive common potential required')
    zero = arb(0)
    values = [(arb(1), potential)]+[(zero, zero)]*total
    binomials = [[comb(k, h) for h in range(k+1)] for k in range(total+1)]
    for table in reversed(regions):
        if check_time:
            check_time()
        following = []
        for k in range(total+1):
            left, right = arb(0), arb(0)
            for h in range(k+1):
                a, b, c, d = table[total-h]
                first, second = values[k-h]
                factor = binomials[k][h]
                left += factor*(a*first+b*second)
                right += factor*(c*first+d*second)
            following.append((left, right))
        values = following
    result = []
    for double in range(total+1):
        left, right = arb(0), arb(0)
        for holes in range(double+1):
            factor = comb(double, holes)*(q-1-len(regions))**(double-holes)
            left += factor*values[holes][0]
            right += factor*values[holes][1]
        result.append((left/(q-1)**double, right/(q-1)**double))
    return tuple(result)


def category_caps_arb(regions, potential, maximum, *, check_time=None):
    caps = {}
    for total in range(maximum+1):
        for double, (first, second) in enumerate(hole_vectors_arb(
                regions, potential, total, check_time=check_time)):
            upper = max(first.upper(), (second/potential).upper())
            if not upper.is_finite() or not upper > 0:
                raise ArithmeticError('Invalid positive category cap')
            caps[total-double, double] = upper
    return caps


def sparse_coefficients_arb(caps, n1, n2, *, q=256):
    if any(type(n) is not int or n < 0 for n in (n1, n2)) or q < 2:
        raise ValueError('Nonnegative cohort sizes required')
    result = [arb(0) for _ in range(n1+2*n2+1)]
    for active_single in range(n1+1):
        for partial in range(n2+1):
            for double in range(n2-partial+1):
                degree = active_single+partial+2*double
                factor = (comb(n1, active_single)*comb(n2, partial)
                          *comb(n2-partial, double)*2**partial*(q-1)**degree)
                result[degree] += factor*caps[active_single+partial, double]
    return tuple(result)


def sparse_tail_arb(caps, n1, n2, *, columns=32, dimension=16, q=256):
    if not 0 < dimension <= columns:
        raise ValueError('Valid constituent dimensions required')
    rows = n1+2*n2
    polynomial = arb_poly(list(sparse_coefficients_arb(caps, n1, n2, q=q)))**columns
    floor = rows*(columns-dimension+1)
    return sum((polynomial[j] for j in range(floor, rows*columns+1)), arb(0))/(q-1)**(rows*(columns-dimension))


def receipt(precision=256, *, max_seconds=600, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('This replay uses exactly256 or384 bits')
    if (type(max_seconds) not in (int, float) or not 0 < max_seconds <= 600):
        raise ValueError('The replay budget must be positive and at most600 seconds')
    started = perf_counter()
    counts = checked_proposal()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Directed replay budget exhausted; no complete receipt produced')
    entries, cohorts, costs = [], [], []
    with ctx.workprec(precision):
        p = arb(2)**rational(LOG_POTENTIAL)
        prefix = rational(MARK)**(-LIMIT)/(255*p)
        raw = arb_pair_regions(MARK, GROUPS)
        tables = regional_envelopes(raw, GROUPS)
        total, cohort_values = arb(0), [arb(0) for _ in range(18)]
        for source, count in sorted(counts.items()):
            check_time()
            source_start = perf_counter()
            regions = [tables[source[i % 2]] for i in range(8)]
            caps = category_caps_arb(regions, p, GROUPS, check_time=check_time)
            source_total = arb(0)
            for n1 in range(18):
                check_time()
                n2 = GROUPS-n1
                value = (prefix*count*comb(GROUPS, n1)*8**n1*28**n2
                         *sparse_tail_arb(caps, n1, n2))
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError('Finite positive cell bounds required')
                entries.append(dict(cell_scope(source, count, n1),
                                    failure_log2_upper=str(ceiling(value.log()/arb(2).log()))))
                source_total += value
                cohort_values[n1] += value
                total += value
            cost = dict(source_matching_counts=list(source), elapsed_seconds=perf_counter()-source_start,
                        failure_log2_upper=str(ceiling(source_total.log()/arb(2).log())))
            costs.append(cost)
            if progress:
                progress(dict(stage='directed-source', completed=len(costs), total=81,
                              **cost, total_elapsed_seconds=perf_counter()-started))
        check_time()
        cohorts = [dict(cohort_counts=[n1, GROUPS-n1, 0],
            failure_log2_upper=str(ceiling(value.log()/arb(2).log())))
            for n1, value in enumerate(cohort_values)]
        aggregate = ceiling(total.log()/arb(2).log())
    return dict(EXPECTED_SCOPE, precision_bits=precision, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS, entries=entries, cohorts=cohorts, source_costs=costs,
        failure_log2_upper=str(aggregate), scoped_class_reaches_40_bits=aggregate < -40,
        compute_budget_seconds=max_seconds, elapsed_seconds=perf_counter()-started)


def checked_rows(data, precision, counts):
    require_fields(data, dict(EXPECTED_SCOPE, precision_bits=precision,
        frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest()))
    found = {}
    for row in data.get('entries', []):
        source = tuple(row.get('source_matching_counts', ()))
        cohort = row.get('cohort_counts', [])
        if (source not in counts or len(cohort) != 3 or type(cohort[0]) is not int
                or not 0 <= cohort[0] <= GROUPS):
            raise ValueError('Unknown source/cohort cell')
        n1 = cohort[0]
        key = (*source, n1)
        if key in found:
            raise ValueError('Duplicate source/cohort cell')
        require_fields(row, cell_scope(source, counts[source], n1))
        if type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Exact rational cell ceiling required')
        found[key] = F(row['failure_log2_upper'])
    if set(found) != {(*source, n1) for source in counts for n1 in range(18)}:
        raise ValueError('Incomplete1458-cell receipt')
    return found


def rounded_union_512(values):
    with ctx.workprec(512):
        total = sum((arb(2)**rational(value) for value in values), arb(0))
        if not total.is_finite() or not total > 0:
            raise ValueError('Finite nonempty union required')
        return ceiling(total.log()/arb(2).log())


def verify(directory=DIRECTORY):
    counts = checked_proposal()
    receipts, rows, inputs = [], [], []
    for precision in (256, 384):
        path = Path(directory)/f'hd_two_sum_affine_sparse_r17_receipt{precision}.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        receipts.append(data)
        rows.append(checked_rows(data, precision, counts))
        inputs.append(dict(filename=path.name, sha256=sha256(raw).hexdigest(),
                           precision_bits=precision, elapsed_seconds=data['elapsed_seconds']))
    if (rows[0] != rows[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']
            or receipts[0]['cohorts'] != receipts[1]['cohorts']):
        raise ValueError('Dual precision cell, cohort or aggregate ceilings disagree')
    rounded = rounded_union_512([rows[0][key] for key in sorted(rows[0])])
    if rounded < F(receipts[0]['failure_log2_upper']):
        raise ArithmeticError('Rounded cell union lies below the raw union ceiling')
    return dict(EXPECTED_SCOPE, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS, receipt_precision_bits=[256, 384],
        independent_sum_precision_bits=512, exact_partition_checked=True,
        all1458_cell_ceilings_agree=True, all18_cohort_ceilings_agree=True,
        aggregate_ceilings_agree=True, raw_aggregate_log2_upper=receipts[0]['failure_log2_upper'],
        rounded_cell_union_log2_upper=str(rounded), scoped_class_reaches_40_bits=rounded < -40,
        inputs=inputs, receipt_arithmetic_trust='observed-fresh-directed-replays-not-proved-by-byte-hashes')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--precision', type=int, choices=(256, 384))
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--max-seconds', type=float, default=600)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen directed artifact')
    result = verify(args.directory) if args.verify else receipt(args.precision, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                     if key not in ('entries', 'cohorts', 'source_costs')}, allow_nan=False), flush=True)
