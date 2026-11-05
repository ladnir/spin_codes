"""Directed six-arity replay for one dense group under the g8 affine mixer.

Scope is r17/source(8,0), sixteen two-row groups and one dense group.
The exact affine support kernel and fixed-row-set Holder moments are new;
the scalar positive-polynomial and checked-root routines are reused.
No old XOR support law, saved FLOAT bound or whole-code claim is imported.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_two_sum_affine_dense import ARITY_CLASSES, affine_class_moments
from hd_two_sum_affine_outer import matrix
from hd_two_sum_affine_sparse_certificate import (
    rational, ceiling, require_fields, regional_envelopes, rounded_union_512)
from hd_xor_dense_holder_certificate import root_upper, holder_branch_caps_arb
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_mask_certificate import column_coefficients_arb, bivariate_tail_arb
from hd_xor_pair_route import arb_pair_regions


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_two_sum_affine_dense_8_0_two_frozen.json'
PROPOSAL_SHA256 = '922d920286a6a513af0ae1a98ba41f51a5379c6fa26de2b9c249d5d64a393b9a'
HELPER_PINS = {
    'hd_two_sum_affine_dense.py': '31b387b8886c0a4cf2a07468812871231bcc965dae1cb26bbb0d64b7d5ab4063',
    'hd_two_sum_affine_outer.py': '52d9c3ad412e79d4ac424270734490d3831fc4bbab02e238fedf7ce24fa3d2af',
    'hd_two_sum_affine_sparse_certificate.py': '3fbf95494210ef7ee0aeaf82f7b9682fb386f5a52a3a034d7607adee3f362522',
    'hd_xor_dense_holder_certificate.py': '0934613d5690c280f81d94ecf41f56f25b1a8bd020f532d42b264458a701428e',
    'hd_xor_mixed_mask_certificate.py': 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b',
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
}
MARK = F('8970988817957597/9007199254740992')
LOG_POTENTIAL = F('-2266954815/33554432')
SOURCE_COUNT = 15708305994709246428578365440
EXPECTED_SCOPE = dict(
    construction='two-sum-affine-gf256-h4-alternating-swapped-pairs-k20-rate-half',
    labels=list(range(8, 16)), field_polynomial='0x11b', field_order=256,
    forbidden_weight=734003, active_groups=17, cohort_counts=[0, 16, 1],
    source_matching_counts=[8, 0], source_classes=1, source_supports=SOURCE_COUNT,
    source_colorings=17, sparse_pair_assignments=28**16, pair_swaps=True,
    total_columns=32, sparse_support_floor=544, dense_nonempty_column_floor=17,
    sparse_normalizer_exponent=512, holder_exponent=32, fixed_row_sets=219,
    deterministic_classes=6, dense_type_partition='all-fixed-row-sets-by-arity',
    dense_rowsets_charged_once_inside_holder=True, extra_dense_beta_factor=False,
    extra_dense_rowset_factor=False, fixed_rows_across_columns=True,
    fixed_source_geometry=True, iid_rowset_resampling_assumed=False,
    affine_mask_kernel_rebuilt=True, old_xor_mask_counts_used=False,
    labelled_dense_masks_retained=True, sparse_hole_histogram_exact=True,
    potential_quotient_after_hole_and_dense_mixtures=True, one_global_mark_and_potential=True,
    regional_matrices_are_upper_bounds=True, ordered_regional_matrices=True,
    pair_multiplicity_strata_used=False, column_occupancy_resource_maximum_used=False,
    projective_normalizations=1, endpoint_potential_factors=1, alpha=1,
    mark=str(MARK), log2_potential=str(LOG_POTENTIAL), witness_label='inherited-old-dense-a8',
    witness_banks=1, fitted_new_witnesses=0, saved_float_estimates_used=False,
    exact_positive_bivariate_tail=True, exact_root_endpoint_checks=True,
    directed=True, old_certificates_transferred=False, whole_code_certificate=False,
    complete_selected_cohort_and_source=True, complete_source_coverage=False)


def checked_proposal():
    for name, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed helper changed: {name}')
    raw = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen affine dense proposal changed')
    data = json.loads(raw)['additional']
    extra = {'saved_float_estimates_used', 'exact_positive_bivariate_tail',
             'exact_root_endpoint_checks', 'column_occupancy_resource_maximum_used'}
    require_fields(data, {key: value for key, value in EXPECTED_SCOPE.items()
                          if key not in extra and key != 'directed'})
    require_fields(data, dict(directed=False, resource_maximum_used=False))
    digest = sha256(bytes(value for row in matrix() for value in row)).hexdigest()
    require_fields(data, dict(affine_matrix_sha256=digest))
    if matching_joint_counts(8192, 17)[8, 0] != SOURCE_COUNT:
        raise ArithmeticError('Exact fixed-source count changed')
    expected = dict(ARITY_CLASSES)
    seen = set()
    for row in data['entries']:
        arity = row.get('active_constituents')
        if type(arity) is not int or arity not in expected or arity in seen:
            raise ValueError('The proposal must contain exactly the six arity classes')
        seen.add(arity)
        require_fields(row, dict(active_row_masks=list(expected[arity]), fixed_row_sets=comb(8, arity)))
    if seen != expected.keys() or sum(len(masks) for masks in expected.values()) != 219:
        raise ValueError('Incomplete219-fixed-row-set partition')
    return dict(affine_matrix_sha256=digest)


def class_weights_arb(masks):
    numerators, denominator = affine_class_moments(masks)
    weights = {mask: root_upper(numerator, denominator) for mask, numerator in enumerate(numerators) if numerator}
    digest = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
    return weights, digest


def offset_hole_suffix_arb(regions, potential, occupied, *, check_time=None):
    """All labelled dense masks, preserving the exact ordered hole recurrence."""
    if (type(occupied) is not int or occupied < 0 or not 1 <= len(regions) <= 8
            or any(occupied+1 >= len(table) for table in regions)
            or not potential.is_finite() or not potential > 0):
        raise ValueError('Complete ordered tables including one dense offset required')
    zero = arb(0)
    values = [[(arb(1), potential)]+[(zero, zero)]*occupied]
    binomials = [[comb(k, h) for h in range(k+1)] for k in range(occupied+1)]
    for table in reversed(regions):
        if check_time:
            check_time()
        following = []
        for suffix in values:
            for bit in (0, 1):
                vector = []
                for k in range(occupied+1):
                    left, right = arb(0), arb(0)
                    for holes in range(k+1):
                        a, b, c, d = table[occupied-holes+bit]
                        first, second = suffix[k-holes]
                        factor = binomials[k][holes]
                        left += factor*(a*first+b*second)
                        right += factor*(c*first+d*second)
                    vector.append((left, right))
                following.append(vector)
        values = following
    return values


def dense_mask_vectors_arb(regions, potential, sparse_groups=16, *, q=256, check_time=None, progress=None):
    if (type(sparse_groups) is not int or not 0 <= sparse_groups <= 16
            or q-1 <= len(regions)):
        raise ValueError('Bounded sparse count and positive no-hole mass required')
    result = {}
    for occupied in range(sparse_groups+1):
        values = offset_hole_suffix_arb(regions, potential, occupied, check_time=check_time)
        for double in range(occupied+1):
            factors = [comb(double, k)*(q-1-len(regions))**(double-k) for k in range(double+1)]
            masks = {}
            for mask, vector in enumerate(values):
                left, right = arb(0), arb(0)
                for k, factor in enumerate(factors):
                    left += factor*vector[k][0]
                    right += factor*vector[k][1]
                masks[mask] = left/(q-1)**double, right/(q-1)**double
            result[occupied-double, double] = masks
        if progress:
            progress(dict(stage='directed-dense-hole-bank', completed=occupied+1, total=sparse_groups+1))
    return result


def receipt(precision=256, *, max_seconds=600, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('This replay uses256 or384 bits')
    if type(max_seconds) not in (int, float) or not 0 < max_seconds <= 600:
        raise ValueError('A positive replay budget at most600 seconds is required')
    started = perf_counter()
    metadata = checked_proposal()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Directed dense replay budget exhausted; no complete receipt produced')
    entries = []
    with ctx.workprec(precision):
        p = arb(2)**rational(LOG_POTENTIAL)
        weights = {arity: class_weights_arb(masks) for arity, masks in ARITY_CLASSES}
        check_time()
        raw = arb_pair_regions(MARK, 17)
        envelopes = regional_envelopes(raw, 17)
        tables = [envelopes[(8, 0)[i % 2]] for i in range(8)]
        vectors = dense_mask_vectors_arb(tables, p, check_time=check_time, progress=progress)
        prefix = SOURCE_COUNT*17*28**16*rational(MARK)**(-734003)/(255*p*255**512)
        total = arb(0)
        for arity, masks in ARITY_CLASSES:
            check_time()
            w, moment_digest = weights[arity]
            empty, nonempty = holder_branch_caps_arb(vectors, w, p)
            tail, _ = bivariate_tail_arb(column_coefficients_arb(empty), column_coefficients_arb(nonempty),
                                       columns=32, x_floor=544, y_floor=17)
            value = prefix*tail
            if not value.is_finite() or not value > 0:
                raise ArithmeticError('Finite positive arity bounds required')
            row = dict(active_constituents=arity, active_row_masks=list(masks), fixed_row_sets=len(masks),
                exact_moments_sha256=moment_digest, root_upper_checks_passed=True,
                positive_mask_weights=len(w), failure_log2_upper=str(ceiling(value.log()/arb(2).log())))
            entries.append(row)
            total += value
            if progress:
                progress(dict(stage='directed-arity', **{key: value for key, value in row.items() if key != 'active_row_masks'}))
        check_time()
        aggregate = ceiling(total.log()/arb(2).log())
    return dict(EXPECTED_SCOPE, **metadata, precision_bits=precision,
        frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        entries=entries, failure_log2_upper=str(aggregate), scoped_class_reaches_40_bits=aggregate < -40,
        compute_budget_seconds=max_seconds, elapsed_seconds=perf_counter()-started)


def checked_rows(data, precision):
    require_fields(data, dict(EXPECTED_SCOPE, precision_bits=precision,
        frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(), **checked_proposal()))
    expected, seen = dict(ARITY_CLASSES), {}
    for row in data.get('entries', []):
        arity = row.get('active_constituents')
        if type(arity) is not int or arity not in expected or arity in seen:
            raise ValueError('Unexpected or repeated arity')
        numerators, denominator = affine_class_moments(expected[arity])
        digest = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
        require_fields(row, dict(active_row_masks=list(expected[arity]), fixed_row_sets=comb(8, arity),
            exact_moments_sha256=digest, root_upper_checks_passed=True,
            positive_mask_weights=sum(value > 0 for value in numerators)))
        if type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Exact rational arity ceiling required')
        seen[arity] = F(row['failure_log2_upper'])
    if seen.keys() != expected.keys():
        raise ValueError('Incomplete six-arity receipt')
    return seen


def verify(directory=DIRECTORY):
    metadata = checked_proposal()
    receipts, rows, inputs = [], [], []
    for precision in (256, 384):
        path = Path(directory)/f'hd_two_sum_affine_dense_h1_receipt{precision}.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        receipts.append(data)
        rows.append(checked_rows(data, precision))
        inputs.append(dict(filename=path.name, sha256=sha256(raw).hexdigest(),
                           precision_bits=precision, elapsed_seconds=data['elapsed_seconds']))
    if rows[0] != rows[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('Dual precision arity or aggregate ceilings disagree')
    rounded = rounded_union_512([rows[0][arity] for arity in sorted(rows[0])])
    if rounded < F(receipts[0]['failure_log2_upper']):
        raise ArithmeticError('Rounded-class union lies below the raw ceiling')
    return dict(EXPECTED_SCOPE, **metadata, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS, receipt_precision_bits=[256, 384],
        independent_sum_precision_bits=512, exact_partition_checked=True,
        all_six_arity_ceilings_agree=True, aggregate_ceilings_agree=True,
        raw_aggregate_log2_upper=receipts[0]['failure_log2_upper'],
        rounded_class_union_log2_upper=str(rounded), scoped_class_reaches_40_bits=rounded < -40,
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
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
