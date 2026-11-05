"""Directed g16 h1 replay using exact cardinality Holder moments.

Scope: r17/source(8,0), sixteen two-row groups and one dense group, whose
fixed row set has arity3..16. No saved FLOAT estimate is consumed. The
positive weight DP retains ordered regional matrices and exact hole law.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS, matrix
from hd_two_sum_affine_sparse_certificate import (
    rational, ceiling, require_fields, regional_envelopes, rounded_union_512)
from hd_xor_dense_holder_certificate import root_upper, holder_branch_caps_arb
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_mask_certificate import column_coefficients_arb, bivariate_tail_arb
from hd_xor_pair_route import arb_pair_regions


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_two_sum_gf16_arity_dense_8_0_all14_float.json'
PROPOSAL_SHA256 = 'c9de91df117fea070552fb4a2b5ab75ef1e7ef265c78b9adf6ec69f77f7c4f99'
HELPER_PINS = {
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
    'hd_two_sum_gf16_sparse.py': 'e225f2fccb2f39d95dbfd7c6896c48957c3fd29898f4f574f406a090f6b95c17',
    'hd_two_sum_gf16_arity_dense.py': '539e8c819edc36a351c21e75ca41e91513c0de271d8c23419503ee68c770142a',
    'hd_two_sum_affine_sparse_certificate.py': '3fbf95494210ef7ee0aeaf82f7b9682fb386f5a52a3a034d7607adee3f362522',
    'hd_xor_dense_holder_certificate.py': '0934613d5690c280f81d94ecf41f56f25b1a8bd020f532d42b264458a701428e',
    'hd_xor_mixed_mask_certificate.py': 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b',
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
    'hd_xor_group_geometry.py': '880b51f316341ed4d77d5020694627805cfd708cdf6278a95e142b8ddf5aff17',
}
MARK = F('8970988817957597/9007199254740992')
LOG_POTENTIAL = F('-2266954815/33554432')
SOURCE_COUNT = 15708305994709246428578365440
ARITIES = tuple(range(3, 17))
FROZEN_SCOPE = dict(
    construction=CONSTRUCTION, labels=list(LABELS), field_order=256, field_polynomial='0x11b',
    active_groups=17, cohort_counts=[0, 16, 1], source_matching_counts=[8, 0],
    source_supports=SOURCE_COUNT, source_colorings=17, sparse_pair_assignments=120**16,
    forbidden_weight=734003, total_columns=16, constituent_dimension=8,
    sparse_support_floor=288, dense_nonempty_column_floor=9, sparse_normalizer_exponent=256,
    holder_exponent=16, arities=list(ARITIES), fixed_dense_row_sets=65399,
    deterministic_arity_classes=14, dense_weight_states=17, enumerated_dense_masks=0,
    dense_rowsets_charged_once_inside_holder=True, extra_dense_rowset_factor=False,
    extra_target_mask_binomial=False, extra_dense_beta_factor=False,
    fixed_rows_across_columns=True, fixed_source_geometry=True, iid_rowset_resampling_assumed=False,
    pair_swaps=True, labelled_region_order_preserved=True, hole_mixture_exact=True,
    potential_quotient_after_hole_and_dense_mixtures=True, projective_normalizations=1,
    endpoint_potential_factors=1, mark=str(MARK), log2_potential=str(LOG_POTENTIAL), alpha=1,
    witness_banks=1, fitted_new_witnesses=0, raw_regional_banks=1,
    complete_source_coverage=False, complete_higher_dense_cohort_coverage=False,
    whole_code_certificate=False, old_certificates_transferred=False,
    complete_selected_cohort_and_source=True)
EXPECTED_SCOPE = dict(FROZEN_SCOPE, directed=True, saved_float_estimates_used=False,
    exact_positive_bivariate_tail=True, exact_root_endpoint_checks=True,
    regional_matching_count_envelopes_used=True, column_occupancy_resource_maximum_used=False,
    no_hole_numerator=239, common_global_mark_and_potential=True)


def exact_moments(arity):
    numerators, denominator = arity_holder_moments(arity, 16, 16, 8)
    digest = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
    return numerators, denominator, digest


def checked_proposal():
    for name, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed helper changed: {name}')
    raw = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen g16 arity proposal changed')
    data = json.loads(raw)
    require_fields(data, dict(FROZEN_SCOPE, directed=False, status='complete-selected-cohort-and-source'))
    if matching_joint_counts(8192, 17)[8, 0] != SOURCE_COUNT:
        raise ArithmeticError('The exact source-support count changed')
    seen = set()
    for row in data['entries']:
        arity = row.get('active_constituents')
        if type(arity) is not int or arity not in ARITIES or arity in seen:
            raise ValueError('Distinct arities3..16 required')
        seen.add(arity)
        require_fields(row, dict(fixed_row_sets=comb(16, arity), exact_moments_sha256=exact_moments(arity)[2]))
    if seen != set(ARITIES) or sum(comb(16, a) for a in seen) != 65399:
        raise ValueError('Incomplete fixed-row-set partition')
    return dict(mixer_matrix_sha256=sha256(bytes(x for row in matrix() for x in row)).hexdigest())


def arity_weights_arb(arity):
    if type(arity) is not int or arity not in ARITIES:
        raise ValueError('This replay covers arities3..16')
    numerators, denominator, digest = exact_moments(arity)
    return ({w: root_upper(numerator, denominator, 16)
             for w, numerator in enumerate(numerators) if numerator}, digest)


def weight_vectors_arb(regions, potential, sparse_groups=16, *, q=256,
                       check_time=None, progress=None):
    """Positive coefficients by dense weight, with labelled sparse holes.

    V_i[w,k] = sum_bit,h C(k,h) R_i(t-h+bit) V_(i+1)[w-bit,k-h].
    V_L[0,0]=v; all other terminal entries vanish. Mix the k-hole vectors
    with C(d,k)(q-1-L)^(d-k)/(q-1)^d for each category(t-d,d).
    """
    if (type(sparse_groups) is not int or not 0 <= sparse_groups <= 16
            or type(q) is not int or q < 3 or not 1 <= len(regions) <= 16
            or len(regions) >= q-1 or any(len(table) <= sparse_groups+1 for table in regions)
            or not potential.is_finite() or not potential > 0):
        raise ValueError('Complete ordered tables, positive potential and bounded sparse count required')
    if any(len(matrix) != 4 or any(not x.is_finite() or not x >= 0 for x in matrix)
           for table in regions for matrix in table):
        raise ValueError('Finite nonnegative regional matrices required')
    result = {}
    for occupied in range(sparse_groups+1):
        values = [[(arb(1), potential)]+[(arb(0), arb(0)) for _ in range(occupied)]]
        binomials = [[comb(k, h) for h in range(k+1)] for k in range(occupied+1)]
        for table in reversed(regions):
            if check_time:
                check_time()
            following = [[[arb(0), arb(0)] for _ in range(occupied+1)] for _ in range(len(values)+1)]
            for weight, suffix in enumerate(values):
                for k in range(occupied+1):
                    for holes in range(k+1):
                        first, second = suffix[k-holes]
                        factor = binomials[k][holes]
                        for bit in (0, 1):
                            a, b, c, d = table[occupied-holes+bit]
                            target = following[weight+bit][k]
                            target[0] += factor*(a*first+b*second)
                            target[1] += factor*(c*first+d*second)
            values = following
        for double in range(occupied+1):
            factors = [comb(double, k)*(q-1-len(regions))**(double-k) for k in range(double+1)]
            by_weight = {}
            for weight, vector in enumerate(values):
                left, right = arb(0), arb(0)
                for k, factor in enumerate(factors):
                    left += factor*vector[k][0]
                    right += factor*vector[k][1]
                by_weight[weight] = left/(q-1)**double, right/(q-1)**double
            result[occupied-double, double] = by_weight
        if progress:
            progress(dict(stage='directed-weight-hole-bank', completed=occupied+1, total=sparse_groups+1))
    if check_time:
        check_time()
    return result


def receipt(precision=256, *, max_seconds=600, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh replay precision must be256 or384')
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 600:
        raise ValueError('A finite positive replay budget at most600 seconds is required')
    started = perf_counter()
    metadata = checked_proposal()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Directed g16 replay budget exceeded; no complete receipt produced')
    entries = []
    with ctx.workprec(precision):
        p = arb(2)**rational(LOG_POTENTIAL)
        weights = {a: arity_weights_arb(a) for a in ARITIES}
        check_time()
        raw = arb_pair_regions(MARK, 17)
        envelopes = regional_envelopes(raw, 17)
        regions = [envelopes[(8, 0)[i % 2]] for i in range(16)]
        vectors = weight_vectors_arb(regions, p, check_time=check_time, progress=progress)
        prefix = SOURCE_COUNT*17*120**16*rational(MARK)**(-734003)/(255*p*255**256)
        total = arb(0)
        for arity in ARITIES:
            check_time()
            w, moment_digest = weights[arity]
            empty, nonempty = holder_branch_caps_arb(vectors, w, p)
            tail, _ = bivariate_tail_arb(column_coefficients_arb(empty), column_coefficients_arb(nonempty),
                                       columns=16, x_floor=288, y_floor=9)
            value = prefix*tail
            if not value.is_finite() or not value > 0:
                raise ArithmeticError('Finite positive arity bounds required')
            row = dict(active_constituents=arity, fixed_row_sets=comb(16, arity),
                exact_moments_sha256=moment_digest, root_upper_checks_passed=True,
                positive_weight_states=len(w), failure_log2_upper=str(ceiling(value.log()/arb(2).log())))
            entries.append(row)
            total += value
            if progress:
                progress(dict(stage='directed-arity', **row))
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
    seen = {}
    for row in data.get('entries', []):
        arity = row.get('active_constituents')
        if type(arity) is not int or arity not in ARITIES or arity in seen:
            raise ValueError('Unexpected or repeated arity')
        numerators, _, digest = exact_moments(arity)
        require_fields(row, dict(fixed_row_sets=comb(16, arity), exact_moments_sha256=digest,
            root_upper_checks_passed=True, positive_weight_states=sum(n > 0 for n in numerators)))
        if type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Exact rational arity ceiling required')
        seen[arity] = F(row['failure_log2_upper'])
    if set(seen) != set(ARITIES):
        raise ValueError('Incomplete14-arity receipt')
    return seen


def verify(directory=DIRECTORY):
    metadata = checked_proposal()
    receipts, rows, inputs = [], [], []
    for precision in (256, 384):
        path = Path(directory)/f'hd_two_sum_gf16_arity_h1_receipt{precision}.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        receipts.append(data)
        rows.append(checked_rows(data, precision))
        inputs.append(dict(filename=path.name, sha256=sha256(raw).hexdigest(),
            precision_bits=precision, elapsed_seconds=data['elapsed_seconds']))
    if rows[0] != rows[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('Dual precision arity or aggregate ceilings disagree')
    rounded = rounded_union_512([rows[0][a] for a in ARITIES])
    if rounded < F(receipts[0]['failure_log2_upper']):
        raise ArithmeticError('Rounded-class union lies below the raw ceiling')
    return dict(EXPECTED_SCOPE, **metadata, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS, receipt_precision_bits=[256, 384],
        independent_sum_precision_bits=512, exact_partition_checked=True,
        all14_arity_ceilings_agree=True, aggregate_ceilings_agree=True,
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
