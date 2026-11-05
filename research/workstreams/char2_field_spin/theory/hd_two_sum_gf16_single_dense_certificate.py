"""Directed r1 dense replay retaining ordered2x2 matrix coefficients.

The fourteen fixed arities use three frozen rational marks. There is no
common-potential quotient or endpoint-potential factor. Positive matrix
polynomials retain at least9 nonempty columns out of16.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_two_sum_gf16_arity_certificate import exact_moments, arity_weights_arb
from hd_two_sum_affine_sparse_certificate import rational, ceiling, require_fields, rounded_union_512
from hd_xor_pair_route import arb_pair_regions


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_two_sum_gf16_single_dense_selected_float.json'
PROPOSAL_SHA256 = 'f0bd63b95bed054500796427cc1a7f9f07c96afce911b137ad85f39a0ee44a15'
SOURCE_PINS = {
    'hd_two_sum_gf16_single_dense_matrix_two_old_marks_float.json': 'ec16d6fe548e8bca9e5486663acbb60013bac5e8edf66888eb24dd09e4e823f8',
    'hd_two_sum_gf16_single_dense_four_rational_marks_float.json': 'c3f33706ddd4896c68676e370ecef9cfb6f5a85d3a71759caddb5753abcbe9ad',
}
HELPER_PINS = {
    'hd_two_sum_gf16_single_dense_matrix.py': '5b395eef84e0760c56fd22bc3d3a7ae9bb2e9d060ed9d767517b2e4c3f13b9a0',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
    'hd_two_sum_gf16_sparse.py': 'e225f2fccb2f39d95dbfd7c6896c48957c3fd29898f4f574f406a090f6b95c17',
    'hd_two_sum_gf16_arity_certificate.py': 'c47e378f680a290c9b8c2c0fa5e852655aabc22acda47ae6ed9ac571c7eb466b',
    'hd_two_sum_affine_sparse_certificate.py': '3fbf95494210ef7ee0aeaf82f7b9682fb386f5a52a3a034d7607adee3f362522',
    'hd_xor_dense_holder_certificate.py': '0934613d5690c280f81d94ecf41f56f25b1a8bd020f532d42b264458a701428e',
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
}
ARITIES = tuple(range(3, 17))
MARKS = (F(1999, 2000), F(1199, 1200), F(499451537, 500000000))
SCHEDULE = {a: MARKS[0 if a <= 5 else 1 if a <= 10 else 2] for a in ARITIES}
FROZEN_SCOPE = dict(construction=CONSTRUCTION, field_order=256, field_polynomial='0x11b',
    labels=list(LABELS), active_groups=1, cohort_counts=[0, 0, 1], all_source_supports=True,
    source_supports=8192, source_colorings=1, forbidden_weight=734003,
    total_columns=16, constituent_dimension=8, dense_nonempty_column_floor=9,
    holder_exponent=16, fixed_dense_row_sets=65399, arities_kept_separate=True,
    dense_rowset_multiplicity_inside_holder=True, extra_dense_rowset_factor=False,
    fixed_rows_across_columns=True, pair_swaps=True, initial_state=0, terminal_vector=[1, 1],
    common_potential_used=False, endpoint_potential_factors=0, projective_normalizations=1,
    alpha=1, mark_selection_per_arity=True, extra_mark_bank_factor=False, whole_code_certificate=False)
EXPECTED_SCOPE = dict(FROZEN_SCOPE, directed=True, arities=list(ARITIES),
    source_full_pairs=0, sparse_groups=0, fixed_source_geometry=True,
    ordered_regional_matrix_products=True, ordered_column_matrix_polynomial=True,
    selected_marks=[str(z) for z in MARKS], rebuilt_regional_banks=3,
    fitted_new_witnesses=0, saved_float_estimates_used=False,
    exact_root_endpoint_checks=True, positive_matrix_coefficient_tail=True,
    complete_selected_cohort=True, complete_all_active_group_coverage=False,
    old_certificates_transferred=False)


def checked_proposal():
    for name, digest in {**SOURCE_PINS, **HELPER_PINS}.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Pinned source changed: {name}')
    raw = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen selected r1 proposal changed')
    data = json.loads(raw)
    require_fields(data, dict(FROZEN_SCOPE, directed=False, proposal_requires_fresh_directed_replay=True,
        selection_used_saved_rows_only=True, additional_numerical_evaluations=0,
        matrix_helper_sha256=HELPER_PINS['hd_two_sum_gf16_single_dense_matrix.py'],
        source_artifacts=[dict(name=name, sha256=digest) for name, digest in SOURCE_PINS.items()]))
    seen = set()
    for row in data['entries']:
        a = row.get('active_constituents')
        if type(a) is not int or a not in ARITIES or a in seen:
            raise ValueError('Distinct arities3..16 required')
        seen.add(a)
        source = list(SOURCE_PINS)[1 if a <= 5 else 0]
        require_fields(row, dict(mark=str(SCHEDULE[a]), fixed_row_sets=comb(16, a),
            exact_moments_sha256=exact_moments(a)[2], source_artifact=source))
    if seen != set(ARITIES) or sum(comb(16, a) for a in seen) != 65399:
        raise ValueError('Incomplete r1 dense partition')
    return dict(SCHEDULE)


def _valid_matrix(matrix):
    return len(matrix) == 4 and all(x.is_finite() and x >= 0 for x in matrix)


def matrix_product(left, right):
    """Row-major2x2 product; no commutation or scalar binomial expansion."""
    a, b, c, d = left
    e, f, g, h = right
    return a*e+b*g, a*f+b*h, c*e+d*g, c*f+d*h


def _add(left, right):
    return left[0]+right[0], left[1]+right[1], left[2]+right[2], left[3]+right[3]


def weight_matrices_arb(regions):
    """Q_w sums every labelled ordered regional product of mask weightw."""
    if (not 1 <= len(regions) <= 16 or any(len(table) != 2 for table in regions)
            or any(not _valid_matrix(matrix) for table in regions for matrix in table)):
        raise ValueError('Complete nonnegative occupancy0/1 regional matrices required')
    zero, identity = (arb(0),)*4, (arb(1), arb(0), arb(0), arb(1))
    values = [identity]
    for table in regions:
        following = [zero]*(len(values)+1)
        for w, matrix in enumerate(values):
            following[w] = _add(following[w], matrix_product(matrix, table[0]))
            following[w+1] = _add(following[w+1], matrix_product(matrix, table[1]))
        values = following
    return tuple(values)


def matrix_branches_arb(matrices, weights):
    if (len(matrices) < 2 or 0 not in weights or not set(weights) <= set(range(len(matrices)))
            or any(not _valid_matrix(matrix) for matrix in matrices)
            or any(not w.is_finite() or not w >= 0 for w in weights.values())):
        raise ValueError('Complete positive matrix coefficients and nonnegative Holder weights required')
    empty = tuple(x*weights[0] for x in matrices[0])
    nonempty = (arb(0),)*4
    for w in sorted(set(weights)-{0}):
        nonempty = _add(nonempty, tuple(x*weights[w] for x in matrices[w]))
    return empty, nonempty


def matrix_tail_arb(branches, columns=16, y_floor=9):
    """Positive matrix polynomial, initialzero row summed over final states."""
    if (len(branches) != 2 or any(not _valid_matrix(matrix) for matrix in branches)
            or type(columns) is not int or columns < 1 or type(y_floor) is not int
            or not 0 <= y_floor <= columns):
        raise ValueError('Nonnegative matrix branches and valid column/floor bounds required')
    zero, identity = (arb(0),)*4, (arb(1), arb(0), arb(0), arb(1))
    values = [identity]
    for _ in range(columns):
        following = [zero]*(len(values)+1)
        for degree, matrix in enumerate(values):
            following[degree] = _add(following[degree], matrix_product(matrix, branches[0]))
            following[degree+1] = _add(following[degree+1], matrix_product(matrix, branches[1]))
        values = following
    by_degree = tuple(matrix[0]+matrix[1] for matrix in values[y_floor:])
    return sum(by_degree, arb(0)), by_degree


def receipt(precision=256, *, max_seconds=120, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh replay precision must be256 or384')
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 120:
        raise ValueError('Finite positive replay budget at most120 seconds required')
    started = perf_counter()
    selected = checked_proposal()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Directed r1 matrix budget exceeded; no complete receipt produced')
    entries, costs = [], []
    with ctx.workprec(precision):
        weights = {a: arity_weights_arb(a) for a in ARITIES}
        total = arb(0)
        for mark in MARKS:
            check_time()
            bank_start = perf_counter()
            raw = arb_pair_regions(mark, 1)
            matrices = weight_matrices_arb([(raw[0, 0], raw[1, 0])]*16)
            prefix = 8192*rational(mark)**(-734003)/255
            selected_arities = [a for a in ARITIES if selected[a] == mark]
            for a in selected_arities:
                check_time()
                w, moment_digest = weights[a]
                tail, _ = matrix_tail_arb(matrix_branches_arb(matrices, w))
                value = prefix*tail
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError('Finite positive arity bound required')
                row = dict(active_constituents=a, mark=str(mark), fixed_row_sets=comb(16, a),
                    exact_moments_sha256=moment_digest, root_upper_checks_passed=True,
                    positive_weight_states=len(w), failure_log2_upper=str(ceiling(value.log()/arb(2).log())))
                entries.append(row)
                total += value
                if progress:
                    progress(dict(stage='directed-r1-arity', **row))
            costs.append(dict(mark=str(mark), arities=selected_arities, elapsed_seconds=perf_counter()-bank_start))
        check_time()
        aggregate = ceiling(total.log()/arb(2).log())
    return dict(EXPECTED_SCOPE, precision_bits=precision, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS, proposal_source_sha256=SOURCE_PINS,
        entries=entries, bank_costs=costs, failure_log2_upper=str(aggregate),
        scoped_class_reaches_40_bits=aggregate < -40, compute_budget_seconds=max_seconds,
        elapsed_seconds=perf_counter()-started)


def checked_rows(data, precision):
    checked_proposal()
    require_fields(data, dict(EXPECTED_SCOPE, precision_bits=precision,
        frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        proposal_source_sha256=SOURCE_PINS, source_sha256=sha256(Path(__file__).read_bytes()).hexdigest()))
    seen = {}
    for row in data.get('entries', []):
        a = row.get('active_constituents')
        if type(a) is not int or a not in ARITIES or a in seen:
            raise ValueError('Unexpected or repeated arity')
        nums, _, digest = exact_moments(a)
        require_fields(row, dict(mark=str(SCHEDULE[a]), fixed_row_sets=comb(16, a),
            exact_moments_sha256=digest, root_upper_checks_passed=True,
            positive_weight_states=sum(n > 0 for n in nums)))
        if type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Exact rational arity ceiling required')
        seen[a] = F(row['failure_log2_upper'])
    if set(seen) != set(ARITIES):
        raise ValueError('Incomplete14-arity receipt')
    return seen


def verify(directory=DIRECTORY):
    receipts, rows, inputs = [], [], []
    for precision in (256, 384):
        path = Path(directory)/f'hd_two_sum_gf16_single_dense_receipt{precision}.json'
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
        raise ArithmeticError('Rounded-arity union lies below raw ceiling')
    return dict(EXPECTED_SCOPE, frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        proposal_source_sha256=SOURCE_PINS, receipt_precision_bits=[256, 384], independent_sum_precision_bits=512,
        all14_arity_ceilings_agree=True, aggregate_ceilings_agree=True,
        raw_aggregate_log2_upper=receipts[0]['failure_log2_upper'],
        rounded_arity_union_log2_upper=str(rounded), scoped_class_reaches_40_bits=rounded < -40,
        inputs=inputs, receipt_arithmetic_trust='observed-fresh-directed-replays-not-proved-by-byte-hashes')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--precision', type=int, choices=(256, 384))
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--max-seconds', type=float, default=120)
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
