"""Directed all-source g16 replay of322 fixed small-cohort cells.

Covers r2..17 with zero/one dense group and the two r1 sparse cells.
The excluded r1 dense cell is never imported as coverage. Frozen witness
selection is per whole cohort cell, not per column or sampled setup.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_two_sum_closed_kernel import arity_holder_moments
from hd_two_sum_gf16_sparse import CONSTRUCTION, LABELS
from hd_two_sum_gf16_arity_certificate import weight_vectors_arb
from hd_two_sum_affine_sparse_certificate import (
    rational, ceiling, require_fields, category_caps_arb, sparse_coefficients_arb,
    sparse_tail_arb, rounded_union_512)
from hd_xor_dense_holder_certificate import root_upper, holder_branch_caps_arb
from hd_xor_mixed_mask_certificate import bivariate_tail_arb
from hd_xor_pair_route import arb_pair_regions


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_two_sum_gf16_small_cover_r1_17_allsource_float.json'
PROPOSAL_SHA256 = '8bbc4c4cb61234fa2cd6f0c327117926fcf7ed112d4bfd8b71687c4a5f16724f'
HELPER_PINS = {
    'hd_two_sum_gf16_small_cover.py': '3f6ade5acce2e005dbf65f8ecfbbb0fe7fe527c09eddfcca2581eb38595c52cd',
    'hd_two_sum_closed_kernel.py': '2ee0fe57de750a19466155317d595589695252de450518b436eab6c5c0fbf979',
    'hd_two_sum_gf16_sparse.py': 'e225f2fccb2f39d95dbfd7c6896c48957c3fd29898f4f574f406a090f6b95c17',
    'hd_two_sum_gf16_arity_certificate.py': 'c47e378f680a290c9b8c2c0fa5e852655aabc22acda47ae6ed9ac571c7eb466b',
    'hd_two_sum_affine_sparse_certificate.py': '3fbf95494210ef7ee0aeaf82f7b9682fb386f5a52a3a034d7607adee3f362522',
    'hd_xor_dense_holder_certificate.py': '0934613d5690c280f81d94ecf41f56f25b1a8bd020f532d42b264458a701428e',
    'hd_xor_mixed_mask_certificate.py': 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b',
    'hd_xor_pair_route.py': 'dec1f2e2ba064f60529c31675a3779ef53dc6021035bcd024366e5614b3861bc',
}
WITNESSES = (
    ('sparse-inherited', F('996248217/1000000000'), F('-127/2')),
    ('dense-inherited', F('8970988817957597/9007199254740992'), F('-2266954815/33554432')),
    ('h2-inherited', F('1116892707587883/1125899906842624'), F(-40)),
)
LITERALS = {name: (z, logp) for name, z, logp in WITNESSES}
ALL_KEYS = tuple((r, n1, r-dense-n1, dense) for r in range(1, 18)
                 for dense in (0, 1) for n1 in range(r-dense+1))
EXCLUDED_KEY = (1, 0, 0, 1)
KEYS = tuple(key for key in ALL_KEYS if key != EXCLUDED_KEY)
EXPECTED_SCOPE = dict(
    construction=CONSTRUCTION, labels=list(LABELS), field_order=256, field_polynomial='0x11b',
    forbidden_weight=734003, total_columns=16, constituent_dimension=8, holder_exponent=16,
    active_group_range=[1, 17], zero_dense_group_range=[1, 17], one_dense_group_range=[2, 17],
    excluded_cells=[dict(active_groups=1, cohort_counts=[0, 0, 1])], selected_cells=322,
    all_source_supports=True, source_geometry_envelope=True, pair_swaps=True,
    fixed_rows_across_columns=True, fixed_source_geometry=True,
    all_dense_rowsets_inside_holder=True, fixed_dense_rowsets=65399, merged_arity_moments=True,
    extra_dense_arity_factor=False, extra_dense_rowset_factor=False, extra_target_mask_binomial=False,
    extra_dense_beta_factor=False, iid_rowset_resampling_assumed=False,
    labelled_region_order_preserved=True, hole_mixture_exact=True,
    potential_quotient_after_hole_and_dense_mixtures=True,
    column_occupancy_resource_maximum_used=False, projective_normalizations=1,
    endpoint_potential_factors=1, alpha=1, fitted_new_witnesses=0,
    witness_banks_available=3, rebuilt_selected_banks=2, global_witness_per_cell=True,
    saved_float_estimates_used=False, exact_root_endpoint_checks=True,
    exact_positive_coefficient_tails=True, directed=True,
    complete_selected_cohort_coverage=True, complete_all_cohort_coverage=False,
    whole_code_certificate=False, old_certificates_transferred=False)


def cell_scope(key):
    r, n1, n2, dense = key
    return dict(active_groups=r, cohort_counts=[n1, n2, dense], source_supports=comb(8192, r),
        source_colorings=comb(r, dense)*comb(r-dense, n1), singleton_row_choices=16**n1,
        pair_row_choices=120**n2, sparse_support_floor=9*(n1+2*n2),
        sparse_normalizer_exponent=8*(n1+2*n2), dense_nonempty_column_floor=9*dense)


def row_key(row):
    r, cohort = row.get('active_groups'), row.get('cohort_counts')
    if type(r) is not int or type(cohort) is not list or len(cohort) != 3 or any(type(n) is not int for n in cohort):
        raise ValueError('Exact integer cohort identity required')
    return (r, *cohort)


def check_proposal_rows(rows, *, selected=False):
    seen = {}
    for row in rows:
        key = row_key(row)
        if key not in ALL_KEYS or key in seen:
            raise ValueError('Unexpected or repeated frozen cell')
        fields = cell_scope(key)
        require_fields(row, {name: value for name, value in fields.items()
            if name not in ('singleton_row_choices', 'pair_row_choices', 'sparse_normalizer_exponent')})
        if selected and row.get('witness') not in LITERALS:
            raise ValueError('Unknown frozen witness')
        seen[key] = row.get('witness') if selected else None
    if set(seen) != set(ALL_KEYS):
        raise ValueError('Incomplete323-cell frozen proposal')
    return seen


def checked_selection():
    for name, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed helper changed: {name}')
    raw = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen small-cover proposal changed')
    data = json.loads(raw)
    require_fields(data, dict(construction=CONSTRUCTION, active_group_range=[1, 17],
        covered_dense_counts=[0, 1], all_source_supports=True, source_geometry_envelope=True,
        fixed_rows_across_columns=True, all_dense_rowsets_inside_holder=True, fixed_dense_rowsets=65399,
        extra_dense_arity_factor=False, pair_swaps=True, directed=False, whole_code_certificate=False,
        fitted_new_witnesses=0, selected_cells=323, source_banks=3))
    banks = set()
    for bank in data['witnesses']:
        name = bank.get('name')
        if name not in LITERALS or name in banks:
            raise ValueError('Distinct approved witness banks required')
        banks.add(name)
        z, logp = LITERALS[name]
        require_fields(bank, dict(mark=str(z), log2_potential=str(logp)))
        check_proposal_rows(bank['entries'])
    if banks != set(LITERALS):
        raise ValueError('Incomplete witness bank')
    selected = check_proposal_rows(data['entries'], selected=True)
    result = {key: selected[key] for key in KEYS}
    if len(result) != 322 or set(result.values()) != {'dense-inherited', 'h2-inherited'}:
        raise ValueError('Unexpected selected322-cell bank assignment')
    return result


def combined_moments():
    denominator = 255**128
    numerators = [0]*17
    for a in range(3, 17):
        nums, den = arity_holder_moments(a, 16, 16, 8)
        scale, remainder = divmod(denominator, den)
        if remainder:
            raise ArithmeticError('Exact moment denominators must divide the common denominator')
        for w, value in enumerate(nums):
            numerators[w] += scale*value
    digest = sha256(json.dumps([numerators, denominator], separators=(',', ':')).encode()).hexdigest()
    return tuple(numerators), denominator, digest


def combined_weights_arb():
    numerators, denominator, digest = combined_moments()
    return {w: root_upper(n, denominator, 16) for w, n in enumerate(numerators) if n}, digest


def all_source_envelope(raw, maximum=17):
    if type(maximum) is not int or not 1 <= maximum <= 17:
        raise ValueError('Maximum1..17 required')
    family = tuple(tuple(max(raw[m, full][i].upper() for full in range(m//2+1))
                         for i in range(4)) for m in range(maximum+1))
    if any(not value.is_finite() or not value >= 0 for row in family for value in row):
        raise ValueError('Nonnegative finite regional upper endpoints required')
    return [family]*16


def cell_value(key, mark, potential, sparse_caps, dense_caps):
    r, n1, n2, dense = key
    if key not in KEYS:
        raise ValueError('This cell is outside the selected322-cell scope')
    fields = cell_scope(key)
    prefix = (fields['source_supports']*fields['source_colorings']*16**n1*120**n2
              *rational(mark)**(-734003)/(255*potential))
    if dense:
        empty, nonempty = dense_caps
        tail, _ = bivariate_tail_arb(list(sparse_coefficients_arb(empty, n1, n2)),
            list(sparse_coefficients_arb(nonempty, n1, n2)), columns=16,
            x_floor=fields['sparse_support_floor'], y_floor=9)
        tail /= 255**fields['sparse_normalizer_exponent']
    else:
        tail = sparse_tail_arb(sparse_caps, n1, n2, columns=16, dimension=8)
    return prefix*tail


def receipt(precision=256, *, max_seconds=600, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh replay precision must be256 or384')
    if type(max_seconds) not in (int, float) or not isfinite(max_seconds) or not 0 < max_seconds <= 600:
        raise ValueError('Finite positive replay budget at most600 seconds required')
    started = perf_counter()
    selected = checked_selection()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('Small-cover directed budget exceeded; no complete receipt produced')
    entries, costs = [], []
    with ctx.workprec(precision):
        weights, moment_digest = combined_weights_arb()
        totals = {r: arb(0) for r in range(1, 18)}
        for name, mark, logp in WITNESSES:
            keys = [key for key in KEYS if selected[key] == name]
            if not keys:
                continue
            check_time()
            bank_start = perf_counter()
            p = arb(2)**rational(logp)
            raw = arb_pair_regions(mark, 17)
            tables = all_source_envelope(raw)
            sparse = category_caps_arb(tables, p, 17, check_time=check_time)
            vectors = weight_vectors_arb(tables, p, 16, check_time=check_time)
            dense = holder_branch_caps_arb(vectors, weights, p)
            for index, key in enumerate(keys):
                check_time()
                value = cell_value(key, mark, p, sparse, dense)
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError('Finite positive cell bound required')
                entries.append(dict(cell_scope(key), witness=name, mark=str(mark), log2_potential=str(logp),
                    failure_log2_upper=str(ceiling(value.log()/arb(2).log()))))
                totals[key[0]] += value
                if progress and ((index+1) % 25 == 0 or index+1 == len(keys)):
                    progress(dict(stage='directed-bank-cells', witness=name,
                        completed=index+1, total=len(keys), elapsed_seconds=perf_counter()-started))
            costs.append(dict(witness=name, cells=len(keys), elapsed_seconds=perf_counter()-bank_start))
        check_time()
        total = sum(totals.values(), arb(0))
        aggregate = ceiling(total.log()/arb(2).log())
        by_weight = [dict(active_groups=r, failure_log2_upper=str(ceiling(totals[r].log()/arb(2).log())))
                     for r in range(1, 18)]
    entries.sort(key=row_key)
    return dict(EXPECTED_SCOPE, precision_bits=precision, frozen_proposal_sha256=PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS, combined_moments_sha256=moment_digest,
        root_upper_checks_passed=True, entries=entries, by_active_groups=by_weight, bank_costs=costs,
        failure_log2_upper=str(aggregate), scoped_class_reaches_40_bits=aggregate < -40,
        compute_budget_seconds=max_seconds, elapsed_seconds=perf_counter()-started)


def checked_rows(data, precision):
    selected = checked_selection()
    require_fields(data, dict(EXPECTED_SCOPE, precision_bits=precision,
        frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        combined_moments_sha256=combined_moments()[2], root_upper_checks_passed=True,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest()))
    seen = {}
    for row in data.get('entries', []):
        key = row_key(row)
        if key not in selected or key in seen:
            raise ValueError('Unexpected, excluded or repeated cell')
        name = selected[key]
        mark, logp = LITERALS[name]
        require_fields(row, dict(cell_scope(key), witness=name, mark=str(mark), log2_potential=str(logp)))
        if type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Exact rational ceiling required')
        seen[key] = F(row['failure_log2_upper'])
    if set(seen) != set(KEYS):
        raise ValueError('Incomplete322-cell directed receipt')
    weights = {}
    for row in data.get('by_active_groups', []):
        r = row.get('active_groups')
        if type(r) is not int or r not in range(1, 18) or r in weights:
            raise ValueError('Repeated or unknown weight aggregate')
        if type(row.get('failure_log2_upper')) is not str:
            raise ValueError('Exact weight aggregate ceiling required')
        weights[r] = F(row['failure_log2_upper'])
    if set(weights) != set(range(1, 18)):
        raise ValueError('Incomplete weight aggregates')
    return seen, weights


def verify(directory=DIRECTORY):
    receipts, rows, inputs = [], [], []
    for precision in (256, 384):
        path = Path(directory)/f'hd_two_sum_gf16_small_receipt{precision}.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        receipts.append(data)
        rows.append(checked_rows(data, precision))
        inputs.append(dict(filename=path.name, sha256=sha256(raw).hexdigest(),
            precision_bits=precision, elapsed_seconds=data['elapsed_seconds']))
    if rows[0] != rows[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('Dual precision cell, weight or aggregate ceilings disagree')
    rounded = rounded_union_512([rows[0][0][key] for key in KEYS])
    if rounded < F(receipts[0]['failure_log2_upper']):
        raise ArithmeticError('Rounded-cell union lies below the raw ceiling')
    return dict(EXPECTED_SCOPE, frozen_proposal_sha256=PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        combined_moments_sha256=combined_moments()[2], root_upper_checks_passed=True,
        receipt_precision_bits=[256, 384], independent_sum_precision_bits=512,
        all322_cell_ceilings_agree=True, all17_weight_ceilings_agree=True, aggregate_ceilings_agree=True,
        raw_aggregate_log2_upper=receipts[0]['failure_log2_upper'],
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
    with args.output.open('x', encoding='utf-8') as output:
        json.dump(result, output, indent=2, allow_nan=False)
        output.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('entries', 'by_active_groups')}, allow_nan=False), flush=True)
