"""Directed replay of twelve frozen witness/arity Holder classes.

Scope: r17, source(8,0), sixteen two-row groups and one dense group, all
28^16 fixed sparse pair assignments, existing h4 DHD inner and fair swaps.
Exact W_B(T)^32=sum_(A in B) N_A(T)^32/255^(16|A|) charges dense row sets
and their comparison factors once. Each root's dyadic upper endpoint is
checked by exact integer exponentiation before it is used. Mask weights
are mixed before the shared potential quotient, and both coefficient
floors are retained. The frozen FLOAT estimates are never used.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_dense_holder_classes import class_moments
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_mixed_all_geometry_certificate import checked_classes as fixed_A_classes
from hd_xor_mixed_mask_certificate import (mixed_vectors_arb, column_coefficients_arb,
    bivariate_tail_arb, _rational, _upper_max)
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import arb_pair_regions
from hd_xor_paired_joint_gate import ordered_tables


DIRECTORY = Path(__file__).resolve().parent
PROPOSAL_NAME = 'hd_xor_dense_holder_8_0_h1_arity12.json'
PROPOSAL_SHA256 = '4577166d2aefbe1912ac6bc096c8f08a3604439c03f72ff0a2dd7ad4fef74a0e'
BASE_PROPOSAL_SHA256 = '44fd35d8bf2e78391572a14d071cd283dc547f7bc204b3802c180132eaaab4e6'
HELPER_PINS = {
    'hd_xor_dense_holder.py': '24ea5d4614229549d7c6f286094b5d0498266da6060907740d5d47aeafecc012',
    'hd_xor_dense_holder_classes.py': '8a03df6021cfd8b68051e53edf6d2d9ac90f11d076524ff788fcbb67968e24c1',
    'hd_xor_mixed_mask_certificate.py': 'b243ecbc00527ef39f151a276b0206323f5b4cfdd14361b380326bbf75dc8f6b',
    'hd_xor_mixed_all_geometry_certificate.py': 'e39f4d8419c8361269896acf21f34b328d3b76987b2bee49690a5d73a0f59e63',
}
SOURCE_COUNT = 15708305994709246428578365440
ALL_MASKS = tuple(mask for mask in range(256) if mask.bit_count() >= 3)
EXPECTED_SCOPE = dict(
    construction='xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half',
    cohort_counts=[0, 16, 1], source_matching_counts=[8, 0], mu_range=[1, 16],
    active_groups=17, source_classes=1, field_order=256, forbidden_weight=734003,
    total_columns=32, sparse_active_rows=32, sparse_support_floor=544,
    dense_nonempty_column_floor=17, sparse_normalizer_exponent=512,
    holder_exponent=32, deterministic_classes=12, fixed_row_sets=219,
    source_supports=SOURCE_COUNT, source_colorings=17, sparse_pair_assignments=28**16,
    pair_swaps=True, alpha=1, dense_rowsets_charged_once_inside_holder=True,
    extra_dense_beta_factor=False, extra_dense_rowset_factor=False,
    dense_row_set_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
    iid_rowset_resampling_assumed=False, geometry_resampled_across_columns=False,
    setup_dependent_selection=False, labelled_dense_masks_retained=True,
    masks_and_cancellations_mixed_before_potential=True,
    common_potential_within_each_class=True, same_potential_both_dense_branches=True,
    exact_positive_bivariate_tail=True, projective_normalizations=1,
    endpoint_potential_factors=1, exact_root_endpoint_checks=True,
    root_enclosure='nonnegative-dyadic-upper-power32-checked-against-exact-rational',
    regional_bank_evaluations=5, fitted_new_witnesses=0, saved_float_estimates_used=False,
    directed=True, scoped_class_complete=True, complete_mixed_class=False,
    whole_code_certificate=False)


def require_fields(data, expected):
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError(f'Wrong frozen scope field: {key}')


def validate_scope(scope=None):
    require_fields(EXPECTED_SCOPE if scope is None else scope, EXPECTED_SCOPE)
    if (Q, REGION_LENGTH, LIMIT) != (256, 8192, 734003):
        raise ValueError('Imported construction constants changed')
    for filename, digest in HELPER_PINS.items():
        if sha256((DIRECTORY/filename).read_bytes()).hexdigest() != digest:
            raise ValueError(f'Reviewed helper changed: {filename}')
    return dict(EXPECTED_SCOPE)


def checked_partition():
    """Read only exact masks and rational witnesses; ignore all FLOAT bounds."""
    raw = (DIRECTORY/PROPOSAL_NAME).read_bytes()
    if sha256(raw).hexdigest() != PROPOSAL_SHA256:
        raise ValueError('Frozen twelve-type proposal changed')
    proposal = json.loads(raw)
    keys = ('construction', 'cohort_counts', 'source_matching_counts', 'mu_range',
            'active_groups', 'field_order', 'forbidden_weight', 'total_columns',
            'sparse_support_floor', 'dense_nonempty_column_floor', 'sparse_normalizer_exponent',
            'holder_exponent', 'deterministic_classes', 'fixed_row_sets', 'sparse_pair_assignments',
            'source_colorings', 'extra_dense_beta_factor', 'extra_dense_rowset_factor',
            'iid_rowset_resampling_assumed', 'setup_dependent_selection',
            'labelled_dense_masks_retained', 'projective_normalizations', 'endpoint_potential_factors')
    require_fields(proposal, {key: EXPECTED_SCOPE[key] for key in keys})
    require_fields(proposal, dict(directed=False, fitted_new_witnesses=0,
        same_frozen_witness_each_intersection=True, witness_classes_refined_only_by_arity=True,
        comparator_proposal_sha256=BASE_PROPOSAL_SHA256))
    original = {mask: (label, z, logp) for label, z, logp, masks in fixed_A_classes() for mask in masks}
    seen, classes, identifiers = set(), [], set()
    for row in proposal['entries']:
        label, z, logp = row['witness_label'], F(row['mark']), F(row['log2_potential'])
        arity, masks = row['active_constituents'], tuple(row['active_row_masks'])
        if (type(arity) is not int or not 3 <= arity <= 8 or not masks
                or any(type(mask) is not int or mask not in ALL_MASKS or mask.bit_count() != arity for mask in masks)
                or tuple(sorted(set(masks))) != masks or seen.intersection(masks)
                or row['fixed_row_sets'] != len(masks) or not 0 < z < 1 or logp > 0):
            raise ValueError('Invalid exact twelve-type partition')
        if any(original[mask] != (label, z, logp) for mask in masks):
            raise ValueError('Partition changed the original fixed-A witness')
        identifier = f'{label}:a{arity}'
        if identifier in identifiers:
            raise ValueError('Duplicate witness/arity type')
        identifiers.add(identifier)
        seen.update(masks)
        classes.append((identifier, label, z, logp, arity, masks))
    if (seen != set(ALL_MASKS) or len(classes) != 12
            or len({(row[1], row[2], row[3]) for row in classes}) != 5):
        raise ValueError('The exact partition must cover219 masks in12 types and5 banks')
    return tuple(sorted(classes))


def dyadic_power_dominates(value, numerator, denominator, exponent=32):
    """Exact integer check, not a floating comparison of overlapping balls."""
    if (type(numerator) is not int or type(denominator) is not int
            or numerator < 0 or denominator <= 0 or type(exponent) is not int or exponent < 1
            or not value.is_finite() or not value >= 0):
        raise ValueError('A nonnegative dyadic value and exact rational moment required')
    mantissa, shift = value.man_exp()  # Reject a non-exact Arb ball.
    left, right = int(mantissa)**exponent*denominator, numerator
    shift = int(shift)*exponent
    return (left << shift) >= right if shift >= 0 else left >= (right << -shift)


def root_upper(numerator, denominator, exponent=32):
    """Positive root enclosure with an independently checked upper endpoint.

    x -> x**exponent is monotone on nonnegative reals. Therefore the exact
    dyadic endpoint u is >=(numerator/denominator)**(1/exponent) whenever
    u**exponent*denominator>=numerator. Arb proposes u; integer arithmetic
    verifies the only root direction consumed by the proof.
    """
    if (type(numerator) is not int or type(denominator) is not int
            or numerator < 0 or denominator <= 0 or type(exponent) is not int or exponent < 1):
        raise ValueError('Nonnegative exact numerator, positive denominator and exponent required')
    upper = arb(0) if numerator == 0 else (arb(numerator)/denominator).root(exponent).upper()
    if not dyadic_power_dominates(upper, numerator, denominator, exponent):
        raise ArithmeticError('The proposed root endpoint is not an exact upper bound')
    return upper


def weights_arb(masks):
    numerators, denominator = class_moments(masks)
    weights = {mask: root_upper(numerator, denominator)
               for mask, numerator in enumerate(numerators) if numerator}
    digest = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
    return weights, digest


def holder_branch_caps_arb(vectors, weights, potential):
    """Cancellation is already mixed; mix whole dense masks before quotient."""
    if (not potential.is_finite() or not potential > 0 or 0 not in weights
            or any(not weight.is_finite() or not weight >= 0 for weight in weights.values())
            or any(not weights.keys() <= table.keys() for table in vectors.values())):
        raise ValueError('Positive potential, nonnegative weights, and complete mask vectors required')
    empty, nonempty = {}, {}
    for kind, table in vectors.items():
        left, right = table[0]
        empty[kind] = _upper_max(weights[0]*left, weights[0]*right/potential)
        first, second = arb(0), arb(0)
        for mask in sorted(weights.keys()-{0}):
            left, right = table[mask]
            first += weights[mask]*left
            second += weights[mask]*right
        nonempty[kind] = _upper_max(first, second/potential)
    return empty, nonempty


def receipt(precision=256, *, scope=None, progress=None):
    metadata = validate_scope(scope)
    if type(precision) is not int or precision < 128:
        raise ValueError('At least128 precision bits are required')
    classes = checked_partition()
    if matching_joint_counts(8192, 17)[8, 0] != SOURCE_COUNT:
        raise ArithmeticError('Exact source-support count changed')
    banks = {}
    for row in classes:
        banks.setdefault(row[1:4], []).append(row)
    started, entries, costs = perf_counter(), [], []
    with ctx.workprec(precision):
        total = arb(0)
        for (label, z, logp), rows in banks.items():
            bank_started = perf_counter()
            potential = arb(2)**_rational(logp)
            weighted_rows = [(row, *weights_arb(row[-1])) for row in rows]
            required = sorted({mask for _, weights, _ in weighted_rows for mask in weights})
            if progress:
                progress(dict(stage='directed-bank-start', witness_label=label, types=len(rows),
                              dense_masks=len(required)))
            raw = arb_pair_regions(z, 17)
            regions = ordered_tables(raw, 8, 0, directed=True)
            vectors = mixed_vectors_arb(regions, potential, required, groups=16, mu=16, progress=progress)
            prefix = SOURCE_COUNT*17*28**16*_rational(z)**(-LIMIT)/(255*potential*255**512)
            for row, weights, moment_digest in weighted_rows:
                identifier, _, _, _, arity, masks = row
                empty, nonempty = holder_branch_caps_arb(vectors, weights, potential)
                tail, _ = bivariate_tail_arb(column_coefficients_arb(empty), column_coefficients_arb(nonempty),
                                             columns=32, x_floor=544, y_floor=17)
                value = prefix*tail
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError('The positive type bound must be finite and positive')
                entry = dict(type_id=identifier, witness_label=label, mark=str(z), log2_potential=str(logp),
                    active_constituents=arity, active_row_masks=list(masks), fixed_row_sets=len(masks),
                    exact_moments_sha256=moment_digest, positive_mask_weights=len(weights),
                    root_upper_checks_passed=True,
                    failure_log2_upper=str(rational_ceiling(value.log()/arb(2).log())))
                entries.append(entry)
                total += value
                if progress:
                    progress(dict(stage='directed-type-complete', **entry))
            cost = dict(witness_label=label, types=len(rows), dense_masks=len(required),
                        elapsed_seconds=perf_counter()-bank_started)
            costs.append(cost)
            if progress:
                progress(dict(stage='directed-bank-complete', **cost))
        entries.sort(key=lambda row: row['type_id'])
        metadata.update(precision_bits=precision, frozen_proposal_sha256=PROPOSAL_SHA256,
            base_proposal_sha256=BASE_PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
            entries=entries, bank_costs=costs, elapsed_seconds=perf_counter()-started,
            failure_log2_upper=str(rational_ceiling(total.log()/arb(2).log())),
            scoped_class_reaches_40_bits=bool(total < arb(2)**-40))
    return metadata


def checked_rows(data, precision):
    expected = dict(EXPECTED_SCOPE, precision_bits=precision,
        frozen_proposal_sha256=PROPOSAL_SHA256, base_proposal_sha256=BASE_PROPOSAL_SHA256,
        reviewed_helper_sha256=HELPER_PINS,
        source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())
    require_fields(data, expected)
    assignment = {row[0]: row for row in checked_partition()}
    found = {}
    for entry in data.get('entries', []):
        identifier = entry.get('type_id')
        if identifier not in assignment or identifier in found:
            raise ValueError('Unexpected or duplicate type')
        _, label, z, logp, arity, masks = assignment[identifier]
        numerators, denominator = class_moments(masks)
        digest = sha256(json.dumps([list(numerators), denominator], separators=(',', ':')).encode()).hexdigest()
        require_fields(entry, dict(witness_label=label, mark=str(z), log2_potential=str(logp),
            active_constituents=arity, active_row_masks=list(masks), fixed_row_sets=len(masks),
            exact_moments_sha256=digest, positive_mask_weights=sum(n > 0 for n in numerators),
            root_upper_checks_passed=True))
        if type(entry.get('failure_log2_upper')) is not str:
            raise ValueError('Exact rational type ceiling required')
        found[identifier] = F(entry['failure_log2_upper'])
    if found.keys() != assignment.keys():
        raise ValueError('Incomplete twelve-type receipt')
    return found


def rounded_union_512(ceilings):
    with ctx.workprec(512):
        total = sum((arb(2)**_rational(value) for value in ceilings), arb(0))
        if not total > 0 or not total.is_finite():
            raise ValueError('A finite nonempty positive union is required')
        return rational_ceiling(total.log()/arb(2).log())


def verify(directory=DIRECTORY):
    validate_scope()
    directory = Path(directory)
    receipts, rows, inputs = [], [], []
    for precision in (256, 384):
        path = directory/f'hd_xor_dense_holder_8_0_h1_receipt{precision}.json'
        raw = path.read_bytes()
        data = json.loads(raw)
        rows.append(checked_rows(data, precision))
        receipts.append(data)
        inputs.append(dict(filename=path.name, sha256=sha256(raw).hexdigest(), precision_bits=precision,
                           elapsed_seconds=data['elapsed_seconds']))
    if rows[0] != rows[1] or receipts[0]['failure_log2_upper'] != receipts[1]['failure_log2_upper']:
        raise ValueError('The two precisions disagree on a type or aggregate rational ceiling')
    rounded = rounded_union_512([rows[0][key] for key in sorted(rows[0])])
    if rounded < F(receipts[0]['failure_log2_upper']):
        raise ArithmeticError('Rounded-cell union lies below the raw aggregate ceiling')
    return dict(EXPECTED_SCOPE, frozen_proposal_sha256=PROPOSAL_SHA256,
        base_proposal_sha256=BASE_PROPOSAL_SHA256, reviewed_helper_sha256=HELPER_PINS,
        receipt_precision_bits=[256, 384], independent_sum_precision_bits=512,
        exact_partition_checked=True, all_per_type_ceilings_agree=True, aggregate_ceilings_agree=True,
        raw_aggregate_log2_upper=receipts[0]['failure_log2_upper'], rounded_cell_union_log2_upper=str(rounded),
        scoped_class_reaches_40_bits=rounded < -40, inputs=inputs,
        receipt_arithmetic_trust='observed-fresh-directed-replays-not-proved-by-byte-hashes')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--precision', type=int, choices=(256, 384))
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--directory', type=Path, default=DIRECTORY)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen directed artifact')
    result = verify(args.directory) if args.verify else receipt(args.precision,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
