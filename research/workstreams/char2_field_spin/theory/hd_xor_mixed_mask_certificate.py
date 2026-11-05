"""Fresh directed replay of the frozen mixed source(8,0), mu16 proposal.

The default replays only the previously worst fixed row mask195. A full
219-mask replay is explicit and uses the unchanged frozen assignment.
Regional coefficients, offset resource bounds, both mixtures, and both
positive support floors are rebuilt in the requested Arb precision.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb
import json
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

from hd_shifted_coverage import rational_ceiling
from hd_xor_fixed_mask_floor import comparison_counts
from hd_xor_geometry_strata import category_resources, geometry_stratum_counts
from hd_xor_group_geometry import matching_joint_counts
from hd_xor_joint_columns import categories
from hd_xor_mixed_mask_witnesses import (checked_classes, CONSTRUCTION,
    SOURCE_MATCHING_COUNTS, MU_RANGE, SPARSE_SUPPORT_FLOOR, DENSE_COLUMN_FLOOR,
    PROPOSAL_SHA256)
from hd_xor_multigroup import Q, REGION_LENGTH, LIMIT
from hd_xor_pair_route import arb_pair_regions
from hd_xor_paired_joint_gate import ordered_tables


EXPECTED_SCOPE = dict(
    construction='xor2-gf256-h4-alternating-swapped-pairs-k20-rate-half',
    cohort_counts=(0, 16, 1), source_matching_counts=(8, 0), mu_range=(16, 16),
    active_groups=17, sparse_active_rows=32, field_order=256,
    forbidden_weight=734003, total_columns=32, sparse_support_floor=544,
    dense_nonempty_column_floor=17, sparse_normalizer_exponent=512,
    pair_swaps=True, alpha=1)


def validate_scope(scope=None):
    scope = EXPECTED_SCOPE if scope is None else scope
    if any(scope.get(key) != value for key, value in EXPECTED_SCOPE.items()):
        raise ValueError('Only the frozen cohort/source/mu/field/length scope is supported')
    if (CONSTRUCTION != EXPECTED_SCOPE['construction'] or Q != 256 or REGION_LENGTH != 8192
            or LIMIT != 734003 or SOURCE_MATCHING_COUNTS != (8, 0) or MU_RANGE != (16, 16)
            or SPARSE_SUPPORT_FLOOR != 544 or DENSE_COLUMN_FLOOR != 17):
        raise ValueError('Imported construction constants disagree with the frozen proof scope')
    return dict(EXPECTED_SCOPE)


def _rational(value):
    value = F(value)
    return arb(value.numerator)/value.denominator


def _upper_max(left, right):
    left, right = left.upper(), right.upper()
    return left if left >= right else right


def _regional_uppers(regions):
    tables = tuple(tuple(tuple(entry.upper() for entry in matrix) for matrix in table)
                   for table in regions)
    if not 1 <= len(tables) <= 8 or any(any(len(matrix) != 4 for matrix in table) for table in tables):
        raise ValueError('One to eight labelled tables of flat-four matrices required')
    if any(not entry.is_finite() or not entry >= 0 for table in tables for matrix in table for entry in matrix):
        raise ValueError('Finite nonnegative regional upper endpoints required')
    return tables


def offset_resource_arb(regions, potential, lower, upper, masks, targets=None):
    """Upper vectors for selected fixed dense masks and sparse sum targets.

    Identical suffix-mask codes share a suffix computation. A partial sum
    j is retained whenever j<=target<=j+remaining*width for any target.
    Thus pruning cannot remove a path to any requested coefficient.
    Maxima compare explicit upper endpoints, never overlapping balls.
    """
    tables = _regional_uppers(regions)
    length = len(tables)
    masks = tuple(sorted(set(masks)))
    if (type(lower) is not int or type(upper) is not int or not 0 <= lower <= upper
            or any(upper+1 >= len(table) for table in tables)
            or not masks or any(type(mask) is not int or not 0 <= mask < 1 << length for mask in masks)
            or not potential > 0 or not potential.is_finite()):
        raise ValueError('Invalid offset table, mask selection, occupancy interval, or potential')
    width = upper-lower
    targets = tuple(range(length*width+1)) if targets is None else tuple(sorted(set(targets)))
    if not targets or any(type(target) is not int or not 0 <= target <= length*width for target in targets):
        raise ValueError('Requested resource targets lie outside the attainable sum range')
    zero = arb(0)
    values = {0: {0: (arb(1), potential.upper())}}
    for position in reversed(range(length)):
        processed = length-position
        degrees = [degree for degree in range(processed*width+1)
                   if any(degree <= target <= degree+position*width for target in targets)]
        following = {}
        for code in sorted({mask >> position for mask in masks}):
            bit, suffix = code & 1, values[code >> 1]
            coefficients = {}
            for degree in degrees:
                first, second = zero, zero
                for delta in range(max(0, degree-(processed-1)*width), min(width, degree)+1):
                    previous = suffix.get(degree-delta)
                    if previous is None:
                        continue
                    a, b, c, d = tables[position][lower+delta+bit]
                    left, right = previous
                    candidate0 = (a*left+b*right).upper()
                    candidate1 = (c*left+d*right).upper()
                    if candidate0 > first:
                        first = candidate0
                    if candidate1 > second:
                        second = candidate1
                coefficients[degree] = first, second
            following[code] = coefficients
        values = following
    return values


def mixed_vectors_arb(regions, potential, masks, groups=16, mu=16, progress=None):
    """Exact cancellation mixture, preserving fixedT and both state coordinates."""
    if type(groups) is not int or not 1 <= groups <= 16 or mu != groups or len(regions) != 8:
        raise ValueError('This directed offset evaluator requires eight regions and mu=groups<=16')
    masks = tuple(sorted(set(masks)))
    requests = {}
    for kind in categories(groups):
        single, double = kind
        for cancelled in range(double+1):
            lo, hi, limit, target = category_resources(single, double, cancelled, mu)
            if limit < 2*hi:
                raise ArithmeticError('An active pair constraint cannot be silently discarded')
            weight = arb(comb(double, cancelled)*254**(double-cancelled))/255**double
            requests.setdefault((lo, hi), []).append((kind, target, weight))
    result = {kind: {mask: (arb(0), arb(0)) for mask in masks} for kind in categories(groups)}
    started = perf_counter()
    for completed, ((lo, hi), uses) in enumerate(sorted(requests.items()), 1):
        vectors = offset_resource_arb(regions, potential, lo, hi, masks,
                                      targets={target for _, target, _ in uses})
        for kind, target, weight in uses:
            for mask in masks:
                left, right = result[kind][mask]
                first, second = vectors[mask][target]
                result[kind][mask] = left+weight*first, right+weight*second
        if progress and (completed % 25 == 0 or completed == len(requests)):
            progress(dict(stage='directed-offset-bank', completed=completed, total=len(requests),
                          dense_masks=len(masks), elapsed_seconds=perf_counter()-started))
    return result


def fixed_branch_caps_arb(vectors, active_mask, potential, group=8, q=256):
    """Mix exactN_A(T)/q^a before the common-potential quotient."""
    counts = comparison_counts(active_mask, group, q)
    normalizer = q**active_mask.bit_count()
    required = {mask for mask, count in enumerate(counts) if count}
    if not potential > 0 or not potential.is_finite() or any(not required <= table.keys() for table in vectors.values()):
        raise ValueError('Positive potential and every positive-probability dense mask are required')
    empty, nonempty = {}, {}
    for kind, table in vectors.items():
        first, second = table[0]
        empty[kind] = _upper_max(first/normalizer, second/(normalizer*potential))
        first, second = arb(0), arb(0)
        for mask in sorted(required-{0}):
            left, right = table[mask]
            first += counts[mask]*left
            second += counts[mask]*right
        nonempty[kind] = _upper_max(first/normalizer, second/(normalizer*potential))
    return empty, nonempty


def column_coefficients_arb(caps, groups=16, q=256):
    if set(caps) != set(categories(groups)) or q < 2:
        raise ValueError('Complete positive category caps are required')
    coefficients = [arb(0) for _ in range(2*groups+1)]
    for (single, double), cap in caps.items():
        cap = cap.upper()
        if not cap >= 0 or not cap.is_finite():
            raise ValueError('Finite nonnegative category upper caps required')
        degree = single+2*double
        coefficients[degree] += comb(groups, single)*comb(groups-single, double)*2**single*(q-1)**degree*cap
    return coefficients


def positive_powers_arb(coefficients, maximum):
    if type(maximum) is not int or maximum < 0 or not coefficients:
        raise ValueError('A nonempty coefficient vector and nonnegative power bound are required')
    if any(not value.is_finite() or not value >= 0 for value in coefficients):
        raise ValueError('Polynomial coefficients must be finite and nonnegative')
    base = arb_poly(coefficients)
    result = [arb_poly([1])]
    for _ in range(maximum):
        result.append(result[-1]*base)
    return result


def bivariate_tail_arb(empty, nonempty, columns=32, x_floor=544, y_floor=17,
                       *, empty_powers=None, empty_scale=None):
    """Positive scalar binomial expansion retaining both exact coefficient floors."""
    if (type(columns) is not int or columns < 0 or type(x_floor) is not int or x_floor < 0
            or type(y_floor) is not int or not 0 <= y_floor <= columns):
        raise ValueError('Invalid column count or coefficient floors')
    if empty_powers is None:
        if empty_scale is not None:
            raise ValueError('A cached power bank is required for an external empty scaling')
        empty_powers = positive_powers_arb(empty, columns-y_floor)
        empty_scale = arb(1)
    if (len(empty_powers) <= columns-y_floor or empty_scale is None
            or not empty_scale > 0 or not empty_scale.is_finite()):
        raise ValueError('Complete cached empty powers and a positive exact scaling are required')
    dense_powers = positive_powers_arb(nonempty, columns)
    maximum_x = columns*max(len(empty)-1, len(nonempty)-1)
    total, by_degree = arb(0), []
    for degree in range(y_floor, columns+1):
        polynomial = empty_powers[columns-degree]*dense_powers[degree]
        # Arb polynomial multiplication encloses each true nonnegative
        # coefficient. No subtraction or FLOAT coefficient cutoff is used.
        value = (sum((polynomial[index] for index in range(x_floor, maximum_x+1)), arb(0))
                 *comb(columns, degree)*empty_scale**(columns-degree))
        by_degree.append(value)
        total += value
    return total, by_degree


def receipt(precision=256, active_masks=(195,), *, scope=None, progress=None):
    metadata = validate_scope(scope)
    if type(precision) is not int or precision < 128:
        raise ValueError('At least128 precision bits are required')
    assigned = checked_classes()
    all_masks = {mask for _, _, _, masks in assigned for mask in masks}
    requested = tuple(active_masks)
    if not requested or any(type(mask) is not int for mask in requested):
        raise ValueError('Fixed row masks must be integers')
    selected = tuple(sorted(requested))
    if len(set(selected)) != len(selected) or not set(selected) <= all_masks:
        raise ValueError('Select distinct fixed row masks from the frozen219-class assignment')
    started, entries = perf_counter(), []
    with ctx.workprec(precision):
        total = arb(0)
        source_count = matching_joint_counts(8192, 17)[8, 0]
        geometry_count = geometry_stratum_counts(16)[16]
        if geometry_count != 28:
            raise ArithmeticError('The common sparse-pair stratum must contain exactly28 assignments')
        for label, z, logp, bank_masks in assigned:
            wanted = sorted(set(selected) & set(bank_masks))
            if not wanted:
                continue
            bank_started = perf_counter()
            potential = arb(2)**_rational(logp)
            raw = arb_pair_regions(z, 17)
            regions = ordered_tables(raw, 8, 0, directed=True)
            dense_masks = sorted({mask for active_mask in wanted
                                  for mask, count in enumerate(comparison_counts(active_mask)) if count})
            vectors = mixed_vectors_arb(regions, potential, dense_masks, progress=progress)
            empty_caps = {kind: _upper_max(table[0][0], table[0][1]/potential)
                          for kind, table in vectors.items()}
            unscaled_empty = column_coefficients_arb(empty_caps)
            empty_powers = positive_powers_arb(unscaled_empty, 15)
            prefix = source_count*17*geometry_count*_rational(z)**(-734003)/(255*potential*255**512)
            for active_mask in wanted:
                active = active_mask.bit_count()
                empty, nonempty = fixed_branch_caps_arb(vectors, active_mask, potential)
                empty_poly, nonempty_poly = column_coefficients_arb(empty), column_coefficients_arb(nonempty)
                tail, by_degree = bivariate_tail_arb(empty_poly, nonempty_poly,
                    empty_powers=empty_powers, empty_scale=arb(1)/256**active)
                beta = arb(256**(32*active))/255**(16*active)
                value = prefix*beta*tail
                if not value.is_finite() or not value > 0:
                    raise ArithmeticError('The directed polynomial tail is not a finite positive enclosure')
                logarithm = value.log()/arb(2).log()
                entry = dict(active_row_mask=active_mask, active_constituents=active,
                    witness_label=label, mark=str(z), log2_potential=str(logp),
                    failure_log2_upper=str(rational_ceiling(logarithm)),
                    largest_upper_dense_degree=17+max(range(len(by_degree)), key=lambda index: by_degree[index].upper()))
                entries.append(entry)
                total += value
                if progress:
                    progress(dict(stage='directed-fixed-A', **entry))
            if progress:
                progress(dict(stage='directed-bank-complete', witness_label=label,
                              assigned_masks=len(wanted), dense_masks=len(dense_masks),
                              elapsed_seconds=perf_counter()-bank_started))
        entries.sort(key=lambda row: row['active_row_mask'])
        if [row['active_row_mask'] for row in entries] != list(selected):
            raise ArithmeticError('The directed replay did not visit each selected fixedA exactly once')
        metadata.update(precision_bits=precision, directed=True, fitted_new_witnesses=0,
            frozen_proposal_sha256=PROPOSAL_SHA256,
            frozen_witness_source_sha256=sha256(Path(__file__).with_name('hd_xor_mixed_mask_witnesses.py').read_bytes()).hexdigest(),
            source_supports=source_count, dense_source_colorings=17, two_row_geometries=geometry_count,
            fixed_row_sets=len(entries), scoped_class_complete=set(selected) == all_masks,
            complete_mixed_class=False, whole_code_certificate=False,
            dense_row_set_fixed_across_columns=True, sparse_geometry_fixed_across_columns=True,
            dense_input_includes_zero=True, exact_uniform_dense_mask_mixture=True,
            masks_and_cancellations_mixed_before_potential=True, same_potential_both_dense_branches=True,
            dense_beta_formula='q^(32a)/(q-1)^(16a)', exact_positive_bivariate_tail=True,
            chernoff_dense_floor=False, product_density_comparison_used=False,
            projective_normalizations=1, endpoint_potential_factors=1,
            failure_log2_upper=str(rational_ceiling(total.log()/arb(2).log())),
            selected_classes_reach_40_bits=bool(total < arb(2)**-40), entries=entries,
            elapsed_seconds=perf_counter()-started)
        return metadata


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, default=256)
    choose = parser.add_mutually_exclusive_group()
    choose.add_argument('--mask', type=int, action='append')
    choose.add_argument('--all', action='store_true', help='Explicit full219-class replay; potentially expensive')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Do not overwrite a frozen directed receipt')
    masks = [mask for _, _, _, items in checked_classes() for mask in items] if args.all else args.mask or [195]
    result = receipt(args.precision, masks,
        progress=lambda row: print(json.dumps(row, allow_nan=False), flush=True))
    result['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    with args.output.open('x', encoding='utf-8') as target:
        json.dump(result, target, indent=2, allow_nan=False)
        target.write('\n')
    print(json.dumps({key: value for key, value in result.items() if key != 'entries'}, allow_nan=False), flush=True)
