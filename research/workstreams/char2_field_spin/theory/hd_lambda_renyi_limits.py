"""Conservative and ideal-method Renyi controls for all Hlambda arities.

Fix A of size a and C=Hlambda F^A. A nonzero word with exact support S
is determined by any |S|-d+1 of its nonzero coordinates when d(C)>=d.
Thus n_C(S)<=255^(|S|-d+1). After the nonzero output diagonal, the
normalized likelihood at nonzero y is at most h=256^(16-a)/255^(d-1).
Its nonzero probability mass is1-256^-a. Therefore its nonzero p-moment
is at most h^(p-1)*(1-256^-a). The zero contribution is exactly
z=256^(p*(16-a)-16).

The distance lower bounds for a=1..15 are16,6,5,4,3,2,...,2. Each Hlambda
column is full. Its nonzero affine correction has weight at least8;
its zero-correction code has minimum weight4. Finally, Hlambda^-1=Hlambda
with all entries nonzero excludes output weight1 unless a=16.
At a=16 the likelihood is exactly1, handled separately.

The IDEAL control instead uses Jensen's lower bound for the nonzero
moment, (1-q^-a)^p/(1-q^-16)^(p-1). It is a lower comparator for THIS
norm-based method, never an upper bound on failure or a security claim.

Both controls retain at least9 nonzero columns using a positive binomial
tail. Each raw group weight includes its rowset count and RS normalization:
B_a=C(16,a)*256^(16a)/255^(8a)*tail^(1/p). The later polynomial has these
RAW weights; no additional beta factor is applied. Original rowsets stay
fixed across columns. No inventory or new construction is required.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

import hd_lambda_renyi_activity as activity

DIRECTORY = Path(__file__).resolve().parent
A, binding = activity.A, activity.binding
ALPHA, MARK = F(1, 5), F(691, 700)
ORDERS, ARITIES, R = (2, 3, 4), tuple(range(1, 17)), 18
MIN_TOTAL, MAX_TOTAL = 18, 245
COARSE, IDEAL = 'coarse-upper', 'ideal-lower-on-this-method'
KINDS = (COARSE, IDEAL)
PINS = {'hd_lambda_renyi_activity.py': '09f5245aecb903097438b65ded8757c6a66da40a13804289c394fef777186c77'}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return {**activity.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def minimum_distance_bound(arity):
    if type(arity) is not int or not 1 <= arity <= 16:
        raise ValueError('Arity must be1..16')
    return 16 if arity == 1 else max(2, 8-arity) if arity < 16 else 1


def coarse_nonzero_moment(arity, order, distance, *, q=256, width=16):
    """General shortened-code bound, given a PROVED distance lower bound."""
    if (type(q) is not int or q not in (4, 256) or type(width) is not int or not 1 <= width <= 16
            or type(arity) is not int or not 1 <= arity <= width
            or order not in ORDERS or type(order) is not int
            or type(distance) is not int or not 1 <= distance <= width-arity+1):
        raise ValueError('Valid dimension, distance, GF4/GF256, and order2/3/4 required')
    h = F(q)**(width-arity)/F(q-1)**(distance-1)
    return h**(order-1)*(1-F(q)**(-arity))


def column_atoms(arity, order, kind, *, q=256):
    """Return exact (zero, nonzero) p-moment atoms for16 coordinates."""
    if kind not in KINDS or type(q) is not int or q not in (4, 256) or type(order) is not int or order not in ORDERS:
        raise ValueError('Specify coarse-upper or ideal-method, GF4/GF256, and order2/3/4')
    distance = minimum_distance_bound(arity)
    zero = F(q)**(order*(16-arity)-16)
    if arity == 16:
        nonzero = 1-F(q)**(-16)
    elif kind == COARSE:
        nonzero = coarse_nonzero_moment(arity, order, distance, q=q)
    else:
        nonzero = (1-F(q)**(-arity))**order/(1-F(q)**(-16))**(order-1)
    if zero <= 0 or nonzero <= 0:
        raise ArithmeticError('Moment atoms must be positive')
    return zero, nonzero


def raw_group_weights(order, kind):
    """Call under an Arb context; returns sixteen RAW B_a and diagnostics."""
    if order not in ORDERS or type(order) is not int or kind not in KINDS:
        raise ValueError('Use order2/3/4 and an explicit comparison kind')
    values, rows = [], []
    for arity in ARITIES:
        zero, nonzero = column_atoms(arity, order, kind)
        tail = activity.activity_tail(zero, (nonzero,)*16, 9)
        raw = (arb(comb(16, arity))*arb(256)**(16*arity)/arb(255)**(8*arity)
               * tail**A.uniform.rational(F(1, order)))
        if not raw > 0 or not raw.is_finite():
            raise ArithmeticError('Invalid raw group comparison weight')
        values.append(raw)
        rows.append(dict(arity=arity, labelled_rowsets=comb(16, arity),
            distance_lower_bound=minimum_distance_bound(arity), zero_moment=str(zero),
            nonzero_moment=str(nonzero), raw_weight_log2_enclosure=str(A.log2(raw))))
    return values, rows


def total_log_bound(raw_coefficient, order, uniform_log, threshold_log):
    """Raw coefficient includes ALL dimensional/RS factors; no beta here."""
    if order not in ORDERS or type(order) is not int or not raw_coefficient > 0:
        raise ValueError('Positive raw coefficient and order2/3/4 required')
    theta = ALPHA*F(order, order-1)
    return (A.uniform.rational(theta)*(A.log2(raw_coefficient)-A.log2(arb(255)))
            + threshold_log + uniform_log)


def _floor(value):
    return -A.uniform.rational_ceiling(-value)


def controls_from_uniform_log(uniform_log):
    """Evaluate both controls under the SAME fixed uniform-provider bound."""
    threshold_log = -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))
    controls = []
    for kind in KINDS:
        for order in ORDERS:
            weights, diagnostics = raw_group_weights(order, kind)
            # Coefficient index0 means total arity18, since each group starts1.
            coefficients = arb_poly(weights)**R
            rows = []
            for total in range(MIN_TOTAL, MAX_TOTAL+1):
                value = total_log_bound(coefficients[total-MIN_TOTAL], order, uniform_log, threshold_log)
                row = dict(total_arity=total)
                if kind == COARSE:
                    row['failure_log2_upper'] = str(A.uniform.rational_ceiling(value))
                else:
                    row['method_expression_log2_lower'] = str(_floor(value))
                    row['method_expression_log2_upper'] = str(A.uniform.rational_ceiling(value))
                rows.append(row)
            controls.append(dict(kind=kind, moment_order=order, outer_theta=str(ALPHA*F(order, order-1)),
                raw_group_weights_include_rowsets_and_RS=True, additional_beta_factor=False,
                weight_diagnostics=diagnostics, totals=rows))
    coarse = [c for c in controls if c['kind'] == COARSE]
    ideal = [c for c in controls if c['kind'] == IDEAL]
    selected, coarse_sum, ideal_sum = [], arb(0), arb(0)
    for index, total in enumerate(range(MIN_TOTAL, MAX_TOTAL+1)):
        upper_choice = min(coarse, key=lambda c: F(c['totals'][index]['failure_log2_upper']))
        lower_choice = min(ideal, key=lambda c: F(c['totals'][index]['method_expression_log2_lower']))
        ideal_upper_choice = min(ideal, key=lambda c: F(c['totals'][index]['method_expression_log2_upper']))
        upper = F(upper_choice['totals'][index]['failure_log2_upper'])
        lower = F(lower_choice['totals'][index]['method_expression_log2_lower'])
        ideal_upper = F(ideal_upper_choice['totals'][index]['method_expression_log2_upper'])
        coarse_sum += arb(2)**A.uniform.rational(upper)
        ideal_sum += arb(2)**A.uniform.rational(lower)
        selected.append(dict(total_arity=total, coarse_selected_order=upper_choice['moment_order'],
            coarse_failure_log2_upper=str(upper), ideal_lower_selected_order=lower_choice['moment_order'],
            ideal_upper_selected_order=ideal_upper_choice['moment_order'],
            ideal_method_expression_log2_lower=str(lower), ideal_method_expression_log2_upper=str(ideal_upper)))
    return dict(controls=controls, selected_totals=selected,
        coarse_all_totals_failure_log2_upper=str(A.uniform.rational_ceiling(A.log2(coarse_sum))),
        ideal_all_totals_method_expression_log2_lower=str(_floor(A.log2(ideal_sum))))


def receipt(precision=256):
    """One fresh uniform bank shared by all six prescribed controls."""
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Use fresh256/384-bit precision')
    pins, started = dependencies(), perf_counter()
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions)
        if len(terms) != 99:
            raise ArithmeticError('Expected all99 indexed source profiles')
        uniform_log = A.log2(sum((term for _, _, _, term in terms), arb(0)))
        result = controls_from_uniform_log(uniform_log)
        uniform_enclosure = str(uniform_log)
    if dependencies() != pins:
        raise ValueError('Reviewed sources changed during limits control')
    return dict(result, construction=binding.SCOPE['construction'], active_groups=R,
        every_group_arity_range=[1, 16], total_arity_range=[MIN_TOTAL, MAX_TOTAL],
        original_rowsets_fixed_across_columns=True, minimum_active_columns_per_group=9,
        restricted_density_not_renormalized=True, source_profiles=99,
        source_count_outside_theta=True, rowset_sum_inside_theta=True,
        projective_normalizations=1, field_order=256, source_dimension=1048576,
        output_length=2097152, forbidden_weight=734003, mark=str(MARK), uniform_alpha=str(ALPHA),
        uniform_source_sum_log2_enclosure=uniform_enclosure,
        coarse_is_failure_upper_bound=True, ideal_is_failure_bound=False,
        ideal_is_lower_comparator_only_for_this_fixed_norm_method=True,
        no_exact_moment_inventory_required=True, other_group_counts_covered=False,
        whole_code_certificate=False, precision_bits=precision, fresh_regional_banks=1,
        arithmetic_enclosures_directed=True, reviewed_source_sha256=pins,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), default=256)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.precision)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('controls', 'selected_totals', 'reviewed_source_sha256')}))
