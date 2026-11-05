"""Directed Renyi comparison retaining at least nine active outer columns.

An active RS[16,8] group has at least nine nonzero columns. For a fixed
original row set A, the zero-column contribution to E_U L_j^p is exactly
z=256^(p*(16-|A|)-16). Let n_j=I_p(A,j)-z. Restricting the auxiliary
density to at least nine active columns gives the positive moment tail
sum_{k>=9} [v^k] product_j(z+n_j*v). Do not renormalize this density.

This controller covers groups of arity11..16, retains every labelled
row set across the actual16-column schedule, and uses a common p before
taking a mixed-arity coefficient. Order2 uses the safe I2<=I3 bound for
arities12..16; orders3/4 and arity11 use exact rational moments.
It is a partial message-class certificate, never a whole-code claim.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

import hd_lambda_renyi_middle as middle
import hd_lambda_renyi_mixed as dense

DIRECTORY = Path(__file__).resolve().parent
A, binding = dense.A, dense.binding
ALPHA, MARK = dense.ALPHA, dense.MARK
ORDERS, ARITIES, R = (2, 3, 4), tuple(range(11, 17)), 18
PINS = {
    'hd_lambda_renyi_mixed.py': '9dc5d2cfba096bab6e3fc3b2ceb39425f62450ea0364a53407b7eff3405c1723',
    'hd_lambda_renyi_middle.py': '1e43e08877c2f23195605c293d810f40063caee10073abf8383a7388f63fb115',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**dense.dependencies(), **middle.checked_dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def checked_middle(data):
    expected = dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=188,
        complement_sizes=[5], arities=[11], moment_orders=list(ORDERS),
        exact_rational_arithmetic=True, complete=True, source_pins_unchanged=True,
        full_parity_includes_identity=True, code='ker(Hlambda[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=48, expected_labelled_complements=4368,
        reviewed_source_sha256=middle.checked_dependencies())
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Wrong exact middle inventory scope: '+key)
    plan = middle.complement_orbits(5)
    expected_rows = {member: (i, row['representative'], image)
        for i, row in enumerate(plan) for member, image in row['members']}
    result, representatives = {}, {}
    for row in data.get('rows', []):
        key = tuple(row['complement'])
        if key not in expected_rows or key in result:
            raise ValueError('Unexpected or duplicate complement')
        index, representative, image = expected_rows[key]
        if (row.get('arity') != 11 or row.get('orbit_index') != index
                or row.get('representative') != list(representative)
                or row.get('representative_to_complement_permutation') != list(image)
                or set(row.get('moments', {})) != {'2', '3', '4'}):
            raise ValueError('Wrong middle coordinate-map witness')
        if any(type(v) is not str for v in row['moments'].values()):
            raise ValueError('Exact rational moments must be strings')
        values = {p: F(row['moments'][str(p)]) for p in ORDERS}
        if any(x < 1 for x in values.values()):
            raise ValueError('Normalized likelihood moment is below one')
        if representatives.setdefault(index, values) != values:
            raise ValueError('Orbit members disagree')
        result[key] = values
    if set(result) != set(expected_rows):
        raise ValueError('Incomplete middle inventory')
    return result


def zero_atom(arity, order, q=256, width=16):
    if (type(arity) is not int or not 1 <= arity <= width
            or type(order) is not int or order < 2):
        raise ValueError('Positive arity and integer order at least two required')
    return F(q)**(order*(width-arity)-width)


def activity_tail(zero, nonzero, minimum=9):
    """Positive polynomial tail: never subtract almost equal large totals."""
    if (type(minimum) is not int or not 0 <= minimum <= len(nonzero)
            or zero < 0 or any(x < 0 for x in nonzero)):
        raise ValueError('Nonnegative atoms and a feasible activity floor required')
    poly = arb_poly([1])
    for atom in nonzero:
        poly *= arb_poly([A.uniform.rational(zero), A.uniform.rational(atom)])
    return sum((poly[k] for k in range(minimum, len(nonzero)+1)), arb(0))


def scheduled_weights(moments, order, minimum=9, arities=ARITIES):
    if order not in ORDERS or not arities or any(a not in ARITIES for a in arities):
        raise ValueError('Use order2/3/4 and arities11..16')
    schedules = [binding.column_permutations('input-inverse-output-xor', j)[0]
                 for j in range(16)]
    weights, diagnostics = [], []
    for arity in arities:
        z = zero_atom(arity, order)
        signatures = Counter()
        for complement in combinations(range(16), 16-arity):
            signature = tuple(sorted(moments[tuple(sorted(before[i] for i in complement))][order]
                                     for before in schedules))
            signatures[signature] += 1
        total = arb(0)
        for signature, multiplicity in signatures.items():
            # This subtraction is exact even when the zero atom dominates I_p.
            nonzero = tuple(value-z for value in signature)
            tail = activity_tail(z, nonzero, minimum)
            total += multiplicity*tail**A.uniform.rational(F(1, order))
        if not total > 0 or not total.is_finite():
            raise ArithmeticError('Nonpositive or nonfinite rowset weight')
        weights.append(total)
        diagnostics.append(dict(arity=arity, labelled_rowsets=sum(signatures.values()),
                                distinct_moment_signatures=len(signatures)))
    return weights, diagnostics


def receipt(small_data, middle_data, precision=256):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh256/384-bit replay required')
    pins = dependencies()
    moments = dense.checked_inventory(small_data)
    for row in moments.values():
        row[2] = row[3]  # E L^2 <= (E L^3)^(1/2) <= E L^3, since E L=1.
    moments.update(checked_middle(middle_data))
    started = perf_counter()
    minimum_total, maximum_total = 11*R, 245
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions)
        if len(terms) != 99:
            raise ArithmeticError('Expected all99 source profiles')
        uniform_log = A.log2(sum((term for _, _, _, term in terms), arb(0)))
        threshold_log = -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))
        log_beta = 128-8*A.log2(arb(255))
        controls = []
        for order in ORDERS:
            theta = ALPHA*F(order, order-1)
            weights, diagnostics = scheduled_weights(moments, order)
            coefficients = arb_poly(weights)**R
            rows = []
            for total in range(minimum_total, maximum_total+1):
                pref = A.uniform.rational(theta)*(total*log_beta
                    + A.log2(coefficients[total-minimum_total])-A.log2(arb(255)))+threshold_log
                rows.append(dict(total_arity=total,
                    failure_log2_upper=str(A.uniform.rational_ceiling(pref+uniform_log))))
            controls.append(dict(moment_order=order, outer_theta=str(theta),
                arity_weight_log2_enclosures=[str(A.log2(x)) for x in weights],
                rowset_diagnostics=diagnostics, totals=rows))
        selected, cumulative = [], arb(0)
        for index, total in enumerate(range(minimum_total, maximum_total+1)):
            choice = min(controls, key=lambda c: F(c['totals'][index]['failure_log2_upper']))
            bound = F(choice['totals'][index]['failure_log2_upper'])
            cumulative += arb(2)**A.uniform.rational(bound)
            selected.append(dict(total_arity=total, selected_moment_order=choice['moment_order'],
                failure_log2_upper=str(bound),
                cumulative_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative)))))
    if dependencies() != pins:
        raise ValueError('Sources changed during replay')
    return dict(construction=binding.SCOPE['construction'], active_groups=R,
        every_group_arity_range=[11, 16], total_arity_range=[minimum_total, maximum_total],
        minimum_active_columns_per_group=9, restricted_density_not_renormalized=True,
        exact_zero_atom_subtraction=True, order2_large_arity_uses_I3_upper=True,
        original_rowsets_fixed_across_columns=True, all6885_labelled_rowsets=True,
        source_profiles=99, source_count_outside_theta=True, rowset_sum_inside_theta=True,
        projective_normalizations=1, field_order=256, source_dimension=1048576,
        output_length=2097152, forbidden_weight=734003, mark=str(MARK), uniform_alpha=str(ALPHA),
        controls=controls, selected_totals=selected, precision_bits=precision,
        fresh_regional_banks=1, directed=True, whole_code_certificate=False,
        groups_below_arity11_covered=False, other_group_counts_covered=False,
        reviewed_source_sha256=pins, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--small-inventory', type=Path, required=True)
    parser.add_argument('--middle-inventory', type=Path, required=True)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    small_raw, middle_raw = args.small_inventory.read_bytes(), args.middle_inventory.read_bytes()
    result = receipt(json.loads(small_raw), json.loads(middle_raw), args.precision)
    result['inventory_sha256'] = dict(small=sha256(small_raw).hexdigest(), middle=sha256(middle_raw).hexdigest())
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('controls', 'reviewed_source_sha256')}))
