"""Directed mixed-arity bridge with scheduled column likelihood moments.

All18 active groups have arities in12..16. Each original row set remains
fixed across16 columns. For common integer moment p, define
d_A=(product_j I_p(P_j(A^c)))^(1/p), w_a=sum_|A|=a d_A, and
beta=256^16/255^8. The raw likelihood comparison coefficient at total
arity k is beta^k [u^(k-216)](sum_a w_a u^(a-12))^18.

Use theta=(1/5)*p/(p-1), as proved by indicator Holder. The uniform
regional provider still has alpha1/5 and mark691/700. Different moment
orders supply alternative bounds for the same class; do not sum them.
Different total arities are disjoint and their failure bounds are summed.
No whole-code claim or timing change results from this proof controller.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

import hd_lambda_high_arity_certificate as binding
import hd_lambda_renyi_orbits as local

DIRECTORY = Path(__file__).resolve().parent
A = binding.arithmetic
R, MIN_TOTAL, MAX_TOTAL = 18, 216, 246
ALPHA, MARK = F(1, 5), F(691, 700)
ORDERS = (3, 4)
PINS = {
    'hd_lambda_high_arity_certificate.py': '156ca16d6f18370eb2484e6a1bf74ef83cbede670b63b2ba03db27f465d2e036',
    'hd_lambda_renyi_orbits.py': '08e828b3d227fca874bce4622266a994550686f9df7997933890baaaae44799a',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**binding.checked_dependencies(), **local.checked_dependencies(), **PINS,
        Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def checked_inventory(data):
    expected = dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=188,
        complement_sizes=[0, 1, 2, 3, 4], arities=[12, 13, 14, 15, 16],
        moment_orders=list(ORDERS), exact_rational_arithmetic=True, complete=True,
        full_parity_includes_identity=True, code='ker(Hlambda[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=48, source_pins_unchanged=True, expected_labelled_complements=2517,
        reviewed_source_sha256=local.checked_dependencies())
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Wrong exact orbit inventory scope: '+key)
    plan = [orbit for b in range(5) for orbit in local.complement_orbits(b)]
    expected_rows = {member: (index, orbit['representative'], image)
        for index, orbit in enumerate(plan) for member, image in orbit['members']}
    result, representative_values = {}, {}
    for row in data.get('rows', []):
        key = tuple(row['complement'])
        if key not in expected_rows or key in result:
            raise ValueError('Unexpected or repeated complement')
        index, representative, image = expected_rows[key]
        if (row.get('arity') != 16-len(key) or row.get('orbit_index') != index
                or row.get('representative') != list(representative)
                or row.get('representative_to_complement_permutation') != list(image)
                or set(row.get('moments', {})) != {'3', '4'}):
            raise ValueError('Wrong local coordinate-equivalence witness')
        values = {p: F(row['moments'][str(p)]) for p in ORDERS}
        if any(v < 1 for v in values.values()):
            raise ValueError('Invalid normalized density moment')
        if representative_values.setdefault(index, values) != values:
            raise ValueError('Orbit members disagree on moment values')
        result[key] = values
    if set(result) != set(expected_rows):
        raise ValueError('Incomplete2517-complement inventory')
    return result


def scheduled_weights(moments, order):
    """Arb weights retain every original labelled row set, without quotienting."""
    if order not in ORDERS:
        raise ValueError('Use the prescribed common order3 or4')
    values = {key: A.uniform.rational(row[order]) for key, row in moments.items()}
    schedules = [binding.column_permutations('input-inverse-output-xor', j)[0]
                 for j in range(16)]
    result = []
    for arity in range(12, 17):
        total = arb(0)
        for complement in combinations(range(16), 16-arity):
            product = arb(1)
            for before in schedules:
                mapped = tuple(sorted(before[i] for i in complement))
                product *= values[mapped]
            total += product**A.uniform.rational(F(1, order))
        result.append(total)
    return result


def total_coefficients(weights, groups=R):
    if len(weights) != 5 or type(groups) is not int or not 1 <= groups <= 18:
        raise ValueError('Five arity weights and1..18 groups required')
    if any(not x.is_finite() or not x > 0 for x in weights):
        raise ValueError('Strictly positive finite arity weights required')
    return arb_poly(weights)**groups


def receipt(data, precision=256):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Use a fresh256/384-bit replay')
    pins = dependencies()
    moments = checked_inventory(data)
    started = perf_counter()
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions)
        uniform_log = A.log2(sum((term for _, _, _, term in terms), arb(0)))
        threshold_log = -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))
        log_beta = 128-8*A.log2(arb(255))
        controls = []
        for order in ORDERS:
            theta = ALPHA*F(order, order-1)
            weights = scheduled_weights(moments, order)
            coefficients = total_coefficients(weights)
            rows = []
            for total in range(MIN_TOTAL, MAX_TOTAL+1):
                pref = A.uniform.rational(theta)*(total*log_beta
                    + A.log2(coefficients[total-MIN_TOTAL])-A.log2(arb(255)))+threshold_log
                value = pref+uniform_log
                rows.append(dict(total_arity=total,
                    failure_log2_upper=str(A.uniform.rational_ceiling(value))))
            controls.append(dict(moment_order=order, outer_theta=str(theta),
                arity_weight_log2_enclosures=[str(A.log2(x)) for x in weights], totals=rows))
        best = []
        cumulative = arb(0)
        for index, total in enumerate(range(MIN_TOTAL, MAX_TOTAL+1)):
            choice = min(controls, key=lambda control: F(control['totals'][index]['failure_log2_upper']))
            bound = F(choice['totals'][index]['failure_log2_upper'])
            cumulative += arb(2)**A.uniform.rational(bound)
            best.append(dict(total_arity=total, selected_moment_order=choice['moment_order'],
                failure_log2_upper=str(bound),
                cumulative_from216_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative)))))
    if dependencies() != pins:
        raise ValueError('Reviewed sources changed during replay')
    return dict(construction=binding.SCOPE['construction'], active_groups=18,
        every_group_arity_range=[12, 16], total_arity_range=[MIN_TOTAL, MAX_TOTAL],
        original_rowsets_fixed_across_columns=True, scheduled_moment_products=True,
        all2517_labelled_rowsets=True, source_profiles=99,
        source_count_outside_theta=True, rowset_sum_inside_theta=True,
        projective_normalizations=1, field_order=256, source_dimension=1048576,
        output_length=2097152, forbidden_weight=734003,
        mark=str(MARK), uniform_alpha=str(ALPHA), controls=controls, selected_totals=best,
        precision_bits=precision, fresh_regional_banks=1, directed=True,
        groups_below_arity12_covered=False, other_group_counts_covered=False,
        whole_code_certificate=False, reviewed_source_sha256=pins,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, required=True)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    raw = args.inventory.read_bytes()
    result = receipt(json.loads(raw), args.precision)
    result['inventory_sha256'] = sha256(raw).hexdigest()
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('controls', 'reviewed_source_sha256')}))
