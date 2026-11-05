"""Fresh directed activity comparison for all group arities10..16 at r18.

The construction is the UNCHANGED scheduled Hlambda, not the new double
transvection. The full exact arity10 inventory is combined with the frozen
arity11 and12..16 inventories. Every original row set remains fixed through
the actual16-column schedule. The nonzero-column floor is nine. Beta is
charged once per active RS row, the projective divisor once inside theta,
and all99 indexed source profiles once outside theta. This is partial
class coverage, never a whole-code certificate.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx
import hd_lambda_renyi_activity as activity

DIRECTORY = Path(__file__).resolve().parent
A, R, ORDERS = activity.A, 18, (2, 3, 4)
ARITIES = tuple(range(10, 17))
INPUTS = {
    'hd_lambda_renyi_orbit_inventory34.json': '227fd5c05f5946698d6fac8a703abe4657ec154f817f48352d0ce346f818e673',
    'hd_lambda_renyi_a11_inventory234.json': 'f08e68c79dea8599f2a63c8afef67dbe10eede78b56a4b3d33f0be03bc1842ba',
    'hd_lambda_renyi_a10_inventory234.json': 'd40ea671e1b1ff8056241cdbe3006efb3b0b92f14ed0c41f8a5e7f95b415784b',
}
ACTIVITY_SHA256 = '09f5245aecb903097438b65ded8757c6a66da40a13804289c394fef777186c77'


def dependencies():
    if sha256((DIRECTORY/'hd_lambda_renyi_activity.py').read_bytes()).hexdigest() != ACTIVITY_SHA256:
        raise ValueError('Frozen activity helper changed')
    return {**activity.dependencies(), Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def checked_a10(data):
    expected = dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=188,
        complement_sizes=[6], arities=[10], moment_orders=list(ORDERS),
        exact_rational_arithmetic=True, complete=True, source_pins_unchanged=True,
        full_parity_includes_identity=True, code='ker(Hlambda[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=48, expected_labelled_complements=8008,
        expected_representatives=272, actual_support_tables=272,
        reviewed_source_sha256=activity.middle.checked_dependencies())
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Wrong arity10 inventory scope: '+key)
    plan = activity.middle.complement_orbits(6)
    expected_rows = {member: (i, row['representative'], image)
        for i, row in enumerate(plan) for member, image in row['members']}
    result, representatives = {}, {}
    for row in data.get('rows', []):
        key = tuple(row['complement'])
        if any(type(i) is not int for i in key) or key not in expected_rows or key in result:
            raise ValueError('Unexpected or duplicate complement')
        index, representative, image = expected_rows[key]
        if (type(row.get('arity')) is not int or row['arity'] != 10
                or row.get('orbit_index') != index or row.get('representative') != list(representative)
                or row.get('representative_to_complement_permutation') != list(image)
                or set(row.get('moments', {})) != {'2', '3', '4'}):
            raise ValueError('Wrong coordinate-map witness')
        if any(type(v) is not str for v in row['moments'].values()):
            raise ValueError('Exact moment strings required')
        values = {p: F(row['moments'][str(p)]) for p in ORDERS}
        if any(values[p] < max(F(1), activity.zero_atom(10, p)) for p in ORDERS):
            raise ValueError('Moment below its zero atom or normalization')
        if representatives.setdefault(index, values) != values:
            raise ValueError('Orbit members disagree')
        result[key] = values
    if set(result) != set(expected_rows):
        raise ValueError('Incomplete labelled arity10 inventory')
    return result


def inputs():
    data = []
    for name, digest in INPUTS.items():
        raw = (DIRECTORY/name).read_bytes()
        if sha256(raw).hexdigest() != digest:
            raise ValueError('Observed inventory changed: '+name)
        data.append(json.loads(raw))
    result = activity.dense.checked_inventory(data[0])
    for row in result.values():
        row[2] = row[3]  # Safe upper, not an assertion that the moments coincide.
    result.update(activity.checked_middle(data[1]))
    result.update(checked_a10(data[2]))
    return result


def scheduled_weights(moments, order):
    if order not in ORDERS:
        raise ValueError('Moment order2,3,4 required')
    schedules = [activity.binding.column_permutations('input-inverse-output-xor', j)[0]
                 for j in range(16)]
    weights, diagnostics = [], []
    for arity in ARITIES:
        signatures = Counter(tuple(sorted(moments[tuple(sorted(before[i] for i in complement))][order]
                                  for before in schedules))
            for complement in combinations(range(16), 16-arity))
        z, total = activity.zero_atom(arity, order), arb(0)
        for signature, multiplicity in signatures.items():
            tail = activity.activity_tail(z, tuple(value-z for value in signature), 9)
            total += multiplicity*tail**A.uniform.rational(F(1, order))
        if not total > 0 or not total.is_finite() or sum(signatures.values()) != comb(16, arity):
            raise ArithmeticError('Complete positive finite rowset sum required')
        weights.append(total)
        diagnostics.append(dict(arity=arity, labelled_rowsets=sum(signatures.values()),
                                distinct_moment_signatures=len(signatures)))
    return weights, diagnostics


def receipt(precision):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh256/384-bit replay required')
    pins, moments, started = dependencies(), inputs(), perf_counter()
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(activity.MARK, activity.ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions)
        if len(terms) != 99 or {(e, o): n for e, o, n, _ in terms} != A.profile_counts(R):
            raise ArithmeticError('All99 indexed profiles required')
        uniform_log = A.log2(sum((value for _, _, _, value in terms), arb(0)))
        threshold = -A.uniform.rational(activity.ALPHA)*734003*A.log2(A.uniform.rational(activity.MARK))
        log_beta = 128-8*A.log2(arb(255))
        controls = []
        for order in ORDERS:
            weights, diagnostics = scheduled_weights(moments, order)
            coefficients = arb_poly(weights)**R
            theta = activity.ALPHA*F(order, order-1)
            rows = []
            for total in range(180, 246):
                bound = (A.uniform.rational(theta)*(total*log_beta
                    + A.log2(coefficients[total-180])-A.log2(arb(255)))+threshold+uniform_log)
                rows.append(dict(total_arity=total,
                    failure_log2_upper=str(A.uniform.rational_ceiling(bound))))
            controls.append(dict(moment_order=order, outer_theta=str(theta),
                diagnostics=diagnostics, arity_weight_log2_enclosures=[str(A.log2(w)) for w in weights],
                totals=rows))
        selected, cumulative = [], arb(0)
        for index, total in enumerate(range(180, 246)):
            chosen = min(controls, key=lambda c: F(c['totals'][index]['failure_log2_upper']))
            bound = F(chosen['totals'][index]['failure_log2_upper'])
            cumulative += arb(2)**A.uniform.rational(bound)
            selected.append(dict(total_arity=total, selected_moment_order=chosen['moment_order'],
                failure_log2_upper=str(bound),
                cumulative_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative)))))
    if dependencies() != pins:
        raise ValueError('Sources changed during replay')
    for name, digest in INPUTS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Inventory changed during replay')
    return dict(construction=activity.binding.SCOPE['construction'], active_groups=18,
        every_group_arity_range=[10, 16], total_arity_range=[180, 245],
        all14893_labelled_rowsets=True, original_rowsets_fixed_across_columns=True,
        minimum_active_columns_per_group=9, restricted_density_not_renormalized=True,
        source_profiles=99, source_count_outside_theta=True, rowset_sum_inside_theta=True,
        projective_normalizations=1, order2_large_arities_use_I3_upper=True,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        field_order=256, mark=str(activity.MARK), uniform_alpha=str(activity.ALPHA),
        controls=controls, selected_totals=selected, inventory_sha256=INPUTS,
        reviewed_source_sha256=pins, precision_bits=precision, directed=True,
        fresh_regional_banks=1, whole_code_certificate=False, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.precision)
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('controls', 'reviewed_source_sha256')}))
