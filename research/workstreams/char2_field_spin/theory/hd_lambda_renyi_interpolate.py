"""Directed interpolation of nonzero column moments before the activity tail.

For a fixed mapped source A and 3 <= p <= 4, let n_k=E_U[L^k;Y!=0].
Holder gives n_p <= n_3^(4-p)*n_4^(p-3). The zero contribution remains
exactly z_p=256^(p*(16-|A|)-16). Apply these bounds before multiplying
the sixteen scalar column polynomials z_p+n_p*v and retaining degree>=9.
This bounds the p-th moment of an UNNORMALIZED restricted density.

Every original A is held fixed through its actual scheduled columns.
Sorting paired (n_3,n_4) signatures only reuses this symmetric scalar
moment calculation; it does not quotient source labels or inner routes.
One common p is used for each mixed-arity coefficient. The uniform
provider keeps alpha=1/5 and Z=691/700; theta=alpha*p/(p-1).
This isolated controller covers only r18, arities11..16, totals198..245.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

import hd_lambda_renyi_activity as activity

DIRECTORY = Path(__file__).resolve().parent
A, binding = activity.A, activity.binding
ALPHA, MARK, R, ARITIES = activity.ALPHA, activity.MARK, 18, tuple(range(11, 17))
ORDERS = tuple(F(k, 8) for k in range(25, 32))
ALLOWED_ORDERS = (F(3),)+ORDERS+(F(4),)
MINIMUM_TOTAL, MAXIMUM_TOTAL = 198, 245
PINS = {'hd_lambda_renyi_activity.py':
        '09f5245aecb903097438b65ded8757c6a66da40a13804289c394fef777186c77'}
INVENTORY_PINS = {
    'small': '227fd5c05f5946698d6fac8a703abe4657ec154f817f48352d0ce346f818e673',
    'middle': 'f08e68c79dea8599f2a63c8afef67dbe10eede78b56a4b3d33f0be03bc1842ba',
}
DEFAULT_SMALL = DIRECTORY/'hd_lambda_renyi_orbit_inventory34.json'
DEFAULT_MIDDLE = DIRECTORY/'hd_lambda_renyi_a11_inventory234.json'


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen activity source changed: '+name)
    return {**activity.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def _orders(orders):
    orders = tuple(orders)
    if (not orders or any(type(p) not in (int, F) for p in orders)
            or len(set(orders)) != len(orders) or any(p not in ALLOWED_ORDERS for p in orders)):
        raise ValueError('Use distinct prescribed eighth-step orders from3 through4')
    return tuple(sorted(F(p) for p in orders))


def checked_inputs(small_path=DEFAULT_SMALL, middle_path=DEFAULT_MIDDLE):
    """Read externally pinned inventories and validate every labelled witness."""
    tables, paths = {}, {'small': Path(small_path), 'middle': Path(middle_path)}
    for name, path in paths.items():
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != INVENTORY_PINS[name]:
            raise ValueError('Unreviewed exact inventory: '+name)
        tables[name] = json.loads(raw)
        for row in tables[name].get('rows', []):
            if any(type(value) is not str for value in row.get('moments', {}).values()):
                raise ValueError('Exact moments must be rational strings')
    moments = activity.dense.checked_inventory(tables['small'])
    moments.update(activity.checked_middle(tables['middle']))
    return moments, paths


def inputs_unchanged(paths):
    return all(sha256(path.read_bytes()).hexdigest() == INVENTORY_PINS[name]
               for name, path in paths.items())


def nonzero_interpolation(n3, n4, order):
    """Arb enclosure of the Holder upper bound, not an exact p-moment."""
    order = F(order)
    if not 3 <= order <= 4 or n3 < 0 or n4 < 0:
        raise ValueError('Nonnegative nonzero moments and3<=p<=4 required')
    if order == 3:
        return A.uniform.rational(n3)
    if order == 4:
        return A.uniform.rational(n4)
    return (A.uniform.rational(n3)**A.uniform.rational(4-order)
            * A.uniform.rational(n4)**A.uniform.rational(order-3))


def zero_contribution(arity, order):
    order = F(order)
    if type(arity) is not int or arity not in ARITIES or not 3 <= order <= 4:
        raise ValueError('Use arities11..16 and3<=p<=4')
    return arb(256)**A.uniform.rational(order*(16-arity)-16)


def positive_tail(zero, nonzero, minimum=9):
    """A positive scalar polynomial tail; inputs may be Arb enclosures."""
    if type(minimum) is not int or not 0 <= minimum <= len(nonzero):
        raise ValueError('Invalid activity floor')
    if zero < 0 or any(value < 0 for value in nonzero):
        raise ValueError('Negative moment atom')
    polynomial = arb_poly([1])
    for value in nonzero:
        polynomial *= arb_poly([zero, value])
    return sum((polynomial[k] for k in range(minimum, len(nonzero)+1)), arb(0))


def prepare_signatures(moments, arities=ARITIES):
    """Pair the two moments of each ACTUAL column before sorting signatures."""
    arities = tuple(arities)
    if not arities or len(set(arities)) != len(arities) or any(a not in ARITIES for a in arities):
        raise ValueError('Use distinct arities11..16')
    schedules = [binding.column_permutations('input-inverse-output-xor', j)[0]
                 for j in range(16)]
    result = {}
    for arity in arities:
        z3, z4 = activity.zero_atom(arity, 3), activity.zero_atom(arity, 4)
        mapped_pairs = {}
        for complement in combinations(range(16), 16-arity):
            values = moments[complement]
            pair = values[3]-z3, values[4]-z4
            if pair[0] < 0 or pair[1] < 0:
                raise ValueError('A supplied moment is below its exact zero contribution')
            mapped_pairs[complement] = pair
        signatures = Counter()
        for complement in combinations(range(16), 16-arity):
            pairs = [mapped_pairs[tuple(sorted(before[i] for i in complement))]
                     for before in schedules]
            signatures[tuple(sorted(pairs))] += 1
        result[arity] = signatures
    return result


def scheduled_weights(prepared, order, minimum=9):
    """Sum complete restricted p-norm bounds over every original row set."""
    order = _orders((order,))[0]
    weights, diagnostics = [], []
    for arity, signatures in prepared.items():
        zero = zero_contribution(arity, order)
        unique_pairs = {pair for signature in signatures for pair in signature}
        atoms = {pair: nonzero_interpolation(*pair, order) for pair in unique_pairs}
        total = arb(0)
        for signature, multiplicity in signatures.items():
            tail = positive_tail(zero, tuple(atoms[pair] for pair in signature), minimum)
            total += multiplicity*tail**A.uniform.rational(1/order)
        if not total > 0 or not total.is_finite():
            raise ArithmeticError('Nonpositive or nonfinite rowset weight')
        weights.append(total)
        diagnostics.append(dict(arity=arity, labelled_rowsets=sum(signatures.values()),
                                distinct_joint_moment_signatures=len(signatures)))
    return weights, diagnostics


def receipt(small_path=DEFAULT_SMALL, middle_path=DEFAULT_MIDDLE, precision=256,
            *, orders=ORDERS, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Use a fresh256/384-bit replay')
    orders, pins = _orders(orders), dependencies()
    moments, paths = checked_inputs(small_path, middle_path)
    prepared = prepare_signatures(moments)
    started = perf_counter()
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions)
        if len(terms) != 99:
            raise ArithmeticError('Expected all99 source profiles')
        uniform_log = A.log2(sum((term for _, _, _, term in terms), arb(0)))
        threshold_log = -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))
        log_beta = 128-8*A.log2(arb(255))
        controls = []
        for order in orders:
            theta = ALPHA*order/(order-1)
            weights, diagnostics = scheduled_weights(prepared, order)
            coefficients = arb_poly(weights)**R
            rows = []
            for total in range(MINIMUM_TOTAL, MAXIMUM_TOTAL+1):
                prefix = A.uniform.rational(theta)*(total*log_beta
                    + A.log2(coefficients[total-MINIMUM_TOTAL])-A.log2(arb(255)))+threshold_log
                rows.append(dict(total_arity=total,
                    failure_log2_upper=str(A.uniform.rational_ceiling(prefix+uniform_log))))
            controls.append(dict(moment_order=str(order), outer_theta=str(theta),
                arity_weight_log2_enclosures=[str(A.log2(value)) for value in weights],
                rowset_diagnostics=diagnostics, totals=rows))
            if progress:
                progress(dict(completed_moment_order=str(order), elapsed_seconds=perf_counter()-started))
        selected, cumulative = [], arb(0)
        for index, total in enumerate(range(MINIMUM_TOTAL, MAXIMUM_TOTAL+1)):
            choice = min(controls, key=lambda row: F(row['totals'][index]['failure_log2_upper']))
            bound = F(choice['totals'][index]['failure_log2_upper'])
            cumulative += arb(2)**A.uniform.rational(bound)
            selected.append(dict(total_arity=total, selected_moment_order=choice['moment_order'],
                failure_log2_upper=str(bound),
                cumulative_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative)))))
    if dependencies() != pins or not inputs_unchanged(paths):
        raise ValueError('Reviewed sources or inventories changed during replay')
    return dict(construction=binding.SCOPE['construction'], active_groups=R,
        every_group_arity_range=[11, 16], total_arity_range=[MINIMUM_TOTAL, MAXIMUM_TOTAL],
        moment_orders=[str(p) for p in orders], minimum_active_columns_per_group=9,
        nonzero_interpolation='n_p<=n_3^(4-p)*n_4^(p-3), before the activity tail',
        zero_contribution='256^(p*(16-arity)-16)', restricted_density_not_renormalized=True,
        original_rowsets_fixed_across_columns=True, paired_moments_preserved_in_signatures=True,
        all6885_labelled_rowsets=True, source_profiles=99, source_count_outside_theta=True,
        rowset_sum_inside_theta=True, projective_normalizations=1, field_order=256,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        mark=str(MARK), uniform_alpha=str(ALPHA), controls=controls, selected_totals=selected,
        precision_bits=precision, fresh_regional_banks=1, directed=True, whole_code_certificate=False,
        groups_below_arity11_covered=False, other_group_counts_covered=False,
        reviewed_source_sha256=pins, inventory_sha256=dict(INVENTORY_PINS),
        source_pins_unchanged=True, inventory_pins_unchanged=True, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--small-inventory', type=Path, default=DEFAULT_SMALL)
    parser.add_argument('--middle-inventory', type=Path, default=DEFAULT_MIDDLE)
    parser.add_argument('--orders', default=','.join(map(str, ORDERS)))
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.small_inventory, args.middle_inventory, args.precision,
        orders=tuple(F(value) for value in args.orders.split(',')),
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('controls', 'reviewed_source_sha256')}))
