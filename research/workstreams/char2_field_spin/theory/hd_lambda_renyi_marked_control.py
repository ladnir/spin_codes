"""Directed RS-floor comparison from the pinned exact marked inventory.

For each original source set A, retain its sixteen actual scheduled images.
At common order p, let T_A be the coefficient tail of product_j J_(P_j A,p)
at degree at least 9*|A|*p. Its p-th root bounds the raw restricted-density
norm. The per-group coefficient is

    B_a(p) = 255^(-8*a) sum_{|A|=a} T_A^(1/p).

These are raw moments: no beta or additional power of 256 belongs here.
Only scalar polynomial multiplication is commutative. Equal sorted tuples
of sixteen local representatives share arithmetic, not source multiplicity.

The tail is computed positively in reversed degree, truncated to deficit
7*a*p. Every coefficient and root is rounded to an Arb upper endpoint.
No midpoint, total-minus-prefix subtraction, or FLOAT coefficient is used.

For r=18 and one common p, the homogeneous class has coefficient B_a^r;
the class allowing every group arity1..5 has coefficient (sum_a B_a)^r.
Indicator Holder gives theta=p/(5*(p-1)). The existing uniform provider
uses mark691/700 and alpha1/5 unchanged. Its indexed-source union remains
outside theta; the single projective factor255 is raised to theta.
Homogeneous bounds are diagnostics contained in the mixed class, not
additional events to union. This is not a whole-code certificate.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_high_arity_certificate as binding
import hd_lambda_renyi_marked as local

DIRECTORY = Path(__file__).resolve().parent
A = binding.arithmetic
ARITIES, ORDERS, R = (1, 2, 3, 4, 5), (2, 3, 4), 18
ALPHA, MARK = F(1, 5), F(691, 700)
INVENTORY_NAME = 'hd_lambda_renyi_marked_inventory234.json'
INVENTORY_SHA256 = 'db3f7d554cf9e6168d29d660e476610c2d0934a359a4a7edd802bbf89873ecfe'
PINS = {
    'hd_lambda_high_arity_certificate.py': '156ca16d6f18370eb2484e6a1bf74ef83cbede670b63b2ba03db27f465d2e036',
    'hd_lambda_renyi_marked.py': '4808a6a1ebdc5ee3535f0c96fad3e0dd9be9b9873c90c43bbbec591b8c9b57d6',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**binding.checked_dependencies(), **local.checked_dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def _same(value, expected):
    if type(value) is not type(expected):
        return False
    if isinstance(expected, list):
        return len(value) == len(expected) and all(_same(v, e) for v, e in zip(value, expected))
    if isinstance(expected, dict):
        return set(value) == set(expected) and all(_same(value[k], e) for k, e in expected.items())
    return value == expected


def _require(data, expected):
    if not isinstance(data, dict):
        raise ValueError('A receipt object is required')
    for key, value in expected.items():
        if not _same(data.get(key), value):
            raise ValueError('Wrong marked inventory field: '+key)


def checked_inventory(data):
    """Validate complete scope, exact coefficients, pins, and every map witness.

    This validates structure without redoing the exact inventory arithmetic.
    Production callers also use load_inventory, which pins the observed file.
    """
    plan = tuple((a, orbit) for a in ARITIES for orbit in local.source_orbits(a))
    expected = dict(field_order=256, field_polynomial='0x11b', lambda_aes_byte=188,
        width=16, arities=list(ARITIES), moment_orders=list(ORDERS),
        coefficients_order='ascending original input degree 0..arity*moment',
        formula='J_p(X)=256^(16*(p-1))*sum_S K_S(X)^p/255^((p-1)*popcount(S))',
        raw_likelihood=True, normalized_column_identity='J_p(1)=256^(arity*p)*I_p',
        exact_rational_arithmetic=True, complete=True, stop_reason=None,
        subgroup_size=48, symmetry='bit permutations and i->i+parity(i)*1111',
        no_larger_symmetry_claim=True, whole_word_chronology_quotiented=False,
        original_source_multiplicity_removed=False, rs_normalization_included=False,
        scheduled_columns_composed=False, input_floor_applied=False,
        source_pins_unchanged=True, regional_provider_runs=0, whole_code_certificate=False,
        expected_labelled_sources=6884, completed_labelled_sources=6884,
        expected_representatives=len(plan), completed_representatives=len(plan),
        maps=[list(image) for image in local.orbits.orthogonal_maps()],
        reviewed_source_sha256=local.checked_dependencies())
    _require(data, expected)
    representatives = data.get('representatives')
    if not isinstance(representatives, list) or len(representatives) != len(plan):
        raise ValueError('Incomplete marked representative inventory')
    polynomials = []
    for index, ((arity, orbit), row) in enumerate(zip(plan, representatives)):
        _require(row, dict(orbit_index=index, arity=arity,
            representative=list(orbit['representative']), multiplicity=len(orbit['members']),
            raw_allowed_input_mass=256**arity))
        if (type(row.get('joint_entries')) is not int or type(row.get('output_masks')) is not int
                or not 2 <= row['output_masks'] <= min(65536, 256**arity)
                or not row['output_masks'] <= row['joint_entries'] <= (arity+1)*row['output_masks']):
            raise ValueError('Wrong representative count diagnostics')
        coefficients = row.get('coefficients')
        if not isinstance(coefficients, dict) or set(coefficients) != {'2', '3', '4'}:
            raise ValueError('Exact coefficients for all three orders are required')
        values = {}
        for p in ORDERS:
            raw = coefficients[str(p)]
            if (not isinstance(raw, list) or len(raw) != arity*p+1
                    or any(type(x) is not str for x in raw)):
                raise ValueError('Wrong exact polynomial degree or coefficient type')
            try:
                value = tuple(F(x) for x in raw)
            except (ValueError, ZeroDivisionError) as error:
                raise ValueError('Malformed exact coefficient') from error
            if (any(c < 0 for c in value) or value[0] != 256**(16*(p-1))
                    or any(value[d] for d in range(1, p)) or value[-1] <= 0
                    or sum(value) < 256**(arity*p)):
                raise ValueError('Invalid raw moment mass or zero/input-degree branch')
            values[p] = value
        polynomials.append(values)
    expected_rows = {member: (a, i, image) for i, (a, orbit) in enumerate(plan)
                     for member, image in orbit['members']}
    rows = data.get('rows')
    if not isinstance(rows, list) or len(rows) != 6884:
        raise ValueError('All 6884 labelled mapped sets are required')
    source_index = {}
    for row in rows:
        source = row.get('mapped_source') if isinstance(row, dict) else None
        if (not isinstance(source, list) or any(type(i) is not int for i in source)):
            raise ValueError('Integer mapped-source labels are required')
        key = tuple(source)
        if key not in expected_rows or key in source_index:
            raise ValueError('Unexpected or repeated mapped source')
        arity, index, image = expected_rows[key]
        _require(row, dict(arity=arity, orbit_index=index,
            representative_to_source_permutation=list(image)))
        source_index[key] = index
    if set(source_index) != set(expected_rows):
        raise ValueError('Incomplete labelled source coverage')
    return dict(polynomials=tuple(polynomials), source_index=source_index,
                inventory_sha256=INVENTORY_SHA256, labelled_sources=6884)


def load_inventory(path=DIRECTORY/INVENTORY_NAME):
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != INVENTORY_SHA256:
        raise ValueError('The observed exact marked inventory hash changed')
    return checked_inventory(json.loads(raw))


def _upper(value):
    if not value.is_finite() or not value >= 0:
        raise ArithmeticError('A finite nonnegative directed value is required')
    return value.upper()


def exact_upper(value):
    """Export a dyadic upper endpoint, never a printed midpoint."""
    endpoint = _upper(value).fmpq()
    return F(int(endpoint.p), int(endpoint.q))


def _convolve_positive(left, right, maximum):
    """Direct nonnegative convolution, truncated only in increasing deficit."""
    out = [arb(0) for _ in range(min(maximum+1, len(left)+len(right)-1))]
    # Explicit sparse terms avoid operations on the exact degree gaps.
    rterms = tuple((j, y) for j, y in enumerate(right) if not y.is_zero())
    for i, x in enumerate(left):
        if not x.is_zero():
            for j, y in rterms:
                if i+j >= len(out):
                    break
                out[i+j] += x*y
    return tuple(_upper(value) for value in out)


def positive_coefficient_tail(polynomials, minimum):
    """Upper bound a positive coefficient tail by reversed-degree convolution.

    Inputs are nonempty exact nonnegative Fraction polynomials. Trailing
    zeros are permitted: their declared degrees still define the reversal.
    The product has no negative deficits, so truncation loses no tail term.
    """
    polynomials = tuple(tuple(F(c) for c in poly) for poly in polynomials)
    if (not polynomials or any(not poly or any(c < 0 for c in poly) for poly in polynomials)
            or type(minimum) is not int or minimum < 0):
        raise ValueError('Nonnegative exact polynomials and a nonnegative degree floor required')
    maximum_degree = sum(len(poly)-1 for poly in polynomials)
    if minimum > maximum_degree:
        return arb(0)
    deficit = maximum_degree-minimum
    value = (arb(1),)
    for poly in polynomials:
        reversed_upper = tuple(_upper(A.uniform.rational(c)) for c in reversed(poly))[:deficit+1]
        value = _convolve_positive(value, reversed_upper, deficit)
    return _upper(sum(value, arb(0)))


def scheduled_signatures(checked, arity):
    """Count every original A once; cache only commutative scalar products."""
    if type(arity) is not int or arity not in ARITIES:
        raise ValueError('Arity1..5 required')
    schedules = tuple(binding.column_permutations('input-inverse-output-xor', j)[0]
                      for j in range(16))
    signatures = Counter()
    for original in combinations(range(16), arity):
        signature = tuple(sorted(checked['source_index'][tuple(sorted(before[i] for i in original))]
                                 for before in schedules))
        signatures[signature] += 1
    if sum(signatures.values()) != comb(16, arity):
        raise ArithmeticError('Original rowset multiplicity was lost')
    return signatures


def scheduled_weights(checked, order, *, arities=ARITIES, check_time=None, progress=None):
    """Return reusable B_a upper endpoints and labelled-signature diagnostics."""
    if (type(order) is not int or order not in ORDERS or not arities
            or len(set(arities)) != len(arities) or any(type(a) is not int or a not in ARITIES for a in arities)):
        raise ValueError('Use common order2/3/4 and distinct arities1..5')
    weights, diagnostics = {}, []
    for arity in arities:
        signatures = scheduled_signatures(checked, arity)
        total = arb(0)
        for index, (signature, multiplicity) in enumerate(signatures.items()):
            if check_time:
                check_time()
            polys = tuple(checked['polynomials'][i][order] for i in signature)
            tail = positive_coefficient_tail(polys, 9*arity*order)
            root = _upper(tail.root(order))
            total += multiplicity*root
            if progress and (index == 0 or (index+1) % 64 == 0 or index+1 == len(signatures)):
                progress(dict(stage='marked-word-tails', order=order, arity=arity,
                    completed_signatures=index+1, total_signatures=len(signatures)))
        weight = _upper(total/arb(255**(8*arity)))
        if not weight > 0:
            raise ArithmeticError('A strictly positive B_a upper bound is required')
        weights[arity] = weight
        diagnostics.append(dict(arity=arity, labelled_rowsets=comb(16, arity),
            distinct_scheduled_signatures=len(signatures), moment_order=order,
            input_degree_floor=9*arity*order, reversed_deficit_maximum=7*arity*order,
            rs_nonzero_exponent=8*arity))
    return weights, diagnostics


def class_prefix(weight, order, groups=R):
    """Raw weights already include rowset counts and RS density factors."""
    if (type(order) is not int or order not in ORDERS or type(groups) is not int or groups < 1
            or not weight.is_finite() or not weight > 0):
        raise ValueError('Positive group weight and prescribed moment order required')
    theta = ALPHA*F(order, order-1)
    return (A.uniform.rational(theta)*(groups*A.log2(weight)-A.log2(arb(255)))
            - A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK)))


def _profile_rows(terms):
    expected = A.profile_counts(R)
    seen, rows = set(), []
    for even, odd, count, term in terms:
        if (type(even) is not int or type(odd) is not int or (even, odd) not in expected
                or (even, odd) in seen or type(count) is not int or count != expected[even, odd]
                or not term.is_finite() or not term > 0):
            raise ValueError('Wrong indexed uniform source profile')
        seen.add((even, odd))
        rows.append(dict(even_full=even, odd_full=odd, indexed_source_count_upper=count,
                         uniform_contribution_log2_upper=str(A.uniform.rational_ceiling(A.log2(term)))))
    if seen != set(expected) or len(rows) != 99:
        raise ValueError('All 99 indexed source profiles are required')
    return rows


def receipt(inventory_path=DIRECTORY/INVENTORY_NAME, precision=256, *,
            weights_only=False, max_seconds=600, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Use a fresh256/384-bit replay')
    if type(weights_only) is not bool:
        raise ValueError('weights_only must be Boolean')
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 1800):
        raise ValueError('Use a positive bounded replay budget, at most1800 seconds')
    pins = dependencies()
    checked = load_inventory(inventory_path)
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('No complete marked control produced before the budget')
    controls, weight_records = [], []
    with ctx.workprec(precision):
        profile_rows, uniform_log = [], None
        if not weights_only:
            if progress:
                progress(dict(stage='fresh-uniform-regional-bank', precision_bits=precision))
            regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
            check_time()
            terms = A.source_profile_terms(R, regions, check_time=check_time)
            profile_rows = _profile_rows(terms)
            uniform_log = A.log2(sum((term for _, _, _, term in terms), arb(0)))
        for order in ORDERS:
            weights, diagnostics = scheduled_weights(checked, order, check_time=check_time, progress=progress)
            weight_records.append(dict(moment_order=order,
                weights=[dict(arity=a, B_upper=str(exact_upper(weights[a])),
                              B_log2_upper=str(A.uniform.rational_ceiling(A.log2(weights[a])))) for a in ARITIES],
                diagnostics=diagnostics))
            if not weights_only:
                homogeneous = []
                for arity in ARITIES:
                    prefix = class_prefix(weights[arity], order)
                    value = prefix+uniform_log
                    homogeneous.append(dict(every_group_arity=arity,
                        prefix_log2_upper=str(A.uniform.rational_ceiling(prefix)),
                        failure_log2_upper=str(A.uniform.rational_ceiling(value))))
                mixed_prefix = class_prefix(_upper(sum(weights.values(), arb(0))), order)
                controls.append(dict(moment_order=order, outer_theta=str(ALPHA*F(order, order-1)),
                    homogeneous=homogeneous,
                    mixed=dict(every_group_arity_range=[1, 5],
                        prefix_log2_upper=str(A.uniform.rational_ceiling(mixed_prefix)),
                        failure_log2_upper=str(A.uniform.rational_ceiling(mixed_prefix+uniform_log)))))
        selected = None
        if controls:
            choice = min(controls, key=lambda row: F(row['mixed']['failure_log2_upper']))
            selected = dict(moment_order=choice['moment_order'],
                            failure_log2_upper=choice['mixed']['failure_log2_upper'])
        uniform_ceiling = None if weights_only else str(A.uniform.rational_ceiling(uniform_log))
        check_time()
    if dependencies() != pins:
        raise ValueError('Reviewed source pins changed during replay')
    # Recheck the observed exact artifact as well as its helper sources.
    if sha256(Path(inventory_path).read_bytes()).hexdigest() != INVENTORY_SHA256:
        raise ValueError('Exact marked inventory changed during replay')
    return dict(construction=binding.SCOPE['construction'], field_order=256,
        lambda_aes_byte=188, active_groups=R, every_group_arity_range=[1, 5],
        original_rowsets_fixed_across_columns=True, all6884_labelled_rowsets=True,
        local_equivalence_only=True, sorted_signatures_are_scalar_product_cache_only=True,
        original_RS_row_floor=9, aggregate_input_floor_per_group='9*a*p',
        per_row_RS_floors_relaxed_to_aggregate=True, raw_marked_likelihood=True,
        no_extra_density_beta=True, rs_factor_outside_each_root='255^(-8*a)',
        positive_deficit_tail=True, coefficient_method='direct-positive-truncated-convolution-upper-endpoints',
        no_midpoints=True, moment_orders=list(ORDERS), weight_bounds=weight_records,
        weights_only=weights_only, mark=str(MARK), uniform_alpha=str(ALPHA),
        projective_normalizations=1, projective_exponent='theta=p/(5*(p-1))',
        rowset_union_inside_theta=True, indexed_source_union_outside_theta=True,
        source_profiles=profile_rows, profile_count=len(profile_rows),
        uniform_profile_union_log2_upper=uniform_ceiling,
        controls=controls, selected_mixed=selected, homogeneous_classes_not_added_to_mixed=True,
        precision_bits=precision, fresh_regional_banks=int(not weights_only), directed=True,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        other_arities_covered=False, other_active_group_weights_covered=False,
        whole_code_certificate=False, inventory_sha256=INVENTORY_SHA256,
        reviewed_source_sha256=pins, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, default=DIRECTORY/INVENTORY_NAME)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--weights-only', action='store_true')
    parser.add_argument('--max-seconds', type=float, default=600)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.inventory, args.precision, weights_only=args.weights_only,
        max_seconds=args.max_seconds, progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                     if k not in ('source_profiles', 'weight_bounds', 'reviewed_source_sha256')}))
