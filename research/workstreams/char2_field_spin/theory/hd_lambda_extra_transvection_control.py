"""Variant-only activity comparison and exact two-source support checks.

The new outer is H'=Hlambda+v*v^T with V=(0,1,2,3), not the measured
Hlambda construction. Its complete complement0..5 moment inventory must
use the actual four-map stabilizer and retain all6885 labelled row sets.
Scheduled moment products keep each original row set fixed across columns.

Only the matrix-independent activity-tail and uniform regional arithmetic
are reused. For common p, the >=9-nonzero-column restricted likelihood norm
gives w_a. At total arity k, the counting coefficient is

    (256^16/255^8)^k [u^(k-198)] (sum_(a=11..16) w_a u^(a-11))^18.

Indicator Holder uses theta=p/(5*(p-1)), uniform alpha1/5, and mark691/700.
The one projective255 factor is inside theta; indexed-source counts are
outside. Totals198..245 are disjoint. Moment orders are alternatives, not
events to union. No old certificate, low-arity table, or timing transfers.

The saved old activity receipt is a comparison only. Its p2 bound used
I3 upper bounds for arities12..16; the new inventory supplies exact I2.
Thus the comparison is not an isolated measurement of the matrix change.
"""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb, isfinite
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

import hd_lambda_extra_transvection as variant
import hd_lambda_renyi_activity as activity

DIRECTORY = Path(__file__).resolve().parent
A = activity.A
R, ARITIES, ORDERS = 18, (11, 12, 13, 14, 15, 16), (2, 3, 4)
TOTALS = tuple(range(198, 246))
ALPHA, MARK = F(1, 5), F(691, 700)
BASELINE_NAME = 'hd_lambda_renyi_activity_a11_256.json'
BASELINE_SHA256 = '77edcd1ea1f70b3557e249198912a267721b37dca82a3bba4d8abac1fe841464'
INVENTORY_NAME = 'hd_lambda_extra_transvection_inventory234.json'
INVENTORY_SHA256 = '2089c87acfb442319465a8e888cb06d5091ca5ebedde26d067dadebe5b8390e5'
PINS = {
    'hd_lambda_extra_transvection.py': 'f5a72f885a58a8ea1455cced05874e5fa613271469d867686bdb09df88ff2f7c',
    'hd_lambda_renyi_activity.py': '09f5245aecb903097438b65ded8757c6a66da40a13804289c394fef777186c77',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Reviewed dependency changed: '+name)
    return {**variant.checked_dependencies(), **activity.dependencies(), **PINS,
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
            raise ValueError('Wrong variant scope or witness: '+key)


def checked_inventory(data):
    """Require the COMPLETE new four-map inventory, not an old48-map table."""
    plan = tuple((b, row) for b in range(6) for row in variant.complement_orbits(b))
    _require(data, dict(construction=variant.CONSTRUCTION, transvection_support=[0, 1, 2, 3],
        field_order=256, field_polynomial='0x11b', lambda_aes_byte=188,
        complement_sizes=list(range(6)), arities=list(ARITIES), moment_orders=list(ORDERS),
        exact_rational_arithmetic=True, complete=True, stop_reason=None,
        full_parity_includes_identity=True, code='ker((Hlambda+v*v^T)[complement,:])',
        output_diagonal='independent uniform nonzero field multipliers',
        normalized_likelihood='256^|complement|*n_C(S)/255^|S|',
        subgroup_size=4, old_subgroup_size=48,
        symmetry='stabilizer of transvection support inside proved orthogonal subgroup',
        no_larger_symmetry_claim=True, all_maps_check_full256_entries=True,
        maps=[list(image) for image in variant.stabilizer_maps()],
        expected_labelled_complements=6885, expected_representatives=len(plan),
        actual_support_tables=len(plan), original_rowsets_not_resampled=True,
        reviewed_source_sha256=variant.checked_dependencies(), source_pins_unchanged=True,
        unchanged_disjoint_complements_by_size={str(b): comb(12, b) for b in range(6)},
        old_receipts_transferred=False, regional_provider_runs=0, native_benchmarks=0,
        native_performance_claim=False, whole_code_certificate=False))
    representatives = data.get('representatives')
    if not isinstance(representatives, list) or len(representatives) != len(plan):
        raise ValueError('Incomplete variant representative table')
    representative_values, maxima = {}, {16-b: {p: F(1) for p in ORDERS} for b in range(6)}
    for index, ((b, orbit), row) in enumerate(zip(plan, representatives)):
        _require(row, dict(orbit_index=index, complement=list(orbit['representative']),
                          arity=16-b, multiplicity=len(orbit['members'])))
        raw = row.get('moments')
        if (not isinstance(raw, dict) or set(raw) != {'2', '3', '4'}
                or any(type(v) is not str for v in raw.values())):
            raise ValueError('Exact moment strings for orders2,3,4 required')
        try:
            values = {p: F(raw[str(p)]) for p in ORDERS}
        except (ValueError, ZeroDivisionError) as error:
            raise ValueError('Malformed moment') from error
        if any(values[p] < max(F(1), activity.zero_atom(16-b, p)) for p in ORDERS):
            raise ValueError('Moment below probability or exact zero-atom mass')
        representative_values[index] = values
        for p in ORDERS:
            maxima[16-b][p] = max(maxima[16-b][p], values[p])
    _require(data, dict(maximum_moments_by_arity={str(a): {str(p): str(v) for p, v in values.items()}
                                               for a, values in maxima.items()}))
    expected = {member: (b, index, row['representative'], image)
                for index, (b, row) in enumerate(plan) for member, image in row['members']}
    rows = data.get('rows')
    if not isinstance(rows, list) or len(rows) != 6885:
        raise ValueError('All6885 labelled complement rows are required')
    moments = {}
    for row in rows:
        complement = row.get('complement') if isinstance(row, dict) else None
        if not isinstance(complement, list) or any(type(i) is not int for i in complement):
            raise ValueError('Integer mapped complement labels required')
        key = tuple(complement)
        if key not in expected or key in moments:
            raise ValueError('Unexpected or repeated complement')
        b, index, representative, image = expected[key]
        _require(row, dict(arity=16-b, orbit_index=index, representative=list(representative),
            representative_to_complement_permutation=list(image),
            moments={str(p): str(v) for p, v in representative_values[index].items()}))
        moments[key] = dict(representative_values[index])
    if set(moments) != set(expected):
        raise ValueError('Incomplete labelled complement coverage')
    return moments


def _sha(value):
    if (type(value) is not str or len(value) != 64
            or any(c not in '0123456789abcdef' for c in value)):
        raise ValueError('An explicit lowercase SHA256 trust anchor is required')
    return value


def load_inventory(path=DIRECTORY/INVENTORY_NAME, expected_sha256=INVENTORY_SHA256):
    raw = Path(path).read_bytes()
    if _sha(expected_sha256) != INVENTORY_SHA256 or sha256(raw).hexdigest() != expected_sha256:
        raise ValueError('Variant inventory hash differs from the observed artifact')
    return checked_inventory(json.loads(raw))


def checked_baseline(data):
    """Validate the pinned old receipt only as a comparison, never as new evidence."""
    _require(data, dict(construction=activity.binding.SCOPE['construction'], active_groups=18,
        every_group_arity_range=[11, 16], total_arity_range=[198, 245],
        minimum_active_columns_per_group=9, restricted_density_not_renormalized=True,
        exact_zero_atom_subtraction=True, order2_large_arity_uses_I3_upper=True,
        original_rowsets_fixed_across_columns=True, all6885_labelled_rowsets=True,
        source_profiles=99, source_count_outside_theta=True, rowset_sum_inside_theta=True,
        projective_normalizations=1, field_order=256, forbidden_weight=734003,
        mark=str(MARK), uniform_alpha=str(ALPHA), precision_bits=256,
        fresh_regional_banks=1, directed=True, whole_code_certificate=False,
        reviewed_source_sha256=activity.dependencies()))
    controls = data.get('controls')
    if not isinstance(controls, list) or len(controls) != 3:
        raise ValueError('Three old comparison controls required')
    values = {}
    for p, control in zip(ORDERS, controls):
        _require(control, dict(moment_order=p, outer_theta=str(ALPHA*F(p, p-1))))
        rows = control.get('totals')
        if not isinstance(rows, list) or len(rows) != len(TOTALS):
            raise ValueError('Incomplete old total-arity comparison')
        for total, row in zip(TOTALS, rows):
            _require(row, dict(total_arity=total))
            if type(row.get('failure_log2_upper')) is not str:
                raise ValueError('Old rational failure ceiling required')
            values[p, total] = F(row['failure_log2_upper'])
    selected = data.get('selected_totals')
    if not isinstance(selected, list) or len(selected) != len(TOTALS):
        raise ValueError('Incomplete old selected-total table')
    result = {}
    for total, row in zip(TOTALS, selected):
        p = min(ORDERS, key=lambda p: values[p, total])
        _require(row, dict(total_arity=total, selected_moment_order=p,
                           failure_log2_upper=str(values[p, total])))
        result[total] = (p, values[p, total])
    return result


def load_baseline(path=DIRECTORY/BASELINE_NAME):
    raw = Path(path).read_bytes()
    if sha256(raw).hexdigest() != BASELINE_SHA256:
        raise ValueError('Pinned old comparison receipt changed')
    return checked_baseline(json.loads(raw))


def pair_branches(source, *, use_variant=True, lambda_value=188, q=256):
    """Seven unmerged raw branches for an ordered pair of input coordinates."""
    source = tuple(source)
    if (len(source) != 2 or len(set(source)) != 2
            or any(type(i) is not int or not 0 <= i < 16 for i in source)
            or type(use_variant) is not bool):
        raise ValueError('Two distinct ordered labels and Boolean variant selector required')
    h = variant.matrix(lambda_value, q) if use_variant else variant.original.matrix(lambda_value, q)
    ratios = variant.original.exceptional_ratios(lambda_value, q)
    first, second = source
    masks = tuple(sum(bool(variant.original._multiply(row[first], ratio, q) ^ row[second]) << i
                      for i, row in enumerate(h)) for ratio in ratios)
    full, d = 65535, q-1
    branches = ((0, 0, 1), (1, full, d), (1, full, d),
                *((2, mask, d) for mask in masks), (2, full, d*(q-4)))
    if sum(count for _, _, count in branches) != q*q:
        raise ArithmeticError('Raw pair mass changed')
    return branches


def two_source_profile():
    """Exact120-pair profile; this does not transfer any old small-arity proof."""
    rows, minimum, collapsed, histograms = [], {'old': 16, 'new': 16}, {'old': [], 'new': []}, {'old': {}, 'new': {}}
    exceptional_histograms, worst = {'old': {}, 'new': {}}, {}
    for source in combinations(range(16), 2):
        data = dict(source=list(source))
        for key, use_variant in (('old', False), ('new', True)):
            branches = pair_branches(source, use_variant=use_variant)
            positive_weights = [mask.bit_count() for weight, mask, count in branches if weight and count]
            least = min(positive_weights)
            data[key+'_branches'] = [list(branch) for branch in branches]
            data[key+'_minimum_output_weight'] = least
            minimum[key] = min(minimum[key], least)
            histograms[key][str(least)] = histograms[key].get(str(least), 0)+1
            for index, (_, mask, count) in enumerate(branches[3:6]):
                weight = mask.bit_count()
                exceptional_histograms[key][str(weight)] = exceptional_histograms[key].get(str(weight), 0)+1
                if key not in worst or weight < worst[key]['output_weight']:
                    worst[key] = dict(source=list(source), exceptional_ratio_index=index,
                                      output_weight=weight, mask=mask, raw_count=count)
            if least == 0:
                collapsed[key].append(list(source))
        rows.append(data)
    return dict(new_construction=variant.CONSTRUCTION, old_construction=activity.binding.SCOPE['construction'],
        q=256, lambda_aes_byte=188, pair_count=120, branch_order=[
            'zero', 'left-single', 'right-single', 'ratio1', 'ratio-lambda-over-lambda-plus1',
            'ratio-lambda-plus1-over-lambda', 'remaining252-ratios'],
        singleton_output_weight=16, raw_pair_mass=256**2, minimum_output_weight=minimum,
        minimum_weight_histograms=histograms, collapsed_pairs=collapsed, pairs=rows,
        exceptional_ratio_weight_histograms=exceptional_histograms,
        first_worst_exceptional_branch=worst,
        exact_integer_counts=True, sampled_amplitudes=False, source_order_retained=True,
        old_lowarity_certificate_transferred=False, whole_code_certificate=False)


def receipt(inventory_path=DIRECTORY/INVENTORY_NAME, inventory_sha256=INVENTORY_SHA256,
            precision=256, *, max_seconds=180, progress=None):
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Use a fresh256/384-bit variant replay')
    if (type(max_seconds) not in (int, float) or not isfinite(max_seconds)
            or not 0 < max_seconds <= 600):
        raise ValueError('Use a bounded positive budget at most600 seconds')
    pins = dependencies()
    moments = load_inventory(inventory_path, inventory_sha256)
    baseline = load_baseline()
    started = perf_counter()
    def check_time():
        if perf_counter()-started >= max_seconds:
            raise TimeoutError('No complete variant control before time budget')
    with ctx.workprec(precision):
        if progress:
            progress(dict(stage='fresh-outer-independent-uniform-bank', precision_bits=precision))
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        terms = A.source_profile_terms(R, regions, check_time=check_time)
        expected_counts = A.profile_counts(R)
        if len(terms) != 99 or {(e, o): n for e, o, n, _ in terms} != expected_counts:
            raise ArithmeticError('All99 indexed source profiles required')
        profile_rows = [dict(even_full=e, odd_full=o, indexed_source_count_upper=n,
            uniform_contribution_log2_upper=str(A.uniform.rational_ceiling(A.log2(term))))
            for e, o, n, term in terms]
        uniform_log = A.log2(sum((term for _, _, _, term in terms), arb(0)))
        threshold = -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))
        log_beta = 128-8*A.log2(arb(255))
        controls = []
        for p in ORDERS:
            check_time()
            if progress:
                progress(dict(stage='actual-scheduled-activity-products', moment_order=p))
            weights, diagnostics = activity.scheduled_weights(moments, p)
            coefficients = arb_poly(weights)**R
            rows = []
            theta = ALPHA*F(p, p-1)
            for total in TOTALS:
                pref = A.uniform.rational(theta)*(total*log_beta
                    + A.log2(coefficients[total-198])-A.log2(arb(255)))+threshold
                rows.append(dict(total_arity=total,
                    failure_log2_upper=str(A.uniform.rational_ceiling(pref+uniform_log))))
            controls.append(dict(moment_order=p, outer_theta=str(theta),
                arity_weight_log2_enclosures=[str(A.log2(x)) for x in weights],
                rowset_diagnostics=diagnostics, totals=rows))
        selected, cumulative, old_cumulative = [], arb(0), arb(0)
        for index, total in enumerate(TOTALS):
            choice = min(controls, key=lambda row: F(row['totals'][index]['failure_log2_upper']))
            bound = F(choice['totals'][index]['failure_log2_upper'])
            old_order, old_bound = baseline[total]
            cumulative += arb(2)**A.uniform.rational(bound)
            old_cumulative += arb(2)**A.uniform.rational(old_bound)
            selected.append(dict(total_arity=total, selected_moment_order=choice['moment_order'],
                failure_log2_upper=str(bound), old_selected_moment_order=old_order,
                old_failure_log2_upper=str(old_bound), old_minus_new_bound_bits=str(old_bound-bound),
                cumulative_new_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative))),
                cumulative_old_log2_upper=str(A.uniform.rational_ceiling(A.log2(old_cumulative)))))
        union = str(A.uniform.rational_ceiling(A.log2(cumulative)))
        uniform_ceiling = str(A.uniform.rational_ceiling(uniform_log))
    check_time()
    if dependencies() != pins:
        raise ValueError('Sources changed during variant replay')
    if sha256(Path(inventory_path).read_bytes()).hexdigest() != inventory_sha256:
        raise ValueError('Variant inventory changed during replay')
    if sha256((DIRECTORY/BASELINE_NAME).read_bytes()).hexdigest() != BASELINE_SHA256:
        raise ValueError('Old comparison receipt changed during replay')
    return dict(construction=variant.CONSTRUCTION, transvection_support=[0, 1, 2, 3],
        input_schedule=activity.binding.SCOPE['input_schedule'], output_schedule=activity.binding.SCOPE['output_schedule'],
        column_map='D2_j Q_j (Hlambda+v*v^T) P_j D1_j', field_order=256, lambda_aes_byte=188,
        active_groups=R, every_group_arity_range=[11, 16], total_arity_range=[198, 245],
        original_rowsets_fixed_across_columns=True, all6885_labelled_rowsets=True,
        actual_stabilizer_size=4, old48_map_quotient_used=False, minimum_active_columns_per_group=9,
        restricted_density_not_renormalized=True, exact_zero_atom_subtraction=True,
        all_orders_use_exact_local_moments=True, source_profiles=profile_rows, profile_count=99,
        source_count_outside_theta=True, rowset_sum_inside_theta=True, projective_normalizations=1,
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        mark=str(MARK), uniform_alpha=str(ALPHA), uniform_profile_union_log2_upper=uniform_ceiling,
        controls=controls, selected_totals=selected, rounded_total_union_log2_upper=union,
        precision_bits=precision, fresh_regional_banks=1, directed=True,
        inventory_sha256=inventory_sha256, reviewed_source_sha256=pins,
        old_comparison=dict(file=BASELINE_NAME, sha256=BASELINE_SHA256,
            construction=activity.binding.SCOPE['construction'], exact_same_classes=True,
            old_p2_large_arities_used_I3_upper=True, isolated_matrix_effect_comparison=False),
        exact_two_source_profile=two_source_profile(), old_receipts_transferred=False,
        old_lowarity_certificates_transferred=False, native_performance_claim=False,
        whole_code_certificate=False, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', type=Path, default=DIRECTORY/INVENTORY_NAME)
    parser.add_argument('--inventory-sha256', default=INVENTORY_SHA256)
    parser.add_argument('--precision', type=int, choices=(256, 384), required=True)
    parser.add_argument('--max-seconds', type=float, default=180)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.inventory, args.inventory_sha256, args.precision, max_seconds=args.max_seconds,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
        if k not in ('controls', 'selected_totals', 'source_profiles', 'reviewed_source_sha256', 'exact_two_source_profile')}))
