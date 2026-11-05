"""Bounded cached-bank Q_eta comparison for mixed arities11..16 at r18.

For each prescribed eta, select the smallest COMPLETED reference upper
bound among nine flag marks. Use that same eta in every column density
moment. No density/reference cross-selection is permitted. Each fixed
original rowset retains its sixteen scheduled moments before the >=9
nonzero-column tail and pth root. Sum rowsets, then take the eighteenth
power of the arity polynomial. For total arity A the bound is

 theta*(A*log2(256^16/255^8)+log2(coefficient)-log2(255))
       -alpha*734003*log2(691/700)+reference_source_log2,

where theta=alpha*p/(p-1), alpha=1/5 and p is3 or4. The density weights
are norm-only, so the beta factor belongs here exactly once. Source
counting already belongs to the reference factor, outside theta. The
projective divisor appears once inside theta. The reference flag floor
and density nonzero-column floor are the two compatible restrictions
proved in hd_lambda_renyi_zero_reference.py.

All controls cover the SAME fixed total-arity class, so choose the least
upper bound per total only after completing each control. Then sum those
disjoint total classes. The optional eta0 control may additionally use
the saved tighter alternating-profile uniform reference. That old
reference is never paired with a nonzero eta density.

Only pinned bank endpoints and exact moment inventories are read.
A384-bit replay requires the separately generated384-bit bank; raising
the arithmetic precision of the256-bit bank is explicitly rejected.
This is a selected mixed-class bound, not a whole-code certificate.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter

from flint import arb, arb_poly, ctx

import hd_lambda_renyi_interpolate as interpolate
import hd_lambda_renyi_zero_reference as reference

DIRECTORY = Path(__file__).resolve().parent
A, binding = reference.A, reference.binding
R, ALPHA, MARK = 18, F(1, 5), F(691, 700)
ARITIES, ORDERS = tuple(range(11, 17)), (3, 4)
ETA_BITS = (96, 104, 112, 120, 128, 136)
ETAS = tuple(F(1, 2**bits) for bits in ETA_BITS)
FLAG_MARKS = tuple(F(2**k) for k in range(9))
MIN_TOTAL, MAX_TOTAL = 198, 245
PINS = {
    'hd_lambda_renyi_zero_reference.py': '1ad60e56a5953ac7cd2bcaf21b49c20d41bc2eab872a79ddab63a207b89157e8',
    'hd_lambda_renyi_interpolate.py': '9ab6575f4f33a3c9d4de78227ab7f2fe1e0128b7d0b6ccb8534f6dbb5fbac4dd',
}
BANK_PINS = {
    256: 'ae289c45b1e51159425e5645001c60c81b71ed3a1cac4a833ccc65baea033d78',
    384: '96af0a052dcd6d2c16e28741e23d3e87b96451534de2535a9b4f81e60a4ec155',
}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen comparison source changed: '+name)
    return {**interpolate.dependencies(), **reference.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def validate_bank(data, precision):
    """Validate scope in addition to the external file hash checked by load_bank."""
    expected = dict(construction=binding.SCOPE['construction'], reference='eta*delta0+(1-eta)*U16',
        active_groups=18, field_order=256, mark=str(MARK), alpha=str(ALPHA),
        regions_per_column=16, columns=16, total_regions=256, whole_column_flags=True,
        active_subset_varies_per_column=True, original_source_fixed=True,
        source_count_outside_alpha=True, source_count_charged_once=True,
        flag_weights_powered_by_alpha=True, flag_mark_outside_alpha=True,
        zero_input_regional_evolution_retained=True, matching_f_maximized_per_entry=True,
        alternating_profile_gain_discarded=True, original_rowsets_not_resampled=True,
        flag_floor_requires_RS_support_event=True, actual_nonzero_column_floor_per_group=9,
        fresh_regional_banks=1, precision_bits=precision, packet_mantissa_bits=60,
        reference_factor_only=True, density_weights_included=False, threshold_included=False,
        projective_factor_included=False, failure_bound=False, whole_code_certificate=False,
        arithmetic_enclosures_directed=True, source_profiles=99,
        regional_max_f_entries_are_exact_rational_upper_bounds=True,
        reviewed_source_sha256=reference.dependencies())
    for key, value in expected.items():
        if type(data.get(key)) is not type(value) or data[key] != value:
            raise ValueError('Incorrect saved reference scope: '+key)
    counts = A.profile_counts(R)
    source_count = sum(counts.values())
    if len(counts) != 99 or data.get('source_count') != str(source_count):
        raise ValueError('Original source count does not match all99 indexed profiles')
    rows = data.get('regional_max_f_upper')
    if (not isinstance(rows, list) or len(rows) != 19
            or any(not isinstance(row, list) or len(row) != 4 for row in rows)
            or any(type(x) is not str for row in rows for x in row)):
        raise ValueError('All19 exact rational max-f matrices are required')
    envelopes = tuple(tuple(F(x) for x in row) for row in rows)
    if any(x < 0 for row in envelopes for x in row):
        raise ValueError('Negative saved regional entry')
    # Initial empty input stays zero. The other diagonal MUST be retained.
    if envelopes[0][0] != 1 or envelopes[0][1] or envelopes[0][2] or not envelopes[0][3] > 0:
        raise ValueError('Incorrect zero-input regional boundary law')
    if type(data.get('uniform_comparison_log2_upper')) is not str:
        raise ValueError('Missing directed uniform comparison upper bound')
    uniform_log = F(data['uniform_comparison_log2_upper'])
    return envelopes, source_count, uniform_log


def load_bank(path, precision):
    if type(precision) is not int or precision not in BANK_PINS:
        raise ValueError('Use the separately generated256/384-bit reference bank')
    path = Path(path)
    raw = path.read_bytes()
    if sha256(raw).hexdigest() != BANK_PINS[precision]:
        raise ValueError('Unreviewed bank or precision-mismatched bank')
    data = json.loads(raw)
    return data, validate_bank(data, precision)


def choose_reference(prepared, eta, source_count, uniform_log):
    """Complete all nine prescribed flags; eta0 alone may use the old U bound."""
    eta = F(eta)
    if eta not in (F(0),)+ETAS:
        raise ValueError('Only the prescribed eta grid is supported')
    rows = []
    for v in FLAG_MARKS:
        value = reference.reference_log(prepared, eta, source_count=source_count,
                                         flag_mark=v, minimum_flags=9*R)
        rows.append(dict(flag_mark=str(v), kind='max-f-column-flags',
                         reference_source_log2_upper=str(A.uniform.rational_ceiling(value))))
    if not eta:
        rows.append(dict(flag_mark='1', kind='saved-uniform-alternating-profiles',
                         reference_source_log2_upper=str(F(uniform_log))))
    selected = min(rows, key=lambda row: F(row['reference_source_log2_upper']))
    return dict(eta=str(eta), all_reference_controls=rows, selected_reference=dict(selected))


def weights_from_signatures(prepared, order, eta):
    """Reuse paired actual-column signatures; transform atoms using THIS eta."""
    eta = F(eta)
    if type(order) is not int or order not in ORDERS or eta not in (F(0),)+ETAS:
        raise ValueError('Prescribed eta and order3/4 required')
    if tuple(prepared) != ARITIES:
        raise ValueError('Complete ordered arities11..16 required')
    weights, diagnostics = [], []
    reference_zero = eta+(1-eta)*F(256)**-16
    scale = F(1-eta)**(1-order)
    for arity, signatures in prepared.items():
        zero = F(256)**(-arity*order)/reference_zero**(order-1)
        atoms = {pair: pair[order-3]*scale for signature in signatures for pair in signature}
        if any(atom < 0 for atom in atoms.values()):
            raise ValueError('Negative nonzero density atom')
        value = arb(0)
        for signature, count in signatures.items():
            tail = reference.activity.activity_tail(zero, tuple(atoms[pair] for pair in signature), 9)
            value += count*tail**A.uniform.rational(F(1, order))
        if not value > 0 or not value.is_finite():
            raise ArithmeticError('Positive finite density norm weight required')
        weights.append(value)
        diagnostics.append(dict(arity=arity, labelled_rowsets=sum(signatures.values()),
            distinct_joint_moment_signatures=len(signatures), zero_moment=str(zero),
            nonzero_moment_scale=str(scale), weight_log2_enclosure=str(A.log2(value))))
    return weights, diagnostics


def total_bound(coefficient, total, order, reference_log):
    """Norm-only coefficient: beta^total and projective divisor each once."""
    if (type(total) is not int or not MIN_TOTAL <= total <= MAX_TOTAL
            or type(order) is not int or order not in ORDERS or not coefficient > 0):
        raise ValueError('Valid positive mixed coefficient, total and order required')
    theta = ALPHA*F(order, order-1)
    log_beta = 128-8*A.log2(arb(255))
    threshold = -A.uniform.rational(ALPHA)*734003*A.log2(A.uniform.rational(MARK))
    return (A.uniform.rational(theta)*(total*log_beta+A.log2(coefficient)-A.log2(arb(255)))
            + threshold + A.uniform.rational(F(reference_log)))


def combined_controls(prepared_moments, prepared_reference, source_count, uniform_log,
                      *, include_uniform=True, progress=None):
    if type(include_uniform) is not bool:
        raise ValueError('Uniform baseline selection must be Boolean')
    etas = ((F(0),) if include_uniform else ()) + ETAS
    references, controls = [], []
    for eta in etas:
        selected = choose_reference(prepared_reference, eta, source_count, uniform_log)
        references.append(selected)
        ref_log = F(selected['selected_reference']['reference_source_log2_upper'])
        for order in ORDERS:
            weights, diagnostics = weights_from_signatures(prepared_moments, order, eta)
            polynomial = arb_poly(weights)**R
            totals = []
            for total in range(MIN_TOTAL, MAX_TOTAL+1):
                bound = total_bound(polynomial[total-MIN_TOTAL], total, order, ref_log)
                totals.append(dict(total_arity=total,
                                   failure_log2_upper=str(A.uniform.rational_ceiling(bound))))
            control = dict(eta=str(eta), moment_order=order,
                outer_theta=str(ALPHA*F(order, order-1)), reference_source_log2_upper=str(ref_log),
                selected_reference=dict(selected['selected_reference']),
                density_reference_eta_match=True, norm_only_weights=True,
                weight_diagnostics=diagnostics, totals=totals)
            controls.append(control)
            if progress:
                progress(dict(completed_eta=str(eta), completed_order=order,
                              reference_log2_upper=str(ref_log)))
    selected_totals, cumulative = [], arb(0)
    for index, total in enumerate(range(MIN_TOTAL, MAX_TOTAL+1)):
        choice = min(controls, key=lambda row: F(row['totals'][index]['failure_log2_upper']))
        bound = F(choice['totals'][index]['failure_log2_upper'])
        cumulative += arb(2)**A.uniform.rational(bound)
        selected_totals.append(dict(total_arity=total, selected_eta=choice['eta'],
            selected_moment_order=choice['moment_order'],
            selected_reference=dict(choice['selected_reference']), failure_log2_upper=str(bound),
            cumulative_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative)))))
    return dict(reference_controls=references, controls=controls, selected_totals=selected_totals,
                union_log2_upper=str(A.uniform.rational_ceiling(A.log2(cumulative))))


def receipt(bank_path, *, precision=256, small_path=interpolate.DEFAULT_SMALL,
            middle_path=interpolate.DEFAULT_MIDDLE, include_uniform=True, progress=None):
    pins, started = dependencies(), perf_counter()
    _, (envelopes, source_count, uniform_log) = load_bank(bank_path, precision)
    moments, inventory_paths = interpolate.checked_inputs(small_path, middle_path)
    moment_signatures = interpolate.prepare_signatures(moments)
    with ctx.workprec(precision):
        bank = reference.reference_from_envelope(envelopes)
        result = combined_controls(moment_signatures, bank, source_count, uniform_log,
                                   include_uniform=include_uniform, progress=progress)
    if (dependencies() != pins or not interpolate.inputs_unchanged(inventory_paths)
            or sha256(Path(bank_path).read_bytes()).hexdigest() != BANK_PINS[precision]):
        raise ValueError('Source, exact inventory, or cached bank changed during control')
    return dict(result, construction=binding.SCOPE['construction'], field_order=256,
        active_groups=18, every_group_arity_range=[11, 16], total_arity_range=[198, 245],
        moment_orders=list(ORDERS), eta_powers=list(ETA_BITS), uniform_baseline=include_uniform,
        flag_marks=list(map(str, FLAG_MARKS)), minimum_active_flags=162,
        minimum_nonzero_columns_per_group=9, reference='eta*delta0+(1-eta)*U16',
        original_rowsets_fixed_across_columns=True, paired_actual_column_moments=True,
        all6885_labelled_rowsets=True, source_profiles=99, source_count=str(source_count),
        source_count_charged_once=True, source_count_outside_theta=True,
        rowset_sum_inside_theta=True, norm_only_weights=True, beta_charged_once=True,
        projective_normalizations=1, mark=str(MARK), uniform_alpha=str(ALPHA),
        source_dimension=1048576, output_length=2097152, forbidden_weight=734003,
        restricted_density_not_renormalized=True, density_reference_eta_match=True,
        precision_bits=precision, cached_bank_precision_bits=precision,
        cached_bank_sha256=BANK_PINS[precision], cached_bank_scope_validated=True,
        inventory_sha256=dict(interpolate.INVENTORY_PINS), reviewed_source_sha256=pins,
        source_pins_unchanged=True, inventory_pins_unchanged=True, bank_pin_unchanged=True,
        fresh_regional_banks=0, new_moment_inventories=0, directed=True,
        whole_code_certificate=False, other_group_counts_covered=False,
        groups_below_arity11_covered=False, elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), default=256)
    parser.add_argument('--bank', type=Path)
    parser.add_argument('--small-inventory', type=Path, default=interpolate.DEFAULT_SMALL)
    parser.add_argument('--middle-inventory', type=Path, default=interpolate.DEFAULT_MIDDLE)
    parser.add_argument('--no-uniform-baseline', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    bank_path = args.bank or DIRECTORY/f'hd_lambda_renyi_zero_reference_final_{args.precision}.json'
    result = receipt(bank_path, precision=args.precision, small_path=args.small_inventory,
        middle_path=args.middle_inventory, include_uniform=not args.no_uniform_baseline,
        progress=lambda row: print(json.dumps(row), flush=True))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ('controls', 'reference_controls', 'reviewed_source_sha256')}))
