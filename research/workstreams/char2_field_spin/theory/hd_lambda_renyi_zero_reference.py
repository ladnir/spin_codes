"""Whole-column zero-mixture reference for the scheduled Hlambda comparison.

Fix the original r-group source G. In each column and group, draw a flag:
with probability eta output zero; otherwise draw a uniform F256^16 vector.
Thus Q_eta=eta*delta_0+(1-eta)*U16. One flag controls all sixteen regions
of its column. It is not sixteen independent site flags.

After expansion over column flags, alpha-subadditivity powers their
probabilities. Conditional on the active subset H, the frozen uniform
regional provider applies with m=|H|. Let Rbar_m be its entrywise maximum
over all matching counts f. The positive column matrix is

 K(v)=sum_m C(r,m)*eta^(alpha*(r-m))*(1-eta)^(alpha*m)
              *v^m*Rbar_m^16.

The reference word bound is v^(-minimum_flags)*e0*K(v)^16*ones.
Only the original G is counted externally. Its indexed-source upper
count is the sum of the frozen profile counts, with no extra binomial
source factor. The max-f relaxation discards alternating-profile gains.
Rbar_0 is the actual zero-input regional matrix, NOT an identity matrix.

The optional minimum_flags=9*r restriction uses the RS support event E:
each group has at least nine nonzero columns. Lift the actual distribution
P to flags using Q_eta(flags|output). The likelihood remains P/Q_eta.
Apply Holder to L*1_E and to the marked event restricted to at least9*r
active flags; E implies that flag event. The density norm below retains
E via its nonzero-column tail. The reference factor relaxes E to the
flag-count event. Consequently v^(-9*r) is outside alpha. This does not
assert independent flag counts after conditioning on E.

For a fixed original rowset of arity a and old U16 moment I_p, the new
zero/nonzero p-moment atoms are exactly
 z_eta=256^(-a*p)/(eta+(1-eta)*256^-16)^(p-1),
 n_eta=(1-eta)^(1-p)*(I_p-256^(p*(16-a)-16)).
The old nonzero subtraction is rational. Retain the actual sixteen
scheduled column moments before taking the activity tail and pth root.
No density renormalization or rowset resampling occurs.

This module supplies a reference factor and reusable density weights.
Its reference-only receipt is NOT a failure bound or a full certificate.
No new provider, encoder, or setup distribution is introduced.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from time import perf_counter

from flint import arb, ctx

import hd_lambda_renyi_activity as activity

DIRECTORY = Path(__file__).resolve().parent
A, binding = activity.A, activity.binding
R, ALPHA, MARK = 18, F(1, 5), F(691, 700)
ETAS = tuple(F(1, 2**k) for k in (64, 72, 80, 88))
ARITIES, ORDERS = tuple(range(11, 17)), (2, 3, 4)
PINS = {'hd_lambda_renyi_activity.py':
        '09f5245aecb903097438b65ded8757c6a66da40a13804289c394fef777186c77'}


def dependencies():
    for name, digest in PINS.items():
        if sha256((DIRECTORY/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Frozen dependency changed: '+name)
    return {**activity.dependencies(), **PINS,
            Path(__file__).name: sha256(Path(__file__).read_bytes()).hexdigest()}


def upper_fraction(value):
    """Exact dyadic upper endpoint; no fixed decimal rounding of tiny entries."""
    if not isinstance(value, arb) or not value.is_finite():
        raise ValueError('Finite Arb entry required')
    endpoint = value if value.is_exact() else value.upper()
    rational = endpoint.fmpq()
    result = F(int(rational.p), int(rational.q))
    if result < 0:
        raise ValueError('Nonnegative regional upper bound required')
    return result


def prepare_reference(regions, *, r=R, regions_per_column=16, columns=16):
    """Cache nineteen positive matrix powers; caller retains Arb precision.

    Smaller dimensions exist for tiny algebra tests, not a production
    construction claim. Production receipt fixes18 groups and16x16 sites.
    """
    if (type(r) is not int or not 1 <= r <= 18
            or type(regions_per_column) is not int or not 1 <= regions_per_column <= 16
            or type(columns) is not int or not 1 <= columns <= 16):
        raise ValueError('Bounded positive group/column dimensions required')
    envelopes, powers = [], []
    for m in range(r+1):
        rows = []
        for f in range(m//2+1):
            row = regions.get((m, f))
            if row is None or len(row) != 4:
                raise ValueError('Missing full regional matrix, including (0,0)')
            if any(not value.is_finite() or not value >= 0 for value in row):
                raise ValueError('Finite nonnegative Arb regional entries required')
            rows.append(tuple(upper_fraction(value) for value in row))
        envelope = tuple(max(row[i] for row in rows) for i in range(4))
        envelopes.append(envelope)
        powers.append(A.uniform.scalar_power(tuple(map(A.uniform.rational, envelope)),
                                             regions_per_column))
    return dict(r=r, regions_per_column=regions_per_column, columns=columns,
                regional_envelope_upper=tuple(envelopes), column_powers=tuple(powers))


def reference_from_envelope(envelopes, *, r=R, regions_per_column=16, columns=16):
    """Rehydrate exact saved upper entries without rebuilding a provider."""
    if len(envelopes) != r+1 or any(len(row) != 4 for row in envelopes):
        raise ValueError('One complete matrix per occupancy required')
    table = {(m, f): tuple(A.uniform.rational(F(x)) for x in envelopes[m])
             for m in range(r+1) for f in range(m//2+1)}
    return prepare_reference(table, r=r, regions_per_column=regions_per_column, columns=columns)


def column_matrix(prepared, eta, *, alpha=ALPHA, flag_mark=F(1)):
    """Fractional flag weights, then one ordered column matrix per flag set."""
    eta, alpha, flag_mark = F(eta), F(alpha), F(flag_mark)
    if not 0 <= eta <= 1 or not 0 < alpha <= 1 or flag_mark < 1:
        raise ValueError('Probability eta, alpha in(0,1], and flag mark>=1 required')
    r = prepared['r']
    result = [arb(0) for _ in range(4)]
    for m, power in enumerate(prepared['column_powers']):
        # Explicit endpoints retain 0^0=1 and exact absent branches.
        if (eta == 0 and m < r) or (eta == 1 and m > 0):
            continue
        coefficient = (arb(comb(r, m))*A.uniform.rational(eta)**A.uniform.rational(alpha*(r-m))
                       * A.uniform.rational(1-eta)**A.uniform.rational(alpha*m)
                       * A.uniform.rational(flag_mark)**m)
        for entry in range(4):
            result[entry] += coefficient*power[entry]
    return tuple(result)


def reference_log(prepared, eta, *, source_count=1, alpha=ALPHA,
                  flag_mark=F(1), minimum_flags=0):
    """Log2 upper expression, including source count ONCE but no threshold.

    minimum_flags>0 requires the support-event lifting described above.
    Set9*r for the production RS comparison; v=1 is the unrestricted
    baseline regardless of the floor. This function performs no fit.
    """
    flag_mark = F(flag_mark)
    if (type(source_count) is not int or source_count < 1
            or type(minimum_flags) is not int
            or not 0 <= minimum_flags <= prepared['r']*prepared['columns']):
        raise ValueError('Positive source count and feasible integer flag floor required')
    matrix = column_matrix(prepared, eta, alpha=alpha, flag_mark=flag_mark)
    word = A.uniform.scalar_power(matrix, prepared['columns'])
    value = word[0]+word[1]
    if not value > 0 or not value.is_finite():
        raise ArithmeticError('Positive finite initial-zero to terminal-ones value required')
    return A.log2(arb(source_count))+A.log2(value)-minimum_flags*A.log2(A.uniform.rational(flag_mark))


def reference_moment_atoms(arity, order, uniform_moment, eta, *, q=256, width=16):
    """Exact rational zero/nonzero atoms; old moment may itself be an upper bound."""
    eta, uniform_moment = F(eta), F(uniform_moment)
    if (type(q) is not int or q not in (4, 256) or type(width) is not int or not 1 <= width <= 16
            or type(arity) is not int or not 1 <= arity <= width
            or type(order) is not int or order not in ORDERS or not 0 <= eta < 1):
        raise ValueError('Valid field, dimension, order2/3/4 and eta<1 required')
    old_zero = F(q)**(order*(width-arity)-width)
    if uniform_moment < 1 or uniform_moment < old_zero:
        raise ValueError('Invalid normalized uniform-reference moment')
    reference_zero = eta+(1-eta)*F(q)**(-width)
    zero = F(q)**(-arity*order)/reference_zero**(order-1)
    nonzero = (1-eta)**(1-order)*(uniform_moment-old_zero)
    return zero, nonzero


def scheduled_weights(moments, order, eta, *, minimum=9, arities=ARITIES):
    """Restricted norm weights ONLY; labelled original rowsets are summed once.

    Input moments[(mapped complement tuple)][p] follows the checked frozen
    inventory API. The caller must validate inventory provenance before
    using these weights in a certificate. This function does not load files.
    The later bridge supplies beta^(total arity); these are not RAW weights.
    """
    eta = F(eta)
    if (type(order) is not int or order not in ORDERS or not 0 <= eta < 1
            or not arities or any(type(a) is not int or a not in ARITIES for a in arities)
            or len(set(arities)) != len(arities)
            or type(minimum) is not int or not 0 <= minimum <= 16):
        raise ValueError('Distinct arities11..16, order2/3/4, eta<1 and valid floor required')
    schedules = [binding.column_permutations('input-inverse-output-xor', j)[0] for j in range(16)]
    weights, diagnostics = [], []
    for arity in arities:
        signatures = Counter()
        for original in combinations(range(16), 16-arity):
            signature = tuple(sorted(F(moments[tuple(sorted(before[i] for i in original))][order])
                                     for before in schedules))
            signatures[signature] += 1
        total = arb(0)
        for signature, count in signatures.items():
            atoms = [reference_moment_atoms(arity, order, value, eta) for value in signature]
            zero = atoms[0][0]
            tail = activity.activity_tail(zero, tuple(n for _, n in atoms), minimum)
            total += count*tail**A.uniform.rational(F(1, order))
        if not total > 0 or not total.is_finite():
            raise ArithmeticError('Positive finite restricted norm required')
        weights.append(total)
        diagnostics.append(dict(arity=arity, labelled_rowsets=sum(signatures.values()),
                                distinct_scheduled_signatures=len(signatures)))
    return weights, diagnostics


def controls_from_regions(regions, *, flag_marks=(F(1),)):
    """Fixed r18 controls; no optimizer, new bank, or moment inventory."""
    flag_marks = tuple(map(F, flag_marks))
    if not flag_marks or any(v < 1 for v in flag_marks) or len(set(flag_marks)) != len(flag_marks):
        raise ValueError('Distinct supplied flag marks>=1 required')
    counts = A.profile_counts(R)
    if len(counts) != 99:
        raise ArithmeticError('Expected all99 indexed source profiles')
    source_count = sum(counts.values())
    prepared = prepare_reference(regions)
    old_terms = A.source_profile_terms(R, regions)
    uniform_log = A.log2(sum((term for _, _, _, term in old_terms), arb(0)))
    controls = []
    for eta in ETAS:
        for flag_mark in flag_marks:
            value = reference_log(prepared, eta, source_count=source_count,
                                  flag_mark=flag_mark, minimum_flags=9*R)
            controls.append(dict(eta=str(eta), flag_mark=str(flag_mark), minimum_active_flags=9*R,
                reference_source_log2_upper=str(A.uniform.rational_ceiling(value)),
                reference_source_log2_enclosure=str(value)))
    return dict(controls=controls, source_count=str(source_count), source_profiles=len(counts),
        uniform_comparison_log2_upper=str(A.uniform.rational_ceiling(uniform_log)),
        uniform_comparison_log2_enclosure=str(uniform_log),
        uniform_max_f_log2_upper=str(A.uniform.rational_ceiling(
            reference_log(prepared, F(0), source_count=source_count))),
        regional_max_f_upper=[[str(value) for value in row] for row in prepared['regional_envelope_upper']],
        regional_max_f_entries_are_exact_rational_upper_bounds=True)


def receipt(precision=256, *, flag_marks=(F(1),)):
    """ONE fresh uniform fractional bank; root owns any production invocation."""
    if type(precision) is not int or precision not in (256, 384):
        raise ValueError('Fresh256/384-bit precision required')
    pins, started = dependencies(), perf_counter()
    with ctx.workprec(precision):
        regions = A.uniform.regional_directed(MARK, ALPHA, R, 60)
        result = controls_from_regions(regions, flag_marks=flag_marks)
    if dependencies() != pins:
        raise ValueError('Frozen sources changed during reference control')
    return dict(result, construction=binding.SCOPE['construction'], reference='eta*delta0+(1-eta)*U16',
        active_groups=R, field_order=256, mark=str(MARK), alpha=str(ALPHA),
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
        arithmetic_enclosures_directed=True, reviewed_source_sha256=pins,
        elapsed_seconds=perf_counter()-started)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precision', type=int, choices=(256, 384), default=256)
    parser.add_argument('--flag-mark', type=F, action='append', dest='flag_marks')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Refuse to overwrite evidence')
    result = receipt(args.precision, flag_marks=args.flag_marks or (F(1),))
    with args.output.open('x', encoding='utf-8') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('reviewed_source_sha256', 'regional_max_f_upper')}))
