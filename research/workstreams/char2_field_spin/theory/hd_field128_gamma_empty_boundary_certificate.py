"""Directed positive-univariate replay for one zero-empty matching profile.

Only the ordinary top-degree slice of the frozen all-pair cap polynomial
is used. Its regional coefficient is [t**f] B(t)**64 / C(4096,f), with
f=r-4096. Ordered raw packet products precede entrywise alpha. Exact
positive tilting and polynomial reversal preserve that coefficient.
A complete receipt covers only(a,b)=(f,f), never all profiles at weight r.
"""

from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from time import monotonic

from flint import arb, arb_poly, ctx

import hd_field128_gamma_all_pair_certificate as provider
import hd_field128_gamma_midpoint_counts as midpoint
import hd_field128_gamma_profile_certificate as profiles


base, tail = provider.base, provider.tail
write_new = profiles.write_new


def source_pins():
    pins = {**provider.source_pins(), **midpoint.source_pins(), **profiles.source_pins()}
    path = Path(__file__).resolve()
    pins[path.name] = sha256(path.read_bytes()).hexdigest()
    return dict(sorted(pins.items()))


def check_source_pins(pins):
    directory = Path(__file__).resolve().parent
    for name, digest in pins.items():
        if sha256((directory/name).read_bytes()).hexdigest() != digest:
            raise RuntimeError('Pinned boundary source changed: '+name)


def top_coefficients(table, alpha, *, slots=32, budget=None):
    """Positive directed coefficients of B(t); generic slots serve tiny tests."""
    alpha = base.exact(alpha)
    if (type(slots) is not int or not 1 <= slots <= 32 or not 0 < alpha <= 1
            or set(table) != set(range(2*slots+1)) or any(len(row) != 4 for row in table.values())):
        raise ValueError('Expected complete packet caps, bounded slots and exact0<alpha<=1')
    budget = budget or base.Budget()
    rows = [tuple(base.upper(x) for x in table[j]) for j in range(slots, 2*slots+1)]
    result = [[arb(0) for unused in range(2*slots+1)] for entry in range(4)]
    exponent = base.rational(alpha)
    for i, left in enumerate(rows):
        budget.check()
        for j, right in enumerate(rows):
            count = comb(slots, i)*comb(slots, j)
            raw = base.matrix_product(left, right)
            for entry, value in enumerate(raw):
                if not value.is_zero():
                    result[entry][i+j] += count*base.upper(value**exponent)
    return tuple(tuple(base.upper(value) for value in row) for row in result)


def polynomial_product(left, right, maximum, budget):
    """Multiply chronologically, truncate exact high degrees, round upward."""
    budget.check()
    result = []
    for i in range(2):
        for j in range(2):
            value = (left[2*i]*right[j]+left[2*i+1]*right[2+j]).truncate(maximum+1)
            result.append(arb_poly([base.upper(value[k]) for k in range(min(len(value), maximum+1))]))
    budget.check()
    return tuple(result)


def regional_boundary(table, alpha, full, *, slots=32, paired_power=64,
                      marker=F(1, 3), reverse=True, budget=None):
    """Generic coefficient extraction; only the public wrapper fixes K20."""
    marker = base.exact(marker)
    if (type(paired_power) is not int or not 1 <= paired_power <= 128
            or type(full) is not int or not 0 <= full <= 2*slots*paired_power
            or marker <= 0 or type(reverse) is not bool):
        raise ValueError('Invalid exact coefficient, positive tilt or polynomial power')
    budget = budget or base.Budget()
    coefficients = top_coefficients(table, alpha, slots=slots, budget=budget)
    pairs = 2*slots*paired_power
    reversed_polynomial = reverse and full > pairs-full
    degree = pairs-full if reversed_polynomial else full
    used_marker = 1/marker if reversed_polynomial else marker
    if reversed_polynomial:
        coefficients = tuple(row[::-1] for row in coefficients)
    tilts = tuple(base.rational(used_marker**k) for k in range(2*slots+1))
    polynomial = tuple(arb_poly([base.upper(value*tilts[k])
                                for k, value in enumerate(row[:degree+1])]) for row in coefficients)
    result = tuple(arb_poly([int(i in (0, 3))]) for i in range(4))
    exponent = paired_power
    while exponent:
        budget.check()
        if exponent & 1:
            result = polynomial_product(result, polynomial, degree, budget)
        exponent //= 2
        if exponent:
            polynomial = polynomial_product(polynomial, polynomial, degree, budget)
    normalizer = comb(pairs, full)*used_marker**degree
    normalized = tuple(base.upper(row[degree]/base.rational(normalizer)) for row in result)
    return normalized, dict(polynomial_reversed=reversed_polynomial, requested_full_pairs=full,
        extracted_degree=degree, regional_slots=pairs, paired_polynomial_power=paired_power,
        original_positive_marker=str(marker), used_positive_marker=str(used_marker),
        exact_coefficient_normalizer=str(normalizer), histogram_normalizer=str(comb(pairs, full)),
        maximum_retained_degree=degree, coefficient_method='positive-Arb-polynomial-with-exact-truncation',
        top_component='ordinary total histogram degree2*slots, not weighted occupancy degree',
        chronological_product_before_alpha=True, counts_outside_alpha=True)


def profile_value(matrix, count, *, repetitions=256):
    if (len(matrix) != 4 or type(count) is not int or count < 0
            or type(repetitions) is not int or repetitions < 2 or repetitions % 2):
        raise ValueError('Invalid profile matrix, count or even repetition count')
    matrix = tuple(base.upper(x) for x in matrix)
    run = base.matrix_power(base.matrix_product(matrix, matrix), repetitions//2)
    transfer = base.upper(run[0]+run[1])
    return base.upper(count*transfer), transfer


def boundary_receipt(z, alpha, weight, *, state=30, distance=F(2, 5),
                     precision=256, bits=80, marker=F(1, 3), reverse=True,
                     max_seconds=180., failure_bits=54, witness_source=None):
    """Fresh directed packet caps and one complete boundary-profile bound."""
    z, alpha, distance, marker = map(base.exact, (z, alpha, distance, marker))
    construction = provider.local.construction(state)
    if (not 0 < z < 1 or not 0 < alpha <= 1 or not 0 < distance < F(1, 2)
            or marker <= 0 or type(reverse) is not bool or type(weight) is not int
            or not base.PAIRS <= weight <= base.GROUPS
            or type(precision) is not int or not 128 <= precision <= 4096
            or type(bits) is not int or not 32 <= bits <= 512
            or type(failure_bits) is not int or not 0 <= failure_bits <= 4096):
        raise ValueError('Invalid boundary construction, exact witness or arithmetic parameters')
    pins, budget = source_pins(), base.Budget(max_seconds)
    source_path = source_raw = None
    if witness_source is not None:
        source_path = Path(witness_source['path']).resolve()
        source_raw = source_path.read_bytes()
        if sha256(source_raw).hexdigest() != witness_source['sha256']:
            raise RuntimeError('Witness source hash mismatch')
        saved = json.loads(source_raw)
        witness = saved['witness']
        if F(witness['z']) != z or F(witness['alpha']) != alpha:
            raise ValueError('Source exact marks differ')
        for key, expected in (('construction', construction), ('state_dimension', state),
                              ('field_order', base.FIELD_ORDER), ('target_distance', str(distance))):
            actual = saved.get(key, witness.get(key, expected))
            if actual != expected:
                raise ValueError('Source construction differs: '+key)
    def check():
        check_source_pins(pins)
        if source_path is not None and source_path.read_bytes() != source_raw:
            raise RuntimeError('Witness source changed')
    check()
    full, forbidden = weight-base.PAIRS, int(distance*base.OUTPUT_LENGTH)
    metadata = dict(schema='gamma-empty-boundary-profile-directed-v1', construction=construction,
        field_order=base.FIELD_ORDER, state_dimension=state, inner_width=64, target_distance=str(distance),
        source_dimension=base.DIMENSION, output_length=base.OUTPUT_LENGTH, active_groups=weight,
        even_full=full, odd_full=full, single_pairs=2*base.PAIRS-weight, empty_pairs=0,
        precision_bits=precision, packet_mantissa_bits=bits, forbidden_weight=forbidden,
        numerical_backend='python-flint Arb/arb_poly with upward coefficient endpoints',
        packet_coefficient_producer='gamma_all_pair_certificate.frozen.sparse.packet_table',
        outer_measure='invertible-toeplitz', gamma_extension_degree=16,
        repetitions=256, paired_region_power=128, requested_failure_bits=failure_bits,
        source_sha256=pins, witness_source=witness_source,
        witness=dict(z=str(z), alpha=str(alpha), positive_coefficient_marker=str(marker)),
        old_numerical_bounds_reused=False, float_coefficients_imported=False,
        full_alias_bank_built=False, complete_selected_weight=False,
        covers_all_supports_at_weight=False, whole_code_certificate=False, max_seconds=max_seconds)
    failure = None
    numerical = {}
    try:
        with ctx.workprec(precision):
            table = provider.frozen.sparse.packet_table(z, 64, state=state, bits=bits, budget=budget)
            matrix, coefficient = regional_boundary(table, alpha, full, marker=marker,
                                                     reverse=reverse, budget=budget)
            counts = midpoint.MidpointCounts(base.PAIRS, weight, max_singletons=64, budget=budget)
            count = counts.count(full, full)
            support, transfer = profile_value(matrix, count)
            prefix = base.rational(alpha)*(weight*tail.log2arb(base.rational(base.outer_density(base.FIELD_ORDER)))
                       -tail.log2arb(arb(base.FIELD_ORDER-1))-forbidden*tail.log2arb(base.rational(z)))
            zero = support.is_zero()
            value = None if zero else prefix+tail.log2arb(support)
            ceiling = None if zero else base.rational_ceiling(value)
            provenance = counts.diagnostics()
            provenance['float_log_count_path_used'] = False
            numerical = dict(regional_matrix_upper=list(map(str, matrix)),
                coefficient_extraction=coefficient, count_upper=count, count_provenance=provenance,
                profile_transfer_upper=str(transfer), profile_support_upper=str(support),
                original_weight_outer_prefix_enclosure=str(prefix), zero=zero,
                log2_upper_rational=None if zero else str(ceiling),
                log2_upper_enclosure=None if zero else str(value),
                meets_requested_profile_target=zero or ceiling <= -failure_bits)
            budget.check()
    except (TimeoutError, ArithmeticError) as error:
        failure = dict(error_type=type(error).__name__, message=str(error))
    finally:
        check()
    complete = failure is None
    if not complete:
        numerical = dict(log2_upper_rational=None, log2_upper_enclosure=None,
                         meets_requested_profile_target=False)
    return dict(metadata, **numerical, status='complete' if complete else 'incomplete',
                directed=complete, complete_boundary_profile=complete,
                profile_pairs_evaluated=int(complete), failure=failure,
                source_pins_unchanged_at_exit=True, elapsed_seconds=monotonic()-budget.started,
                scope='Only the zero-empty matching profile(a,b)=(r-4096,r-4096); no whole-weight or code certificate')
